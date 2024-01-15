
from lookups.constants import AMINO_ACIDS, WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, INTEGER_ORDERED_AMINO_ACIDS, PROTON_MASS, AMMONIUM
import numpy as np

def get_b_ions(sequence: str, charge: int = None): 
    masses = []
    length = len(sequence)
    if charge is None or charge == 1:
        total = SINGLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total)
    if charge is None or charge == 2:
        total = DOUBLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total/2)
    return masses

def get_y_ions(sequence: str, charge: int = None): 
    masses = []
    length = len(sequence)
    if charge is None or charge == 1:
        total = SINGLY_CHARGED_Y_BASE
        for i in range (0,length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total)
    if charge is None or charge == 2:
        total = DOUBLY_CHARGED_Y_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total/2)
            
    return masses

def calc_masses(sequence: str, charge: int =None, ion: str = None):
    masses = []
    length = len(sequence)
    total = WATER_MASS
    for i in range(length):
        total +=  AMINO_ACIDS[sequence[i]]
    pre_mz_charge = 2 if charge is None else charge
    pre_mz = (total+pre_mz_charge*PROTON_MASS)/pre_mz_charge
    if ion is None or ion == 'b': 
        masses += get_b_ions(sequence, charge=charge)
    if ion is None or ion == 'y': 
        masses += get_y_ions(sequence, charge=charge)
    return masses, pre_mz

def calc_masses_no_water(sequence: str, charge: int =None, ion: str = None):
    masses = []
    length = len(sequence)
    total = 0
    for i in range(length):
        total +=  AMINO_ACIDS[sequence[i]]
    pre_mz_charge = 2 if charge is None else charge
    pre_mz = (total+pre_mz_charge*PROTON_MASS)/pre_mz_charge
    if ion is None or ion == 'b': 
        masses += get_b_ions(sequence, charge=charge)
    if ion is None or ion == 'y': 
        masses += get_y_ions(sequence, charge=charge)
    return masses, pre_mz

def calc_masses_no_ammonium(sequence: str, charge: int =None, ion: str = None):
    masses = []
    length = len(sequence)
    total = WATER_MASS
    for i in range(length):
        total +=  AMINO_ACIDS[sequence[i]]
    total -= AMMONIUM
    pre_mz_charge = 2 if charge is None else charge
    pre_mz = (total+pre_mz_charge*PROTON_MASS)/pre_mz_charge
    if ion is None or ion == 'b': 
        masses += get_b_ions(sequence, charge=charge)
    if ion is None or ion == 'y': 
        masses += get_y_ions(sequence, charge=charge)
    return masses, pre_mz

def get_max_mass(seqeunce: str, ion: str, charge: int):
    if ion == 'y':
        total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])
        mz = total / charge
        return mz
    else:
        total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])
    mz = total / charge
    return mz

def get_total_sum(precursor_mass,precursor_charge):
    total_sum = (precursor_mass * precursor_charge) - (precursor_charge * PROTON_MASS) - WATER_MASS
    return total_sum    

def convert_precursor_to_ion(precursor_mass, precursor_charge):
    total_sum = get_total_sum(precursor_mass, precursor_charge)
    normalized_b = (DOUBLY_CHARGED_B_BASE + total_sum) / 2
    normalized_y = (DOUBLY_CHARGED_Y_BASE + total_sum) / 2
    return normalized_b, normalized_y

def normalize_to_two(b_mass, y_mass):
    normalized_b = (DOUBLY_CHARGED_B_BASE + b_mass) / 2
    normalized_y = (DOUBLY_CHARGED_Y_BASE + y_mass) / 2
    return normalized_b, normalized_y
    
def get_sum_from_ion(b_mass, y_mass, b_charge, y_charge):
    b_sum, y_sum = get_raw_mass(b_mass, 0, b_charge), get_raw_mass(y_mass, 1, y_charge)
    return b_sum, y_sum
    
def get_converted_missing_mass(precursor_mass, precursor_charge, b_mass, y_mass, b_charge, y_charge):
    precursor_total_sum = get_total_sum(precursor_mass, precursor_charge)
    b_sum, y_sum = get_sum_from_ion(b_mass, y_mass, b_charge, y_charge)
    b_missing_sum = precursor_total_sum - b_sum
    y_missing_sum = precursor_total_sum - y_sum
    normalized_b, normalized_y = normalize_to_two(b_missing_sum, y_missing_sum)
    return normalized_b, normalized_y
    
def get_precursor(sequence: str, charge: int = 1):
    total = WATER_MASS
    for aa in sequence:
        total +=  AMINO_ACIDS[aa]
    return (total + charge * PROTON_MASS) / charge

def get_raw_mass(mass, ion, charge):
    if ion == 0:
        if charge == 1:
            raw_mass = (mass * charge) - SINGLY_CHARGED_B_BASE
        else:
            raw_mass = (mass * charge) - DOUBLY_CHARGED_B_BASE
    else:
        if charge == 1:
            raw_mass = (mass * charge) - SINGLY_CHARGED_Y_BASE
        else:
            raw_mass = (mass * charge) - DOUBLY_CHARGED_Y_BASE
        
    return raw_mass    

def convert_ion_to_precursor(mass, ion, charge, prec_charge):
    total = get_raw_mass(mass, ion, charge)
    return (total + WATER_MASS + (prec_charge * PROTON_MASS)) / prec_charge

def convert_raw_to_precursor(total, charge):
    return (total + WATER_MASS + (charge * PROTON_MASS)) / charge

def calc_precursor_as_disjoint(b_mass, y_mass, b_charge, y_charge, precursor_charge):
    b_sequence_mass = get_raw_mass(b_mass, 0, b_charge)
    y_sequence_mass = get_raw_mass(y_mass, 1, y_charge)
    total_sequence_mass = b_sequence_mass + y_sequence_mass
    precursor = convert_raw_to_precursor(total_sequence_mass, precursor_charge)
    return precursor

def generate_spectrum(sequence: str, charge: int = None, ion: str = None):
    this_entry = {}
    masses, pre_mz = calc_masses(sequence, charge=charge, ion=ion)
    this_entry['spectrum'] = masses
    this_entry['precursor_mass'] = pre_mz
    return this_entry

def generate_spectra(sequences: list, charge=None, ion=None):
    return [generate_spectrum(seq, charge=charge, ion=ion) for seq in sequences]

def generate_min_ordering(sequence: str):
    middle = [np.int8(INTEGER_ORDERED_AMINO_ACIDS[aa]) for aa in sequence if aa in INTEGER_ORDERED_AMINO_ACIDS]
    if len(middle) == 0:
        return []
    return middle

def calc_combined_mass(total_raw, target_ion):
    if target_ion == 0:
        converted = (DOUBLY_CHARGED_B_BASE + total_raw) / 2
    else:
        converted = (DOUBLY_CHARGED_Y_BASE + total_raw) / 2
    return converted