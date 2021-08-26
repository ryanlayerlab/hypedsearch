from constants import AMINO_ACIDS, WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, INTEGER_ORDERED_AMINO_ACIDS, PROTON_MASS
import numpy as np

def b_ions(sequence: str, charge: int = None): 
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

def y_ions(sequence: str, charge: int = None): 
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
        masses += b_ions(sequence, charge=charge)
    if ion is None or ion == 'y': 
        masses += y_ions(sequence, charge=charge)
    return masses, pre_mz

def max_mass(seqeunce: str, ion: str, charge: int):
    if ion == 'y':
        total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])
        mz = total / charge
        return mz
    if ion == 'b':
        total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])
    mz = total / charge
    return mz

def get_precursor(sequence: str, charge: int = 1):
    total = WATER_MASS
    for aa in sequence:
        total +=  AMINO_ACIDS[aa]
    return (total + charge * PROTON_MASS) / charge 

def gen_spectrum(sequence: str, charge: int = None, ion: str = None):
    this_entry = {}
    masses, pre_mz = calc_masses(sequence, charge=charge, ion=ion)
    this_entry['spectrum'] = masses
    this_entry['precursor_mass'] = pre_mz
    return this_entry

def gen_spectra(sequences: list, charge=None, ion=None):
    return [gen_spectrum(seq, charge=charge, ion=ion) for seq in sequences]

def gen_min_ordering(sequence: str):
    middle = [np.int8(INTEGER_ORDERED_AMINO_ACIDS[aa]) for aa in sequence if aa in INTEGER_ORDERED_AMINO_ACIDS]
    if len(middle) == 0:
        return []
    return middle

   