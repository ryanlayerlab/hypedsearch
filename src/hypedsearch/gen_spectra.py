from constants import AMINO_ACIDS, WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, INTEGER_ORDERED_AMINO_ACIDS, PROTON_MASS

import numpy as np

def b_ions(sequence: str, charge: int = None): 
    '''calculate b ion masses for a sequence for a given charge(s)
    
    :param sequence: amino acid sequence to calculate the b ion masses for
    :type sequnce: str
    :param charge: the charge of the b ions to calculate the mass for. If left
        as None, both singly and doubly charged b ions will be calculated. 
        (default is None)
    :type charge: int

    :returns: the b ion masses for the input sequence
    :rtype: list
    '''

    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #b+
        total = SINGLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total)
            #Since z (the charge) is equal to one, the total here is the m/z
            
    if charge is None or charge == 2:
        #b++
        total = DOUBLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total/2)
            
    return masses

def y_ions(sequence: str, charge: int = None): 
    '''calculate y ion masses for a sequence for a given charge(s)
    
    :param sequence: amino acid sequence to calculate the y ion masses for
    :type sequnce: str
    :param charge: the charge of the y ions to calculate the mass for. If left
        as None, both singly and doubly charged y ions will be calculated. 
        (default is None)
    :type charge: int

    :returns: the y ion masses for the input sequence
    :rtype: list
    '''

    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #y+
        total = SINGLY_CHARGED_Y_BASE
        for i in range (0,length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total)
            
    if charge is None or charge == 2:
        #y++
        total = DOUBLY_CHARGED_Y_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total/2)
            
    return masses

def calc_masses(sequence: str, charge: int =None, ion: str = None) -> (list, float):
    '''Calculate the molecular weight (Daltons) of an Amino Acid sequence
    
    :param sequence: amino acid sequence to calculate ion masses for
    :type sequence: str
    :param charge: the charge of the ions to calculate the mass for. If left
        as None, both singly and doubly charged ions will be calculated. 
        (default is None)
    :type charge: int
    :param ion: ion type to calculate masses for. Values are ['b', 'y']. If set 
        to None, both 'b' and 'y' ions are calculated.
        (default is None)
    :type ion: str

    :returns: the first return value is the list of masses calculated for the 
        input sequence in no order. the second return value is the calculated
        precursor mass of the sequence
    :rtype: list, float
    '''

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

def max_mass(seqeunce: str, ion: str, charge: int) -> float:
    '''Calculate the maximum mass of a sequence of an ion type and charge

    :param sequence: the sequence to generate the max mass for
    :type sequence: str
    :param ion: the ion type for which we calculate the mass. Options: 'b', 'y'
    :type ion: str
    :param charge: the charge to calculate the mass for. Options are: [1, 2]
    :type charge: int

    :returns: the maximum mass
    :rtype: float
    '''

    # all y ions
    if ion == 'y':
        total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])

        # divide by charge
        mz = total / charge
        return mz

    # otherwise do the b
    if ion == 'b':
        total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])

    # divide by charge
    mz = total / charge
    return mz

def get_precursor(sequence: str, charge: int = 1) -> float:
    '''Calculate JUST the precursor mass of the input sequence at the charge provided.

    :param sequence: the amino acid sequence to calculate the precursor of
    :type sequence: str
    :param charge: the charge for which to calculate the precursor mass. 
        (default is 1)
    :type charge: int

    :returns: the percursor mass of the sequence
    :rtype: float
    '''

    total = WATER_MASS
    for aa in sequence:
        total +=  AMINO_ACIDS[aa]
    # proton mass is 1.00...
    return (total + charge * PROTON_MASS) / charge  


def gen_spectrum(sequence: str, charge: int = None, ion: str = None) -> dict:
    '''Generate a spectrum for a single sequence
    
    :param sequence: amino acid sequence to calculate ion masses for
    :type sequence: str
    :param charge: the charge of the ions to calculate the mass for. If left
        as None, both singly and doubly charged ions will be calculated. 
        (default is None)
    :type charge: int
    :param ion: ion type to calculate masses for. Values are ['b', 'y']. If set 
        to None, both 'b' and 'y' ions are calculated.
        (default is None)
    :type ion: str

    :returns: a dictionary with the spectrum and precursor mass in the form 
        {'spectrum': list, 'precursor_mass: float}
    :rtype: dict
    '''
    
    this_entry = {}
    masses, pre_mz = calc_masses(sequence, charge=charge, ion=ion)
    this_entry['spectrum'] = masses
    this_entry['precursor_mass'] = pre_mz
    return this_entry

def gen_spectra(sequences: list, charge=None, ion=None) -> list:
    '''Generates mass spectra for a list of sequences

    :param sequences: amino acid sequences to calculate ion masses for
    :type sequences: list
    :param charge: the charge of the ions to calculate the mass for. If left
        as None, both singly and doubly charged ions will be calculated. 
        (default is None)
    :type charge: int
    :param ion: ion type to calculate masses for. Values are ['b', 'y']. If set 
        to None, both 'b' and 'y' ions are calculated.
        (default is None)
    :type ion: str

    :returns: dictionaries of {'spectrum': list, 'precursor_mass': float} for 
        all sequences in the order they were passed in
    :rtype: list
    '''

    return [gen_spectrum(seq, charge=charge, ion=ion) for seq in sequences]

def gen_min_ordering(sequence: str) -> list:
    '''Generates an np array the length of the sequence that is the minimal representation
    of a spectrum (for ordering purposes). Each amino acid is represented as an integer
    (an 8 bit integer by NumPy). The integer is the sorted value (lowest to highest) by mass. 
    For example, G has the lowest mass, so its index is 0. W, the heaviest, has an index of 19.
    This is done for the smallest (memory) representation for sorting a list of k-mers
    (needed for simple DAWG construction)

    :param sequence: the sequence of amino acids to convert
    :type sequence: str
    
    :returns: a list of 16 bit integers
    :rtype: list
    '''

    middle = [np.int8(INTEGER_ORDERED_AMINO_ACIDS[aa]) for aa in sequence if aa in INTEGER_ORDERED_AMINO_ACIDS]
    if len(middle) == 0:
        return []
    return middle

   