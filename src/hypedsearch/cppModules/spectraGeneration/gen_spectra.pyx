# distutils: language = c++

from libcpp.string cimport string 
from libcpp.vector cimport vector
from libcpp cimport bool

cimport genSpectra

def gen_spectrum(sequence: str, ion='', charge=-1, sort=True) -> list:
    '''
    Generate the theoretical spectrum for a sequence

    Inputs:
        sequence:   (str) the peptide (AA sequence) to convert
    kwargs:
        ion:        (str) the ion type to generate. Leave blank to generate for both 
                            b and y ions
        charge:     (int) the charge to generate the spectrum for. Leave blank for both
                            singly and doubly
        sort:       (bool) to sort the resulting spectrum. Default=True
    Outputs:
        (list) the floats of the spectrum
    '''
    if len(sequence) == 0:
        return []

    seq = str.encode(sequence)
    i = str.encode(ion)

    spec = genSpectra.genSpectrum(seq, i, charge, sort)
    return [x for x in spec]


