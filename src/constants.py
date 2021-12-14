AMINO_ACIDS={
    "A": 71.037114,
    "R": 156.101111,
    "N": 114.042927,
    "D": 115.026943,
    "C": 103.009185,
    "E": 129.042593,
    "Q": 128.058578,
    "G": 57.021464,
    "H": 137.058912,
    "I": 113.084064,
    "L": 113.084064,
    "K": 128.094963,
    "M": 131.040485,
    "F": 147.068414,
    "P": 97.052764,
    "S": 87.032028,
    "T": 101.047679,
    "U": 150.95363,
    "W": 186.079313,
    "Y": 163.06332,
    "V": 99.068414,
    "X": 0, # added to ignore. TODO: figure out what to do with it
    "B": 113.084064, # added to ignore. TODO: figure out what to do with it
    "Z": 0, # added to ignore. TODO: figure out what to do with it
}

# used for the minimum correct ordering of a sequence. Order is by mass
INTEGER_ORDERED_AMINO_ACIDS = {
    'X': 0,
    'Z': 0, 
    'G': 1,
    'A': 2, 
    'S': 3, 
    'P': 4, 
    'V': 5,
    'T': 6,
    'C': 7,
    'I': 8,
    'L': 8,
    'B': 8,
    'N': 9,
    'D': 10,
    'Q': 11,
    'K': 12,
    'E': 13,
    'M': 14,
    'H': 15,
    'F': 16,
    'U': 17,
    'R': 18,
    'Y': 19,
    'W': 20
}

#This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
WATER_MASS = 2 * 1.007825035 + 15.99491463 
PROTON_MASS = 1.007276466879

SINGLY_CHARGED_Y_BASE = 3 * 1.007825035 + 15.99491463 - 0.0005486 #for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
DOUBLY_CHARGED_Y_BASE = 4 * 1.007825035 + 15.99491463 - 2 * 0.0005486 #another proton to make doubly charged

SINGLY_CHARGED_B_BASE = 1.007825035 - 0.0005486 #for the H to turn the residue NH on the N-terminus into NH2
DOUBLY_CHARGED_B_BASE = 2 * 1.007825035 - 2 * 0.0005486 #adding one more proton this time to make it doubly charged
