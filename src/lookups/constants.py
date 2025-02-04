# Amino acid masses
# Source for dictionary below is the code here https://www.rapidnovor.com/mxie-customized/massCal.js
# which is used by this website to compute masses: https://www.rapidnovor.com/mass-calculator/
AMINO_ACID_MASSES = {
    "A": 71.0371138,
    "R": 156.1011110,
    "N": 114.0429274,
    "D": 115.0269430,
    "C": 103.0091845,
    "E": 129.0425931,
    "Q": 128.0585775,
    "G": 57.0214637,
    "H": 137.0589119,
    "I": 113.0840640,
    "L": 113.0840640,
    "K": 128.0949630,
    "M": 131.0404846,
    "F": 147.0684139,
    "P": 97.0527639,
    "S": 87.0320284,
    "T": 101.0476785,
    "W": 186.0793130,
    "Y": 163.0633285,
    "U": 168.9641990,
    "V": 99.0684139,
}

# Sources for dictionary below are:
# - https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
# - https://proteomicsresource.washington.edu/protocols06/masses.php
# AMINO_ACID_MASSES = {
#     "A": 71.037114,
#     "R": 156.101111,
#     "N": 114.042927,
#     "D": 115.026943,
#     "C": 103.009185,
#     "E": 129.042593,
#     "Q": 128.058578,
#     "G": 57.021464,
#     "H": 137.058912,
#     "I": 113.084064,
#     "L": 113.084064,
#     "K": 128.094963,
#     "M": 131.040485,
#     "F": 147.068414,
#     "P": 97.052764,
#     "S": 87.032028,
#     "T": 101.047679,
#     "U": 150.95363,
#     "W": 186.079313,
#     "Y": 163.06332,
#     "V": 99.068414,
#     "O": 237.147727,
# }

# Other masses
HYDROGEN_MASS = 1.007825035
WATER_MASS = 18.010564686  # from https://www.rapidnovor.com/mxie-customized/massCal.js
OXYGEN_MASS = 15.99491463
PROTON_MASS = 1.00727646688  # from Scott
# This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
WATER_MASS = (2 * HYDROGEN_MASS) + OXYGEN_MASS


# What are these doing
SINGLY_CHARGED_Y_BASE = (
    3 * 1.007825035 + 15.99491463 - 0.0005486
)  # for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
DOUBLY_CHARGED_Y_BASE = (
    4 * 1.007825035 + 15.99491463 - 2 * 0.0005486
)  # another proton to make doubly charged
SINGLY_CHARGED_B_BASE = (
    1.007825035 - 0.0005486
)  # for the H to turn the residue NH on the N-terminus into NH2
DOUBLY_CHARGED_B_BASE = (
    2 * 1.007825035 - 2 * 0.0005486
)  # adding one more proton this time to make it doubly charged
AMMONIUM = 3 * 1.00782503517 + 14.003074  # NH3 weight
