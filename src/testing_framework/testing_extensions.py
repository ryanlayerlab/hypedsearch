import os
import sys
from typing import Tuple
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)


import database
import testing_utils

from constants import WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE
from utils import ppm_to_da
from preprocessing import merge_search, preprocessing_utils
from collections import defaultdict

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20


datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)

# print(input_spectra[167][0])
# print(input_spectra[390][0])
# print(input_spectra[894][0])
# print(input_spectra[958][0])
# for spectrum in input_spectra:
#     for mz in spectrum[0]:
#         tolerance = ppm_to_da(mz, ppm_tolerance)
#         if (57.021464 + WATER_MASS + SINGLY_CHARGED_Y_BASE > mz - tolerance) and (57.021464 + WATER_MASS + SINGLY_CHARGED_Y_BASE < mz + tolerance):
#             print('Found a G at: ', mz)
#         if (57.021464 + WATER_MASS + SINGLY_CHARGED_B_BASE > mz - tolerance) and (57.021464 + WATER_MASS + SINGLY_CHARGED_B_BASE < mz + tolerance):
#             print('Found a G at: ', mz)


correct_sequences = testing_utils.generate_truth_set(datasets[0])
correct_sequence = correct_sequences[0]

# for i, sequence in enumerate(correct_sequences):
#     if sequence[0] == 'G':
#         print(i, sequence)

#Generate all subsequences of each protein
db = database.build(dataset[2])
matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)

print('Starting mapping...')
mz_mapping = defaultdict(set)
for spectrum in input_spectra:
    for i, mz in enumerate(spectrum[0]):
        boundaries = preprocessing_utils.make_boundaries(mz, ppm_tolerance)
        for boundary in matched_masses_b.keys():
            if boundary == str(boundaries[0]) + '-' + str(boundaries[1]):
                for matching in matched_masses_b[boundary]:
                    mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))
        for boundary in matched_masses_y.keys():
            if boundary == str(boundaries[0]) + '-' + str(boundaries[1]):
                for matching in matched_masses_b[boundary]:
                    mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))

print('Done')
mz_mapping = sorted(mz_mapping)
print(mz_mapping)




print('I want dinner...')








#Length is 3867 for b and 3793 for y after changes
#Length is 4000 for b and 3927 for y after changes

#Map to (P_y, S_i, P_j, seq, b/y)
# Where P_y is protein this was found in, S_i is m/z number, P_j is location within that protein, seq and b/y are straightforward

# Find initial hits

# for i, mz in enumerate(input_spectrum):
#     #Map to (P_y, S_i, P_j, seq, b/y)
#     # Where P_y is protein this was found in, S_i is m/z number, P_j is location within that protein, seq and b/y are straightforward
#     for prot in proteins
