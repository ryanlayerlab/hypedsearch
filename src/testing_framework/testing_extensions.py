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

input_spectrum = input_spectra[0]

# for i, sequence in enumerate(correct_sequences):
#     if sequence[0] == 'G':
#         print(i, sequence)

#Generate all subsequences of each protein
db = database.build(dataset[2])
matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)
matched_b_list = []
matched_y_list = []
for key in matched_masses_b.keys():
    matched_b_list.append((key, matched_masses_b[key]))
for key in matched_masses_y.keys():
    matched_y_list.append((key, matched_masses_y[key]))
with open('matched_masses_b.txt', 'w') as b:
    [b.write(str(x[0]) + ': ' + str(x[1]) + '\n') for x in matched_b_list]
with open('matched_masses_y.txt', 'w') as y:
    [y.write(str(x[0]) + ': ' + str(x[1]) + '\n') for x in matched_y_list]
print('Starting mapping...')

mz_mapping = defaultdict(set)
testing_utils.find_matches_in_spectrum(input_spectrum, mz_mapping, ppm_tolerance, matched_masses_b, matched_masses_y)
print('Done')

#Sorting dict
print('Writing data...')
mapping_list = []
for key in sorted(mz_mapping.keys()):
    mapping_list.append(mz_mapping[key])

with open('metadata.txt', 'w') as m:
    [m.write(str(x) + '\n') for x in mapping_list]
print('Done')