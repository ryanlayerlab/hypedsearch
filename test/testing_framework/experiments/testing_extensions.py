import os
import sys
from typing import Tuple
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)


module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)


import database, gen_spectra
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


correct_sequences = testing_utils.generate_truth_set(datasets[0])

db = database.build(dataset[2])

matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)
print('Finished matching masses')

print('Starting mapping...')
b_hits =[]
y_hits = []
b_set = set()
y_set = set()
correct_hits = []
mz_miss_set = set()
for i, spectrum in enumerate(input_spectra):
    #Remember to add in abundance if it is helpful
    input_num = i+1
    correct_sequence = correct_sequences[i]
    testing_utils.find_hits(mz_mapping, boundaries, spectrum, input_num, matched_masses_b, matched_masses_y, b_hits, y_hits, b_set, y_set, mz_miss_set)
    testing_utils.append_correct_hits(correct_hits, correct_sequence, spectrum, ppm_tolerance)
print('Done')

print('Writing data...')
all_hits = b_set | y_set
with open('all_hits.txt', 'w') as w:
    [w.write(str(x) + '\n') for x in all_hits]
print('Done')

print(len(all_hits))
print(len(mz_miss_set))