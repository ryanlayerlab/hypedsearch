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
from preprocessing import merge_search, preprocessing_utils
from collections import defaultdict
from objects import Spectrum

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20

datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


correct_sequences = testing_utils.generate_truth_set(datasets[0])

path = dataset[2]
db = database.build(path)

matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)
print('Finished matching masses')

b_hits =[]
y_hits = []
b_set = set()
y_set = set()
mz_miss_set = set()

for spectrum_num in range(0,1086):
    correct_sequence = correct_sequences[spectrum_num]

    input_spectrum = input_spectra[spectrum_num]

    #Remember to add in abundance if it is helpful
    testing_utils.find_hits(mz_mapping, boundaries, input_spectrum, spectrum_num, matched_masses_b, matched_masses_y, b_hits, y_hits, b_set, y_set, mz_miss_set)
    print('Done')

print(len(b_hits))
print(len(y_hits))
with open("b_hits.txt", 'w') as b:
    for x in b_hits:
        pep_id = x[0]
        w = x[1][0]
        prot_id = x[1][1]
        seq = x[1][2]
        loc = x[1][3]
        ion = x[1][4]
        charge = x[1][5]
        out = [pep_id, w, prot_id, seq, loc, ion, charge]
        b.write('\t'.join([str(i) for i in out]) + '\n')
with open("y_hits.txt", 'w') as y_file:
    for y in y_hits:
        pep_id = y[0]
        w = y[1][0]
        prot_id = y[1][1]
        seq = y[1][2]
        loc = y[1][3]
        ion = y[1][4]
        charge = y[1][5]
        out = [pep_id, w, prot_id, seq, loc, ion, charge]
        y_file.write('\t'.join([str(i) for i in out]) + '\n')
print('Done')