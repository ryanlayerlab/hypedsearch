import os
import sys
from typing import Tuple

module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)


import database, gen_spectra
import testing_utils

from constants import WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE
from utils import ppm_to_da
from preprocessing import merge_search, preprocessing_utils
from collections import defaultdict
from objects import Spectrum

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20

import matplotlib.pyplot as plt

def collect_hit_or_miss_abundances(input_spectrum, mz_hit_set, correct_hits):
    hit_tuples = []
    miss_tuples = []
    correct_hit_tuples = []
    hit_list = []
    
    [hit_list.append(x[1]) for x in mz_hit_set]
    for i, x in enumerate(input_spectrum[0]):
        if x in hit_list:
            hit_tuples.append((x, input_spectrum[1][i]))
            if x in correct_hits:
                correct_hit_tuples.append((x, input_spectrum[1][i]))
        else:
            miss_tuples.append((x, input_spectrum[1][i]))

    # pick apart tuples
    hit_list = []
    hit_abundances = []
    miss_list = []
    miss_abundances = []
    correct_hit_abundances = []
    c_hits = []
    
    for x in hit_tuples:
        hit_list.append(x[0])
        hit_abundances.append(x[1])
        
    for x in correct_hit_tuples:
        c_hits.append(x[0])
        correct_hit_abundances.append(x[1])
        
    for x in miss_tuples:
        miss_list.append(x[0])
        miss_abundances.append(x[1])
        
    return hit_list, hit_abundances, miss_list, miss_abundances, correct_hit_abundances, c_hits

datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


correct_sequences = testing_utils.generate_truth_set(datasets[0])

# path = '/home/ncol107453/jaime_hypedsearch/hypedsearch/data/database/single_prot.fasta'
path = dataset[2]
db = database.build(path)

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
    correct_sequence = correct_sequences[i]
    testing_utils.find_hits(mz_mapping, boundaries, spectrum, i, matched_masses_b, matched_masses_y, b_hits, y_hits, b_set, y_set, mz_miss_set)
    testing_utils.append_correct_hits(correct_hits, correct_sequence, spectrum, ppm_tolerance)
print('Done')

#Writing b and y hits
# print('Writing data...')
# with open("b_hits.txt", 'w') as b:
#     for x in b_hits:
#         pep_id = x[0]
#         w = x[1][0]
#         prot_id = x[1][1]
#         seq = x[1][2]
#         loc = x[1][3]
#         ion = x[1][4]
#         charge = x[1][5]
#         out = [pep_id, w, prot_id, seq, loc, ion, charge]
#         b.write('\t'.join([str(i) for i in out]) + '\n')
# with open("y_hits.txt", 'w') as y_file:
#     for y in y_hits:
#         pep_id = y[0]
#         w = y[1][0]
#         prot_id = y[1][1]
#         seq = y[1][2]
#         loc = y[1][3]
#         ion = y[1][4]
#         charge = y[1][5]
#         out = [pep_id, w, prot_id, seq, loc, ion, charge]
#         y_file.write('\t'.join([str(i) for i in out]) + '\n')
# print('Done')

# Calculating total length
total_len = 0
for spectrum in input_spectra:
    total_len = total_len + len(spectrum[0])
mz_hit_set = b_set | y_set
    
print('Total number of hits:', len(mz_hit_set), 'out of', total_len, 'total mz values', '(' + str(round((len(mz_hit_set) / total_len * 100))) + '%)')
print('Total number of misses:', len(mz_miss_set), 'out of', total_len, 'total mz values', '(' + str(round((len(mz_miss_set) / total_len * 100))) + '%)')
print('Num b hits:', len(b_set))
print('Num y hits:', len(y_set))
print(len(mz_hit_set), len(mz_miss_set), len(mz_hit_set) + len(mz_miss_set))