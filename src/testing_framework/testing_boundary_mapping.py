import random
import string

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
from utils import ppm_to_da, overlap_intervals
from preprocessing import merge_search, preprocessing_utils
from collections import defaultdict
from objects import Spectrum

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20




def define_single_spectrum(mz_list, ppm_tol):
    boundaries = [make_boundaries(mz, ppm_tol) for mz in mz_list]

    # make overlapped boundaries larger boundaries
    boundaries = overlap_intervals(boundaries)

    # make a mapping for mz -> boundaries
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(mz_list):
        
        # if the current boundary encapsulates s_i, add to list
        if boundaries[b_i][0] <= mz_list[s_i] <= boundaries[b_i][1]:
            mz_mapping[mz_list[s_i]] = b_i 
            s_i += 1

        # else if the s_i < boundary, increment s_i
        elif mz_list[s_i] < boundaries[b_i][0]:
            s_i += 1

        # else if s_i > boundary, incrment b_i
        elif mz_list[s_i] > boundaries[b_i][1]:
            b_i += 1
    return boundaries, mz_mapping

def make_boundaries(mz, ppm_tol):
    da_tol = ppm_to_da(mz, ppm_tol)
    return [mz - da_tol, mz + da_tol]

datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


correct_sequences = testing_utils.generate_truth_set(datasets[0])

# path = '/home/ncol107453/jaime_hypedsearch/hypedsearch/data/database/single_prot.fasta'
path = dataset[2]
db = database.build(path)






#Determined sequence
correct_sequence = 'MSSP'
ion = 'b'
charge = 1
mz_array = gen_spectra.gen_spectrum(correct_sequence, 1, ion)['spectrum']

#Random sequence
# correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,6))
# print(correct_sequence)
# mz_array = gen_spectra.gen_spectrum(correct_sequence, 1, ion)['spectrum']

abundances = []
[abundances.append(random.random() * 10 * random.randrange(0,10000)) for x in range(0,len(mz_array))]

input_spectrum = Spectrum(mz_array, abundances)

mz_array.sort()


#Real data
# spectrum_num = 0
# input_spectrum = input_spectra[spectrum_num]
# correct_sequence = correct_sequences[spectrum_num]

boundaries, mz_mapping = define_single_spectrum(input_spectrum[0], ppm_tolerance)

matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)









print('Starting mapping...')

b_hits =[]
y_hits = []
b_set = set()
y_set = set()
correct_hits = []
mz_miss_set = set()
#Remember to add in abundance if it is helpful
input_num = 1
testing_utils.find_hits(mz_mapping, boundaries, input_spectrum, input_num, matched_masses_b, matched_masses_y, b_hits, y_hits, b_set, y_set, mz_miss_set)
testing_utils.append_correct_hits(correct_hits, correct_sequence, input_spectrum, ppm_tolerance)
print('Done')

#Writing b and y hits
print('Writing data...')
with open("b_hits.txt", 'w') as b:
    for x in b_hits:
        pep_id = x[0]
        w = x[1][0]
        prot_id = x[1][1]
        seq = x[1][2]
        loc = x[1][3]
        write_ion = x[1][4]
        charge = x[1][5]
        out = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        b.write('\t'.join([str(i) for i in out]) + '\n')
with open("y_hits.txt", 'w') as y_file:
    for y in y_hits:
        pep_id = y[0]
        w = y[1][0]
        prot_id = y[1][1]
        seq = y[1][2]
        loc = y[1][3]
        write_ion = y[1][4]
        charge = y[1][5]
        out = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        y_file.write('\t'.join([str(i) for i in out]) + '\n')
print('Done')

# with open("mz_ab.txt", 'w') as m:
#     mz_list = input_spectrum[0]
#     abundance_list = input_spectrum[1]
#     for i,mz in enumerate(input_spectrum[0]):
#         abundance = abundance_list[i]
#         out = [mz, abundance]
#         m.write('\t'.join([str(i) for i in out]) + '\n')


# What we expect: [4, 0, 1, 4, MSSP]
intervals = testing_utils.map_hits_to_intervals(ion)
grouped_intervals = testing_utils.group_intervals(intervals)
scored_intervals = testing_utils.merge_intervals(grouped_intervals)
scored_intervals.sort(key=lambda x: x[3], reverse=True)
with open("scored_intervals.txt", 'w') as s:
    for interval in scored_intervals:
        # parent_prot = interval[0]
        # start_pos = interval[1]
        # end_pos = interval[2]
        # score = interval[3]
        # seq = interval[4]
        # out = [score, parent_prot, start_pos, end_pos, seq]
        out = []
        for entry in interval:
            out.append(entry)
        s.write('\t'.join([str(i) for i in out]) + '\n')
