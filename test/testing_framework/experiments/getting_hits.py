import os
import sys
import random
import string
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
from src.constants.objects import Spectrum

#Assumptions:
max_peptide_length = 20
ppm_tolerance = 20


datasets = testing_utils.define_data()

dataset = datasets[0]

input_spectra_path = dataset[0]
# path = '/home/ncol107453/jaime_hypedsearch/hypedsearch/data/database/single_prot.fasta'
path = dataset[2]
db = database.build(path)

input_spectra, boundaries, mz_mapping = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


correct_sequences = testing_utils.generate_truth_set(datasets[0])



#Determined sequence
# correct_sequence = 'MSSP'
# ion = None
# charge = 1
# mz_array = gen_spectra.gen_spectrum(correct_sequence, charge, ion)['spectrum']

# #Random sequence
# correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,15))
# while testing_utils.wrong_char(correct_sequence):
#     correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,6))
# print(correct_sequence)
# ion = None
# charge = 1
# mz_array = gen_spectra.gen_spectrum(correct_sequence, charge, ion)['spectrum']

# abundances = []
# [abundances.append(random.random() * 10 * random.randrange(0,10000)) for x in range(0,len(mz_array))]

# input_spectrum = Spectrum(mz_array, abundances)

# mz_array.sort()

#Real data
spectrum_num = 0
input_spectrum = input_spectra[spectrum_num]
correct_sequence = correct_sequences[spectrum_num]
print(correct_sequence)
ion = None
print(len(input_spectrum.spectrum))

boundaries, mz_mapping = testing_utils.define_single_spectrum(input_spectrum[0], ppm_tolerance)
matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)

print('Getting hits...')

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
with open("Spectrum_b_hits.txt", 'w') as b:
    for i, x in enumerate(b_hits):
        pep_id = x[0]
        w = x[1][0]
        prot_id = x[1][1]
        seq = x[1][2]
        loc = x[1][3]
        write_ion = x[1][4]
        charge = x[1][5]
        x = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        b_hits[i] = x
        b.write('\t'.join([str(i) for i in x]) + '\n')
with open("Spectrum_y_hits.txt", 'w') as y_file:
    for i, y in enumerate(y_hits):
        pep_id = y[0]
        w = y[1][0]
        prot_id = y[1][1]
        seq = y[1][2]
        loc = y[1][3]
        write_ion = y[1][4]
        charge = y[1][5]
        y = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        y_hits[i] = y
        y_file.write('\t'.join([str(i) for i in y]) + '\n')
print('Done')