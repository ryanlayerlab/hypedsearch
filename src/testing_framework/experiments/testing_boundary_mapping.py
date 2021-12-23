import random

import os
import sys
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..'))
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
input_spectra, boundaries = testing_utils.preprocess_input_spectra(input_spectra_path, ppm_tolerance)


correct_sequences = testing_utils.generate_truth_set(datasets[0])

# path = '/home/ncol107453/jaime_hypedsearch/hypedsearch/data/database/prot_0_268.fasta'
path = dataset[2]
db = database.build(path)






#Determined sequence
correct_sequence = 'MSSP'
ion = 'b'
charge = 1
mz_array = gen_spectra.gen_spectrum(correct_sequence, charge, ion)['spectrum']

# #Random sequence
# correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,15))
# while testing_utils.wrong_char(correct_sequence):
#     correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,6))
# ion = None
# charge = 1
# mz_array = gen_spectra.gen_spectrum(correct_sequence, charge, ion)['spectrum']

abundances = []
[abundances.append(random.random() * 10 * random.randrange(0,10000)) for x in range(0,len(mz_array))]

input_spectrum = Spectrum(mz_array, abundances)

mz_array.sort()

#Real data
# spectrum_num = 0
# input_spectrum = input_spectra[spectrum_num]
# correct_sequence = correct_sequences[spectrum_num]
# print(correct_sequence)
# ion = None






boundaries, mz_mapping = testing_utils.define_single_spectrum(input_spectrum[0], ppm_tolerance)
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
with open("y_hits.txt", 'w') as y_file:
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

# with open("mz_ab.txt", 'w') as m:
#     mz_list = input_spectrum[0]
#     abundance_list = input_spectrum[1]
#     for i,mz in enumerate(input_spectrum[0]):
#         abundance = abundance_list[i]
#         out = [mz, abundance]
#         m.write('\t'.join([str(i) for i in out]) + '\n')

# Generate scored kmers
scored_b, scored_y = testing_utils.score_hits(b_hits, y_hits, input_spectrum)
top_x = 50
filtered_b, filtered_y = testing_utils.filter_hits(scored_b, scored_y, top_x)

#Check if top 50 b and y have good kmers
# good_kmers = testing_utils.check_scores_for_kmer(filtered_b, filtered_y, correct_sequence)


# Generate overlapping boundaries
intervals = testing_utils.map_hits_to_intervals(ion)

# Writing data
with open("MSSP" + "_scored_intervals.txt", 'w') as s:
    s.write(str(correct_sequence) + '\n')
    [s.write(str(x) + ', ') for x in input_spectrum.spectrum]
    s.write('\n')
    [s.write(str(x) + ', ') for x in input_spectrum.abundance]
    s.write('\n')
    s.write('\n')
    for interval in scored_intervals:
        out = []
        for entry in interval:
            out.append(entry)

        s.write('\t'.join([str(i) for i in out]) + '\n')

with open("MSSP" + "_scored_kmers.txt", 'w') as s:
    s.write(str(correct_sequence) + '\n')
    [s.write(str(x) + ', ') for x in input_spectrum.spectrum]
    s.write('\n')
    [s.write(str(x) + ', ') for x in input_spectrum.abundance]
    s.write('\n')
    s.write('\n')
    s.write('Filtered_b: \n')
    for kmer in filtered_b:
        out = []
        for entry in filtered_b:
            out.append(entry)

        s.write('\t'.join([str(i) for i in out]) + '\n')
    s.write('\n \n \n \n \n \n \n \n \n')
    s.write("Filtered_y")
    for kmer in filtered_y:
        out = []
        for entry in filtered_y:
            out.append(entry)

        s.write('\t'.join([str(i) for i in out]) + '\n')