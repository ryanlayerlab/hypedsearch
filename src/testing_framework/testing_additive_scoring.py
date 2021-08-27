import random
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

from scoring import mass_comparisons
from objects import Spectrum
from utils import hashable_boundaries


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






#Determined sequence
# correct_sequence = 'MSSP'
# ion = 'b'
# charge = 1
# spectrum = gen_spectra.gen_spectrum(correct_sequence, charge, ion)
# mz_array = spectrum["spectrum"]
# precursor_mass = spectrum["precursor_mass"]

# #Random sequence
# correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,15))
# while testing_utils.wrong_char(correct_sequence):
#     correct_sequence = ''.join(random.choice(string.ascii_uppercase) for _ in range(0,6))
# ion = None
# charge = 1
# mz_array = gen_spectra.gen_spectrum(correct_sequence, charge, ion)['spectrum']

# abundances = []
# [abundances.append(random.random() * 10 * random.randrange(0,10000)) for x in range(0,len(mz_array))]

# input_spectrum = Spectrum(mz_array, abundances, 0, 0, -1, precursor_mass, charge)

# mz_array.sort()

#Real data
def full_analysis(spectrum_num): 
    input_spectrum = input_spectra[spectrum_num]
    correct_sequence = correct_sequences[spectrum_num]





    # boundaries, mz_mapping = testing_utils.define_single_spectrum(input_spectrum[0], ppm_tolerance)
    # matched_masses_b, matched_masses_y, db = testing_utils.modified_match_masses(boundaries, db, max_peptide_length)




    # print('Grabbing hits...')

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
    # print('Done')

    b_hits, y_hits = testing_utils.parse_hits(b_hits, y_hits)
    b_hits, y_hits = testing_utils.filter_hits(b_hits, y_hits, input_spectrum.precursor_mass, input_spectrum.precursor_charge)

    scored_b, scored_y = testing_utils.score_hits(b_hits, y_hits, input_spectrum)
    TOP_X = 50
    top_50_scored_b, top_50_scored_y = testing_utils.get_top_X(scored_b, scored_y, TOP_X)

    good_kmers, bad_kmers = testing_utils.get_good_kmers(top_50_scored_b, top_50_scored_y, correct_sequence, TOP_X)
    return good_kmers, bad_kmers

print("Starting full analysis")
for i in range(0,1085):
    g_k = []
    b_k = []
    good_kmers, bad_kmers = full_analysis(i)
    [g_k.append(x) for x in good_kmers]
    [b_k.append(x) for x in bad_kmers]


with open("good_kmers.txt", 'w') as gk:
    [gk.write(x + '\n') for x in g_k]
with open("bad_kmers.txt", 'w') as bk:
    [bk.write(x + '\n') for x in b_k]