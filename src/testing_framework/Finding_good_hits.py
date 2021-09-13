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
from utils import ppm_to_da, overlap_intervals
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


print('Grabbing hits...')
for spectrum_num, input_spectrum in enumerate(input_spectra):
    #spectrum_num = 0
    correct_sequence = correct_sequences[spectrum_num]
    print(spectrum_num, correct_sequence)

    # input_spectrum = input_spectra[spectrum_num]

    b_hits =[]
    y_hits = []
    b_set = set()
    y_set = set()
    correct_hits = []
    mz_miss_set = set()
    #Remember to add in abundance if it is helpful
    input_num = 0
    testing_utils.find_hits(mz_mapping, boundaries, input_spectrum, input_num, matched_masses_b, matched_masses_y, b_hits, y_hits, b_set, y_set, mz_miss_set)
    testing_utils.append_correct_hits(correct_hits, correct_sequence, input_spectrum, ppm_tolerance)
    # print('Done')

    testing_utils.write_data(b_hits, y_hits)

    ion = 'b'
    testing_utils.cluster_hits(ion)
    ion = 'y'
    testing_utils.cluster_hits(ion)

    b_sorted_clusters, y_sorted_clusters = testing_utils.sort_clusters()
    # testing_utils.print_clusters(b_sorted_clusters, y_sorted_clusters)

    b_hit_arr, y_hit_arr = testing_utils.get_hits_from_cluster(b_sorted_clusters, y_sorted_clusters)
    b_good_hits = []
    y_good_hits = []
    for i, hit in enumerate(b_hit_arr):
        assessment, score = testing_utils.is_good_hit(hit, 'b', correct_sequence)
        if assessment:
            tuple = (i, hit, score)
            b_good_hits.append(tuple)
    for i, hit in enumerate(y_hit_arr):
        assessment, score = testing_utils.is_good_hit(hit, 'y', correct_sequence)
        if assessment:
            tuple = (i, hit, score)
            y_good_hits.append(tuple)
    with open('b_good_hits.txt', 'w') as b:
        b.write("Spectrum_num: " + str(spectrum_num) + '\n')
        b.write(correct_sequence + ": \n\n\n\n")
        [b.write(str(x) + '\n') for x in b_good_hits]
    with open('y_good_hits.txt', 'w') as y:
        y.write("Spectrum_num: " + str(spectrum_num) + '\n\n\n\n')
        [y.write(str(x) + '\n') for x in y_good_hits]
        y.write(correct_sequence + ": \n\n\n\n")