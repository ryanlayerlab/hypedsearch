import os
import sys

module_path = os.path.abspath(os.path.join('hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
# module_path = os.path.abspath(os.path.join('../../..'))
# if module_path not in sys.path:
#     sys.path.append(module_path) 

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search
import database

import time

ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
peak_filter = 25
relative_abundance_filter = 0.1
size_lim = 10
cap_size = 20
get = False
no_k = True

def get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter):
    datasets = testing_utils.define_data()
    dataset = datasets[0]
    input_spectra_path = [os.path.join(dataset[0], 'NOD2_E3.mzML')]
    input_spectra, boundaries = preprocessing_utils.load_spectra(input_spectra_path, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
    correct_sequences = testing_utils.generate_truth_set(datasets[0])
    path = dataset[2]
    db = database.build(path)

    return input_spectra, boundaries, correct_sequences, db

input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)

write_path = os.path.abspath(os.path.join(module_path, 'intermediate_files'))
if get:
    print("Grabbing masses from file...")
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_from_file(os.path.join(write_path, 'matched_masses_b.txt'), os.path.join(write_path, 'matched_masses_y.txt'), os.path.join(write_path, 'kmer_set.txt'), no_k)
    print('Done')
else:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_modified_match_masses(boundaries, db, max_peptide_length, True, write_path)
print('Finished matching masses')

kmers = dict()
for mass in matched_masses_b:
    locations = []
    for kmer in matched_masses_b[mass]:
        seq = kmer[2]
        data = (mass, kmer[1], kmer[3])
        if seq not in kmers.keys():
            kmers[seq] = [data]
        else:
            kmers[seq].append(data)
for mass in matched_masses_y:
    locations = []
    for kmer in matched_masses_y[mass]:
        seq = kmer[2]
        data = (mass, kmer[1], kmer[3])
        if seq not in kmers.keys():
            kmers[seq] = [data]
        else:
            kmers[seq].append(data)

with open("Matched_kmers.txt", "w") as m:
    [m.write(str(x) + ": " + str(kmers[x]) + "\n") for x in kmers]