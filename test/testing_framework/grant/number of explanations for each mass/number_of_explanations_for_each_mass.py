from heapq import merge
import os
import sys

# module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src')) #For running from vscode
# if module_path not in sys.path:
#     sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../../..')) #For running in terminal
if module_path not in sys.path:
    sys.path.append(module_path)

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search, clustering
from gen_spectra import gen_spectrum
import identification
import database
import utils
from gen_spectra import get_precursor

import matplotlib.pyplot as plt
 
ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
peak_filter = 25
relative_abundance_filter = 0.1
size_lim = 10
cap_size = 20
get = True
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
    matched_masses_b, matched_masses_y, kmer_set = merge_search.modified_match_masses(boundaries, db, max_peptide_length, True, write_path)
print('Finished matching masses')

unique_b,unique_y = testing_utils.get_unique_matched_masses(boundaries, matched_masses_b, matched_masses_y)
unique_m = {**unique_b, **unique_y}

unique_mass_lengths = []
for mass_key in unique_m.keys():
    unique_mass_lengths.append(len(unique_m[mass_key]))

plt.bar(unique_m.keys(), unique_mass_lengths)
plt.title('Number of explanations for each mass')
plt.xlabel('Mass size')
plt.ylabel('Number of explanations')
plt.legend()
plt.savefig("Number_of_explanations_for_each_mass")