import os
import sys
# module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search, clustering
import identification
import database

ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
peak_filter = 25
relative_abundance_filter = 0.1
size_lim = 10
cap_size = 20
get = True


import matplotlib.pyplot as plt

def get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter):
    datasets = testing_utils.define_data()
    dataset = datasets[0]
    input_spectra_path = [os.path.join(dataset[0], 'NOD2_E3.mzML')]
    input_spectra, boundaries = preprocessing_utils.load_spectra(input_spectra_path, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
    correct_sequences = testing_utils.generate_truth_set(datasets[0])
    path = dataset[2]
    db = database.build(path)

    return input_spectra, boundaries, correct_sequences, db

def filter_spectra_by_size(correct_answers, spectra, size_lim):
    new_spectra = []
    for i, spectrum in enumerate(spectra):
        if len(correct_answers[i]) <= size_lim:
            new_spectra.append(spectrum)
    return new_spectra

def filter_for_correct(clusters, ion, correct_seq):
    correct_clusters = []
    for c in clusters:
        score = c[0]
        pid = c[1]
        seq = c[2]
        mz = c[3]
        start = c[4]
        end = c[5]
        indices = c[6:]
        assessment, _ = testing_utils.is_good_hit(seq, ion, correct_seq)
        if assessment:
            correct_clusters.append(c)
    
    return correct_clusters

def remove_duplicates(merged_seqs):
    new_seqs = set()
    [new_seqs.add(x) for x in merged_seqs]
    return new_seqs

def tally_seqs(b_sorted_clusters, y_sorted_clusters):
    mz_list = []
    for c in b_sorted_clusters:
        for index in c.indices:
            mz = index[3]
            mz_list.append(mz)
    for c in y_sorted_clusters:
        for index in c.indices:
            mz = index[3]
            mz_list.append(mz)
    return mz_list

def isclose(hit1, hit2, dist):
    if abs(hit1-hit2) >= dist:
        return True
    else:
        return False

def bundle(hit_seqs, bin_width):
    hit_seqs = sorted(hit_seqs)
    bundled_hits = []
    last_bundle = 0
    current_count = 0
    for j, hit in enumerate(hit_seqs):
        if isclose(hit, hit_seqs[last_bundle], bin_width):
            current_count = current_count + 1
        else:
            bundled_hits.append((hit_seqs[last_bundle], current_count))
            current_count = 0
    return bundled_hits


input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)
input_spectra = filter_spectra_by_size(correct_sequences, input_spectra, size_lim)

write_path = os.path.abspath(os.path.join(module_path, 'intermediate_files'))
if get:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_from_file(os.path.join(write_path, 'matched_masses_b.txt'), os.path.join(write_path, 'matched_masses_y.txt'), os.path.join(write_path, 'kmer_set.txt'), False)
else:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.modified_match_masses(boundaries, db, max_peptide_length, True, write_path)
print('Finished matching masses')
print('Getting unique matched masses...')
unique_b, unique_y = testing_utils.get_unique_matched_masses(boundaries, matched_masses_b, matched_masses_y)
print('Done')

hit_seqs = []
for i,input_spectrum in enumerate(input_spectra):
    b_hits,y_hits = identification.create_hits(i, input_spectrum, matched_masses_b, matched_masses_y, False, write_path)
    for ion in 'by':
        clusters = testing_utils.create_clusters(ion, b_hits, y_hits)
        correct_sequence = correct_sequences[i]
        clusters = filter_for_correct(clusters, ion, correct_sequence)
        if ion == 'b':
            b_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_b)
        else:
            y_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_y)
    [hit_seqs.append(x) for x in tally_seqs(b_sorted_clusters, y_sorted_clusters)]

bundled_hits = bundle(hit_seqs, cap_size)
new_path = os.path.join(module_path, 'testing_framework/data/plots')
h, k = [], []
[h.append(x[0]) for x in bundled_hits]
[k.append(x[1]) for x in bundled_hits]
ax = plt.bar(h, k)
plt.xlabel("Size of hit")
plt.ylabel("Frequency")
plt.title("Number of times correct hits come from mz based on size")
plt.savefig('Frequency by Size')