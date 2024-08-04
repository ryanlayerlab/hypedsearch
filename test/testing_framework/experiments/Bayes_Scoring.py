import os
import sys
import operator
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

def remove_duplicates(merged_seqs):
    new_seqs = set()
    [new_seqs.add(x) for x in merged_seqs]
    return new_seqs

def get_top_data(b_c, y_c, correct_sequence):
    top, top10, top50 = False, False, False

    c = b_c[0]
    seq = c.seq
    assessment, _ = testing_utils.is_good_hit(seq, 'b', correct_sequence)
    if assessment == True:
        c = y_c[0]
        seq = c.seq
        assessment, _ = testing_utils.is_good_hit(seq, 'y', correct_sequence)
        if assessment == True:
            top = True
            top10 = True
            top50 = True
            return top, top10, top50
    for i in range(1,10):
        c = b_c[i]
        seq = c.seq
        assessment, _ = testing_utils.is_good_hit(seq, 'b', correct_sequence)
        if assessment == True:
            c = y_c[i]
            seq = c.seq
            assessment, _ = testing_utils.is_good_hit(seq, 'y', correct_sequence)
            if assessment == True:
                top10 = True
                top50 = True
                return top, top10, top50
    for i in range(10,50):
        c = b_c[i]
        seq = c.seq
        assessment, _ = testing_utils.is_good_hit(seq, 'b', correct_sequence)
        if assessment == True:
            c = y_c[i]
            seq = c.seq
            assessment, _ = testing_utils.is_good_hit(seq, 'y', correct_sequence)
            if assessment == True:
                top50 = True
                return top, top10, top50

    return top,top10,top50



input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)

write_path = os.path.abspath(os.path.join(module_path, 'intermediate_files'))
if get:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_from_file(os.path.join(write_path, 'matched_masses_b.txt'), os.path.join(write_path, 'matched_masses_y.txt'), os.path.join(write_path, 'kmer_set.txt'), False)
else:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_all_matched_rows(boundaries, db, max_peptide_length, True, write_path)
print('Finished matching masses')
print('Getting unique matched masses...')
unique_b, unique_y = testing_utils.get_unique_matched_masses(boundaries, matched_masses_b, matched_masses_y)
print('Done')

hit_seqs = []
top_c, top_10c, top_50c = [],[],[]
for spectrum_num,input_spectrum in enumerate(input_spectra):
    correct_sequence = correct_sequences[spectrum_num]
    b_hits,y_hits = identification.create_hits(spectrum_num, input_spectrum, matched_masses_b, matched_masses_y, False, write_path)
    for ion in 'by':
        clusters = testing_utils.create_clusters(ion, b_hits, y_hits)
        if ion == 'b':
            b_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_b)
        else:
            y_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_y)
    b_sorted_clusters = sorted(b_sorted_clusters, key=operator.attrgetter('prob','score','pid'))
    y_sorted_clusters = sorted(y_sorted_clusters, key=operator.attrgetter('prob','score','pid'))
    top, top_10, top_50 = get_top_data(b_sorted_clusters, y_sorted_clusters, correct_sequence)
    if top == True:
        top_c.append(spectrum_num)
        top_10c.append(spectrum_num)
        top_50c.append(spectrum_num)
    elif top_10 == True:
        top_10c.append(spectrum_num)
        top_50c.append(spectrum_num)
    elif top_50 == True:
        top_50c.append(spectrum_num)

print("Hit was the bottom hit ", len(top_c), " times")
print("Hit was in the bottom 10 ", len(top_10c), " times")
print("Hit was in the bottom 50 ", len(top_50c), " times")