import os
import sys
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
get = False


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
            correct_clusters.append(seq)
    
    return correct_clusters


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

new_path = os.path.join(module_path, 'testing_framework/data')
for spectrum_num,input_spectrum in enumerate(input_spectra):
    correct_sequence = correct_sequences[spectrum_num]
    b_hits,y_hits = identification.create_hits(spectrum_num, input_spectrum, matched_masses_b, matched_masses_y, False, write_path)
    for ion in 'by':
        clusters = testing_utils.create_clusters(ion, b_hits, y_hits)
        if ion == 'b':
            b_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_b)
        else:
            y_sorted_clusters = testing_utils.Bayes_clusters(ion, clusters, write_path, kmer_set, unique_y)
        
    for i, cluster in enumerate(b_sorted_clusters):
        score = cluster.score
        post_prob = cluster.prob
        seq = cluster.seq
        cluster_num = i
        ion = 'b'
        assessment, _ = testing_utils.is_good_hit(cluster.seq, ion, correct_sequence)


        with open(os.path.join(new_path, str(spectrum_num)+ "_data.txt"), 'a') as d:
            d.write(str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')
        
        with open(os.path.join(new_path, "total_data.txt"), 'a') as d:
            d.write(str(i) + '\t' + str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')

    for i, cluster in enumerate(y_sorted_clusters):
        score = cluster.score
        post_prob = cluster.prob
        seq = cluster.seq
        cluster_num = i
        ion = 'y'
        assessment, _ = testing_utils.is_good_hit(cluster.seq, ion, correct_sequence)
    
        with open(os.path.join(new_path, str(spectrum_num)+ "_data.txt"), 'a') as d:
            d.write(str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')
        
        with open(os.path.join(new_path, "total_data.txt"), 'a') as d:
            d.write(str(i) + '\t' + str(score) + '\t' + str(post_prob) + '\t' + seq + '\t' + str(cluster_num) + '\t' + str(assessment) + '\t' + ion + '\n')
print('Done')
