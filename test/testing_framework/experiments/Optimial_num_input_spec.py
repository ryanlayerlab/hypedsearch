import os
import sys
# module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search
import identification
import database
import gen_spectra
import utils

ppm_tolerance = 20
precursor_tolerance = 10
max_peptide_length = 23
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

def get_top_data(merged_seqs, correct_sequence):
    top, top10, top50 = False, False, False

    m_seq = merged_seqs[0]
    b_seq = m_seq[3][4]
    y_seq = m_seq[4][4]
    assessment, _ = testing_utils.is_good_hit(b_seq, 'b', correct_sequence)
    if assessment == True:
        assessment, _ = testing_utils.is_good_hit(y_seq, 'y', correct_sequence)
        if assessment == True:
            top = True
            top10 = True
            top50 = True
            return top, top10, top50
    for i in range(1,10):
        m_seq = merged_seqs[i]
        b_seq = m_seq[3][4]
        y_seq = m_seq[4][4]
        assessment, _ = testing_utils.is_good_hit(b_seq, 'b', correct_sequence)
        if assessment == True:
            assessment, _ = testing_utils.is_good_hit(y_seq, 'y', correct_sequence)
            if assessment == True:
                top10 = True
                top50 = True
                return top, top10, top50
    for i in range(10,50):
        m_seq = merged_seqs[i]
        b_seq = m_seq[3][4]
        y_seq = m_seq[4][4]
        assessment, _ = testing_utils.is_good_hit(b_seq, 'b', correct_sequence)
        if assessment == True:
            assessment, _ = testing_utils.is_good_hit(y_seq, 'y', correct_sequence)
            if assessment == True:
                top50 = True
                return top, top10, top50

    return False, False, False

def filter_by_precursor(merged_sequences, ppm_tolerance, target_precursor, correct_sequence):
    new_m = []
    for m_seq in merged_sequences:
        b_seq = m_seq[3][4]
        y_seq = m_seq[4][4]
        b_end = m_seq[3][2]
        b_start = m_seq[3][1]
        y_start = m_seq[4][1]
        y_end = m_seq[4][2]
        if b_seq == y_seq:
            new_seq = b_seq
        elif (y_end >= b_end) and (y_start <= b_start):
            new_seq = y_seq
        elif (b_end >= y_end) and (b_start <= y_start):
            new_seq = b_seq
        elif (b_end >= y_start) and (b_start <= y_start):
            new_b_end = b_end - y_start
            new_seq = b_seq[0:len(b_seq)-new_b_end-1] + y_seq
        else:
            new_seq = b_seq + y_seq
        tol = utils.ppm_to_da(target_precursor, ppm_tolerance)
        if gen_spectra.get_precursor(new_seq, 2) > target_precursor + tol:
            continue
        else:
            new_m.append(m_seq)
    return new_m

current_max = 0
for peak_filter in range(25,200):
    print("On spec num", peak_filter)
    input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)

    write_path = os.path.abspath(os.path.join(module_path, 'intermediate_files'))
    matched_masses_b, matched_masses_y, kmer_set = merge_search.modified_match_masses(boundaries, db, max_peptide_length, True, write_path)
    unique_b, unique_y = testing_utils.get_unique_matched_masses(boundaries, matched_masses_b, matched_masses_y)

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
        m = testing_utils.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
        m.sort(key = lambda x: x[0]) 
        # m = filter_by_precursor(m, precursor_tolerance, input_spectrum.precursor_mass, correct_sequence)
        top, top_10, top_50 = get_top_data(m, correct_sequence)
        if top_10 == True:
            top_10c.append(spectrum_num)

    if current_max < len(top_10c):
        best_num = peak_filter
        current_max = len(top_10c)

print("Best peak filter is:", peak_filter, "at", current_max, "number of spectra in the top 10")