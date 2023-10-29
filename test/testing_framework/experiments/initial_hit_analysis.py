import os
import sys

module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
# module_path = os.path.abspath(os.path.join('../..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search, clustering
import identification
import database
import utils
from gen_spectra import get_precursor


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

def remove_duplicates(merged_seqs):
    new_seqs = set()
    [new_seqs.add(x) for x in merged_seqs]
    return new_seqs

def filter_matched_masses(masses, matched_masses_b, matched_masses_y):
    filtered_b, filtered_y = dict(), dict()
    for mass in masses:
        if mass in matched_masses_b.keys():
            filtered_b[mass] = matched_masses_b[mass]
        if mass in matched_masses_y.keys():
            filtered_y[mass] = matched_masses_y[mass]
    return filtered_b, filtered_y

def get_top_comb(pure_seqs, hybrid_seqs):
    merged_top = []
    pure_index, hybrid_index = 0,0
    while len(merged_top) < 50:
        pure = pure_seqs[pure_index]
        hybrid = hybrid_seqs[hybrid_index]
        if pure[0] >= hybrid[0]: #We give ties to the non-hybrid sequences
            merged_top.append(pure)
            pure_index = pure_index + 1
        else:
            merged_top.append(hybrid)
            hybrid_index = hybrid_index + 1
    return merged_top

def filter_by_precursor(mseqs, obs_prec, precursor_tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if b_seq != y_seq:
            new_seq = b_seq + y_seq
        else:
            new_seq = b_seq
        tol = utils.ppm_to_da(obs_prec, precursor_tol)
        if not (get_precursor(new_seq, charge) > obs_prec + tol):
            filtered_seqs.append(comb_seq)
    return filtered_seqs

def to_percent(index, total):
    return int(100 * (index)/total)

def get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start):
    seq = ''
    if y_start > b_end:
        return b_seq + y_seq
    else:
        for i in range(b_start, y_start):
            seq = seq + b_seq[i]
        return seq
def overlap(comb_seq):
    b_seq = comb_seq[3][4]
    y_seq = comb_seq[4][4]
    b_pid = comb_seq[3][0]
    y_pid = comb_seq[4][0]
    if b_pid == y_pid:
        y_start = comb_seq[4][1]
        b_end = comb_seq[3][2]
        if (y_start - b_end > 0) & (y_start - b_end < 10):
            b_start = comb_seq[3][1]
            return get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start)
        else:
            return b_seq + y_seq
    else:
        return b_seq + y_seq

def find_next_mass(comb_seq, ion):
    if ion == 'b':
        b_tup = comb_seq[3]
        target_index = b_tup[2]
        target_prot = b_tup[0]
        for i, prot_name in enumerate(db.proteins):
            if i == target_prot:
                protein = db.proteins[prot_name]
                prot_seq = protein[0][1]
                to_add = prot_seq[target_index] if target_index < len(prot_seq) else ''
                    
    else:
        y_tup = comb_seq[4]
        target_index = y_tup[1]
        target_prot = y_tup[0]
        for i, prot_name in enumerate(db.proteins):
            if i == target_prot:
                protein = db.proteins[prot_name]
                prot_seq = protein[0][1]
                to_add = prot_seq[target_index] if target_index < len(prot_seq) else ''
    
    return to_add

def filter_by_missing_mass(mseqs, obs_prec, precursor_tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        new_seq = overlap(comb_seq)
        tol = utils.ppm_to_da(obs_prec, precursor_tol)
        if abs(obs_prec - get_precursor(new_seq, charge)) <= tol:
            filtered_seqs.append(comb_seq)
        else:
            next_b = find_next_mass(comb_seq, 'b')
            b_seq = comb_seq[3][4]
            y_seq = comb_seq[4][4]
            b_dif = obs_prec + tol - get_precursor(b_seq + next_b + y_seq, charge)
            next_y = find_next_mass(comb_seq, 'y')
            y_dif = obs_prec + tol - get_precursor(b_seq + next_y + y_seq, charge)
            if b_dif >= 0 and y_dif >= 0:
                filtered_seqs.append(comb_seq)
                
    return filtered_seqs

# def evaluate_top_merges()


input_spectra, boundaries, correct_sequences, db = get_spectra_and_db(ppm_tolerance, peak_filter, relative_abundance_filter)

write_path = os.path.abspath(os.path.join(module_path, 'intermediate_files'))
if get:
    print("Grabbing masses from file...")
    matched_masses_b, matched_masses_y, kmer_set = merge_search.get_from_file(os.path.join(write_path, 'matched_masses_b.txt'), os.path.join(write_path, 'matched_masses_y.txt'), os.path.join(write_path, 'kmer_set.txt'), no_k)
    print('Done')
else:
    matched_masses_b, matched_masses_y, kmer_set = merge_search.modified_match_masses(boundaries, db, max_peptide_length, True, write_path)
print('Finished matching masses')
# print('Getting unique matched masses...')
# unique_b, unique_y = testing_utils.get_unique_matched_masses(boundaries, matched_masses_b, matched_masses_y)
# print('Done')

top_count, top_10_count, top_50_count = False, False, False
for spectrum_num,input_spectrum in enumerate(input_spectra):
    print(f'Getting seeds for {spectrum_num+1}/{len(input_spectra)} [{to_percent(spectrum_num+1, len(input_spectra))}%]', end='\r')
    correct_sequence = correct_sequences[spectrum_num]
    filtered_mm_b, filtered_mm_y = filter_matched_masses(input_spectrum.mz_values, matched_masses_b, matched_masses_y)
    # unique_b, unique_y = testing_utils.get_unique_matched_masses(input_spectrum.mz_values, filtered_mm_b, filtered_mm_y)
    b_hits,y_hits = identification.create_hits(spectrum_num, input_spectrum, filtered_mm_b, filtered_mm_y, False, write_path)
    for ion in "by":
            clusters = clustering.create_clusters(ion, b_hits, y_hits)
            if ion ==  'b':
                b_sorted_clusters = clustering.Score_clusters(ion, clusters)
            else:
                y_sorted_clusters = clustering.Score_clusters(ion, clusters)
    merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
    merged_seqs.sort(key = lambda x: x[0], reverse = True)
    merged_seqs = filter_by_precursor(merged_seqs, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)
    hybrid_merged = clustering.get_hybrid_matches(input_spectrum, b_sorted_clusters, y_sorted_clusters, precursor_tolerance)
    hybrid_merged = filter_by_precursor(hybrid_merged, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)
    hybrid_merged = filter_by_missing_mass(hybrid_merged, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)  

    merged_top = get_top_comb(merged_seqs, hybrid_merged)
    # top, top_10, top_50 = evaluate_top_merges(merged_top, correct_sequence)

    # if top == True:
    #     top_count = top_count + 1
    # if top_10 == True:
    #     top_10_count = top_10_count + 1
    # if top_50 == True:
    #     top_50_count = top_50_count + 1