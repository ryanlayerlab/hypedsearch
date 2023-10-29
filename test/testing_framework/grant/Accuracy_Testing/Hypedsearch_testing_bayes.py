import os
import sys

# module_path = os.path.abspath(os.path.join('hypedsearch', 'src'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../../..'))
if module_path not in sys.path:
    sys.path.append(module_path) 

from testing_framework import testing_utils
from preprocessing import preprocessing_utils, merge_search, clustering
from alignment import alignment
from scoring import scoring
from gen_spectra import gen_spectrum
import identification
import database
import utils
from gen_spectra import get_precursor

import time

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

def combine_merges(pure_seqs, hybrid_seqs):
    merged_top = []
    pure_index, hybrid_index = 0,0
    if len(hybrid_seqs) == 0:
        return pure_seqs
    if len(pure_seqs) == 0:
        return hybrid_seqs
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

def filter_by_precursor(mseqs, obs_prec, tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if b_seq != y_seq:
            new_seq = b_seq + y_seq
        else:
            new_seq = b_seq
        if not (get_precursor(new_seq, charge) > obs_prec + tol):
            filtered_seqs.append(comb_seq)
    return filtered_seqs

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

def comb_seq_precursor(comb_seq, obs_prec, precursor_tol, charge): #Currently unused. Created for testing how close initial hits are to the precursor
    new_seq = overlap(comb_seq)
    tol = utils.ppm_to_da(obs_prec, precursor_tol)
    return obs_prec + tol - get_precursor(new_seq, charge)

def modified_find_next_mass(cluster, ion):
    if ion == 'b':
        target_index = cluster[2] + 1
    else:
        target_index = cluster[1]-1
    target_prot = cluster[0]
    for i, prot_name in enumerate(db.proteins):
        if i == target_prot:
            protein = db.proteins[prot_name]
            prot_seq = protein[0][1]
            to_add = prot_seq[target_index] if (target_index < len(prot_seq) and target_index > 0) else ''
            break
    
    return to_add

# def in_tol(mass, tolerance, queried_mass):
#     if queried_mass >= mass - tol:
#         if queried_mass <= mass + tol:
#             return True
#     return False
def filter_by_missing_mass(mseqs, obs_prec, tol, charge):
    filtered_seqs = []
    for comb_seq in mseqs:
        new_seq = overlap(comb_seq)
        dif = obs_prec + tol - get_precursor(new_seq, charge)
        if dif <= 1: #tol can vary but i'm not sure how much. Tol is .05 for spec 4 Other hacks are 2*tol
            filtered_seqs.append(comb_seq)
        else:
            next_b = modified_find_next_mass(comb_seq[3], 'b')
            b_seq = comb_seq[3][4]
            y_seq = comb_seq[4][4]
            b_dif = obs_prec + tol - get_precursor(b_seq + next_b + y_seq, charge)
            next_y = modified_find_next_mass(comb_seq[4], 'y')
            y_dif = obs_prec + tol - get_precursor(b_seq + next_y + y_seq, charge)
            if b_dif >= 0 or y_dif >= 0:
                filtered_seqs.append(comb_seq)
                
    return filtered_seqs

def make_merge(b, y, b_seq, y_seq):
    new_b = (b[0], b[1], b[2], b[3], b_seq)
    new_y = (y[0], y[1], y[2], y[3], y_seq)
    return (b[3] + y[3], b[1] - y[2], y[2]-b[1], new_b, new_y)    
    
            
def add_amino_acids(alignment_list, missing_mass, b_c, y_c, comb_seq, b_seq, y_seq, precursor_charge, prec_mass, tol, stop_b):
    #This function recursively adds in amino acids    
    if abs(get_precursor(b_seq + y_seq, precursor_charge) - prec_mass) <= tol:
        alignment_list.append(make_merge(b_c, y_c, b_seq, y_seq))
        return
    
    if get_precursor(b_seq + y_seq, precursor_charge) > prec_mass + tol:
        return
    
    next_b = modified_find_next_mass(b_c, 'b')
    next_y = modified_find_next_mass(y_c, 'y')
    
    if get_precursor(b_seq + y_seq, precursor_charge) < prec_mass - tol and (next_b != "") and stop_b == False:
        mod_b = b_seq + next_b
        mod_b_c = (b_c[0], b_c[1], b_c[2]+1, b_c[3], mod_b)
        add_amino_acids(alignment_list, missing_mass, mod_b_c, y_c, comb_seq, mod_b, y_seq, precursor_charge, prec_mass, tol, stop_b)
    stop_b = True
    if get_precursor(b_seq + y_seq, precursor_charge) < prec_mass - tol and (next_y != ""):
        mod_y = next_y + y_seq
        mod_y_c = (y_c[0], y_c[1]-1, y_c[2], y_c[3], mod_y)
        add_amino_acids(alignment_list, missing_mass, b_c, mod_y_c, comb_seq, b_seq, mod_y, precursor_charge, prec_mass, tol, stop_b)
        
    return

        
def find_alignments(merged_seqs, obs_prec, prec_charge, tol):
    alignments = []
    for comb_seq in merged_seqs:
        b_cluster = comb_seq[3]
        y_cluster = comb_seq[4]
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if b_seq != y_seq:
            new_seq = b_seq + y_seq
            missing_mass = obs_prec - get_precursor(new_seq, prec_charge)
            add_amino_acids(alignments, missing_mass, b_cluster, y_cluster, comb_seq, b_seq, y_seq, prec_charge, obs_prec, tol, False)           
        else:
            new_seq = b_seq
            if (abs(get_precursor(new_seq, prec_charge) - obs_prec) <= tol):
                alignments.append(comb_seq)
            
    return alignments

def write_clusters(clusters, filepath):
    with open(os.path.join(filepath, "clusters.txt"), 'w') as c:
        [c.write(str(x) + '\n') for x in clusters]

def check_merged(merged_top, correct_sequence):
    for comb_seq in merged_top:
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if (correct_sequence[:len(b_seq)] == b_seq) and (correct_sequence[len(correct_sequence) - len(y_seq):] == y_seq):
            return True

    return False

def score_by_dist(comb_seq, obs_prec, charge):
    b_seq = comb_seq[3][4]
    y_seq = comb_seq[4][4]
    if b_seq != y_seq:
        new_seq = b_seq + y_seq
    else:
        new_seq = b_seq
    dist = abs(get_precursor(new_seq, charge) - obs_prec)
    return dist

def rescore(comb_seq, input_masses, tol):
    total_score = 0
    b_seq = comb_seq[3][4]
    y_seq = comb_seq[4][4]
    if b_seq == y_seq:
        sequence = b_seq
    else:
        sequence = b_seq + y_seq
    spectrum = gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical:
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif observed - tol <= theoretical and observed + tol >= theoretical:
            total_score = total_score + 1
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
        
    return(total_score)

def second_scoring(alignments, input_spectrum, tol):
    new_merges = []
    for comb_seq in alignments:
        dist = score_by_dist(comb_seq, input_spectrum.precursor_mass, input_spectrum.precursor_charge)
        score = rescore(comb_seq, input_spectrum.mz_values, tol)
        new_merges.append((score, 1/dist, comb_seq))
    return new_merges

def evaluate_top_merges(aligned_merges, correct_sequence):
    for i, comb_seq in enumerate(aligned_merges):
        b_seq = comb_seq[2][3][4]
        y_seq = comb_seq[2][4][4]
        if b_seq != y_seq:
            new_seq = b_seq+"-"+y_seq
        else:
            new_seq = b_seq
        if new_seq == correct_sequence:
            if i == 0:
                return True, True, True
            elif i <= 10:
                return False, True, True
            else:
                return False, False, False
    return False, False, False

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

with open(os.path.join(write_path, "Spec_data.txt"), 'w') as s:
    s.write("")

top_count, top_10_count, top_50_count, miss_count, missed_but_top_count = 0, 0, 0, 0, 0
top_array, top_10_array, top_50_array, missed_array, missed_but_top_array, better_hybrids_array = [], [], [], [], [], []
found_hybrids = []

hybrid_spectra = [input_spectra[4], input_spectra[5], input_spectra[6], input_spectra[7], input_spectra[8], input_spectra[9], input_spectra[10], input_spectra[11]]
hybrid_correct_seqs = [correct_sequences[4], correct_sequences[5], correct_sequences[6], correct_sequences[7], correct_sequences[8], correct_sequences[9], correct_sequences[10], correct_sequences[11]]
# input_spectra = hybrid_spectra
# correct_sequences = hybrid_correct_seqs

# input_spectra = [hybrid_spectra[7]]
# correct_sequences = [hybrid_correct_seqs[7]]

for spectrum_num,input_spectrum in enumerate(input_spectra):
    input_spectrum = input_spectra[spectrum_num]
    # print(f'Getting seeds for {spectrum_num+1}/{len(input_spectra)} [{to_percent(spectrum_num+1, len(input_spectra))}%]', end='\r')
    print("Spectrum num:", spectrum_num)
    correct_sequence = correct_sequences[spectrum_num]
    filtered_mm_b, filtered_mm_y = filter_matched_masses(input_spectrum.mz_values, matched_masses_b, matched_masses_y)
    # unique_b, unique_y = testing_utils.get_unique_matched_masses(input_spectrum.mz_values, filtered_mm_b, filtered_mm_y)
    hit_time = time.time()
    b_hits,y_hits = identification.create_hits(spectrum_num, input_spectrum, filtered_mm_b, filtered_mm_y, False, write_path)
    hit_time = time.time() - hit_time
    cluster_time = time.time()
    for ion in "by":
            clusters = clustering.create_clusters(ion, b_hits, y_hits)
            write_clusters(clusters, write_path)
            if ion ==  'b':
                b_sorted_clusters = clustering.Score_clusters(ion, clusters)
            else:
                y_sorted_clusters = clustering.Score_clusters(ion, clusters)
    cluster_time = time.time() - cluster_time
    
    Ryan_merge_time = time.time()
    merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
    merged_seqs.sort(key = lambda x: x[0], reverse = True)
    Ryan_merge_time = time.time() - Ryan_merge_time
    
    prec_tol = utils.ppm_to_da(input_spectrum.precursor_mass, precursor_tolerance)
    
    prec_filter_time = time.time()
    merged_seqs = clustering.filter_by_precursor(merged_seqs, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    prec_filter_time = time.time() - prec_filter_time
    
    mm_filter_time = time.time()
    merged_seqs = clustering.filter_by_missing_mass(db, merged_seqs, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    mm_filter_time = time.time() - mm_filter_time
    
    Hybrid_merge_time = time.time()
    hybrid_merged = clustering.get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    Hybrid_merge_time = time.time() - Hybrid_merge_time
    
    hybrid_prec_filter_time = time.time()
    hybrid_merged = clustering.filter_by_precursor(hybrid_merged, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    hybrid_prec_filter_time = time.time() - hybrid_prec_filter_time
    
    hybrid_mm_filter_time = time.time()
    hybrid_merged = clustering.filter_by_missing_mass(db, hybrid_merged, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    hybrid_mm_filter_time = time.time() - mm_filter_time
    
    top_50_time = time.time()
    merged_top = clustering.combine_merges(merged_seqs, hybrid_merged, 50)
    top_50_time = time.time() - top_50_time

    merge_check_time = time.time()
    in_merged = check_merged(hybrid_merged, correct_sequence)
    merge_check_time = time.time() - merge_check_time

    alignment_time = time.time()
    alignments = alignment.find_alignments(merged_top, input_spectrum.precursor_mass, input_spectrum.precursor_charge, prec_tol, db)
    alignment_time = time.time() - alignment_time
    
    ss_time = time.time()
    rescored_alignments = scoring.second_scoring(alignments, input_spectrum, ppm_tolerance)
    rescored_alignments = sorted(rescored_alignments, key = lambda x: (x[0], x[1]), reverse=True)
    ss_time = time.time() - ss_time

    # eval_time = time.time()
    # top, top_10, top_50 = evaluate_top_merges(rescored_alignments, correct_sequence)
    # eval_time = time.time() - eval_time
    
    # if top == True:
    #     top_count = top_count + 1
    #     top_array.append(spectrum_num)
    # if top_10 == True and top == False:
    #     top_10_count = top_10_count + 1
    #     top_10_array.append(spectrum_num)
    #     if spectrum_num in range(4,11):
    #         found_hybrids.append(spectrum_num)
    # if in_merged == True and top_10 == False and top == False:
    #     missed_but_top_array.append(spectrum_num)
    #     missed_but_top_count = missed_but_top_count + 1
    # if top_50 == True and top_10 == False and top == False:
    #     top_50_count = top_50_count + 1
    #     top_50_array.append(spectrum_num)
    # if (top == False) and (top_10 == False) and (top_50 == False):
    #     miss_count = miss_count + 1
    #     missed_array.append(spectrum_num)
        
    with open(os.path.join(write_path, "Spec_data.txt"), 'a') as s:
        s.write("Spectrum num: " + str(spectrum_num) + '\n')
        s.write("Sizes:\n")
        s.write("\t b_hits: " + str(len(b_hits)) + '\n')
        s.write("\t y_hits: " + str(len(y_hits)) + '\n')
        s.write("\t b_clusters: " + str(len(b_sorted_clusters)) + '\n')
        s.write("\t y_clusters: " + str(len(y_sorted_clusters)) + '\n')
        s.write("\t Non-hybrid merges: " + str(len(merged_seqs)) + '\n')
        s.write("\t Hybrid merges: " + str(len(hybrid_merged)) + '\n')
        s.write("\t Alignments: " + str(len(alignments)) + '\n')
        s.write("Timing data:\n")
        s.write("\t Hits: " + str(hit_time) + '\n')
        s.write("\t Getting Clusters: " + str(cluster_time) + '\n')
        s.write("\t Ryan merge: " + str(Ryan_merge_time) + '\n')
        s.write("\t Filter by Prec: " + str(prec_filter_time) + '\n')
        s.write("\t Getting hybrids: " + str(Hybrid_merge_time) + '\n')
        s.write("\t Hybrid Filter by Prec: " + str(hybrid_prec_filter_time) + '\n')
        s.write("\t Filtering by missing mass: " + str(mm_filter_time) + '\n')
        # s.write("\t Getting top 50: " + str(top_50_time) + '\n')
        # s.write("\t Checking if in merged top: " + str(merge_check_time) + '\n')
        s.write("\t Alignments time: " + str(alignment_time) + '\n')
        s.write("\t Second scoring: " + str(ss_time) + '\n')
        # s.write("\t Evaluating hits: " + str(eval_time) + '\n')
        s.write("\n")

print("Number of correct top_alignments:", top_count)
print("Number of correct top_10 alignments:", top_10_count)
print("Number of correct Hybrids", str(len(found_hybrids))+'/'+str(len(hybrid_spectra)))
print("Number of missed alignments with a correct hit in the top 50:", missed_but_top_count)
print("Number of correct top_50 alignments:", top_50_count)
print("Number of missed alignments:", miss_count)

with open(os.path.join(write_path, "Spectra_with_correct_top.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in top_array]
with open(os.path.join(write_path, "Spectra_with_correct_top_10.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in top_10_array]
with open(os.path.join(write_path, "Spectra_with_missed_but_top.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in missed_but_top_array]
with open(os.path.join(write_path, "Spectra_with_correct_top_50.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in top_50_array]
with open(os.path.join(write_path, "Spectra_with_missed.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in missed_array]
with open(os.path.join(write_path, "Spectra_with_better_hybrids.txt"), 'w') as t:
    [t.write(str(x) + '\n') for x in better_hybrids_array]