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
    if len(hybrid_seqs) == 0:
        [merged_top.append(pure_seqs[i]) for i in range(0,10)]
        return merged_top
    if len(pure_seqs) == 0:
        [merged_top.append(hybrid_seqs[i]) for i in range(0,10)]
        return merged_top
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
                        
    return to_add

def make_merge(b, y, b_seq, y_seq):
    new_b = (b[0], b[1], b[2], b[3], b_seq)
    new_y = (y[0], y[1], y[2], y[3], y_seq)
    return (b[3] + y[3], b[1] - y[2], y[2]-b[1], new_b, new_y)

def add_amino_acids(alignment_list, missing_mass, b_c, y_c, comb_seq, b_seq, y_seq, precursor_charge, prec_mass, tol):
    #This function recursively adds in amino acids    
    if abs(get_precursor(b_seq + y_seq, precursor_charge) - prec_mass) <= tol:
        alignment_list.add(make_merge(b_c, y_c, b_seq, y_seq))
        return
    
    if get_precursor(b_seq + y_seq, precursor_charge) > prec_mass + tol:
        return
    
    next_b = modified_find_next_mass(b_c, 'b')
    next_y = modified_find_next_mass(y_c, 'y')
    
    if get_precursor(b_seq + y_seq, precursor_charge) < prec_mass - tol and (next_b != ""):
        mod_b = b_seq + next_b
        mod_b_c = (b_c[0], b_c[1], b_c[2]+1, b_c[3], mod_b)
        add_amino_acids(alignment_list, missing_mass, mod_b_c, y_c, comb_seq, mod_b, y_seq, precursor_charge, prec_mass, tol)
    if get_precursor(b_seq + y_seq, precursor_charge) < prec_mass - tol and (next_y != ""):
        mod_y = next_y + y_seq
        mod_y_c = (y_c[0], y_c[1]-1, y_c[2], y_c[3], mod_y)
        add_amino_acids(alignment_list, missing_mass, b_c, mod_y_c, comb_seq, b_seq, mod_y, precursor_charge, prec_mass, tol)
        
    return

        
def find_alignments(merged_seqs, obs_prec, prec_charge, tol):
    alignments = set()
    for comb_seq in merged_seqs:
        b_cluster = comb_seq[3]
        y_cluster = comb_seq[4]
        b_seq = comb_seq[3][4]
        y_seq = comb_seq[4][4]
        if b_seq != y_seq:
            new_seq = b_seq + y_seq
            missing_mass = obs_prec - get_precursor(new_seq, prec_charge)
            add_amino_acids(alignments, missing_mass, b_cluster, y_cluster, comb_seq, b_seq, y_seq, prec_charge, obs_prec, tol)            
        else:
            new_seq = b_seq
            if (abs(get_precursor(new_seq, prec_charge) - obs_prec) <= tol):
                alignments.add(comb_seq)
            
    return alignments

def write_clusters(clusters, filepath):
    with open(os.path.join(filepath, "clusters.txt"), 'w') as c:
        [c.write(str(x) + '\n') for x in clusters]
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
    spectrum_num = 5
    input_spectrum = input_spectra[spectrum_num]
    print(f'Getting seeds for {spectrum_num+1}/{len(input_spectra)} [{to_percent(spectrum_num+1, len(input_spectra))}%]', end='\r')
    correct_sequence = correct_sequences[spectrum_num]
    filtered_mm_b, filtered_mm_y = filter_matched_masses(input_spectrum.mz_values, matched_masses_b, matched_masses_y)
    # unique_b, unique_y = testing_utils.get_unique_matched_masses(input_spectrum.mz_values, filtered_mm_b, filtered_mm_y)
    b_hits,y_hits = identification.create_hits(spectrum_num, input_spectrum, filtered_mm_b, filtered_mm_y, False, write_path)
    for ion in "by":
            clusters = clustering.create_clusters(ion, b_hits, y_hits)
            write_clusters(clusters, write_path)
            if ion ==  'b':
                b_sorted_clusters = clustering.Score_clusters(ion, clusters)
            else:
                y_sorted_clusters = clustering.Score_clusters(ion, clusters)
    merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
    merged_seqs.sort(key = lambda x: x[0], reverse = True)
    merged_seqs = filter_by_precursor(merged_seqs, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)
    hybrid_merged = clustering.get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)
    hybrid_merged = filter_by_precursor(hybrid_merged, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)
    hybrid_merged = filter_by_missing_mass(hybrid_merged, input_spectrum.precursor_mass, precursor_tolerance, input_spectrum.precursor_charge)  

    merged_top = get_top_comb(merged_seqs, hybrid_merged)

    tol = utils.ppm_to_da(input_spectrum.precursor_mass, precursor_tolerance)
    sample_merged_seqs = [(8, 0, 0, (0, 59, 62, 2, 'DPQV'), (0, 59, 74, 6, 'DPQVAQLELGGGPGAG'))]
    alignments = find_alignments(sample_merged_seqs, input_spectrum.precursor_mass, input_spectrum.precursor_charge, tol)
    print(alignments)
    # top, top_10, top_50 = evaluate_top_merges(merged_top, correct_sequence)

    # if top == True:
    #     top_count = top_count + 1
    # if top_10 == True:
    #     top_10_count = top_10_count + 1
    # if top_50 == True:
    #     top_50_count = top_50_count + 1