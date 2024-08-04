from multiprocessing import Pool, set_start_method
import operator, time, json, os
from typing import Any
from collections import ChainMap
from postprocessing.postprocessing_utils import postprocessing
from lookups.objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry, KMer, ConvertedPrecursors,AlignedSpectrumsParams,AlignedSpectrumParams,AlignedSpectrum
from alignment import alignment
from lookups.utils import ppm_to_da, to_percent, is_json, is_file
from preprocessing import merge_search, preprocessing_utils, clustering, evaluation
from file_io import JSON
import multiprocessing as mp
from computational_pipeline.gen_spectra import convert_precursor_to_ion, calculate_masses
from scoring.scoring import second_scoring, rescore_merges
from alignment.alignment import find_alignments
import computational_pipeline.finding_seqs


ID_SPECTRUM = 0
MULTIPROCESSING = 0

global alignment_times
alignment_times = []
global b_scoring_times
b_scoring_times = []
global y_scsoring_times
y_scoring_times = []
global filter_times
filter_times = []

TOP_X = 50

def handle_DEV_truth(filtered_b,filtered_y,b_results,keep_b_count,y_results,keep_y_count,fall_off,_id,is_hybrid,truth_seq,spectrum):
    metadata = {
        'top_x_b_hits': filtered_b, 
        'top_x_y_hits': filtered_y, 
        'excluded_b_hits': [x[0] for x in b_results[keep_b_count:]],
        'excluded_y_hits': [x[0] for x in y_results[keep_y_count:]], 
        'cut_off_b_score': b_results[keep_b_count - 1][1], 
        'cut_off_y_score': y_results[keep_y_count - 1][1]
    }
    fall_off[_id] = DEVFallOffEntry(
        is_hybrid, 
        truth_seq, 
        'top_x_filtering', 
        metadata
    )

def write_hits(b_hits, y_hits, location):
    with open(os.path.join(location, "b_hits.txt"), 'w+') as b:
        for x in b_hits:
            pep_id = x[0]
            w = x[1]
            prot_id = x[2][1]
            loc = x[2][2]
            ion = x[2][3]
            charge = x[2][4]
            out = [pep_id, w, prot_id, loc, ion, charge]
            b.write('\t'.join([str(i) for i in out]) + '\n')
    with open(os.path.join(location, "y_hits.txt"), 'w+') as b:
        for y in y_hits:
            pep_id = y[0]
            w = y[1]
            prot_id = y[2][1]
            loc = y[2][2]
            ion = y[2][3]
            charge = y[2][4]
            out = [pep_id, w, prot_id, loc, ion, charge]
            b.write('\t'.join([str(i) for i in out]) + '\n')

def get_hits_from_file(bf, yf):
    b_hits, y_hits = [], []
    with open(bf, 'r') as b:
        for line in b:
            A = line.rstrip().split('\t')
            pep_id = int(A[0])
            w = float(A[1])
            prot_id = int(A[2])
            seq = A[3]
            loc = A[4]
            ion = A[5]
            charge = int(A[6])
            out = [pep_id, w, prot_id, seq, loc, ion, charge]
            b_hits.append(out)
    with open(yf, 'r') as b:
        for line in b:
            A = line.rstrip().split('\t')
            pep_id = int(A[0])
            w = float(A[1])
            prot_id = int(A[2])
            seq = A[3]
            loc = A[4]
            ion = A[5]
            charge = int(A[6])
            out = [pep_id, w, prot_id, seq, loc, ion, charge]
            y_hits.append(out)
    return b_hits, y_hits

def create_hits(spec_num,spectrum,matched_masses_b,matched_masses_y,b_prec,y_prec):
    b_hits, y_hits = [], []
    for mz in spectrum.mz_values:
        if mz in matched_masses_b:
            for tuple in matched_masses_b[mz]:
                tup = (spec_num, mz, tuple)
                b_hits.append(tup)
        if mz in matched_masses_y:
            for tuple in matched_masses_y[mz]:
                tup = (spec_num, mz, tuple)
                y_hits.append(tup)
    
    for tuple in matched_masses_b[b_prec]:
        tup = (spec_num, b_prec, tuple)
        b_hits.append(tup)
    for tuple in matched_masses_y[y_prec]:
        tup = (spec_num, y_prec, tuple)
        y_hits.append(tup)
    return b_hits, y_hits

def handle_DEV_setup(truth):
    truth = mp.Manager().dict(truth)

def handle_DEV_result(output_dir,fall_off,cores):
    output_dir = output_dir + '/' if output_dir[-1] != '/' else output_dir
    safe_write_fall_off = {}
    for k, v in fall_off.items():
        safe_write_fall_off[k] = v._asdict()
    JSON.save_dict(output_dir + 'fall_off.json', safe_write_fall_off)
    if cores == 1:
        identification_instrumentation = objects.Identification_Instrumentation(
        average_b_scoring_time = sum(b_scoring_times)/len(b_scoring_times),
        average_y_scoring_time = sum(y_scoring_times)/len(y_scoring_times),
        time_to_filter_out_top_50_kmers = sum(filter_times)/len(filter_times),
        average_extension_time = sum(alignment.extension_times)/len(alignment.extension_times),
        average_non_hybrid_refinement_time = sum(alignment.Non_hybrid_refine_time)/len(alignment.Non_hybrid_refine_time),
        average_non_hybrid_scoring_time = sum(alignment.non_hybrid_scoring_times)/len(alignment.non_hybrid_scoring_times),
        average_hybrid_refinement_time = sum(alignment.Hybrid_refine_times)/len(alignment.Hybrid_refine_times),
        average_hybrid_scoring_time = sum(alignment.hybrid_scoring_times)/len(alignment.hybrid_scoring_times),
        average_alignment_time = sum(alignment_times)/len(alignment_times)
        )
            
def check_top_location(top_naturals, top_hybrids, natural_seqs, hybrid_seqs):
    top_natural, top_hybrid = top_naturals[0], top_hybrids[0]
    top_natural_cluster, top_hybrid_cluster = top_natural[2], top_hybrid[2]
    top_nat_location, top_hyb_location = -1, -1
    for i,seq in enumerate(natural_seqs):
        b_cluster, target_b_cluster = seq[3], top_natural_cluster[0]
        y_cluster, target_y_cluster = seq[4], top_natural_cluster[1]
        if (b_cluster[0] == target_b_cluster[5]) and (b_cluster[1] == target_b_cluster[1]) and (y_cluster[2] == target_y_cluster[2]) and (b_cluster[2] <= target_b_cluster[2]) and (y_cluster[1] >= target_y_cluster[1]):
            top_nat_location = i
            break
    
    for i,seq in enumerate(hybrid_seqs):
        b_cluster, target_b_cluster = seq[3], top_hybrid_cluster[0]
        y_cluster, target_y_cluster = seq[4], top_hybrid_cluster[1]
        if (b_cluster[0] == target_b_cluster[5]) and (y_cluster[0] == target_y_cluster[5]):
            if (y_cluster[2] == target_y_cluster[2]) and (b_cluster[1] == target_b_cluster[1]):
                if (b_cluster[2] <= target_b_cluster[2]) and (y_cluster[1] >= target_y_cluster[1]):
                    top_hyb_location = i
                    break
    
    if top_hyb_location == -1:
        for i,seq in enumerate(natural_seqs):
            b_cluster, target_b_cluster = seq[3], top_hybrid_cluster[0]
            y_cluster, target_y_cluster = seq[4], top_hybrid_cluster[1]
            if (b_cluster[0] == target_b_cluster[5]) and (b_cluster[1] == target_b_cluster[1]) and (y_cluster[2] == target_y_cluster[2]) and (b_cluster[2] <= target_b_cluster[2]) and (y_cluster[1] >= target_y_cluster[1]):
                top_hyb_location = i
                break
    
    with open("locations.txt", 'a') as l:
        l.write(str(top_nat_location) + '\t' + str(top_hyb_location) + '\n')        
    return
    
def find_sequence(b_sequence, y_sequence, b_pid, y_pid, protein_list):
    b_prot_sequence = protein_list[b_pid][1]
    b_target_starts, y_target_ends = [],[]
    for i in range(0,len(b_prot_sequence)-len(b_sequence)+1):
        testing_b = b_prot_sequence[i:i+len(b_sequence)]
        if testing_b == b_sequence:
            b_target_starts.append(i)
    y_prot_sequence = protein_list[y_pid][1]
    for i in range(0, len(y_prot_sequence)-len(y_sequence)+1):
        testing_y = y_prot_sequence[i:i+len(y_sequence)]
        if testing_y == y_sequence:
            y_target_ends.append(i+len(y_sequence))
            break
            
    return b_target_starts, y_target_ends
    
def find_target_clusters(b_sorted_clusters, y_sorted_clusters, b_sequence, y_sequence, b_pid, y_pid, protein_list):
    b_target_starts, y_target_ends = find_sequence(b_sequence, y_sequence, b_pid, y_pid, protein_list)
    print("\n")
    
    print("For b:")
    for i, cluster in enumerate(b_sorted_clusters):
        if cluster[1] == b_pid and cluster[2] in b_target_starts:
            print(i, cluster)
    print("\n For y:")
    for i, cluster in enumerate(y_sorted_clusters):
        if cluster[1] == y_pid and cluster[3] in y_target_ends: 
            print(i, cluster)
            
def group_by_uniqueness(natives, hybrids):
    unique_merges = dict()
    for merge in hybrids:
        left_seq, right_seq = merge[0][6], merge[1][6]
        full_seq = left_seq + right_seq
        if (full_seq, 1) not in unique_merges.keys():
            unique_merges[(full_seq, 1)] = []
        unique_merges[(full_seq, 1)].append(merge)
    
    for merge in natives:
        full_seq = merge[0][6]
        if (full_seq, 0) not in unique_merges.keys():
            unique_merges[(full_seq, 0)] = []
        unique_merges[(full_seq, 0)].append(merge)
    return unique_merges

def get_aligned_spectrums_from_postprocessed_alignments(postprocessed_alignments):
    aligned_spectrums = []
    for postprocessed_alignment in postprocessed_alignments:
        hybrid = postprocessed_alignment[0]
        left_proteins,right_proteins = postprocessed_alignment[1], postprocessed_alignment[2]
        sequence = postprocessed_alignment[3]
        b_scores = postprocessed_alignment[4]
        y_scores = postprocessed_alignment[5]
        total_score = postprocessed_alignment[6]
        total_gaussian_score = str(postprocessed_alignment[7])
        extensions = postprocessed_alignment[8]
        precursor_mass, precursor_charge = str(postprocessed_alignment[9]), str(postprocessed_alignment[10])
        total_mass_error = str(postprocessed_alignment[11])
        total_count = str(0)
        aligned_spectrum = AlignedSpectrum(hybrid,left_proteins,right_proteins,sequence,b_scores,y_scores,
                                    total_score, total_gaussian_score,extensions,precursor_mass,precursor_charge,
                                    total_mass_error,total_count)
        aligned_spectrums.append(aligned_spectrum)
    return aligned_spectrums

# def get_fragment_hits(aligned_spectrum_params, converted_precursors, matched_masses):
#     spectrum = aligned_spectrum_params.spectrum
#     converted_b = converted_precursors.converted_precursor_b
#     converted_y = converted_precursors.converted_precursor_y
#     matched_masses_b = matched_masses.matched_masses_b
#     matched_masses_y = matched_masses.matched_masses_y
#     b_hits,y_hits = create_hits(spectrum.num,spectrum,matched_masses_b,matched_masses_y,converted_b, converted_y)
#     hits = objects.FragmentHits(b_hits=b_hits,y_hits=y_hits)
#     return hits

def get_clusters(base_alignment_params, b_kmers,y_kmers):
    sqllite_database = base_alignment_params.sqllite_database
    ppm_tolerance = base_alignment_params.ppm_tolerance
    b_clusters = clustering.create_b_clusters(b_kmers,sqllite_database)
    y_clusters = clustering.create_y_clusters(y_kmers,sqllite_database) 
    return (b_clusters,y_clusters)

def create_rescored_alignments(rescored_naturals, rescored_hybrids):
    rescored_alignments = sorted(rescored_naturals + rescored_hybrids, key = lambda x: (x[0], x[1]), reverse = True)
    rescored_alignments = [x for x in rescored_alignments if x[0] > 6]
    return rescored_alignments

def create_postprocessed_alignments(aligned_spectrum_params, rescored_alignments):
    spectrum = aligned_spectrum_params.spectrum
    base_alignment_params = aligned_spectrum_params.base_alignment_params
    sqllite_database = base_alignment_params.sqllite_database
    number_hybrids = base_alignment_params.number_hybrids
    number_natives = base_alignment_params.number_natives
    postprocessed_alignments = postprocessing(rescored_alignments, sqllite_database, spectrum, number_hybrids, number_natives)
    return postprocessed_alignments

def get_converted_precursors(aligned_spectrum_params):
    spectrum = aligned_spectrum_params.spectrum
    converted_precursor_b, converted_precursor_y = convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    converted_precursors =  ConvertedPrecursors(converted_precursor_b=converted_precursor_b,converted_precursor_y=converted_precursor_y)
    return converted_precursors

def get_matched_kmers(aligned_spectrum_params):
    spectrum = aligned_spectrum_params.spectrum
    sqllite_database = aligned_spectrum_params.base_alignment_params.sqllite_database
    ppm_tolerance = aligned_spectrum_params.base_alignment_params.ppm_tolerance
    input_masses = spectrum.mz_values        
    b_rows, y_rows = merge_search.get_all_matched_rows(input_masses, sqllite_database, ppm_tolerance)
    b_kmers = [KMer(*row) for row in b_rows]
    y_kmers = [KMer(*row) for row in y_rows]
    return (b_kmers,y_kmers)

def print_named_tuple(tuple_description, named_tuple):
    file_path = tuple_description + '.txt'
    with open(file_path, 'w') as file:
        file.write(str(named_tuple))

def create_aligned_spectrum_with_target(aligned_spectrum_params):
    return None
#     spectrum = aligned_spectrum_params.spectrum
#     base_alignment_params = aligned_spectrum_params.base_alignment_params
#     precursor_charge = spectrum.precursor_charge
#     precursor_mass = spectrum.precursor_mass
#     precursor_tolerance = base_alignment_params.precursor_tolerance
#     target_data = computational_pipeline.finding_seqs.get_target_data(aligned_spectrum_params)
#     converted_precursors = get_converted_precursors(aligned_spectrum_params)
#     matched_masses = get_matched_kmers(aligned_spectrum_params,converted_precursors, converted_precursors)
#     target_alignment_data = computational_pipeline.finding_seqs.check_in_matched_masses(converted_precursors,matched_masses,target_data)
#     if len(target_alignment_data.matched_masses_b) > 0 or len(target_alignment_data.matched_masses_y) > 0:
#         precursor_hit_result = alignment.get_percursor_hit(aligned_spectrum_params, target_alignment_data)
#         hits = get_fragment_hits(aligned_spectrum_params, target_alignment_data) 
#         sorted_clusters = get_clusters(base_alignment_params, precursor_charge, hits, converted_precursors)
#         target_sorted_clusters = computational_pipeline.finding_seqs.check_in_sorted_clusters(target_alignment_data,sorted_clusters)
#         search_space = clustering.get_search_space(sorted_clusters,precursor_charge)
#         good_searches = computational_pipeline.finding_seqs.check_in_searches(aligned_spectrum_params,target_data,target_alignment_data)
#         unique_native_merged_seqs = alignment.pair_natives(search_space, precursor_mass, precursor_tolerance)
#         score_filter = precursor_hit_result.score_filter
#         unique_hybrid_merged_seqs = alignment.pair_indices(aligned_spectrum_params, search_space, score_filter)
#         unknown = ([],[])
#         good_merged_seqs = computational_pipeline.finding_seqs.check_in_merges(good_searches,unknown)
#         unique_merges = ChainMap(unique_hybrid_merged_seqs, unique_native_merged_seqs)
#         rescored_alignments = rescore_merges(aligned_spectrum_params,unique_merges)
#         rescored_merges = []
#         good_merges = []
#         good_rescored_alignments = computational_pipeline.finding_seqs.check_in_rescored_merges(rescored_merges,good_merges)
#         postprocessed_alignments = create_postprocessed_alignments(aligned_spectrum_params, good_rescored_alignments)
#         aligned_spectrums = get_aligned_spectrums_from_postprocessed_alignments(postprocessed_alignments)
#         return aligned_spectrums
#     else:
#         return None

def create_aligned_spectrum(aligned_spectrum_params):
    spectrum = aligned_spectrum_params.spectrum
    base_alignment_params = aligned_spectrum_params.base_alignment_params
    (b_kmers,y_kmers) = get_matched_kmers(aligned_spectrum_params)
    (b_clusters,y_clusters) = get_clusters(base_alignment_params, b_kmers,y_kmers)
    print(b_clusters)
    # search_space = clustering.get_search_space(clusters,precursor_charge)
    # unique_native_merged_seqs = alignment.pair_natives(search_space, precursor_mass, precursor_tolerance)
    # unique_hybrid_merged_seqs = alignment.pair_indices(aligned_spectrum_params, search_space, score_filter)
    # unique_merges = ChainMap(unique_hybrid_merged_seqs, unique_native_merged_seqs)
    # rescored_alignments = rescore_merges(aligned_spectrum_params,unique_merges)
    # postprocessed_alignments = create_postprocessed_alignments(aligned_spectrum_params, rescored_alignments)
    # aligned_spectrums = get_aligned_spectrums_from_postprocessed_alignments(postprocessed_alignments)
    # return aligned_spectrums
    return None

def get_aligned_spectrums(aligned_spectrums_params:AlignedSpectrumsParams):
    aligned_spectrums = []
    spectrums = aligned_spectrums_params.spectrums
    base_alignment_params = aligned_spectrums_params.base_alignment_params
    target_seq = base_alignment_params.target_seq
    if len(target_seq) > 0:
        for spectrum in spectrums:
            aligned_spectrum_params = AlignedSpectrumParams(spectrum=spectrum,base_alignment_params=base_alignment_params)
            aligned_spectrum = create_aligned_spectrum_with_target(aligned_spectrum_params)
            if aligned_spectrum is not None:
                aligned_spectrums.extend(aligned_spectrum)
    else:
        for spectrum in spectrums:
            aligned_spectrum_params = AlignedSpectrumParams(spectrum=spectrum,base_alignment_params=base_alignment_params)
            aligned_spectrum = create_aligned_spectrum(aligned_spectrum_params)
            if aligned_spectrum is not None:
                aligned_spectrums.extend(aligned_spectrum)
    return aligned_spectrums
