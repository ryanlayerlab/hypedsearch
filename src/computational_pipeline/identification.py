from multiprocessing import Pool, set_start_method
import operator
from typing import Any
from collections import ChainMap
from postprocessing.postprocessing_utils import postprocessing
from lookups.objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
from alignment import alignment
from lookups.utils import ppm_to_da, to_percent, is_json, is_file
from preprocessing import merge_search, preprocessing_utils, clustering, evaluation
from file_io import JSON
import lookups.objects as objects
import time
import multiprocessing as mp
import json
import os
from computational_pipeline.gen_spectra import convert_precursor_to_ion, calculate_masses
from scoring.scoring import second_scoring, rescore_merges
from alignment.alignment import find_alignments
import computational_pipeline.finding_seqs
import lookups.objects

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
        aligned_spectrum = objects.Aligned_Spectrum(hybrid,left_proteins,right_proteins,sequence,b_scores,y_scores,
                                    total_score, total_gaussian_score,extensions,precursor_mass,precursor_charge,
                                    total_mass_error,total_count)
        aligned_spectrums.append(aligned_spectrum)
    return aligned_spectrums

def create_b_and_y_hits(create_hits_params):
    spectrum = create_hits_params.spectrum
    converted_b = create_hits_params.converted_b
    converted_y = create_hits_params.converted_y
    matched_masses_b = create_hits_params.matched_masses_b
    matched_masses_y = create_hits_params.matched_masses_y
    b_hits,y_hits = create_hits(spectrum.num,spectrum,matched_masses_b,matched_masses_y,converted_b, converted_y)
    return b_hits,y_hits

def create_b_and_y_sorted_clusters(sorted_clusters_params):
    built_database = sorted_clusters_params.built_database
    converted_b = sorted_clusters_params.converted_b
    converted_y = sorted_clusters_params.converted_y
    b_hits = sorted_clusters_params.b_hits
    y_hits = sorted_clusters_params.y_hits
    prec_charge = sorted_clusters_params.prec_charge
    ppm_tol = sorted_clusters_params.ppm_tol
    b_clusters = clustering.create_clusters('b', b_hits, y_hits)
    b_sorted_clusters = clustering.old_score_clusters(0, b_clusters, converted_b, built_database.proteins, prec_charge, ppm_tol)
    y_clusters = clustering.create_clusters('y', y_hits, y_hits)
    y_sorted_clusters = clustering.old_score_clusters(1, y_clusters, converted_y, built_database.proteins, prec_charge, ppm_tol)
    return b_sorted_clusters,y_sorted_clusters

def create_rescored_alignments(rescored_naturals, rescored_hybrids):
    rescored_alignments = sorted(rescored_naturals + rescored_hybrids, key = lambda x: (x[0], x[1]), reverse = True)
    rescored_alignments = [x for x in rescored_alignments if x[0] > 6]
    return rescored_alignments

def create_postprocessed_alignments(post_processed_alignments_params):
    db = post_processed_alignments_params.db
    rescored_alignments = post_processed_alignments_params.rescored_alignments
    spectrum = post_processed_alignments_params.spectrum
    num_hybrids = post_processed_alignments_params.num_hybrids
    num_natives = post_processed_alignments_params.num_natives
    postprocessed_alignments = postprocessing(rescored_alignments, db, spectrum, num_hybrids, num_natives)
    return postprocessed_alignments

def prep_data_structures_for_alignment(aligned_spectrum_params):
    spectrum = aligned_spectrum_params.spectrum
    sqllite_database = aligned_spectrum_params.sqllite_database
    ppm_tolerance = aligned_spectrum_params.ppm_tolerance    
    
    input_list = spectrum.mz_values        
    b_precursor, y_precursor = convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
    matched_masses_b, matched_masses_y = merge_search.get_modified_match_masses(input_list, sqllite_database, ppm_tolerance, b_precursor, y_precursor)
    
    alignment_data = lookups.objects.AlignmentData(
    b_precursor=b_precursor,
    y_precursor=y_precursor,
    matched_masses_b=matched_masses_b,
    matched_masses_y=matched_masses_y
)
    return alignment_data

def create_create_aligned_spectrum_params(spectrum,sqllite_database,ppm_tolerance):
    create_aligned_spectrum_params = lookups.objects.CreateAlignedSpectrumParams(
        spectrum = spectrum,
        sqllite_database = sqllite_database,
        ppm_tolerance = ppm_tolerance
    )
    return create_aligned_spectrum_params

def create_target_data_params(aligned_params):
    target_data_params = lookups.objects.TargetDataParams(
        original_target_seq=aligned_params.original_target_seq,
        mz_values=aligned_params.spectrum.mz_values, 
        ppm_tolerance=aligned_params.ppm_tolerance,
        precursor_mass=aligned_params.spectrum.precursor_mass, 
        prec_tol=aligned_params.prec_tol,
        precursor_charge=aligned_params.spectrum.precursor_charge  
    )  
    return target_data_params  

def create_check_matched_masses_params(alignment_data):
    matched_masses_b = alignment_data.matched_masses_b
    matched_masses_y = alignment_data.matched_masses_y
    target_left_pids = alignment_data.converted_b
    target_right_pids = alignment_data.converted_y
    target_left_indices = [i for i in range(len(matched_masses_b))]
    target_right_indices = [i for i in range(len(matched_masses_y))]
    
    check_matched_masses_params = lookups.objects.CheckMatchedMassesParams(
        matched_masses_b=matched_masses_b,
        matched_masses_y=matched_masses_y,
        target_left_pids=target_left_pids,
        target_left_indices=target_left_indices,
        target_right_pids=target_right_pids,
        target_right_indices=target_right_indices
    )
    return check_matched_masses_params

def create_find_from_precursor_params(good_entries,aligned_spectrum_params):
    converted_b = good_entries.good_b_entries
    matched_masses_b = [entry.mass for entry in good_entries.good_b_entries]
    spectrum = aligned_spectrum_params.spectrum
    ppm_tolerance = aligned_spectrum_params.ppm_tolerance
    proteins_database = aligned_spectrum_params.proteins_database

    find_from_precursor_params =  lookups.objects.FindFromPrecursorParams(
        converted_b=converted_b,
        matched_masses_b=matched_masses_b,
        spectrum=spectrum,
        ppm_tolerance=ppm_tolerance,
        proteins_database=proteins_database
    )
    return find_from_precursor_params

def create_create_hits_params(precursor_hit_result):
    best_precursor_hit = precursor_hit_result.best_precursor_hit
    score_filter = precursor_hit_result.score_filter

    num = 0
    converted_b = [f'B{num}_converted']  
    converted_y = [f'Y{num}_converted']
    
    matched_masses_b = [101.1, 102.2, 103.3]  # Example logic, adjust as needed
    matched_masses_y = [201.1, 202.2, 203.3]  # Example logic, adjust as needed
    
    create_hits_params =  lookups.objects.CreateHitsParams(
        spectrum=best_precursor_hit.spectrum,
        converted_b=converted_b,
        converted_y=converted_y,
        matched_masses_b=matched_masses_b,
        matched_masses_y=matched_masses_y,
        num=num
    )

    return create_hits_params


def create_sorted_clusters_params(hits):
    # b_sorted_clusters = sorted(hits.b_hits, key=lambda x: x.cluster_id)
    # y_sorted_clusters = sorted(hits.y_hits, key=lambda x: x.cluster_id)
    # #TODO FIX THIS
    # sorted_clusters_params = objects.SortedClustersParams(
    #     b_sorted_clusters=b_sorted_clusters,
    #     y_sorted_clusters=y_sorted_clusters,
    #     good_b_entries=good_entries.good_b_entries,
    #     good_y_entries=good_entries.good_y_entries
    # )
    # return sorted_clusters_params
    return None

def create_check_sorted_clusters_params(sorted_clusters, good_entries):
    b_sorted_clusters = sorted_clusters.b_sorted_clusters
    y_sorted_clusters = sorted_clusters.y_sorted_clusters
    
    check_sorted_clusters_params =  objects.CheckSortedClustersParams(
        b_sorted_clusters=b_sorted_clusters,
        y_sorted_clusters=y_sorted_clusters,
        good_b_entries=good_entries.good_b_entries,
        good_y_entries=good_entries.good_y_entries
    )

    return check_sorted_clusters_params

def create_search_space_params(good_entries):
    b_sorted_clusters = sorted(good_entries.good_b_entries)
    y_sorted_clusters = sorted(good_entries.good_y_entries)
    precursor_charge = len(b_sorted_clusters) + len(y_sorted_clusters)
    search_space_params = objects.SearchSpaceParams(b_sorted_clusters, y_sorted_clusters, precursor_charge)
    return search_space_params

def create_check_searches_params(search_space):
    b_search_space = search_space.b_search_space
    y_search_space = search_space.y_search_space
    target_left_pids = [1] * len(b_search_space)  # Example data, adjust as needed
    target_right_pids = [2] * len(y_search_space)  # Example data, adjust as needed
    target_left_indices = list(range(len(b_search_space)))  # Example data, adjust as needed
    target_right_indices = list(range(len(y_search_space)))  # Example data, adjust as needed
    target_seq = "TARGET_SEQUENCE"  # Example sequence, replace with actual logic
    precursor_charge = len(b_search_space) + len(y_search_space)  # Example calculation
    ppm_tolerance = 10.0  # Example tolerance, replace with actual logic
    check_searches_params = objects.CheckSearchesParams(b_search_space, y_search_space, target_left_pids, target_right_pids, target_left_indices, target_right_indices, target_seq, precursor_charge, ppm_tolerance)
    return check_searches_params

def create_pair_natives_params(good_searches):
    # Extracting b_search_space and y_search_space from GoodSearches
    b_search_space = good_searches.good_b_searches
    y_search_space = good_searches.good_y_searches
    precursor_mass = sum(b_search_space) + sum(y_search_space)  # Example calculation
    prec_tol = 20.0  # Example tolerance, replace with actual logic
    pair_natives_params = objects.PairNativesParams(b_search_space, y_search_space, precursor_mass, prec_tol)
    return pair_natives_params

def create_pair_indices_params(aligned_spectrum_params):
    # Example calculations/assignments based on the given parameters
    b_search_space = aligned_spectrum_params.built_database.get('b_search_space', [])
    y_search_space = aligned_spectrum_params.built_database.get('y_search_space', [])
    
    precursor_mass = sum(b_search_space) + sum(y_search_space)  # Example calculation, replace with actual logic
    prec_tol = aligned_spectrum_params.prec_tol
    precursor_charge = aligned_spectrum_params.num_hybrids + aligned_spectrum_params.num_natives  # Example calculation, replace with actual logic
    score_filter = max(aligned_spectrum_params.ppm_tolerance, 0)  # Example filter, replace with actual logic
    pair_indices_params =  objects.PairIndicesParams(b_search_space, y_search_space, precursor_mass, prec_tol, precursor_charge, score_filter)
    return pair_indices_params

def create_check_merges_params(unique_native_merged_seqs, unique_hybrid_merged_seqs):
    good_b_searches = []
    good_y_searches = []

    check_merges_params = objects.CheckMergesParams(
        unique_hybrid_merged_seqs=unique_hybrid_merged_seqs,
        unique_native_merged_seqs=unique_native_merged_seqs,
        good_b_searches=good_b_searches,
        good_y_searches=good_y_searches
    )
    return check_merges_params

def create_rescore_merges_params(unique_merges):
    spectrum = None
    ppm_tolerance = None

    rescore_merges_params = objects.RescoreMergesParams(
        unique_merges=unique_merges,
        spectrum=spectrum,
        ppm_tolerance=ppm_tolerance
    )
    return rescore_merges_params

def create_check_rescored_merges_params(unique_rescored):
    good_merged_seqs = []

    check_rescored_merges_params = objects.CheckRescoredMergesParams(
        unique_rescored=unique_rescored,
        good_merged_seqs=good_merged_seqs
    )  

    return check_rescored_merges_params  

def create_post_processed_alignments_params(good_rescored):
    built_database = None
    spectrum = None
    num_hybrids = 0
    num_natives = 0

    post_processed_alignments_params = objects.PostprocessedAlignmentsParams(
        built_database=built_database,
        unique_rescored=good_rescored,
        spectrum=spectrum,
        num_hybrids=num_hybrids,
        num_natives=num_natives
    )
    
    return post_processed_alignments_params    


def create_aligned_spectrum_with_target(create_aligned_spectrum_params):
    target_data_params = create_target_data_params(create_aligned_spectrum_params)
    target_data = computational_pipeline.finding_seqs.get_target_data(target_data_params)
    alignment_data = prep_data_structures_for_alignment(create_aligned_spectrum_params)
    check_matched_masses_params = create_check_matched_masses_params(alignment_data)
    good_entries = computational_pipeline.finding_seqs.check_in_matched_masses(check_matched_masses_params)
    find_from_precursor_params = create_find_from_precursor_params(good_entries, create_aligned_spectrum_params)
    precursor_hit_result = alignment.find_from_precursor(find_from_precursor_params)
    create_hits_params = create_create_hits_params(precursor_hit_result)
    hits = create_b_and_y_hits(create_hits_params)
    sorted_clusters_params = create_sorted_clusters_params(hits)
    sorted_clusters = create_b_and_y_sorted_clusters(sorted_clusters_params)
    check_sorted_clusters_params = create_check_sorted_clusters_params(sorted_clusters)
    good_entries = computational_pipeline.finding_seqs.check_in_sorted_clusters(check_sorted_clusters_params)
    search_space_params = create_search_space_params(good_entries)
    search_space = clustering.get_search_space(search_space_params)
    check_searches_params = create_check_searches_params(search_space)
    good_searches = computational_pipeline.finding_seqs.check_in_searches(check_searches_params)
    pair_natives_params = create_pair_natives_params(good_searches)
    unique_native_merged_seqs = alignment.pair_natives(pair_natives_params)
    pair_indices_params = create_pair_indices_params(unique_native_merged_seqs)
    unique_hybrid_merged_seqs = alignment.pair_indices(pair_indices_params)
    check_merges_params = create_check_merges_params(unique_native_merged_seqs,unique_hybrid_merged_seqs)
    good_merged_seqs = computational_pipeline.finding_seqs.check_in_merges(check_merges_params)
    unique_merges = ChainMap(unique_hybrid_merged_seqs, unique_native_merged_seqs)
    rescore_merges_params = create_rescore_merges_params(unique_merges)
    unique_rescored = rescore_merges(rescore_merges_params)
    check_rescored_merges_params = create_check_rescored_merges_params(unique_rescored)
    good_rescored = computational_pipeline.finding_seqs.check_in_rescored_merges(check_rescored_merges_params)
    post_processed_alignments_params = create_post_processed_alignments_params(good_rescored)
    postprocessed_alignments = create_postprocessed_alignments(post_processed_alignments_params)
    aligned_spectrums = get_aligned_spectrums_from_postprocessed_alignments(postprocessed_alignments)
    return aligned_spectrums

def create_aligned_spectrum(create_aligned_spectrum_params):
    alignment_data = prep_data_structures_for_alignment(create_aligned_spectrum_params)
    # precursor_hit_result = alignment.find_from_precursor(alignment_data)
    # create_hits_params = create_create_hits_params(precursor_hit_result)
    # hits = create_b_and_y_hits(create_hits_params)
    # sorted_clusters_params = create_sorted_clusters_params(hits)
    # sorted_clusters = create_b_and_y_sorted_clusters(sorted_clusters_params)
    # search_space_params = create_search_space_params(sorted_clusters)
    # search_space = clustering.get_search_space(search_space_params)
    # pair_natives_params = create_pair_natives_params(search_space)
    # unique_native_merged_seqs = alignment.pair_natives(pair_natives_params)
    # pair_indices_params = create_pair_indices_params(unique_native_merged_seqs)
    # unique_hybrid_merged_seqs = alignment.pair_indices(pair_indices_params)
    # unique_merges = ChainMap(unique_hybrid_merged_seqs, unique_native_merged_seqs)
    # check_rescored_merges_params = create_check_rescored_merges_params(unique_rescored)
    # unique_rescored = rescore_merges(check_rescored_merges_params)
    # post_processed_alignments_params = create_post_processed_alignments_params(unique_rescored)
    # postprocessed_alignments = create_postprocessed_alignments(post_processed_alignments_params)
    # aligned_spectrums = get_aligned_spectrums_from_postprocessed_alignments(postprocessed_alignments)
    # return aligned_spectrums

def get_aligned_spectrums(aligned_spectrums_params):
    aligned_spectrums = []
    target_seq = aligned_spectrums_params.target_seq
    spectrums = aligned_spectrums_params.spectrums
    sqllite_database = aligned_spectrums_params.sqllite_database
    ppm_tolerance = aligned_spectrums_params.ppm_tolerance
    if len(target_seq) > 0:
        for spectrum in spectrums:
            create_aligned_spectrum_params = create_create_aligned_spectrum_params(spectrum,sqllite_database,ppm_tolerance)
            aligned_spectrum = create_aligned_spectrum_with_target(create_aligned_spectrum_params)
            aligned_spectrums.extend(aligned_spectrum)
    else:
        for spectrum in spectrums:
            create_aligned_spectrum_params = create_create_aligned_spectrum_params(spectrum,sqllite_database,ppm_tolerance)
            aligned_spectrum = create_aligned_spectrum(create_aligned_spectrum_params)
            aligned_spectrums.extend(aligned_spectrum)
    return aligned_spectrums
