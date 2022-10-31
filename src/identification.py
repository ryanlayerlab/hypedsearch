from collections import defaultdict
from hashlib import new
from heapq import merge
from pickle import TRUE
from typing import Dict
from multiprocessing import Pool, set_start_method

from sqlalchemy import false
from postprocessing.postprocessing_utils import postprocessing
from objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
from gen_spectra import get_precursor
from alignment import alignment
from utils import ppm_to_da, to_percent, hashable_boundaries, is_json, is_file
import utils
from scoring import scoring
from preprocessing import digestion, merge_search, preprocessing_utils, clustering, evaluation
import database
from sqlite import database_file
from file_io import JSON
import objects
import time
import multiprocessing as mp
import copy
import json
import os
from scoring.scoring import second_scoring
from alignment.alignment import find_alignments



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
    
def adjust_for_truth_and_fall_off(spectrum,truth,filtered_b,filtered_y,b_results,keep_b_count,y_results,keep_y_count,fall_off):
    _id = spectrum.id
    truth_seq = truth[_id]['sequence']
    is_hybrid = truth[_id]['hybrid']
    if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, filtered_b, filtered_y):
        handle_DEV_truth(filtered_b,filtered_y,b_results,keep_b_count,y_results,keep_y_count,fall_off,_id,is_hybrid,truth_seq,spectrum)
        return Alignments(spectrum, [])

def id_spectrum(spec_num, input_spectrum: Spectrum, db: Database, matched_masses_b, matched_masses_y, ppm_tolerance, precursor_tolerance, location):
    b_hits,y_hits = create_hits(spec_num,input_spectrum,matched_masses_b,matched_masses_y,False,location)
    for ion in "by":
        clusters = clustering.create_clusters(ion, b_hits, y_hits)
        if ion ==  'b':
            b_sorted_clusters = clustering.Score_clusters(ion, clusters)
        else:
            y_sorted_clusters = clustering.Score_clusters(ion, clusters)

    merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
    merged_seqs.sort(key = lambda x: x[0], reverse = True)
    prec_tol = utils.ppm_to_da(input_spectrum.precursor_mass, precursor_tolerance)
    merged_seqs = clustering.filter_by_precursor(merged_seqs, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    merged_seqs = clustering.filter_by_missing_mass(db, merged_seqs, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    
    hybrid_merged = clustering.get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    hybrid_merged = clustering.filter_by_precursor(hybrid_merged, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)
    hybrid_merged = clustering.filter_by_missing_mass(db, hybrid_merged, input_spectrum.precursor_mass, prec_tol, input_spectrum.precursor_charge)

    merged_top = clustering.combine_merges(merged_seqs, hybrid_merged, 50)
    
    alignments = alignment.find_alignments(merged_top, input_spectrum.precursor_mass, input_spectrum.precursor_charge, prec_tol, db)

    rescored_alignments = scoring.second_scoring(alignments, input_spectrum, ppm_tolerance)
    rescored_alignments = sorted(rescored_alignments, key = lambda x: (x[0], x[1]), reverse=True)
    
    postprocessed_alignments = postprocessing(rescored_alignments, db)
    # raw_results = id_spectrum(spectrum, db, b_hits, y_hits, ppm_tolerance, precursor_tolerance,n,digest_type=digest,truth=truth, fall_off=fall_off)
    return postprocessed_alignments

def mp_id_spectrum(input_q: mp.Queue, db_copy: Database, results: dict, fall_off: dict = None, truth: dict = None):
    while True:
        next_entry = input_q.get(True)
        if next_entry == 'exit':
            return 
        if truth is not None and fall_off is not None:
            _id = next_entry.spectrum.id 
            truth_seq = truth[_id]['sequence']
            is_hybrid = truth[_id]['hybrid']
            if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, next_entry.b_hits, next_entry.y_hits):
                metadata = {
                    'initial_b_candidates': next_entry.b_hits, 
                    'initial_y_candidates': next_entry.y_hits
                }
                fall_off[_id] = DEVFallOffEntry(
                    is_hybrid, truth_seq, 'mass_matching', metadata
                )
                results[_id] = Alignments(next_entry.spectrum, [])
                continue
        results[next_entry.spectrum.id] = id_spectrum(next_entry.spectrum, 
            db_copy, next_entry.b_hits, next_entry.y_hits, 
            next_entry.ppm_tolerance, next_entry.precursor_tolerance)
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
def create_hits(spec_num,spectrum,matched_masses_b,matched_masses_y,DEBUG,location):

    if DEBUG:
        filename = "spec_" + str(spec_num) + "_"
        if utils.find_dir(filename + 'b_hits.txt', location) and utils.find_dir(filename + 'y_hits.txt', location):
            b_hits, y_hits = get_hits_from_file(os.path.join(location, 'b_hits.txt'), os.path.join(location, 'y_hits.txt'))
    if not DEBUG or not (utils.find_dir(filename + 'b_hits.txt', location) and utils.find_dir(filename + 'y_hits.txt', location)):
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
        # write_hits(b_hits, y_hits, location)
    return b_hits, y_hits

def align_on_single_core(spectra,matched_masses_b,matched_masses_y,db,ppm_tolerance,precursor_tolerance,results,DEBUG,location):
    for i, input_spectrum in enumerate(spectra):
        print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')
        results[input_spectrum.id]=id_spectrum(i,input_spectrum, db, matched_masses_b, matched_masses_y, ppm_tolerance, precursor_tolerance, location)

def align_on_multi_core(cores,mp_id_spectrum,db,spectra,boundaries,matched_masses_b,matched_masses_y,ppm_tolerance,precursor_tolerance,results,DEBUG,location):
    print('Initializing other processors...')
    results = mp.Manager().dict()
    q = mp.Manager().Queue()
    num_processes = cores
    ps = [ #This code makes the subprocesses and is what takes a while with spinning up the cores
        mp.Process(
            target=mp_id_spectrum, 
            args=(q, copy.deepcopy(db), results)
        ) for _ in range(num_processes) 
    ]

    for p in ps: #starts process and runs first line of mp_id_spectrum
        p.start()    
    print('start each of the process Done.')
    for i, input_spectrum in enumerate(spectra):
        print(f'\rStarting job for {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='')
        b_hits,y_hits = create_hits(i,input_spectrum,matched_masses_b,matched_masses_y,DEBUG,location)

        o = MPSpectrumID(
            b_hits, 
            y_hits, 
            input_spectrum, 
            ppm_tolerance, 
            precursor_tolerance
        )
        
        q.put(o)

    while len(results) < len(spectra):
        print(f'\rCreating an alignment for {len(results)}/{len(spectra)} [{to_percent(len(results), len(spectra))}%]', end='')
        time.sleep(1)

    [q.put('exit') for _ in range(num_processes)]

    for p in ps:
        p.join()
        
    return results

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

class alignment_info:
    def __init__(self, max_peptide_len, location, precursor_tolerance, database, ppm_tolerance, results_len, new): #This is like the named tuple
        self.max_pep_len = max_peptide_len
        self.write_path = location
        self.prec_tol = precursor_tolerance
        self.ppm_tol = ppm_tolerance
        self.db = database
        self.results_len = results_len
        self.make_new = new
    
    def __call__(self, spectrum):
        print(f'\rCreating an alignment for {spectrum.num}/{self.results_len} [{to_percent(spectrum.num, self.results_len)}%]', end='')
        input_list = spectrum.mz_values
        matched_masses_b, matched_masses_y = merge_search.modified_match_masses(input_list, self.db, self.max_pep_len, self.ppm_tol, self.make_new)
        
        #Matched masses data is of form (mass, start, end, ion_int, charge, protein_num)
        # hit_time = time.time()
        b_hits,y_hits = create_hits(spectrum.num,spectrum,matched_masses_b,matched_masses_y,True,self.write_path)
        # hit_time = time.time()-hit_time
        print("hits took:", hit_time)
        for ion in "by":
            clusters = clustering.create_clusters(ion, b_hits, y_hits)
            if ion ==  'b':
                b_sorted_clusters = clustering.Score_clusters(ion, clusters, self.db.proteins)
            else:
                y_sorted_clusters = clustering.Score_clusters(ion, clusters, self.db.proteins)

        merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
        merged_seqs.sort(key = lambda x: x[0], reverse = True)
        prec_tol = ppm_to_da(spectrum.precursor_mass, self.prec_tol)
        merged_seqs = clustering.filter_by_precursor(merged_seqs, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        merged_seqs = clustering.filter_by_missing_mass(self.db, merged_seqs, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        
        hybrid_merged = clustering.get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        hybrid_merged = clustering.filter_by_precursor(hybrid_merged, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        hybrid_merged = clustering.filter_by_missing_mass(self.db, hybrid_merged, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)

        merged_top = clustering.combine_merges(merged_seqs, hybrid_merged, 50)
        
        alignments = find_alignments(merged_top, spectrum.precursor_mass, spectrum.precursor_charge, prec_tol, self.db)

        rescored_alignments = second_scoring(alignments, spectrum, self.ppm_tol)
        rescored_alignments = sorted(rescored_alignments, key = lambda x: (x[0], x[1]), reverse=True)
        
        postprocessed_alignments = postprocessing(rescored_alignments, self.db)
        # raw_results = id_spectrum(spectrum, db, b_hits, y_hits, ppm_tolerance, precursor_tolerance,n,digest_type=digest,truth=truth, fall_off=fall_off)
        return postprocessed_alignments

def align(numcores, spectra, location, precursor_tolerance, db, ppm_tolerance, max_peptide_len,new):
    set_start_method('forkserver')
    p = Pool(numcores)
    y, spec_nums = [], []
    now = time.time()
    [spec_nums.append(i) for i in range(0, len(spectra))]
    x = alignment_info(max_peptide_len,location,precursor_tolerance,db,ppm_tolerance,len(spectra),new)
    y = p.map(x, spectra) #function can only take 1 input so make object
    p.close()
    p.join()
    print("On 16 cores", time.time() - now)
    return y
    
def id_spectra(spectra_files: list, db: database, verbose: bool = True, 
    min_peptide_len: int = 5, max_peptide_len: int = 10, peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0,ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10, digest: str = '',cores: int = 1,
    n: int = 5,DEBUG: bool = False, truth_set: str = "", output_dir: str = ''):
    truth = None
    if is_json(truth_set) and is_file(truth_set):
        DEV = True
        truth = json.load(open(truth_set, 'r'))
    fall_off = None
    make_new = False
    if DEBUG:
        filepath = os.path.abspath(os.path.join("data", "NOD2_E3_results.ssv"))
        correct_sequences = evaluation.generate_truth_set(filepath)
        correct_sequences = [correct_sequences[0], correct_sequences[100], correct_sequences[702]]
    else:
        correct_sequences = []
    verbose and print('Loading spectra...')
    spectra, boundaries = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
    verbose and print('Loading spectra Done')
    dirname = os.path.dirname(os.path.abspath(__file__))
    location = os.path.join(dirname, 'intermediate_files')
    # if utils.find_dir('matched_masses_b.txt', location) and utils.find_dir('matched_masses_y.txt', location):
    #     print("getting matched_masses from file...")
    #     matched_masses_b, matched_masses_y = merge_search.get_from_file(os.path.join(location, 'matched_masses_b.txt'), os.path.join(location, 'matched_masses_y.txt'), True)
    #     print("getting matched_masses from file Done")
    # If needed:
        # merge_search.make_database_file #No file with that name yet    
        
    # build whole database of all mass -> list of tuples of (info)
    # query for input mass?
    # input mass -> list of tuples of (info)
    #input_mass -> [(relevent info), ...]
    
    # TODO
    #matched_masses_b, matched_masses_y = merge_search.match_masses_using_webservice(boundaries, ppm_tolerance)
    results = {}
    DEV = False
    if DEV:
        handle_DEV_setup(truth)
    results = align(cores,spectra,location,precursor_tolerance,db,ppm_tolerance,max_peptide_len,make_new)
    # if cores == 1:
    #     align_on_single_core(spectra,matched_masses_b,matched_masses_y,db,ppm_tolerance,precursor_tolerance,results,DEBUG,location)
    # else:
    #     results = align_on_multi_core(cores,mp_id_spectrum,db,spectra,boundaries,matched_masses_b,matched_masses_y,ppm_tolerance,precursor_tolerance,results,DEBUG,location)
    if DEV:
        handle_DEV_result(output_dir,fall_off,cores)
    return results

