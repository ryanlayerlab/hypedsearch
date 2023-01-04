from multiprocessing import Pool, set_start_method

from postprocessing.postprocessing_utils import postprocessing
from objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
from alignment import alignment
from utils import ppm_to_da, to_percent, is_json, is_file
import utils
from preprocessing import merge_search, preprocessing_utils, clustering, evaluation
import database
from file_io import JSON
import objects
import time
import multiprocessing as mp
import json
import os
from gen_spectra import convert_precursor_to_ion, calc_masses
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
        
def check_duplicates(merged_seqs, hybrid_seqs):
    merged_set = set(merged_seqs)
    hybrid_set = set(hybrid_seqs)
    print(len(merged_set), len(merged_seqs))
    print(len(hybrid_set), len(hybrid_seqs))

def check_duplicates_cross(merged_seqs, hybrid_seqs):
    overlapped_seqs = []
    for seq in merged_seqs:
        if seq in hybrid_seqs:
            overlapped_seqs.append(seq)
    [print(x) for x in overlapped_seqs]
    
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
    #want code that takes in sequences and tells me where this sequence lives in the sorted clusters
    b_target_starts, y_target_ends = find_sequence(b_sequence, y_sequence, b_pid, y_pid, protein_list)
    print("\n")
    
    print("For b:")
    for i, cluster in enumerate(b_sorted_clusters):
        # if cluster[2] in b_target_starts: #for debugging
        #     print(i, cluster)
        if cluster[1] == b_pid and cluster[2] in b_target_starts: #need to check which dictates start position
            print(i, cluster)
    print("\n For y:")
    for i, cluster in enumerate(y_sorted_clusters):
        # if cluster[2] in b_target_starts: #for debugging
        #     print(i, cluster)
        if cluster[1] == y_pid and cluster[3] in y_target_ends: 
            print(i, cluster)
    
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
        total_time = time.time()
        input_list = spectrum.mz_values
        matched_masses_b, matched_masses_y = merge_search.modified_match_masses(input_list, self.db, self.max_pep_len, self.ppm_tol, self.make_new)
                
        #Matched masses data is of form (mass, start, end, ion_int, charge, protein_num)
        hit_time = time.time()
        b_hits,y_hits = create_hits(spectrum.num,spectrum,matched_masses_b,matched_masses_y,True,self.write_path)
        hit_time = time.time()-hit_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Hits took:" + '\t' + str(hit_time) + "\n")
        converted_b, converted_y = convert_precursor_to_ion(spectrum.precursor_mass, spectrum.precursor_charge)
        for ion in "by":
            cluster_time = time.time()
            clusters = clustering.create_clusters(ion, b_hits, y_hits)
            cluster_time = time.time() - cluster_time
            with open('Timing_data.txt', 'a') as t:
                t.write("Clusters took:" + '\t' + str(cluster_time) + "\n")
            if ion ==  'b':
                cluster_time = time.time()
                b_sorted_clusters = clustering.Score_clusters(ion, clusters, self.max_pep_len, converted_b)
                cluster_time = time.time() - cluster_time
                with open('Timing_data.txt', 'a') as t:
                    t.write("Scoring b clusters took:" + '\t' + str(cluster_time) + "\n")
            else:
                cluster_time = time.time()
                y_sorted_clusters = clustering.Score_clusters(ion, clusters, self.max_pep_len, converted_y)
                cluster_time = time.time() - cluster_time
                with open('Timing_data.txt', 'a') as t:
                    t.write("Scoring y clusters took:" + '\t' + str(cluster_time) + "\n")


        # find_target_clusters(b_sorted_clusters, y_sorted_clusters, "DLKIIWNKTKH", "DLKIIWNKTKH", 140, 140, self.db.proteins)


        start_time = time.time()
        merged_seqs = clustering.Ryan_merge(b_sorted_clusters, y_sorted_clusters)
        merged_seqs.sort(key = lambda x: x[0], reverse = True)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Ryan merging and sorting took:" + '\t' + str(end_time) + "\n")
        prec_tol = ppm_to_da(spectrum.precursor_mass, self.prec_tol)
        start_time = time.time()
        merged_seqs = clustering.filter_by_precursor(merged_seqs, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge, self.max_pep_len)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Precursor filtering took:" + '\t' + str(end_time) + "\n")
        # start_time = time.time()
        # merged_seqs = clustering.filter_by_missing_mass(self.db, merged_seqs, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        # end_time = time.time() - start_time
        # with open('Timing_data.txt', 'a') as t:
        #     t.write("Missing mass filtering took:" + '\t' + str(end_time) + "\n")
        
        start_time = time.time()
        hybrid_merged = clustering.get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Finding hybrid merges took:" + '\t' + str(end_time) + "\n")
        start_time = time.time()
        hybrid_merged = clustering.filter_by_precursor(hybrid_merged, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge, self.max_pep_len)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Filtering hybrids by precursor masss took:" + '\t' + str(end_time) + "\n")


        # check_duplicates(merged_seqs, hybrid_merged)
        # check_duplicates_cross(merged_seqs, hybrid_merged)
        # print("Naturals:")
        # [print(x) for x in merged_seqs[-10:]]
        # [print(x) for x in merged_seqs[:10]]
        # print("Hybrids:")
        # [print(x) for x in hybrid_merged[-10:]]
        # [print(x) for x in hybrid_merged[:10]]
        
        start_time = time.time()
        natural_alignments, hybrid_alignments = find_alignments(merged_seqs, hybrid_merged, spectrum.precursor_mass, spectrum.precursor_charge, prec_tol, self.db, self.prec_tol)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Making alignments took:" + '\t' + str(end_time) + "\n")

        start_time = time.time()
        rescored_naturals, rescored_hybrids = second_scoring(natural_alignments, hybrid_alignments, spectrum, self.ppm_tol, self.db.proteins, self.max_pep_len)
        rescored_naturals = sorted(rescored_naturals, key = lambda x: (x[0], x[1]), reverse=True)
        rescored_hybrids = sorted(rescored_hybrids, key = lambda x: (x[0], x[1]), reverse=True)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Second round of scoring and sorting took:" + '\t' + str(end_time) + "\n")
            
        # check_top_location(rescored_naturals[0], rescored_hybrids[0], merged_seqs, hybrid_merged)
        
        rescored_alignments = sorted(rescored_naturals + rescored_hybrids, key = lambda x: (x[0], x[1]), reverse = True)
        
        start_time = time.time()
        postprocessed_alignments = postprocessing(rescored_alignments, self.db)
        end_time = time.time() - start_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Postprocessing took:" + '\t' + str(end_time) + "\n")
        # raw_results = id_spectrum(spectrum, db, b_hits, y_hits, ppm_tolerance, precursor_tolerance,n,digest_type=digest,truth=truth, fall_off=fall_off)
        total_time = time.time() - total_time
        with open('Timing_data.txt', 'a') as t:
            t.write("Analysis of spectrum " + str(spectrum.num) +  " took:" + "\t" + str(total_time) + "\n")
        return postprocessed_alignments

def align(numcores, spectra, location, precursor_tolerance, db, ppm_tolerance, max_peptide_len,new):
    set_start_method('forkserver')
    p = Pool(numcores)
    y, spec_nums = [], []
    now = time.time()
    [spec_nums.append(i) for i in range(0, len(spectra))]
    x = alignment_info(max_peptide_len,location,precursor_tolerance,db,ppm_tolerance,len(spectra),new)
    y = p.map(x, spectra)
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

