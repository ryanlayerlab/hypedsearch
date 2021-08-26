from objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
import gen_spectra
from alignment import alignment
from utils import ppm_to_da, to_percent, overlap_intervals, hashable_boundaries, is_json, is_file
import utils
from scoring import scoring, mass_comparisons
from preprocessing import digestion, merge_search, preprocessing_utils
import database
from file_io import JSON
import objects
import time
import multiprocessing as mp
import copy
import json

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



def id_spectrum(spectrum: Spectrum, db: Database, b_hits: dict, y_hits: dict,
    ppm_tolerance: int, precursor_tolerance: int,n: int,digest_type: str = '',
    truth: dict = None, fall_off: dict = None, is_last: bool = False):
    precursor_tolerance = utils.ppm_to_da(spectrum.precursor_mass, precursor_tolerance)
    score_b_start = time.time()
    b_intermediate = [(kmer, mass_comparisons.optimized_compare_masses(spectrum.mz_values, gen_spectra.gen_spectrum(kmer, ion='b'))) for kmer in b_hits]
    b_results = sorted(b_intermediate, key=lambda x: (x[1], 1/len(x[0])), reverse=True)
    b_scoring_times.append(time.time() - score_b_start)
    score_y_start = time.time()
    y_intermediate = [(kmer, mass_comparisons.optimized_compare_masses(spectrum.mz_values, gen_spectra.gen_spectrum(kmer, ion='y'))) for kmer in y_hits]
    y_results = sorted(y_intermediate, key=lambda x: (x[1], 1/len(x[0])), reverse=True)
    y_scoring_times.append(time.time() - score_y_start)
    filter_start = time.time()
    filtered_b, filtered_y = [], []
    max_b_score = max([x[1] for x in b_results])
    max_y_score = max([x[1] for x in y_results])
    num_max_b = sum([1 for x in b_results if x[1] == max_b_score])
    num_max_y = sum([1 for x in y_results if x[1] == max_y_score])
    keep_b_count = max(TOP_X, num_max_b)
    keep_y_count = max(TOP_X, num_max_y)
    filtered_b = [x[0] for x in b_results[:keep_b_count] if x[1] > 0]
    filtered_y = [x[0] for x in y_results[:keep_y_count] if x[1] > 0]
    filter_times.append(time.time() - filter_start)
    if truth is not None and fall_off is not None:
        _id = spectrum.id
        truth_seq = truth[_id]['sequence']
        is_hybrid = truth[_id]['hybrid']
        if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, filtered_b, filtered_y):
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
            return Alignments(spectrum, [])

    align_start = time.time()
    alignments = alignment.attempt_alignment(spectrum, db, filtered_b, filtered_y, 
        ppm_tolerance=ppm_tolerance, precursor_tolerance=precursor_tolerance,
        n=n, truth=truth, fall_off=fall_off, is_last=is_last)
    TOT_ALIGNMENT = time.time() - align_start
    alignment_times.append(TOT_ALIGNMENT)
    return alignments


def match_on_single_core(spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,db,ppm_tolerance,precursor_tolerance,n,digest,truth,fall_off,results,DEBUG):
    for i, spectrum in enumerate(spectra):
        print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')
        b_hits, y_hits = [], []
        for mz in spectrum.mz_values:
            mapped = mz_mapping[mz]
            b = boundaries[mapped]
            b = hashable_boundaries(b)
            if b in matched_masses_b:
                b_hits += matched_masses_b[b]
            if b in matched_masses_y:
                y_hits += matched_masses_y[b]
        is_last = DEBUG and i == len(spectra) - 1
        raw_results = id_spectrum(
            spectrum, db, 
            b_hits, y_hits, 
            ppm_tolerance, precursor_tolerance,
            n,digest_type=digest,
            truth=truth, fall_off=fall_off, 
            is_last=is_last)
        results[spectrum.id]=raw_results

def match_on_multi_core(DEV,truth,cores,mp_id_spectrum,db,spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,ppm_tolerance,precursor_tolerance,n,digest):
    multiprocessing_start = time.time()
    print('Initializing other processors...')
    results = mp.Manager().dict()
    if DEV:
        fall_off = mp.Manager().dict()
        truth = mp.Manager().dict(truth)
    q = mp.Manager().Queue()
    num_processes = cores
    ps = [
        mp.Process(
            target=mp_id_spectrum, 
            args=(q, copy.deepcopy(db), results, fall_off, truth)
        ) for _ in range(num_processes) 
    ]
    for p in ps:
        p.start()
    time_to_spin_up_cores = time.time() - multiprocessing_start          
    print('start each of the process Done.')
    for i, spectrum in enumerate(spectra):
        print(f'\rStarting job for {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='')
        b_hits, y_hits = [], []
        for mz in spectrum.spectrum:
            mapped = mz_mapping[mz]
            b = boundaries[mapped]
            b = hashable_boundaries(b)

            if b in matched_masses_b:
                b_hits += matched_masses_b[b]

            if b in matched_masses_y:
                y_hits += matched_masses_y[b]

        o = MPSpectrumID(
            b_hits, 
            y_hits, 
            spectrum, 
            ppm_tolerance, 
            precursor_tolerance, 
            n, 
            digest
        )
        
        q.put(o)

    while len(results) < len(spectra):
        print(f'\rCreating an alignment for {len(results)}/{len(spectra)} [{to_percent(len(results), len(spectra))}%]', end='')
        time.sleep(1)

    [q.put('exit') for _ in range(num_processes)]

    for p in ps:
        p.join()

def id_spectra(spectra_files: list, database: database, verbose: bool = True, 
    min_peptide_len: int = 5, max_peptide_len: int = 20, peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0,ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10, digest: str = '',cores: int = 1,
    n: int = 5,DEBUG: bool = False, truth_set: str = '', output_dir: str = ''):
    DEV = False
    truth = None

    if is_json(truth_set) and is_file(truth_set):
        DEV = True
        truth = json.load(open(truth_set, 'r'))

    fall_off = None
    database_start = time.time()
    verbose and print('Loading database...')
    db = database
    verbose and print('Loading database Done')
    time_to_build_database = time.time() - database_start
    spectra_start = time.time()
    verbose and print('Loading spectra...')
    spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Loading spectra Done')
    time_to_load_in_spectra = time.time() - spectra_start
    mapping_start = time.time()
    matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, db, max_peptide_len)
    time_to_map_boundaries_to_kmers = time.time() - mapping_start
    results = {}
    if DEV:
        fall_off = {}
        fall_off = mp.Manager().dict()
        truth = mp.Manager().dict(truth)
    if cores == 1:
        match_on_single_core(spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,db,ppm_tolerance,precursor_tolerance,n,digest,truth,fall_off,results,DEBUG)
    else:
        match_on_multi_core(DEV,truth,cores,mp_id_spectrum,db,spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,ppm_tolerance,precursor_tolerance,n,digest)
    if DEV:
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
    return results

def mp_id_spectrum(
    input_q: mp.Queue, 
    db_copy: Database, 
    results: dict, 
    fall_off: dict = None, 
    truth: dict = None
    ) -> None:
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
        results[next_entry.spectrum.id] = id_spectrum(
            next_entry.spectrum, 
            db_copy, 
            next_entry.b_hits, 
            next_entry.y_hits, 
            next_entry.ppm_tolerance, 
            next_entry.precursor_tolerance,
            next_entry.n,
            truth, 
            fall_off
        )