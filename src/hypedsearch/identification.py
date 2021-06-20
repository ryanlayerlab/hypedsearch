from objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
#from cppModules import gen_spectra
import gen_spectra
from alignment import alignment
from utils import ppm_to_da, to_percent, overlap_intervals, hashable_boundaries, is_json, is_file
import utils
from scoring import scoring, mass_comparisons
from preprocessing import digestion, merge_search, preprocessing_utils
from file_io import JSON
import os
from dataclasses import dataclass
import time
import multiprocessing as mp
import copy
import json
import database

TOP_X = 50

@dataclass
class Id_Spectra_Arguments:
    spectra_files: list = [], 
    database_file: database = None,     
    verbose: bool = True, 
    min_peptide_len: int = 5, 
    max_peptide_len: int = 20, 
    peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0,
    ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10, 
    digest: str = '',
    cores: int = 1,
    n: int = 5,
    is_debug: bool = False, 
    truth_set: str = '', 
    output_dir: str = ''

@dataclass
class Id_Spectrum_Arguments:
    spectrum: Spectrum = None, 
    db: Database = Database,
    b_hits: dict = [], 
    y_hits: dict = [],
    ppm_tolerance: int = 0, 
    precursor_tolerance: int = 0, 
    n: int = 0,
    digest_type: str = '',
    truth: dict = None, 
    fall_off: dict = None, 
    is_last: bool = False

def id_spectrum(idsa: Id_Spectrum_Arguments) -> Alignments:
    precursor_tolerance = utils.ppm_to_da(idsa.spectrum.precursor_mass, idsa.precursor_tolerance)
    b_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(idsa.spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='b'))
        ) for kmer in idsa.b_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    y_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(idsa.spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='y'))
        ) for kmer in idsa.y_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    filtered_b, filtered_y = [], []
    max_b_score = max([x[1] for x in b_results])
    max_y_score = max([x[1] for x in y_results])
    num_max_b = sum([1 for x in b_results if x[1] == max_b_score])
    num_max_y = sum([1 for x in y_results if x[1] == max_y_score])
    keep_b_count = max(TOP_X, num_max_b)
    keep_y_count = max(TOP_X, num_max_y)
    filtered_b = [x[0] for x in b_results[:keep_b_count] if x[1] > 0]
    filtered_y = [x[0] for x in y_results[:keep_y_count] if x[1] > 0]
    if idsa.truth is not None and idsa.fall_off is not None:
        _id = idsa.spectrum.id
        truth_seq = idsa.truth[_id]['sequence']
        is_hybrid = idsa.truth[_id]['hybrid']

        if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, filtered_b, filtered_y):
            metadata = {
                'top_x_b_hits': filtered_b, 
                'top_x_y_hits': filtered_y, 
                'excluded_b_hits': [x[0] for x in b_results[keep_b_count:]],
                'excluded_y_hits': [x[0] for x in y_results[keep_y_count:]], 
                'cut_off_b_score': b_results[keep_b_count - 1][1], 
                'cut_off_y_score': y_results[keep_y_count - 1][1]
            }
            idsa.fall_off[_id] = DEVFallOffEntry(
                is_hybrid, 
                truth_seq, 
                'top_x_filtering', 
                metadata
            )
            return Alignments(idsa.spectrum, [])

    aaa = alignment.Attempt_Alignment_Arguments(spectrum=idsa.spectrum,db=idsa.db,b_hits=idsa.b_hits,y_hits=idsa.y_hits,
        n=idsa.n,ppm_tolerance=idsa.ppm_tolerance, precursor_tolerance=precursor_tolerance,digest_type=idsa.digest_type,
        DEBUG=idsa.DEBUG,is_last=idsa.is_last,truth=idsa.truth,fall_off=idsa.fall_off)            
    alignment_results = alignment.attempt_alignment(aaa)
    return alignment_results

def load_all_spectra(spectra_files,ppm_tolerance,peak_filter,relative_abundance_filter,verbose):
    verbose and print('Loading spectra...')
    lsa = preprocessing_utils.Load_Spectra_Arguments(spectra_files=spectra_files,ppm_tolerance=ppm_tolerance,peak_filter=peak_filter,relative_abundance_filter=relative_abundance_filter)
    lsr = preprocessing_utils.load_spectra(lsa)
    verbose and print('Loading spectra Done')
    return lsr.all_spectra, lsr.boundaries, lsr.mz_mapping

def create_alignment_single_core(spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,is_debug,results,db,ppm_tolerance,precursor_tolerance,n,digest,truth,fall_off):
    for i, spectrum in enumerate(spectra):
        print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')
        b_hits, y_hits = [], []
        for mz in spectrum.spectrum:
            mapped = mz_mapping[mz]
            b = boundaries[mapped]
            b = hashable_boundaries(b)
            if b in matched_masses_b:
                b_hits += matched_masses_b[b]
            if b in matched_masses_y:
                y_hits += matched_masses_y[b]
        is_last = is_debug and i == len(spectra) - 1
        results[spectrum.id] = id_spectrum(
            spectrum, 
            db, 
            b_hits, 
            y_hits, 
            ppm_tolerance, 
            precursor_tolerance,
            n,
            digest_type=digest,
            truth=truth, 
            fall_off=fall_off, 
            is_last=is_last
        )
def create_alignment_multi_core(is_dev,truth,cores,mp_id_spectrum,db,spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,ppm_tolerance,precursor_tolerance,n,digest):
    print('Initializing other processors...')
    results = mp.Manager().dict()
    if is_dev:
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
    print('Done.')
    for i, spectrum in enumerate(spectra):
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

def set_up_for_dev(truth_set):
    if is_json(truth_set) and is_file(truth_set):
        truth= json.load(open(truth_set, 'r'))
        return True, truth
    else:       
        return False, None

def update_truth(is_dev,truth):
    if is_dev:
        fall_off = {}
        fall_off = mp.Manager().dict()
        truth = mp.Manager().dict(truth)

def output_for_dev(is_dev,output_dir,fall_off):
    if is_dev:
        output_dir = output_dir + '/' if output_dir[-1] != '/' else output_dir
        safe_write_fall_off = {}
        for k, v in fall_off.items():
            safe_write_fall_off[k] = v._asdict()
        JSON.save_dict(output_dir + 'fall_off.json', safe_write_fall_off)

def id_spectra(idsa:utils.Id_Spectra_Arguments) -> dict:
    is_dev, truth = set_up_for_dev(idsa.truth_set)
    fall_off = None
    spectra, boundaries, mz_mapping  = load_all_spectra(idsa.spectra_files,idsa.ppm_tolerance,idsa.peak_filter,idsa.relative_abundance_filter,idsa.verbose)
    matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, idsa.database_file, idsa.max_peptide_len)
    results = {}
    update_truth(is_dev,truth)
    if idsa.cores == 1:
        create_alignment_single_core(spectra,mz_mapping,boundaries,matched_masses_b,matched_masses_y,
            idsa.is_debug,results,db,idsa.ppm_tolerance,idsa.precursor_tolerance,idsa.n,idsa.digest,truth,fall_off)
    else:
        create_alignment_multi_core(is_dev,truth,idsa.cores,idsa.mp_id_spectrum,db,spectra,mz_mapping,
            boundaries,matched_masses_b,matched_masses_y,idsa.ppm_tolerance,idsa.precursor_tolerance,idsa.n,idsa.digest)
    output_for_dev(is_dev,idsa.output_dir,fall_off)
    return results

def mp_id_spectrum(
    input_q: mp.Queue, 
    db_copy: Database, 
    results: dict, 
    fall_off: dict = None, 
    truth: dict = None
    ) -> None:
    while True:

        # wait to get something from the input queue
        next_entry = input_q.get(True)

        # if it says 'exit', quit
        if next_entry == 'exit':
            return 

        # if fall off is not none, see if we have the correct value in here
        if truth is not None and fall_off is not None:

            # pull out the id to make it easier
            _id = next_entry.spectrum.id 

            # pull out the truth sequence and the hybrid bool 
            truth_seq = truth[_id]['sequence']
            is_hybrid = truth[_id]['hybrid']

            # see if we still have the correct results
            if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, next_entry.b_hits, next_entry.y_hits):

                # add some metadata. Add the b and y hits we DID have
                metadata = {
                    'initial_b_candidates': next_entry.b_hits, 
                    'initial_y_candidates': next_entry.y_hits
                }
                
                # create the fall off dev object
                fall_off[_id] = DEVFallOffEntry(
                    is_hybrid, truth_seq, 'mass_matching', metadata
                )
            
                # add an empty entry to the results
                results[_id] = Alignments(next_entry.spectrum, [])
                continue

        # otherwise run id spectrum 
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

def build_load_database(database_file,verbose):
    verbose and print('Loading database...')
    db = database.build(database_file)
    verbose and print('Loading database Done')
    return db

def build_load_spectra_files(spectra_folder_path):
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder_path):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    return spectra_files

def create_id_spectra_arguments(args: dict):
    verbose=True,
    speactra_folder_path = args['spectra_folder']
    database_file_path = args['database_file']
    database_file = build_load_database(database_file_path,verbose)
    spectra_files = build_load_spectra_files(speactra_folder_path)
    min_peptide_len=args['min_peptide_len']
    max_peptide_len=args['max_peptide_len']
    peak_filter=args['peak_filter']
    relative_abundance_filter=args['relative_abundance_filter']
    ppm_tolerance=args['tolerance']
    precursor_tolerance=args['precursor_tolerance']
    digest=args['digest']
    max_cores = max(1, args['cores'])
    cores = min(max_cores, mp.cpu_count() - 1)
    n=args['n'] * 10
    is_debug=args['DEBUG']
    truth_set=args['truth_set']
    output_dir=args['output_dir']
    id_spectra_arguments = Id_Spectra_Arguments(spectra_files=spectra_files,
        database_file=database_file,verbose=verbose,
        min_peptide_len=min_peptide_len,max_peptide_len=max_peptide_len,peak_filter=peak_filter,
        relative_abundance_filter=relative_abundance_filter,ppm_tolerance=ppm_tolerance,
        precursor_tolerance=precursor_tolerance,digest=digest,cores=cores,
        n=n,is_debug=is_debug,truth_set=truth_set,output_dir=output_dir)
    return id_spectra_arguments        