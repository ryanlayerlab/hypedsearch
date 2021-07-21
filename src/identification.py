from objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
#from cppModules import gen_spectra
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
# top results to keep for creating an alignment
global alignment_times
alignment_times = []
global b_scoring_times
b_scoring_times = []
global y_scsoring_times
y_scoring_times = []
global filter_times
filter_times = []

TOP_X = 50

def id_spectrum(
    spectrum: Spectrum, 
    db: Database,
    b_hits: dict, 
    y_hits: dict,
    ppm_tolerance: int, 
    precursor_tolerance: int,
    n: int,
    digest_type: str = '',
    truth: dict = None, 
    fall_off: dict = None, 
    is_last: bool = False
    ) -> Alignments:
    '''Given the spectrum and initial hits, start the alignment process for 
    the input spectrum

    :param spectrum: observed spectrum in question
    :type spectrum: Spectrum
    :param db: Holds all the source sequences
    :type db: Database
    :param b_hits: all k-mers found from the b-ion search
    :type b_hits: list
    :param y_hits: all k-mers found from the y-ion search
    :type y_hits: list
    :param ppm_tolerance: the parts per million error allowed when trying to match masses
    :type ppm_tolerance: int
    :param precursor_tolerance: the parts per million error allowed when trying to match
        precursor masses
    :type percursor_tolerance: int
    :param n: the number of alignments to save
    :type n: int
    :param digest_type: the digest performed on the sample
        (default is '')
    :type digest_type: str
    :param truth: a set of id keyed spectra with the desired spectra. A better description of what this looks like can be 
        seen in the param.py file. If left None, the program will continue normally
        (default is None)
    :type truth: dict
    :param fall_off: only works if the truth param is set to a dictionary. This is a dictionary (if using multiprocessing, 
        needs to be process safe) where, if a sequence loses the desired sequence, a key value pair of spectrum id, 
        DevFallOffEntry object are added to it. 
        (default is None)
    :type fall_off: dict
    :param is_last: Only works if DEV is set to true in params. If set to true, timing evaluations are done. 
        (default is False)
    :type is_last: bool

    :returns: Alignments for the spectrum. If no alignment can be created, and empty Alignments object is inserted
    :rtype: Alignments
    '''
    start_time = time.time()
    # convert the ppm tolerance of the precursor to an int for the rest of the time
    precursor_tolerance = utils.ppm_to_da(spectrum.precursor_mass, precursor_tolerance)

    # score and sort these results
    score_b_start = time.time()
    b_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='b'))
        ) for kmer in b_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    b_scoring_times.append(time.time() - score_b_start)
    score_y_start = time.time()
    y_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='y'))
        ) for kmer in y_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    y_scoring_times.append(time.time() - score_y_start)

    filter_start = time.time()
    # filter out the results
    # 1. take all non-zero values 
    # 2. either take the TOP_X or if > TOP_X have the same score, all of those values
    filtered_b, filtered_y = [], []

    # find the highest b and y scores
    max_b_score = max([x[1] for x in b_results])
    max_y_score = max([x[1] for x in y_results])

    # count the number of kmers that have the highest value
    num_max_b = sum([1 for x in b_results if x[1] == max_b_score])
    num_max_y = sum([1 for x in y_results if x[1] == max_y_score])

    # if we have more than TOP_X number of the highest score, take all of them
    keep_b_count = max(TOP_X, num_max_b)
    keep_y_count = max(TOP_X, num_max_y)

    # take the afformentioned number of results that > than zero
    filtered_b = [x[0] for x in b_results[:keep_b_count] if x[1] > 0]
    filtered_y = [x[0] for x in y_results[:keep_y_count] if x[1] > 0]

    filter_times.append(time.time() - filter_start)

    # if fall off and truth are not none, check to see that we can still make the truth seq
    if truth is not None and fall_off is not None:

        # pull out id, hybrid, and truth seq to make it easier
        _id = spectrum.id
        truth_seq = truth[_id]['sequence']
        is_hybrid = truth[_id]['hybrid']

        if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, filtered_b, filtered_y):

            # add some metadata about what we kept and what fell off
            metadata = {
                'top_x_b_hits': filtered_b, 
                'top_x_y_hits': filtered_y, 
                'excluded_b_hits': [x[0] for x in b_results[keep_b_count:]],
                'excluded_y_hits': [x[0] for x in y_results[keep_y_count:]], 
                'cut_off_b_score': b_results[keep_b_count - 1][1], 
                'cut_off_y_score': y_results[keep_y_count - 1][1]
            }

            # make dev fall off object and add to fall off
            fall_off[_id] = DEVFallOffEntry(
                is_hybrid, 
                truth_seq, 
                'top_x_filtering', 
                metadata
            )

            # skip this entry all together
            return Alignments(spectrum, [])

    # create an alignment for the spectrum
    align_start = time.time()
    alignments = alignment.attempt_alignment(
        spectrum, 
        db, 
        filtered_b, 
        filtered_y, 
        ppm_tolerance=ppm_tolerance, 
        precursor_tolerance=precursor_tolerance,
        n=n, 
        truth=truth, 
        fall_off=fall_off, 
        is_last=is_last
    )
    TOT_ALIGNMENT = time.time() - align_start
    alignment_times.append(TOT_ALIGNMENT)
    return alignments


def id_spectra(
    spectra_files: list, 
    database_file: database, 
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
    DEBUG: bool = False, 
    truth_set: str = '', 
    output_dir: str = ''
) -> dict:
    DEV = False
    truth = None

    if is_json(truth_set) and is_file(truth_set):
        DEV = True
        truth = json.load(open(truth_set, 'r'))

    fall_off = None
    database_start = time.time()
    # build/load the database
    verbose and print('Loading database...')
    db = database_file
    verbose and print('Loading database Done')
    #instrumentation
    time_to_build_database = time.time() - database_start
    
    # load all of the spectra
    spectra_start = time.time()
    verbose and print('Loading spectra...')
    spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Loading spectra Done')
    #instrumentation
    time_to_load_in_spectra = time.time() - spectra_start


    # get the boundary -> kmer mappings for b and y ions
    mapping_start = time.time()
    matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, db, max_peptide_len)
    #instrumentation
    time_to_map_boundaries_to_kmers = time.time() - mapping_start

    # keep track of the alingment made for every spectrum
    results = {}

    if DEV:
        fall_off = {}
        fall_off = mp.Manager().dict()
        truth = mp.Manager().dict(truth)

    # if we only get 1 core, don't do the multiprocessing bit
    if cores == 1:
        # go through and id all spectra
        for i, spectrum in enumerate(spectra):

            print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')

            # get b and y hits
            b_hits, y_hits = [], []
            for mz in spectrum.spectrum:

                # get the correct boundary
                mapped = mz_mapping[mz]
                b = boundaries[mapped]
                b = hashable_boundaries(b)

                if b in matched_masses_b:
                    b_hits += matched_masses_b[b]

                if b in matched_masses_y:
                    y_hits += matched_masses_y[b]

            is_last = DEBUG and i == len(spectra) - 1

            # pass it into id_spectrum
            #TODO: HERE
            raw_results = id_spectrum(
                spectrum, db, 
                b_hits, y_hits, 
                ppm_tolerance, precursor_tolerance,
                n,digest_type=digest,
                truth=truth, fall_off=fall_off, 
                is_last=is_last)
            results[spectrum.id]=raw_results.alignments
    else:
        
        multiprocessing_start = time.time()
        print('Initializing other processors...')
        results = mp.Manager().dict()

        if DEV:
            fall_off = mp.Manager().dict()
            truth = mp.Manager().dict(truth)

        # start up processes and queue for parallelizing things
        q = mp.Manager().Queue()
        num_processes = cores
        ps = [
            mp.Process(
                target=mp_id_spectrum, 
                args=(q, copy.deepcopy(db), results, fall_off, truth)
            ) for _ in range(num_processes) 
        ]

        # start each of the process
        for p in ps:
            p.start()
        #instrumentation
        time_to_spin_up_cores = time.time() - multiprocessing_start          
        print('start each of the process Done.')

        # go through and id all spectra
        for i, spectrum in enumerate(spectra):
            print(f'\rStarting job for {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='')
            # get b and y hits
            b_hits, y_hits = [], []
            for mz in spectrum.spectrum:

                # get the correct boundary
                mapped = mz_mapping[mz]
                b = boundaries[mapped]
                b = hashable_boundaries(b)

                if b in matched_masses_b:
                    b_hits += matched_masses_b[b]

                if b in matched_masses_y:
                    y_hits += matched_masses_y[b]

            # create a named tuple to put in the database
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

        # now send 'exit' message to all our processes
        [q.put('exit') for _ in range(num_processes)]

        # join them
        for p in ps:
            p.join()

    # if we have set DEV, we need to dump this to a json
    if DEV:
        output_dir = output_dir + '/' if output_dir[-1] != '/' else output_dir

        safe_write_fall_off = {}

        # we need to convert all our DEVFallOffEntries to dicts
        for k, v in fall_off.items():
            safe_write_fall_off[k] = v._asdict()

        JSON.save_dict(output_dir + 'fall_off.json', safe_write_fall_off)
        #instrumentation
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
    '''Multiprocessing function for to identify a spectrum. Each entry in the 
    input_q must be a MPSpectrumID object

    :param input_q: a queue to pull MPSpectrumID objects from for analysis
    :type input_q: mp.Queue
    :param db_copy: a copy of the original database for alignments
    :type db_copy: Database
    :param results: a multiprocesses safe dictionary to save the alignments in
    :type results: dict
    :param truth_set: dictionary containing all the desired alignments to make. 
        The format of the file is {spectrum_id: {'sequence': str, 'hybrid': bool, 'parent': str}}. 
        If left as None, the program will continue as normal
        (default is None)
    :type truth_set: dict
    :param fall_off: only used if the truth_set param is set to a valid json. Must be a multiprocess
        safe dictionary to store the fall off information to
    :type fall_off: dict

    :returns: None
    :rtype: None
    '''
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
         #TODO: HERE
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