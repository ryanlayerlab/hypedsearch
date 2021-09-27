from typing import Tuple
import string
import collections
import operator
from pyteomics import fasta
import os
import pandas as pd
from collections import defaultdict

from collections import namedtuple
import sys
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)
from preprocessing import preprocessing_utils
from scoring import mass_comparisons
from objects import Database
import gen_spectra, utils

from collections import defaultdict
from typing import Iterable
from utils import overlap_intervals, ppm_to_da
import array as arr
from math import ceil, trunc
from posixpath import split

BATCH_SIZE = 300

db_dict_b = defaultdict(set)
db_dict_y = defaultdict(set)
kmer_set = defaultdict(list)

def generate_truth_set(Dataset):
    correct_sequences = []
    abundances = []
    with open(Dataset[1], 'r') as truth_set:
        for q, line in enumerate(truth_set):
            if q != 0:
                split_line = line.split(';')
                correct_sequences.append(split_line[9])

    return correct_sequences

def db_filter(db_file: str, results_file: str, output_fasta: str) -> None:
    '''
    Create the subset of proteins needed for the database search
    
    Inputs:
        db_file:        (str)  the original fasta file
        results_file:   (str)  the results ssv file from spectrumMill
        output_fasta:   (str)  the fasta file to write to
    '''
    
    # load all protiens into a dictionary
    db = {}
    for entry in fasta.read(db_file):
        name = entry.description.split('|')[2]
        name = name[:name.index('OS=')-1]
        name = ' '.join(name.split(' ')[1:])
        db[name.lower()] = entry

    # load the results ssv into a dataframe 
    res_df = pd.read_csv(results_file, sep=';')
        
    print(f'Number of results: {len(res_df.index)}')

    # keep track of those we want
    filtered = []
    for idx, row in res_df.iterrows():
        key = row['entry_name'].lower()
        
        if key not in db:
            continue
            
        filtered.append(db[key])

    filtered = list(set(filtered))
    
    print(f'Number of proteins in database was reduced from {len(db)} to {len(filtered)}')
    
    fasta.write(filtered, output_fasta, file_mode='w')

def root_path():
    return os.path.abspath(os.sep)

def define_data():

    root = root_path()

    # define the list of datasets
    # make it a list of tuples of (mzml, spectrum mill *sv, database, prefix dir)
    Dataset = namedtuple(
        'Dataset', 
        ['spectra_dir', 'spectrumMill_results', 'full_database', 'highest_dir', 'filtered_fasta']
    )


    raw_prefix = os.path.join(root, 'home', 'naco3124', 'jaime_hypedsearch', 'hypedsearch', 'data')


    NOD2_data = Dataset(
        os.path.join(raw_prefix, 'unused_spectra'),
        os.path.join(raw_prefix, 'NOD2_E3_results.ssv'),
        os.path.join(raw_prefix, 'database', 'sample_database.fasta'),
        os.path.join(raw_prefix) + os.path.sep,
        ''
    )

    # BALB3_data = Dataset(
    #     os.path.join(raw_prefix, BALB3_top_dir, 'mzml') + os.path.sep,
    #     os.path.join(raw_prefix, BALB3_top_dir, 'BALB3_E3_results.ssv'),
    #     os.path.join(raw_prefix, 'mouse_database.fasta'),
    #     os.path.join(raw_prefix, BALB3_top_dir) + os.path.sep,
    #     ''
    # )

    # datasets = [NOD2_data, BALB3_data]
    datasets = [NOD2_data]

    # FILTER THE DATA

    updated_datasets = []

    for dataset in datasets:
            
        # make a file name for the output for the filtered fasta file
        output_fasta = os.path.join(dataset.highest_dir, 'filtered_' + dataset.highest_dir.split(os.path.sep)[-1].replace(os.path.sep, '') + '_database.fasta')
            
        # check to see if we've created it before
        if not os.path.isfile(output_fasta):
            db_filter(dataset[2], dataset[1], output_fasta)
        
        updated_datasets.append(dataset._replace(filtered_fasta=output_fasta))

    datasets = updated_datasets
    return datasets

def preprocess_input_spectra(spectra_folder, ppm_tolerance, peak_filter: int = 25, relative_abundance_filter: float = 0.01):
    # returns a list of spectrum objects
    # Pass in a spectra path and output a list of Spectrum objects

    # get all the spectra file names
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    
    # load all of the spectra
    print('Loading spectra...')
    spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter = peak_filter, 
        relative_abundance_filter = relative_abundance_filter
    )
    print('Done')

    return spectra, boundaries, mz_mapping

def isintolerance(val1, val2, tolerance):
    if abs(val1-val2) <= tolerance:
        return True
    else:
        return False

def get_total_length(solution_set) -> int:
    total_length = 0
    for x in solution_set:
        total_length = total_length + len(x)
    return total_length

def get_average_from_set(input_set) -> float:
    set_total = 0
    for x in input_set:
        set_total = set_total + x

    set_average = set_total / len(input_set)
    return set_average

def modified_match_masses(
    spectra_boundaries: list, 
    db: Database, 
    max_pep_len: int = 30
    ) -> Tuple[dict, dict, Database]:
    '''Take in a list of boundaries from observed spectra and return a b and y
    dictionary that maps boundaries -> kmers

    :param spectra_boundaries: boundaries as lists as [lower_bound, upper_bound]
    :type spectra_boundaries: list
    :param db: source proteins
    :type db: Database
    :param max_pep_len: maximum peptide length in k-mer prefetching
    :type max_pep_len: int

    :returns: mapping of b ion masses to k-mers, mapping of y ion masses to 
        k-mers, updated database
    :rtype: (dict, dict, Database)
    '''

    # keep track of all of the good mass matches and kmers
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)

    # estimate the max len
    estimated_max_len = ceil(spectra_boundaries[-1][1] / 57.021464)
    max_len = min(estimated_max_len, max_pep_len)

    # calc the number of batches needed
    num_batches = ceil(len(db.proteins) / BATCH_SIZE)

    # create batches of proteins in the form of (prot name, prot entry)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

    # go through each batch set and create the list representation, merge, and keep good prots
    for batch_num, batch_set in enumerate(batched_prots):

        print(f'On batch {batch_num + 1}/{num_batches}\n', end='')

        extended_batch_set = [(k, entry) for (k, v) in batch_set for entry in v]

        # create our list representation
        batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = make_database_set(extended_batch_set, max_len)

        # find tha batch matched masses for both b and y ions
        matched_masses_b_batch = merge(batch_b_list, index_list_b, batch_kmer_b, spectra_boundaries)
        matched_masses_y_batch = merge(batch_y_list, index_list_y, batch_kmer_y, spectra_boundaries)

        # add these these hits to our function scoped set of masses and to the kmer set
        for k, v in matched_masses_b_batch.items():
            matched_masses_b[k] += v 
            
            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

        for k, v in matched_masses_y_batch.items():
            matched_masses_y[k] += v 

            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

    #update kmers in db to the kmer set
    db = db._replace(kmers=kmer_set)

    return (matched_masses_b, matched_masses_y, db)

def make_database_set(
    proteins: list, 
    max_len: int
    ) -> Tuple[arr.array, arr.array, list, arr.array, arr.array, list, dict]:
    '''Create parallel lists of (masses, index_maps, kmers) for the merge sort operation
    where index_maps map the massses to a range of positions in the kmers list 

    :param proteins: protein entries of shape (name, entry) where entry has a 
        *.sequence* attribute
    :type proteins: list
    :param max_len: max k-mer length
    :type max_len: int
    
    :returns: b ion masses created from the protein set, mapping from b ion masses 
        to the kmers associated (same length as b ion masses), kmers associated with 
        b ion masses, y ion masses created from the protein set, mapping from y ion masses 
        to the kmers associated (same length as y ion masses), kmers associated with 
        y ion masses, mapping of k-mers to source proteins
    :rtype: (array, array, list, array, array, list, dict)
    '''
    db_dict_b = defaultdict(set)
    db_dict_y = defaultdict(set)
    kmer_set = defaultdict(list)

    # function to add all masses of b+, b++, y+, y++ ions to the correct dict for sorting
    # and keep track of the source proteins
    def add_all(kmer, prot_name, start_location, end_location, protein_number):
        for ion in 'by':
            for charge in [1, 2]:

                #hack
                #spec = pre_spec.spectrum
                #the_type = type(spec) #python = 'dict' #cpp = 'list'
                pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
                spec = pre_spec
                if isinstance(pre_spec,dict):
                    spec = pre_spec.get('spectrum')

                for i, mz in enumerate(spec):
                    start_position = start_location if ion == 'b' else end_location
                    end_position = start_position + i if ion == 'b' else end_location - i
                    kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:] 
                    r_d = db_dict_b if ion == 'b' else db_dict_y
                    # r_d[mz].add(kmer_to_add)
                    if ion == 'b':
                        r_d[mz].add((mz, protein_number, kmer_to_add, str(start_position) + '-' + str(end_position), ion, charge))
                    else:
                        r_d[mz].add((mz, protein_number, kmer_to_add, str(end_position) + '-' + str(start_position), ion, charge))
                    kmer_set[kmer_to_add].append(prot_name)

    plen = len(proteins)

    # go through each protein and add all kmers to the correct dictionary for later sorting
    for i, (prot_name, prot_entry) in enumerate(proteins):
        
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        
        for j in range(1, max_len):
            kmer = prot_entry.sequence[:j]
            start_position = 1
            end_position = j
            add_all(kmer, prot_name, start_position, end_position, i)

        for j in range(len(prot_entry.sequence) - max_len):
            kmer = prot_entry.sequence[j:j+max_len]
            start_position = j + 1
            end_position = j + max_len
            add_all(kmer, prot_name, start_position, end_position, i)

        for j in range(len(prot_entry.sequence) - max_len, len(prot_entry.sequence)):
            kmer = prot_entry.sequence[j:]
            start_position = j+1
            end_position = len(prot_entry.sequence)
            add_all(kmer, prot_name, start_position, end_position, i)

        # for kmer_len in range(1,max_len):
        #     for j in range(len(prot_entry.sequence) - kmer_len + 1):
        #         kmer = prot_entry.sequence[j:j+kmer_len]
        #         add_all(kmer, prot_name, 1, 1, 1)

    print('\nSorting the set of protein masses...')
    
    # arrays take less space than lists
    db_list_b, index_list_b, kmer_list_b = arr.array('f'), arr.array('i'), []
    db_list_y, index_list_y, kmer_list_y = arr.array('f'), arr.array('i'), []

    # sort the keys to make sure that we are going through masses the correct way
    sorted_keys = sorted(db_dict_b.keys())
    for mz in sorted_keys:

        # get all kmers associated with this mass, append it the b ion list, and keep track 
        # of the kmers in the kmer list
        kmers = db_dict_b[mz]
        db_list_b.append(mz)
        offset = 0 if not len(index_list_b) else index_list_b[-1]
        index_list_b.append(len(kmers) + offset)
        kmer_list_b += kmers

    sorted_keys = sorted(db_dict_y.keys())
    for mz in sorted_keys:

        # get all kmers associated with this mass, append it the y ion list, and keep track 
        # of the kmers in the kmer list
        kmers = db_dict_y[mz]
        db_list_y.append(mz)
        offset = 0 if not len(index_list_y) else index_list_y[-1]
        index_list_y.append(len(kmers) + offset)
        kmer_list_y += kmers

    print('Done')

    return db_list_b, index_list_b, kmer_list_b, db_list_y, index_list_y, kmer_list_y, kmer_set

def merge(
    mz_s: Iterable, 
    indices: Iterable, 
    kmers: Iterable, 
    boundaries: Iterable
    ) -> defaultdict:
    '''Perform a linear search of observed mz values that fit into mappings 
    to get a mapping from mz boundaries [*lower_bound*, *upper_bound*] to a set 
    of k-mers

    :param mz_s: mz values to look through
    :type mz_s: Iterable
    :param indices: index mappings from mz values to kmers. The steps to get 
        from an *m/z* value to a set of k-mers is as follows: if we have a m/z
        value at index *i*, we will get the values in range of *indices*[*i*-1] 
        to *indices*[*i*], call it j in J. Then, the k-mers we want are all kmers
        at *kmers*[*j*] for each j in J.
    :type indices: Iterable
    :param kmers: the k-mers associated with the mz_s in the range of indices
        as described in the *indices* param
    :type kmers: Iterable
    :param boundaries: lower upper bounds for a mass in a list like 
        [*lower_bound*, *upper_bound*]
    :type boundaries: Iterable

    :returns: mapping from a string (using src.utils.hashable_boundaries) 
        of <lower_bound>-<upper_bound> to a list of k-mers
    :rtype: defaultdict
    '''

    b_i, mz_i = 0, 0

    matched_masses = defaultdict(list)

    while b_i < len(boundaries) and mz_i < len(mz_s):

        # if mz_s[mz_i] is in the boundary, keep track of it increment mz_i
        if boundaries[b_i][0] <= mz_s[mz_i] <= boundaries[b_i][1]:
            matched_masses[utils.hashable_boundaries(boundaries[b_i])] += kmers[indices[mz_i - 1]:indices[mz_i]]
            mz_i += 1

        # if the upper boundary is less than the mz_i, increment b_i
        elif mz_s[mz_i] > boundaries[b_i][1]:
            b_i += 1

        # if the lower boundary is greater than the mz_i, increment mz_i
        elif mz_s[mz_i] < boundaries[b_i][0]:
            mz_i += 1

    return matched_masses

def find_hits(mz_mapping, boundaries, spectrum, spec_num, matched_masses_b, matched_masses_y,):
    b_hits =[]
    y_hits = []
    b_hit_set = set()
    y_hit_set = set()
    miss_set = set()
    hit_list = []
    for k, mz in enumerate(spectrum.mz_values):
        mapped = mz_mapping[mz]
        b = boundaries[mapped]
        b = utils.hashable_boundaries(b)

        if b in matched_masses_b:
            mz_hit_tuple = (spec_num, k, mz)
            b_hit_set.add(mz_hit_tuple)
            for tuple in matched_masses_b[b]:
                extended_tuple = (spec_num, mz, tuple)
                b_hits.append(extended_tuple)
                hit_list.append(mz)
            # b_hits += matched_masses_b[b]

        if b in matched_masses_y:
            mz_hit_tuple = (spec_num, k, mz)
            y_hit_set.add(mz_hit_tuple)
            for tuple in matched_masses_y[b]:
                extended_tuple = (spec_num, mz, tuple)
                y_hits.append(extended_tuple)
                hit_list.append(mz)
            # y_hits += matched_masses_y[b]
        
        if mz not in hit_list:
            miss_tuple = (spec_num, k,mz)
            miss_set.add(miss_tuple)
   
    return b_hits, y_hits, b_hit_set, y_hit_set, miss_set

def find_misses(input_spectrum, hit_set, miss_set):
    for k, mz in enumerate(input_spectrum.spectrum):
        if (k,mz) not in hit_set:
            miss_tuple = (k,mz)
            miss_set.add(miss_tuple)

def map_mz(input_spectra, ppm_tolerance, matched_masses_b, matched_masses_y):
    #Map to (P_y, S_i, P_j, seq, b/y)
    # Where P_y is protein this was found in, S_i is m/z number, P_j is location within that protein, seq and b/y are straightforward 
    mz_mapping = defaultdict(set)
    for spectrum in input_spectra:
        find_matches_in_spectrum(spectrum, mz_mapping, ppm_tolerance, matched_masses_b, matched_masses_y)

def find_matches_in_spectrum(spectrum, mz_mapping, ppm_tolerance, matched_masses_b, matched_masses_y):
    for i, mz in enumerate(spectrum[0]):
        boundary_set = preprocessing_utils.make_boundaries(mz, ppm_tolerance)
        match_b(matched_masses_b, matched_masses_y, boundary_set, mz_mapping, i, mz)
        match_y(matched_masses_b, matched_masses_y, boundary_set, mz_mapping, i, mz)
    return mz_mapping

def match_b(matched_masses_b, matched_masses_y, boundary_set, mz_mapping, i, mz):
    for boundary in matched_masses_b.keys():
        if boundary == str(boundary_set[0]) + '-' + str(boundary_set[1]):
            add_match_to_dict('b', matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz)

def match_y(matched_masses_b, matched_masses_y, boundary_set, mz_mapping, i, mz):
    for boundary in matched_masses_y.keys():
        if boundary == str(boundary_set[0]) + '-' + str(boundary_set[1]):
            add_match_to_dict('y', matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz)

def add_match_to_dict(ion, matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz):
    if ion == 'b':
        for matching in matched_masses_b[boundary]:
            mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))
    else:
        for matching in matched_masses_y[boundary]:
            mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))

def check_if_max(val1, current_max, location):
    if (val1 > current_max):
        max = val1
        max_location = location
    
    return max, max_location
def collect_metadata(input_spectra, correct_sequences, ppm_tolerance, all_hits, all_misses, hit_abundances, miss_abundances, misleading_abundances):
    for i, spectrum in enumerate(input_spectra):
        max_abundance = 0
        found = False
        initial_hits = []
        input_spectrum = spectrum[0]
        tot_measured_spec_length = tot_measured_spec_length + len(input_spectrum)
        input_abundance_set = spectrum[1]
        precursor_mass = spectrum[5]
        precursor_charge = spectrum[6]
        ideal_spectrum = gen_spectra.gen_spectrum(correct_sequences[i])
        tot_ideal_spec_length = tot_ideal_spec_length + len(ideal_spectrum['spectrum'])
        # Checking input_spectrum for hits
        for a, j in enumerate(input_spectrum):
            #Finding max abundance
            max_abundance, max_abundance_location = check_if_max(input_abundance_set[a], max_abundance, a)

        # Checking precursor mass
        tolerance = utils.ppm_to_da(precursor_mass, ppm_tolerance)
        if (not isintolerance(precursor_mass, gen_spectra.get_precursor(correct_sequences[i], precursor_charge), tolerance)):
            count = count + 1

    # Collecting average abundance
    avg_hit_abundance = get_average_from_set(hit_abundances)
    avg_miss_abundance = get_average_from_set(miss_abundances)
    # Collecting total length
    total_length = get_total_length(correct_sequences)

def append_correct_hits(correct_hits, correct_sequence, input_spectrum, ppm_tolerance):
    correct_hits = []
    ideal_spectrum = gen_spectra.gen_spectrum(correct_sequence)
    for val in ideal_spectrum['spectrum']:
        boundary = preprocessing_utils.make_boundaries(float(val), ppm_tolerance)
        for mz in input_spectrum[0]:
            if mz > boundary[0] and mz < boundary[1]:
                correct_hits.append(mz)

def transform(line):
    A = line.rstrip().split()
    spectrum_num = int(A[0])
    mz = float(A[1])
    parent_prot = int(A[2])
    seq = A[3]
    interval = A[4]
    ion = A[5]
    charge = int(A[6])
    initial_score = 1

    start_end = interval.split('-')
    interval_start = int(start_end[0])
    interval_end = int(start_end[1])
    interval = [parent_prot, interval_start, interval_end, initial_score, seq, mz, ion]
    return interval


def map_to_interval(hit_file):
    interval_list = []
    with open(hit_file, 'r') as h:
        for line in h:
            interval = transform(line)
            interval_list.append(interval)
        return interval_list

def group_intervals(sorted_intervals):
    interval_dict = dict()
    for interval in sorted_intervals:
        if int(interval[0]) not in interval_dict.keys():
            interval_dict[int(interval[0])] = []
            interval_dict[int(interval[0])].append([interval[0], interval[1], interval[2], interval[3], interval[4], interval[5], interval[6]])
        else:
            interval_dict[int(interval[0])].append([interval[0], interval[1], interval[2], interval[3], interval[4], interval[5], interval[6]])
    return interval_dict

def isempty(list):
    if len(list) == 0:
        return True
    else:
        return False

def wrong_char(sequence: str):
    if (('B' in sequence) or ('J' in sequence) or ('O' in sequence) or ('U' in sequence) or ('Z' in sequence) or ('X' in sequence)):
        return True
    else:
        return False
def merge_intervals(interval_dict):
    updated_interval_list = []
    for parent_prot in interval_dict.keys():
        interval_list = interval_dict[parent_prot]
        interval_list.sort(key=lambda x: x[1])
        interval_list = merge_overlapping_interval(interval_list)
        [updated_interval_list.append(x) for x in interval_list]
    return updated_interval_list

def merge_interval(interval1, interval2, stack, mz_hit_list):
    interval1_len = interval1[2] - interval1[1]
    interval2_len = interval2[2] - interval2[1]
    if interval1_len >= interval2_len:
        interval1.append(interval2)
        stack[-1] = interval1
        mz_hit_list.append(interval1[5])
    else:
        interval2.append(interval1)
        stack[-1] = interval2
        mz_hit_list.append(interval2[5])
    return stack, mz_hit_list

def compare_intervals(interval, comparison_interval):
    if (interval != comparison_interval):
        #Case 1: interval encompasses comparison_interval
        if (interval[1] <= comparison_interval[1]) and (interval[2] >= comparison_interval[2]):
            interval[3] = interval[3] + comparison_interval[3]
            return True
        #Case 2: comparison_interval encompasses interval
        elif (interval[1] >= comparison_interval[1]) and (interval[2] <= comparison_interval[2]):
            comparison_interval[3] = comparison_interval[3] + interval[3]
            return True
        #Case 3: Neither interval encompasses the other
        else:
            return False
    return True

def merge_overlapping_interval(interval_list):
    # input: [[170, 175], [170, 180], [1, 2]]
    # output: [[0, 170, 180, 2], [0, 1, 2, 1]]
    stack = []
    stack.append(interval_list[0])
    mz_hit_list = []
    for interval in interval_list:
        if ((compare_intervals(stack[-1], interval)) and (interval[5] not in mz_hit_list)):
            stack, mz_hit_list = merge_interval(stack[-1], interval, stack, mz_hit_list)
        else:
            stack.append(interval)

    return stack

def score_hits(b_hits, y_hits, spectrum):
    b_results = sorted([
        (
            hit, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(hit[3], ion='b'))
        ) for hit in b_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    y_results = sorted([
        (
            hit, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(hit[3], ion='y'))
        ) for hit in y_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )

    return b_results, y_results

def filter_hits(b_hits, y_hits, precursor_mass, precursor_charge):

    filtered_b, filtered_y = [], []
    count = 0
    for hit in b_hits:
        if not (gen_spectra.get_precursor(hit[3]) * hit[6]) > (precursor_mass * precursor_charge):
            filtered_b.append(hit)
        else:
            if hit[3] == "MSSP":
                print(hit[3], hit[6], hit[3] * hit[6], gen_spectra.get_precursor(hit[3]), precursor_mass, precursor_charge, precursor_mass * precursor_charge, )
            count = count + 1
    for hit in y_hits:
        if not (gen_spectra.get_precursor(hit[3]) * hit[6]) > (precursor_mass * precursor_charge):
            filtered_y.append(hit)
        else:
            count = count + 1

    print(str(count), "kmers were filtered out")
    return filtered_b, filtered_y

def parse_hits(b_hits, y_hits):
    for i, x in enumerate(b_hits):
        pep_id = x[0]
        w = x[1][0]
        prot_id = x[1][1]
        seq = x[1][2]
        loc = x[1][3]
        write_ion = x[1][4]
        charge = x[1][5]
        x = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        b_hits[i] = x
    for i, y in enumerate(y_hits):
        pep_id = y[0]
        w = y[1][0]
        prot_id = y[1][1]
        seq = y[1][2]
        loc = y[1][3]
        write_ion = y[1][4]
        charge = y[1][5]
        y = [pep_id, w, prot_id, seq, loc, write_ion, charge]
        y_hits[i] = y
    return b_hits, y_hits

def isSubSequence(string1, string2, m, n):
    # Base Cases
    if m == 0:
        return True
    if n == 0:
        return False
 
    # If last characters of two
    # strings are matching
    if string1[m-1] == string2[n-1]:
        return isSubSequence(string1, string2, m-1, n-1)
 
    # If last characters are not matching
    return isSubSequence(string1, string2, m, n-1)

def get_top_X(b_results, y_results, TOP_X):
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

    return filtered_b, filtered_y

def is_good_hit(kmer: string, ion, correct_sequence):
    #take in a kmer and ion and determine if sequence is good subsequence. Also return score based on length of kmer
    if ion == 'b':
        if (kmer[0] == correct_sequence[0]):
            if (isSubSequence(kmer, correct_sequence, len(kmer), len(correct_sequence))):
                return (True, len(kmer))
            else:
                return (False, 0)
        else:
            return (False, 0)
    else:
        if (kmer[-1] == correct_sequence[-1]):
            if (isSubSequence(kmer, correct_sequence, len(kmer), len(correct_sequence))):
                return (True, len(kmer))
            else:
                return (False, 0)
        else:
            return (False, 0)

def get_good_kmers(b_results, y_results, correct_sequence, TOP_X):
    good_kmers = []
    bad_kmers = []
    for i, b_entry in enumerate(b_results):
        seq = b_entry[3]
        if (seq[0] == correct_sequence[0]):
            if (isSubSequence(seq, correct_sequence, len(seq), len(correct_sequence))):
                b_entry.append(TOP_X - i)
                good_kmers.append(b_entry)
            else:
                b_entry.append(TOP_X - i)
                bad_kmers.append(b_entry)
        else:
            b_entry.append(TOP_X - i)
            bad_kmers.append(b_entry)
    for i, y_entry in enumerate(y_results):
        seq = y_entry[3]
        if (seq[-1] == correct_sequence[-1]):
            if (isSubSequence(seq, correct_sequence, len(seq), len(correct_sequence))):
                y_entry.append(TOP_X - i)
                good_kmers.append(y_entry)
            else:
                y_entry.append(TOP_X - i)
                bad_kmers.append(y_entry)
        else:
            y_entry.append(TOP_X - i)
            bad_kmers.append(y_entry)
    
    return good_kmers, bad_kmers

def make_boundaries(mz, ppm_tol):
    da_tol = ppm_to_da(mz, ppm_tol)
    return [mz - da_tol, mz + da_tol]

def define_single_spectrum(mz_list, ppm_tol):
    boundaries = [make_boundaries(mz, ppm_tol) for mz in mz_list]

    # make overlapped boundaries larger boundaries
    boundaries = overlap_intervals(boundaries)

    # make a mapping for mz -> boundaries
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(mz_list):
        
        # if the current boundary encapsulates s_i, add to list
        if boundaries[b_i][0] <= mz_list[s_i] <= boundaries[b_i][1]:
            mz_mapping[mz_list[s_i]] = b_i 
            s_i += 1

        # else if the s_i < boundary, increment s_i
        elif mz_list[s_i] < boundaries[b_i][0]:
            s_i += 1

        # else if s_i > boundary, incrment b_i
        elif mz_list[s_i] > boundaries[b_i][1]:
            b_i += 1
    return boundaries, mz_mapping

def write_cluster(cluster):
    if len(cluster) == 0 : return None
    O = []

    O.append(len(cluster))
    O.append(cluster[0].pid)

    max_len = 0
    max_hit = None

    for hit in cluster:
        l = hit.end - hit.start + 1
        if l > max_len:
            max_len = l
            max_hit = hit

    O.append(max_hit.seq)
    O.append(max_hit.mz)
    
    for hit in cluster:
        O.append( (hit.start, hit.end, hit.seq, hit.mz) )    # b_hits

    with open('clusters.txt', 'a') as c:
        c.write( '\t'.join( [str(o) for o in O] ) )
        c.write('\n')

    
def print_clusters(b_sorted_clusters, y_sorted_clusters):
    for cluster in b_sorted_clusters:
        non_indices = str(cluster.score) + '\t' + str(cluster.pid) + '\t' + cluster.seq
        print(non_indices + '\t'+ '\t'.join([str(o) for o in cluster.indices]))
    for cluster in y_sorted_clusters:
        non_indices = str(cluster.score) + '\t' + str(cluster.pid) + '\t' + cluster.seq
        print(non_indices + '\t'+ '\t'.join([str(o) for o in cluster.indices]))

def get_hits_from_cluster(b_sorted_clusters, y_sorted_clusters):
    b_hit_arr = []
    y_hit_arr = []
    for cluster in b_sorted_clusters:
        b_hit_arr.append(cluster.seq)
    for cluster in y_sorted_clusters:
        y_hit_arr.append(cluster.seq)
    
    return b_hit_arr, y_hit_arr

def write_hits(b_hits, y_hits):
    with open("b_hits.txt", 'w') as b:
        for x in b_hits:
            pep_id = x[0]
            w = x[1]
            prot_id = x[2][1]
            seq = x[2][2]
            loc = x[2][3]
            ion = x[2][4]
            charge = x[2][5]
            out = [pep_id, w, prot_id, seq, loc, ion, charge]
            b.write('\t'.join([str(i) for i in out]) + '\n')
    with open("y_hits.txt", 'w') as y_file:
        for y in y_hits:
            pep_id = y[0]
            w = y[1]
            prot_id = y[2][1]
            seq = y[2][2]
            loc = y[2][3]
            ion = y[2][4]
            charge = y[2][5]
            out = [pep_id, w, prot_id, seq, loc, ion, charge]
            y_file.write('\t'.join([str(i) for i in out]) + '\n')

def parse_hits(Hit, file_name):
    hits = []
    for l in open(file_name):
            A = l.rstrip().split('\t')
            pid = int(A[2])
            start = int(A[4].split('-')[0])
            end = int(A[4].split('-')[1])
            seq = A[3]
            mz = A[1]

            hits.append( Hit(pid=pid, start=start, end=end, seq=seq, mz=mz) )
    return hits
def create_clusters(ion):
    Hit = collections.namedtuple('Hit', 'pid start end seq mz')

    if ion == 'b':
        file_name = 'b_hits.txt'
        # file_name = sys.argv[1]

        hits = parse_hits(Hit, file_name)

        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'start', 'end'))

        last_pid = None
        last_start = None

        cluster = []

        with open('clusters.txt', 'w') as c:
            c.write('')

        if ion == 'b':
            for hit in sorted_hits:
                if last_pid == hit.pid and last_start == hit.start:
                    cluster.append(hit)
                else:
                    write_cluster(cluster)
                    cluster = [hit]
                last_pid = hit.pid
                last_start = hit.start
    
    else:
        file_name = 'y_hits.txt'
        # file_name = sys.argv[1]

        hits = parse_hits(Hit, file_name)

        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'start', 'end'))

        last_pid = None
        last_start = None

        cluster = []

        with open('clusters.txt', 'w') as c:
            c.write('')


        if ion == 'y':
            for hit in sorted_hits:
                if last_pid == hit.pid and last_end == hit.end:
                    cluster.append(hit)
                else:
                    write_cluster(cluster)
                    cluster = [hit]
                last_pid = hit.pid
                last_end = hit.end

def set_prior(mz, ion, mz_mapping, boundaries, matched_masses_b, matched_masses_y):
    # for i in range (0, len(b_sorted_clusters)):
    # prior = 1 / # of occurances
    mapped = mz_mapping[mz]
    b = boundaries[mapped]
    b = utils.hashable_boundaries(b)
    if ion == 'b':
        if len(matched_masses_b[b]) == 0:
            print(mz)
            print(mapped)
            print(b)
            with open('mz_mapping.txt', 'w') as m:
                [m.write(str(x) + '\n') for x in mz_mapping]
        P_A = 1/len(matched_masses_b[b]) # if (len(matched_masses_b[b]) !=0) else 1
    else:
        if len(matched_masses_y[b]) == 0:
            print(mz)
            print(mapped)
            print(b)
            with open('mz_mapping.txt', 'w') as m:
                [m.write(str(x) + '\n') for x in mz_mapping]
        P_A = 1/len(matched_masses_y[b]) # if (len(matched_masses_b[b]) !=0) else 1
    
    return P_A

def calc_post_prob(prior, indices):
    post = 0
    current_prob = prior
    for element in indices:
        A = element.rstrip().split(',')
        string = A[2]
#         mz = (A[3])
#         mz = mz[:-1]
#         mz = mz[2:]
        post = post + 1/len(string) + prior
    return post

def sort_clusters_by_post_prob(ion, mz_mapping, boundaries, matched_masses_b, matched_masses_y):
    cluster = collections.namedtuple('cluster', 'post_prob prior score pid seq mz indices')
    if ion == 'b':
        b_cluster_array = []
        with open('clusters.txt', 'r') as c:
            for line in c:
                A = line.rstrip().split('\t')
                score = int(A[0])
                pid = int(A[1])
                seq = A[2]
                mz = float(A[3])
                indices = []
                [indices.append(A[x]) for x in range(4,len(A))]
                prior = set_prior(mz, ion, mz_mapping, boundaries, matched_masses_b, matched_masses_y)
                post_prob = calc_post_prob(prior, indices)

                b_cluster_array.append(cluster(post_prob=post_prob, prior=prior, score=score, pid=pid, seq=seq, mz=mz, indices=indices) )

        b_sorted_clusters = sorted(b_cluster_array, key=operator.attrgetter('post_prob', 'score', 'pid', 'prior'), reverse = True)
        return b_sorted_clusters
    else:
        # y_hits
        y_cluster_array = []
        with open('clusters.txt', 'r') as c:
            for line in c:
                A = line.rstrip().split('\t')
                score = int(A[0])
                pid = int(A[1])
                seq = A[2]
                mz = float(A[3])
                indices = []
                [indices.append(A[x]) for x in range(4,len(A))]
                prior = set_prior(mz, 'y', mz_mapping, boundaries, matched_masses_b, matched_masses_y)
                post_prob = calc_post_prob(prior, indices)

                y_cluster_array.append(cluster(post_prob=post_prob, prior=prior, score=score, pid=pid, seq=seq, mz=mz, indices=indices) )

        y_sorted_clusters = sorted(y_cluster_array, key=operator.attrgetter('post_prob', 'score', 'pid', 'prior'), reverse = True)
        return y_sorted_clusters


def write_b_sorted_cluster(b_sorted_clusters):
    # for cluster in b_sorted_clusters:
    #     non_indices = str(cluster.score) + '\t' + str(cluster.post_prob) + '\t' + str(cluster.prior) + '\t' + str(cluster.pid) + '\t' + cluster.seq + '\t' + str(cluster.mz)
    #     print(non_indices) #+ '\t'+ '\t'.join([str(o) for o in cluster.indices]))
    with open("b_sorted_clusters.txt", 'w') as b:
        [b.write(str(x) + '\n') for x in b_sorted_clusters]
def write_y_sorted_cluster(y_sorted_clusters):
    # for cluster in y_sorted_clusters:
        # non_indices = str(cluster.score) + '\t' + str(cluster.post_prob) + '\t' + str(cluster.prior) + '\t' + str(cluster.pid) + '\t' + cluster.seq + '\t' + str(cluster.mz)
        # print(non_indices) #+ '\t'+ '\t'.join([str(o) for o in cluster.indices]))
        with open("y_sorted_clusters.txt", 'w') as y:
            [y.write(str(x) + '\n') for x in y_sorted_clusters]

def test_optimized_compare_masses(
    mz_mapping,
    bbby,
    observed: list, 
    reference: list, 
    ppm_tolerance: int = 20, 
    needs_sorted: bool = False
    ) -> float:
    '''Score two spectra against eachother. Simple additive scoring of ions found

    :param observed: observed set of m/z values
    :type observed: list
    :param reference: reference set of m/z values
    :type reference: list
    :param ppm_tolerance: parts per million mass error allowed when matching masses. 
        (default is 20)
    :type ppm_tolerance: int
    :param needs_sorted: Set to true if either the observed or reference need to 
        be sorted. 
        (default is False)
    :type needs_sorted: bool

    :returns: the number of matched ions
    :rtype: int

    :Example:

    >>> optimized_compare_masses([1, 2, 4], [1, 3, 4], 1, False)
    >>> 2
    '''
    if len(observed) == 0 or len(reference) == 0:
        return 0.0

    if needs_sorted:
        observed.sort()
        reference.sort()
        
    def boundaries(mass):
        tol = ppm_to_da(mass, ppm_tolerance)
        return [mass - tol, mass + tol]
                
    # calculate the boundaries for each of the reference masses for binary search
    observed_boundaries = []
    for obs in observed:
        observed_boundaries += boundaries(obs)

    obs_boundaries = []
    for mz in observed:
        mapped = mz_mapping[mz]
        b = bbby[mapped]
        obs_boundaries.append(b[0])
        obs_boundaries.append(b[1])
    
    if obs_boundaries != observed_boundaries:
        print('WRONG!')
        
    #hack
    #the_type = type(reference) #python = 'dict' #cpp = 'list'
    updated_reference = reference      
    if isinstance(updated_reference,dict):
        reference = updated_reference.get('spectrum')
    
    # local variables for score
    return_value = 0
    for ref in reference:
        if bisect(observed_boundaries, ref) % 2:
            return_value = return_value + 1
    # return_value = sum([1 for ref in reference if bisect(observed_boundaries, ref) % 2])
    return return_value