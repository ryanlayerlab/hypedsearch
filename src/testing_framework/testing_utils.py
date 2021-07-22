from typing import Tuple
from pyteomics import fasta
import os
import pandas as pd
from collections import defaultdict

from collections import namedtuple
import sys
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src', 'hypedsearch'))
if module_path not in sys.path:
    sys.path.append(module_path)
from preprocessing import preprocessing_utils
from objects import Database
import gen_spectra

from utils import hashable_boundaries, predicted_len

from file_io import spectra
from collections import defaultdict
from typing import Iterable
from utils import overlap_intervals, ppm_to_da
import array as arr
from math import ceil
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


    raw_prefix = os.path.join(root, 'mnt', 'c', 'Users', 'Maxim', 'Documents', "Layer_Lab", 'Database', 'raw_inputs')
    NOD2_top_dir = 'NOD2_E3'
    BALB3_top_dir = 'BALB3_E3'


    NOD2_data = Dataset(
        # os.path.join(root, 'mnt', 'c', 'Users', 'Maxim', 'Documents', 'Layer_Lab', 'Database', 'Hybrid_inputs') + os.path.sep,
        os.path.join(raw_prefix, NOD2_top_dir, 'mzml') + os.path.sep,
        os.path.join(raw_prefix, NOD2_top_dir, 'NOD2_E3_results.ssv'),
        os.path.join(root, 'home', 'ncol107453', 'jaime_hypedsearch', 'hypedsearch', 'data', 'database', 'sample_database.fasta'),
        os.path.join(raw_prefix, NOD2_top_dir) + os.path.sep,
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

# def get_proteins(db: Database):
#     return db.proteins
# def get_kmers(db: Database, proteins: list, ):

#     for i, (prot_name, prot_entry) in enumerate(proteins):
        
#         print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        
#         #Create all kmers
#         for kmer_len in range(1,max_len):
#             for j in range(len(prot_entry.sequence) - kmer_len + 1):
#                 kmer = prot_entry.sequence[j:j+kmer_len]
#                 add_all_masses(kmer, prot_name, )

#     b = db_dict_b.keys().sort()
#     y = db_dict_y.keys().sort()

# def add_all_masses(kmer, prot_name, start_location):
#     for ion in 'by':
#         for charge in [1, 2]:
#             pre_spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)
#             spec = pre_spec
#             if isinstance(pre_spec,dict):
#                 spec = pre_spec.get('spectrum')

#             for i, mz in enumerate(spec):
#                 kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
#                 r_d = db_dict_b if ion == 'b' else db_dict_y
#                 r_d[mz].add(kmer_to_add)
#                 kmer_set[kmer_to_add].append((prot_name, start_location))

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
                        r_d[mz].add((protein_number, kmer_to_add, str(start_position) + '-' + str(end_position), ion, charge))
                    else:
                        r_d[mz].add((protein_number, kmer_to_add, str(end_position) + '-' + str(start_position), ion, charge))
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
            matched_masses[hashable_boundaries(boundaries[b_i])] += kmers[indices[mz_i - 1]:indices[mz_i]]
            mz_i += 1

        # if the upper boundary is less than the mz_i, increment b_i
        elif mz_s[mz_i] > boundaries[b_i][1]:
            b_i += 1

        # if the lower boundary is greater than the mz_i, increment mz_i
        elif mz_s[mz_i] < boundaries[b_i][0]:
            mz_i += 1

    return matched_masses

def map_mz(input_spectra, ppm_tolerance, matched_masses_b, matched_masses_y):
    #Map to (P_y, S_i, P_j, seq, b/y)
    # Where P_y is protein this was found in, S_i is m/z number, P_j is location within that protein, seq and b/y are straightforward 
    mz_mapping = defaultdict(set)
    for spectrum in input_spectra:
        find_matches_in_spectrum(spectrum, mz_mapping, ppm_tolerance, matched_masses_b, matched_masses_y)

def find_matches_in_spectrum(spectrum, mz_mapping, ppm_tolerance, matched_masses_b, matched_masses_y):
    for i, mz in enumerate(spectrum[0]):
        boundaries = preprocessing_utils.make_boundaries(mz, ppm_tolerance)
        match_b(matched_masses_b, matched_masses_y, boundaries, mz_mapping, i, mz)
        match_y(matched_masses_b, matched_masses_y, boundaries, mz_mapping, i, mz)
    return mz_mapping

def match_b(matched_masses_b, matched_masses_y, boundaries, mz_mapping, i, mz):
    for boundary in matched_masses_b.keys():
        if boundary == str(boundaries[0]) + '-' + str(boundaries[1]):
            add_match_to_dict('b', matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz)

def match_y(matched_masses_b, matched_masses_y, boundaries, mz_mapping, i, mz):
    for boundary in matched_masses_y.keys():
        if boundary == str(boundaries[0]) + '-' + str(boundaries[1]):
            add_match_to_dict('y', matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz)

def add_match_to_dict(ion, matched_masses_b, matched_masses_y, boundary, mz_mapping, i, mz):
    if ion == 'b':
        for matching in matched_masses_b[boundary]:
            mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))
    else:
        for matching in matched_masses_y[boundary]:
            mz_mapping[mz].add((mz, matching[0], i, matching[2], matching[1], matching[3], matching[4]))