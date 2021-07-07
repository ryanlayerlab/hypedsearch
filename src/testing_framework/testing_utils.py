from pyteomics import fasta
import os
import pandas as pd
from posixpath import split

from collections import namedtuple
import sys
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src', 'hypedsearch'))
if module_path not in sys.path:
    sys.path.append(module_path)
from preprocessing import preprocessing_utils

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


    raw_prefix = os.path.join(root, )
    NOD2_top_dir = 'NOD2_E3'
    BALB3_top_dir = 'BALB3_E3'


    NOD2_data = Dataset(
        # os.path.join(root, 'mnt', 'c', 'Users', 'Maxim', 'Documents', 'Layer_Lab', 'Database', 'Hybrid_inputs') + os.path.sep,
        os.path.join(raw_prefix, NOD2_top_dir, 'mzml') + os.path.sep,
        os.path.join(raw_prefix, NOD2_top_dir, 'NOD2_E3_results.ssv'),
        os.path.join(raw_prefix, 'mouse_database.fasta'),
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

def preprocess_input_spectra(spectra_folder, ppm_tolerance, peak_filter: int = 0, relative_abundance_filter: float = 0.0):
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
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    print('Done')

    return spectra

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
