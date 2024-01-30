# This is an exact clone of identification.py with functions renamed for clarity and all code relating to creating an 
# alignment removed

from typing import Tuple
import sys
import os
path_to_src = (os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(path_to_src)
from src.lookups.objects import Database, Spectrum, MPSpectrumID, DEVFallOffEntry
from src.preprocessing import merge_search, preprocessing_utils
from src import database
from src.file_io import JSON

import time
import os
import copy
import json

# top results to keep for creating an alignment
TOP_X = 50

def database_and_spectra_preprocessing(
    spectra_files: str, 
    database_file: str, 
    write_path: str,
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

    # build/load the database
    verbose and print('Loading database...')
    db = database.build(database_file)
    verbose and print('Done')
    
    # load all of the spectra
    verbose and print('Loading spectra...')
    spectra, boundaries = preprocessing_utils.load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Done')

    # get the boundary -> kmer mappings for b and y ions

    matched_masses_b, matched_masses_y, db = merge_search.modified_match_masses(boundaries, db, max_peptide_len, DEBUG, write_path)
    return db