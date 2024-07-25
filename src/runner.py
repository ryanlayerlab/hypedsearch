import os
import sys
from pyteomics import fasta
from lookups.objects import Database
import lookups.objects
import lookups.utils
import computational_pipeline
import computational_pipeline.identification as cp_id
import multiprocessing as mp
from preprocessing import preprocessing_utils
from preprocessing.sqlite_database import Sqllite_Database
from datetime import datetime
from postprocessing.summary import write_aligned_spectrums_to_disk
from collections import namedtuple

def get_built_database(database_file_path):
    return 

def get_number_of_cores(number_of_cores):
    revised_number_of_cores = number_of_cores = max(1, number_of_cores)
    min_number_of_cores = min(revised_number_of_cores, mp.cpu_count() - 1)
    return min_number_of_cores

def get_spectrums(spectra_file_paths,number_peaks,relative_abundance):
    spectrums = []
    for spectra_file_path in spectra_file_paths:
        spectra = preprocessing_utils.load_spectra(spectra_file_path, number_peaks, relative_abundance)
        spectrums.extend(spectra)
    return spectrums

def create_fasta_database(database_file_path):
    fasta_database = Database(database_file_path,{},{})
    prots = []
    for entry in fasta.read(database_file_path):
        prots.append(entry)
    fasta_database = fasta_database._replace(proteins=prots)
    return fasta_database  

def create_sqllite_database(database_file_path,max_peptide_length,digest_left,digest_right):
    fasta_database = create_fasta_database(database_file_path)
    sqllite_database = Sqllite_Database(max_peptide_length, True)
    kv_proteins = [(k, v) for k, v in fasta_database.proteins]  
    sqllite_database.populate_database(kv_proteins, max_peptide_length, digest_left, digest_right)
    return sqllite_database

def get_existing_sqllite_database(max_peptide_length):
    sqllite_database = Sqllite_Database(max_peptide_length, False)
    return sqllite_database

def get_output_file_name(spectra_file_paths):
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y%m%d%H%M%S')
    first_spectra_file_path = spectra_file_paths[0]
    last_token_with_extension = os.path.basename(first_spectra_file_path)
    filename, extension = os.path.splitext(last_token_with_extension)
    return_value = filename + "_" + formatted_datetime
    return return_value
 
def get_aligned_spectrums_params(args: dict):
    spectrums = get_spectrums(args['spectra_file_paths'],args['number_peaks'],args['relative_abundance'])
    sqllite_database = None
    if args['create_sqllite_database']:
        sqllite_database = create_sqllite_database(args['database_file_path'],args['max_peptide_length'], args['digest_left'], args['digest_right'])
    else:
        sqllite_database = get_existing_sqllite_database(args['max_peptide_length'])
    lookups.utils.make_dir(args['output_folder_path'])
    max_peptide_length=args['max_peptide_length']
    ppm_tolerance=args['ppm_tolerance']
    precursor_tolerance=args['precursor_tolerance']
    number_hybrids=args['number_hybrids']
    number_natives=args['number_natives']
    target_seq = args['target_seq']
    aligned_spectrums_params = lookups.objects.AlignedSpectrumParams(spectrums=spectrums,sqllite_database=sqllite_database, 
        max_peptide_length=max_peptide_length, ppm_tolerance=ppm_tolerance,precursor_tolerance=precursor_tolerance,
        number_hybrids=number_hybrids,number_natives=number_natives,target_seq=target_seq)
    return aligned_spectrums_params

def run(args: dict):
    aligned_spectrums_params = get_aligned_spectrums_params(args)
    output_file_name = get_output_file_name(args['spectra_file_paths']) 
    output_folder_path=args['output_folder_path']
    aligned_spectrums = cp_id.get_aligned_spectrums(aligned_spectrums_params)  
    write_aligned_spectrums_to_disk(aligned_spectrums, output_folder_path, output_file_name)
