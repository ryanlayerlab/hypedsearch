import os
import lookups.objects
import lookups.utils
import computational_pipeline
import computational_pipeline.identification as cp_id
import multiprocessing as mp
from preprocessing import preprocessing_utils
from preprocessing.sqlite_database import Sqllite_Database
from preprocessing import kmer_database
from datetime import datetime
from postprocessing.summary import write_aligned_spectrums_to_disk
from collections import namedtuple

def get_built_database(database_file_path):
    return computational_pipeline.database_generator.build_database(database_file_path)

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

def do_create_kmer_database(built_database, max_peptide_length, digest_left, digest_right):
    dbf = Sqllite_Database(max_peptide_length, True)
    kv_prots = [(k, v) for k, v in built_database.proteins]    
    kmer_database.modified_make_database_set(kv_prots, max_peptide_length, dbf, (digest_left, digest_right))
    
def get_output_file_name(spectra_file_paths):
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y%m%d%H%M%S')
    first_spectra_file_path = spectra_file_paths[0]
    last_token_with_extension = os.path.basename(first_spectra_file_path)
    filename, extension = os.path.splitext(last_token_with_extension)
    return_value = filename + "_" + formatted_datetime
    return return_value

def run(args: dict) -> dict:
    spectrums = get_spectrums(args['spectra_file_paths'],args['number_peaks'],args['relative_abundance'])
    built_database = get_built_database(args['database_file_path'])
    lookups.utils.make_dir(args['output_folder_path'])
    if args['create_kmer_database']:
        do_create_kmer_database(built_database, args['number_peaks'], args['digest_left'], args['digest_right'])
    max_peptide_length=args['max_peptide_length']
    ppm_tolerance=args['ppm_tolerance']
    precursor_tolerance=args['precursor_tolerance']
    number_hybrids=args['number_hybrids']
    number_natives=args['number_natives']
    target_seq = args['target_seq']
    output_file_name = get_output_file_name(args['spectra_file_paths']) 
    output_folder_path=args['output_folder_path']
    aligned_spectrums = cp_id.get_aligned_spectrums(spectrums,built_database,max_peptide_length,ppm_tolerance,precursor_tolerance,number_hybrids,number_natives,target_seq)  
    write_aligned_spectrums_to_disk(aligned_spectrums, output_folder_path, output_file_name)
    
