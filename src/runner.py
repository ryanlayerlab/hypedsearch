import computational_pipeline.identification
import multiprocessing as mp
import lookups.utils
from preprocessing import preprocessing_utils
from computational_pipeline.sqlite import database_file
from preprocessing import merge_search
from datetime import datetime
import os
from postprocessing.summary import write_aligned_spectras_to_disk

def get_built_database(database_file_path):
    return computational_pipeline.database.build(database_file_path)

def get_number_of_cores(number_of_cores):
    revised_number_of_cores = number_of_cores = max(1, number_of_cores)
    min_number_of_cores = min(revised_number_of_cores, mp.cpu_count() - 1)
    return min_number_of_cores

def get_spectras(spectra_file_paths,number_peaks,relative_abundance):
    spectras = []
    for spectra_file_path in spectra_file_paths:
        spectra = preprocessing_utils.load_spectra(spectra_file_path, number_peaks, relative_abundance)
        spectras.append(spectra)
    return spectras

def do_create_kmer_database(built_database, max_peptide_length, digest_left, digest_right):
    dbf = database_file(max_peptide_length, True)
    kv_prots = [(k, v) for k, v in built_database.proteins]    
    merge_search.modified_make_database_set(kv_prots, max_peptide_length, dbf, (digest_left, digest_right))

def get_output_file_name(spectra_file_paths):
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y%m%d%H%M%S')
    first_spectra_file_path = spectra_file_paths[0]
    last_token_with_extension = os.path.basename(first_spectra_file_path)
    filename, extension = os.path.splitext(last_token_with_extension)
    return_value = filename + "_" + formatted_datetime
    return return_value

def run(args: dict) -> dict:
    built_database = get_built_database(args['database_file_path'])
    spectras = get_spectras(args['spectra_file_paths'],args['number_peaks'],args['relative_abundance'])
    lookups.utils.make_dir(args['output_folder_path'])
    output_file_name = get_output_file_name(args['spectra_file_paths'])
    if args['create_kmer_database']:
        do_create_kmer_database(built_database, args['number_peaks'], args['digest_left'], args['digest_right'])

    max_peptide_length=args['max_peptide_length']
    ppm_tolerance=args['ppm_tolerance']
    precursor_tolerance=args['precursor_tolerance']
    number_hybrids=args['number_hybrids']
    number_natives=args['number_natives']
    output_folder_path=args['output_folder_path']
    aligned_spectras = computational_pipeline.identification.get_aligned_spectras(
        spectras,built_database,max_peptide_length,ppm_tolerance,precursor_tolerance,
        number_hybrids,number_natives)    
    write_aligned_spectras_to_disk(aligned_spectras, output_folder_path, output_file_name)
    
