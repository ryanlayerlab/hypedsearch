import computational_pipeline.identification
import multiprocessing as mp
import lookups.utils
from preprocessing import preprocessing_utils

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

def run(args: dict) -> dict:
    number_of_cores = get_number_of_cores(args['number_cores'])
    built_database = get_built_database(args['database_file_path'])
    spectras = get_spectras(args['spectra_file_paths'],args['number_peaks'],args['relative_abundance'])
    lookups.utils.make_dir(args['output_folder_path'])

    matched_spectras = computational_pipeline.identification.get_matched_spectras(
        spectras = spectras, 
        built_database = built_database, 
        max_peptide_length=args['max_peptide_length'], 
        ppm_tolerance=args['ppm_tolerance'], 
        precursor_tolerance=args['precursor_tolerance'],
        number_peaks=args['number_peaks'],
        relative_abundance=args['relative_abundance'],
        digest_left=args['digest_left'], 
        digest_right=args['digest_right'],
        number_hybrids=args['number_hybrids'], 
        number_natives=args['number_natives'], 
        number_of_cores=number_of_cores, 
        create_kmer_database = args["create_kmer_database"],
        verbose=args['verbose'], 
        output_folder_path=args['output_folder_path'])
    return matched_spectras
