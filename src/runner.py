import computational_pipeline.identification
import multiprocessing as mp

def run(args: dict) -> dict:
    number_of_cores = max(1, args['number_cores'])
    number_of_cores = min(number_of_cores, mp.cpu_count() - 1)
    matched_spectra = computational_pipeline.identification.get_matched_spectrums(
        spectra_file_paths = args['spectra_file_paths'], 
        built_database = args['built_database'], 
        max_peptide_length=args['max_peptide_length'], 
        ppm_tolerance=args['ppm_tolerance'], 
        precursor_tolerance=args['precursor_tolerance'],
        number_peaks=args['number_peaks'],
        relative_abundance_filter=args['relative_abundance'],
        digest_left=args['digest_left'], 
        digest_right=args['digest_right'],
        number_hybrids=args['number_hybrids'], 
        number_natives=args['number_natives'], 
        number_of_cores=number_of_cores, 
        create_kmer_database = args["create_kmer_database"],
        verbose=args['verbose'], debug=args['debug'],
        output_folder_path=args['output_folder_path'])
    return matched_spectra
