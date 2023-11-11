import computational_pipeline.identification
import multiprocessing as mp

def run(args: dict) -> dict:
    cores = max(1, args['number_cores'])
    cores = min(cores, mp.cpu_count() - 1)
    matched_spectra = computational_pipeline.identification.id_spectra(args['spectra_file_paths'], args['built_database'], 
        min_peptide_len=args['min_peptide_length'], max_peptide_len=args['max_peptide_length'], 
        ppm_tolerance=args['ppm_tolerance'], precursor_tolerance=args['precursor_tolerance'],
        peak_filter=args['number_peaks'],relative_abundance_filter=args['relative_abundance'],
        digest_left=args['digest_left'], digest_right=args['digest_right'], num_hybrids=args['number_hybrids'], 
        num_natives=args['number_natives'],verbose=args['verbose'], DEBUG=args['debug'], 
        cores=cores, make_new = args["create_kmer_database"], output_dir=args['output_folder_path'])
    return matched_spectra
