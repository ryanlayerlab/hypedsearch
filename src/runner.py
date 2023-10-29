import computational_pipeline.identification
import multiprocessing as mp

def run(args: dict) -> dict:
    cores = max(1, args['cores'])
    cores = min(cores, mp.cpu_count() - 1)
    matched_spectra = computational_pipeline.identification.id_spectra(args['spectra_files'], args['database_file'], min_peptide_len=args['min_peptide_len'], 
        max_peptide_len=args['max_peptide_len'], ppm_tolerance=args['tolerance'], precursor_tolerance=args['precursor_tolerance'],
        peak_filter=args['peak_filter'],relative_abundance_filter=args['relative_abundance_filter'],
        digest_left=args['digest_left'], digest_right=args['digest_right'], num_hybrids=args['num_hybrids'], num_natives=args['num_natives'],verbose=True, DEBUG=args['DEBUG'], cores=cores, make_new = args["new_db"], output_dir=args['output_dir'])
    return matched_spectra
