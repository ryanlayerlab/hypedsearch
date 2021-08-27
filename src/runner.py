import os
import identification
from postprocessing import summary, review
import multiprocessing as mp


def run(args: dict) -> dict:
    cores = max(1, args['cores'])
    cores = min(cores, mp.cpu_count() - 1)
    matched_spectra = identification.id_spectra(args['spectra_files'], args['database_file'], min_peptide_len=args['min_peptide_len'], 
        max_peptide_len=args['max_peptide_len'], ppm_tolerance=args['tolerance'], precursor_tolerance=args['precursor_tolerance'],
        peak_filter=args['peak_filter'],relative_abundance_filter=args['relative_abundance_filter'],
        digest=args['digest'], n=args['n'] * 10,verbose=True, DEBUG=args['DEBUG'], cores=cores,truth_set=args['truth_set'], output_dir=args['output_dir'])
    return matched_spectra
    
    