'''
exec.py

Author: Zachary McGrath 
Date: 6 April 2020

Executor for the program
In charge of the flow of the program
'''
import os
import identification
from postprocessing import summary, review
import multiprocessing as mp

def run(args: dict) -> None:
    '''
    Executing function for the program

    Inputs:
        args:   object arguments from main. Should be validated in main. Attributes of args:
            spectra_folder:             (str) full path the the directory containing all spectra files
            database_file:              (str) full path to the .fasta database file
            output_dir:                 (str) full path the the directory to save output to
            min_peptide_len:            (int) minimum peptide length to consider
            max_peptide_len:            (int) maximum peptide length to consider
            tolerance:                  (int) the ppm tolerance to allow in search
            precursor_tolerance:        (int) the ppm tolerance to allow when matching precursors
            peak_filter:                (int) the number of peaks to filter by 
            relative_abundance_filter:  (float) the percentage of the total abundance a peak must
                                            be to pass the filter
            digest:                     (str) the digest performed
            missed_cleavages:           (int) the number of missed cleavages allowed in digest
            verbose:                    (bool) extra printing
            cores:                      (int) the number of cores allowed to use
            n:                          (int) the number of alignments to keep per spectrum
            DEBUG:                      (bool) debuging print messages. Default=False
    Outputs:
        None
    '''
    # get all the spectra file names
    spectra_files = []
    for (root, _, filenames) in os.walk(args['spectra_folder']):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))

    # make sure cores is: 1 <= cores <= cpu cores
    cores = max(1, args['cores'])
    cores = min(cores, mp.cpu_count() - 1)

    matched_spectra = identification.id_spectra(
        spectra_files, args['database_file'], 
        min_peptide_len=args['min_peptide_len'], 
        max_peptide_len=args['max_peptide_len'], 
        ppm_tolerance=args['tolerance'], 
        precursor_tolerance=args['precursor_tolerance'],
        peak_filter=args['peak_filter'],
        relative_abundance_filter=args['relative_abundance_filter'],
        digest=args['digest'], 
        n=args['n'] * 10,
        verbose=True, 
        DEBUG=args['DEBUG'], 
        cores=cores,
        truth_set=args['truth_set'], 
        output_dir=args['output_dir']
    )
    print('\nFinished search. Writting results to {}...'.format(args['output_dir']))

    # matched_spectra = review.tie_breaker(matched_spectra, '', args['n'])

    summary.generate(matched_spectra, args['output_dir'])
    
    
