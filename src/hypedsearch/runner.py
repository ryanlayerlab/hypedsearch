import os
import utils
import identification
from postprocessing import summary, review

def run(args: dict) -> None:
    id_spectra_arguments = identification.create_id_spectra_arguments(args)
    matched_spectra = identification.id_spectra(id_spectra_arguments)
    print('\nFinished search. Writting results to {}...'.format(args['output_dir']))
    #output_dir = args['output_dir']
    #summary.generate(matched_spectra, output_dir)
    
    