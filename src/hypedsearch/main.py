
import argparse
import sys

import utils, runner
from config_loader import Config

def stringtobool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def set_args(args) -> dict:
    use_params = stringtobool(args.config)
    config = Config()
    spectra_folder = args.spectra_folder if not use_params else config['spectra_dir']
    database_file = args.database_file if not use_params else config['database_file']
    output_dir = args.output_dir if not use_params else config['output_dir']
    min_peptide_len = args.min_peptide_len if not use_params else config['min_peptide_len']
    max_peptide_len = args.max_peptide_len if not use_params else config['max_peptide_len']
    ppm_tolerance = args.tolerance if not use_params else config['ppm_tolerance']
    precursor_tolerance = args.precursor_tolerance if not use_params else config['precursor_tolerance']
    verbose = stringtobool(args.verbose) if not use_params else config['verbose']
    peak_filter = args.peak_filter if not use_params else config['num_peaks']
    relative_abundance_filter = args.rel_abund_filter if not use_params else config['relative_abundance']
    digest = args.digest if not use_params else config['digest']
    cores = args.cores if not use_params else config['cores']
    n = args.n if not use_params else config['top_results']
    debug = config['debug']
    truth_set = config['truth_set']

    if not utils.is_dir(spectra_folder):
        print(f'Error: {spectra_folder} is not a real path. Path to directory with spectra files is necessary.')
        sys.exit(0)
    if not utils.is_fasta(database_file) or not utils.is_file(database_file):
        print(f'Error: {database_file} is not a valid .fasta file. .fasta file needed.')
        sys.exit(0)

    output_dir = utils.make_valid_dir_string(output_dir)
    utils.make_dir(output_dir)

    return {
        'spectra_folder': spectra_folder,
        'database_file': database_file,
        'output_dir': output_dir,
        'min_peptide_len': min_peptide_len,
        'max_peptide_len': max_peptide_len,
        'tolerance': ppm_tolerance,
        'precursor_tolerance': precursor_tolerance,
        'verbose': verbose, 
        'peak_filter': peak_filter, 
        'relative_abundance_filter': relative_abundance_filter,
        'digest': digest, 
        'DEBUG': debug, 
        'cores': cores,
        'n': n,
        'truth_set': truth_set
    }

def main(args: object) -> None:
    arguments = set_args(args)
    runner.run(arguments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')
    parser.add_argument('--spectra-folder', dest='spectra_folder', type=str, default='./', help='Path to folder containing spectra files.')
    parser.add_argument('--database-file', dest='database_file', type=str, default='./', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir', dest='output_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--config', dest='config', type=bool, default=True, help='Use the config.yaml file adjacent to main.py instead of using command line arguments. Default=True')
    parser.add_argument('--min-peptide-len', dest='min_peptide_len', type=int, default=5, help='Minimum peptide length to consider. Default=5')
    parser.add_argument('--max-peptide-len', dest='max_peptide_len', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    parser.add_argument('--tolerance', dest='tolerance', type=int, default=20, help='ppm tolerance to allow in search. Deafult=20')
    parser.add_argument('--precursor-tolerance', dest='precursor_tolerance', type=float, default=1, help='ppm tolerance to accept when matching precursor masses. Default=10')
    parser.add_argument('--peak-filter', dest='peak_filter', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
    parser.add_argument('--abundance-filter', dest='rel_abund_filter', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')
    parser.add_argument('--digest', dest='digest', type=str, default='', help='The digest performed. Default=None')
    parser.add_argument('--verbose', dest='verbose', type=bool, default=True, help='Extra printing to console during run. Default=True')
    parser.add_argument('--cores', dest='cores', type=int, default=1, help='The number of cores allowed to use when searching. Uses at least 1 and at most the number of available cores. Default=1')
    parser.add_argument('--n', dest='n', type=int, default=5, help='The number of alignments to keep per spectrum. Default=5')
    args = parser.parse_args()
    main(args)