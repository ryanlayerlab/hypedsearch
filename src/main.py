import os
import argparse
import sys
import distutils
from src.constants.objects import Database
import src.constants.utils, runner, database
from config_loader import Config
from postprocessing import summary, review

def string_to_bool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def get_spectra_files(spectra_folder):
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    return spectra_files

def get_database_file(database_file_path):
    return database.build(database_file_path)

def set_args(args) -> dict:
    use_params = string_to_bool(args.config)
    if args.config:
        config = Config()
    spectra_folder = args.spectra_folder if not use_params else config['spectra_dir']
    database_file_path = args.database_file if not use_params else config['database_file']
    output_dir = args.output_dir if not use_params else config['output_dir']
    min_peptide_len = args.min_peptide_len if not use_params else config['min_peptide_len']
    max_peptide_len = args.max_peptide_len if not use_params else config['max_peptide_len']
    ppm_tolerance = args.tolerance if not use_params else config['ppm_tolerance']
    precursor_tolerance = args.precursor_tolerance if not use_params else config['precursor_tolerance']
    verbose = string_to_bool(args.verbose) if not use_params else config['verbose']
    peak_filter = args.peak_filter if not use_params else config['num_peaks']
    relative_abundance_filter = args.rel_abund_filter if not use_params else config['relative_abundance']
    digest_left = args.digest_left if not use_params else config['digest_left']
    digest_right = args.digest_right if not use_params else config['digest_right']
    cores = args.cores if not use_params else config['cores']
    make_new = args.new_db if not use_params else config['new_db']
    num_hybrids = args.num_hybrids if not use_params else config['top_hybrids']
    num_natives = args.num_natives if not use_params else config['top_natives']
    debug = args.debug if not use_params else config['debug']

    if not utils.is_dir(spectra_folder):
        print(f'Error: {spectra_folder} is not a real path. Path to directory with spectra files is necessary.')
        sys.exit(0)
    if not utils.is_fasta(database_file_path) or not utils.is_file(database_file_path):
        print(f'Error: {database_file_path} is not a valid .fasta file. .fasta file needed.')
        sys.exit(0)

    spectra_files = get_spectra_files(spectra_folder)
    database_file = get_database_file(database_file_path)
    output_dir = utils.make_valid_dir_string(output_dir)
    utils.make_dir(output_dir)

    return {'spectra_files': spectra_files,'database_file': database_file,'output_dir': output_dir,
        'min_peptide_len': min_peptide_len,'max_peptide_len': max_peptide_len,'tolerance': ppm_tolerance,
        'precursor_tolerance': precursor_tolerance,'verbose': verbose, 'peak_filter': peak_filter, 
        'relative_abundance_filter': relative_abundance_filter, 'digest_left': digest_left, 'digest_right': digest_right, 'DEBUG': debug, 
        'cores': cores,'num_hybrids': num_hybrids, 'num_natives': num_natives, 'new_db': make_new}
    
def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'False'

def main(args: object) -> None:
    arguments = set_args(args)
    runner.run(arguments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')
    parser.add_argument('--spectra-folder', dest='spectra_folder', type=str, default='./', help='Path to folder containing spectra files.')
    parser.add_argument('--database-file', dest='database_file', type=str, default='./', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir', dest='output_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--no-config', dest='config', action='store_false')
    parser.set_defaults(config=True)
    parser.add_argument('--min-peptide-len', dest='min_peptide_len', type=int, default=5, help='Minimum peptide length to consider. Default=5')
    parser.add_argument('--max-peptide-len', dest='max_peptide_len', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    parser.add_argument('--tolerance', dest='tolerance', type=int, default=20, help='ppm tolerance to allow in search. Deafult=20')
    parser.add_argument('--precursor-tolerance', dest='precursor_tolerance', type=float, default=1, help='ppm tolerance to accept when matching precursor masses. Default=10')
    parser.add_argument('--peak-filter', dest='peak_filter', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
    parser.add_argument('--abundance-filter', dest='rel_abund_filter', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')

    parser.add_argument('--digest-left', dest='digest_left', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts left of. Default=None')
    parser.add_argument('--digest-right', dest='digest_right', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts right of. Default=None')

    parser.add_argument('--verbose', dest='verbose', type=lambda x:bool(distutils.util.strtobool(x)))
    parser.add_argument('--cores', dest='cores', type=int, default=1, help='The number of cores allowed to use when searching. Uses at least 1 and at most the number of available cores. Default=1')
    parser.add_argument('--new-database', dest='new_db', type=lambda x:bool(distutils.util.strtobool(x)))
    parser.add_argument('--num-hybrids', dest='num_hybrids', type=int, default=5, help='The number of hybrid alignments to keep per spectrum. Default=5')
    parser.add_argument('--num-natives', dest='num_natives', type=int, default=5, help='The number of native alignments to keep per spectrum. Default=5')

    parser.add_argument('--debug', dest='debug', type=bool, default=False, help='The number of alignments to keep per spectrum. Default=5')
    args = parser.parse_args()

    main(args)
