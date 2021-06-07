import sys
import argparse
import dataclasses
import utils, runner
from config_loader import Config
        
def set_args(args) -> dict:
    config = Config()
    min_peptide_len = args.min_peptide_len if not args.use_config_file else config['min_peptide_len']
    max_peptide_len = args.max_peptide_len if not args.use_config_file else config['max_peptide_len']
    ppm_tolerance = args.tolerance if not args.use_config_file else config['ppm_tolerance']
    precursor_tolerance = args.precursor_tolerance if not args.use_config_file else config['precursor_tolerance']
    verbose = args.verbose if not args.use_config_file else config['verbose']
    peak_filter = args.peak_filter if not args.use_config_file else config['num_peaks']
    relative_abundance_filter = args.rel_abund_filter if not args.use_config_file else config['relative_abundance']
    digest = args.digest if not args.use_config_file else config['digest']
    cores = args.cores if not args.use_config_file else config['cores']
    n = args.n if not args.use_config_file else config['top_results']
    debug = args.cores if not args.use_config_file else config['cores']
    truth_set = args.truth_set if not args.use_config_file else config['truth_set']
    spectra_file_path = args.spectra_file_path if not args.use_config_file else config['spectra_file_path']
    database_file_path = args.database_file_path if not args.database_file_path else config['database_file_path']
    output_dir = args.output_dir if not args.use_config_file else config['output_dir']
    #spectra_file = args.spectra_file is not args.use_io args.spectra_file else utils.load_spectra_file(spectra_file_path)
    spectra_file = None
    #database_file = args.database_file is not args.use_io args.database_file else utils.load_database_file(database_file_path)
    database_file = None
    
    return {'min_peptide_len': min_peptide_len, 'max_peptide_len': max_peptide_len, 'tolerance': ppm_tolerance,
        'precursor_tolerance': precursor_tolerance, 'verbose': verbose, 'peak_filter': peak_filter, 
        'relative_abundance_filter': relative_abundance_filter,'digest': digest, 'DEBUG': debug, 
        'cores': cores,'n': n,'truth_set': truth_set,'spectra_file': spectra_file,'database_file' : database_file}        
        
def main(args: object) -> None:
    arguments = set_args(args)
    runner.run(arguments)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')
    #CLI hypedsearch
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
    parser.add_argument('--truth_set', dest='truth_set', type=str, default='', help='Not sure what this does. Default=None')
    #CLI In Memory Datasets
    parser.add_argument('--use_in_memory', dest='use_in_memory', type=bool, default=True, help='Use In Memory Datasets using file system. Default=True')
    parser.add_argument('--spectra_file', dest='spectra_file', type=list, default=None, help='In Memory Spectra File. Default=None')
    parser.add_argument('--database_file', dest='database_file', type=list, default=None, help='In Memory Database File. Default=None')
    #CLI I/O Datasets
    parser.add_argument('--use_io', dest='use_io', type=bool, default=False, help='Use file system datasets instead of in-memory. Default=False')
    parser.add_argument('--spectra-folder-path', dest='spectra_file_path', type=str, default='./', help='Path to .mzml spectra file.')
    parser.add_argument('--database-file-path', dest='database_file_path', type=str, default='./', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir-path', dest='output_dir_path', type=str, default='~/', help='Directory to save all figures. Default=~/')
    #Config File (always I/O Datasets)
    parser.add_argument('--use_config_file', dest='use_config_file', type=bool, default=False, help='Use .yaml config file instead of command line arguments. Default=False')
    parser.add_argument('--config-file-path', dest='config_file_path', type=bool, default=True, help='Use the config.yaml file adjacent to main.py instead of using command line arguments. Default=True')
    args = parser.parse_args()    
    main(args)
