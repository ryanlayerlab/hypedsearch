import os
import argparse
import sys
import distutils
from lookups.objects import Database
import lookups.utils, runner, computational_pipeline.database
from config_loader import Config
from postprocessing import summary, review

def get_spectra_file_paths(spectra_folder_path):
    spectra_file_paths = []
    for (root, _, file_paths) in os.walk(spectra_folder_path):
        for file_path in file_paths:
            spectra_file_paths.append(os.path.join(root, file_path))
    return spectra_file_paths

def populate_arguments(args) -> dict:
    use_params = lookups.utils.string_to_bool(args.config)
    if args.config: config = Config()
    spectra_file_path = args.spectra_file_path if not use_params else config['spectra_file_path']
    spectra_folder_path = args.spectra_folder_path if not use_params else config['spectra_folder_path']
    database_file_path = args.database_file_path if not use_params else config['database_file_path']
    output_folder_path = args.output_folder_path if not use_params else config['output_folder_path']
    max_peptide_length = args.max_peptide_length if not use_params else config['max_peptide_length']
    ppm_tolerance = args.tolerance if not use_params else config['ppm_tolerance']
    precursor_tolerance = args.precursor_tolerance if not use_params else config['precursor_tolerance']
    number_peaks = args.number_peaks if not use_params else config['number_peaks']
    relative_abundance = args.rel_abund_filter if not use_params else config['relative_abundance']
    digest_left = args.digest_left if not use_params else config['digest_left']
    digest_right = args.digest_right if not use_params else config['digest_right']
    number_cores = args.number_cores if not use_params else config['number_cores']
    create_kmer_database = args.create_kmer_database if not use_params else config['create_kmer_database']
    number_hybrids = args.number_hybrids if not use_params else config['number_hybrids']
    number_natives = args.number_natives if not use_params else config['number_natives']
    verbose = lookups.utils.string_to_bool(args.verbose) if not use_params else config['verbose']

    if spectra_file_path != '':
        spectra_file_paths = [spectra_file_path]
    else:
        spectra_file_paths = get_spectra_file_paths(spectra_folder_path)
        
    output_folder_path = lookups.utils.make_valid_dir_string(output_folder_path)
    

    return {'spectra_file_paths': spectra_file_paths,'database_file_path': database_file_path,'output_folder_path': output_folder_path,
        'max_peptide_length': max_peptide_length,'ppm_tolerance': ppm_tolerance,
        'precursor_tolerance': precursor_tolerance,'number_peaks': number_peaks, 
        'relative_abundance': relative_abundance, 'digest_left': digest_left, 'digest_right': digest_right, 
        'number_cores': number_cores,'number_hybrids': number_hybrids, 'number_natives': number_natives, 
        'create_kmer_database': create_kmer_database,'verbose': verbose}
    
def check_arguments(arguments):
    are_valid = True
    database_file_path = arguments['database_file_path']
    if not lookups.utils.is_fasta(database_file_path) or not lookups.utils.is_file(database_file_path): are_valid = False
    return are_valid

def main(args: object) -> None:
    arguments = populate_arguments(args)
    arguments_valid = check_arguments(arguments)
    if arguments_valid:
        runner.run(arguments)
    else:
        sys.exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')
    parser.add_argument('--spectra-file-paths', dest='spectra_file_paths', type=str, default='./', help='Path to folder containing spectra files.')
    parser.add_argument('--database-file-path', dest='database-file-path', type=str, default='./', help='protein database path')
    parser.add_argument('--output-folder-path', dest='output_folder_path', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--max-peptide-length', dest='max_peptide_length', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    parser.add_argument('--ppm-tolerance', dest='ppm_tolerance', type=int, default=20, help='ppm tolerance to allow in search. Deafult=20')
    parser.add_argument('--precursor-tolerance', dest='precursor_tolerance', type=float, default=1, help='ppm tolerance to accept when matching precursor masses. Default=10')
    parser.add_argument('--number-peaks', dest='number_peaks', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
    parser.add_argument('--relative-abundance', dest='relative_abundance', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')
    parser.add_argument('--digest-left', dest='digest_left', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts left of. Default=None')
    parser.add_argument('--digest-right', dest='digest_right', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts right of. Default=None')
    parser.add_argument('--number-cores', dest='number_cores', type=int, default=1, help='The number of cores allowed to use when searching. Uses at least 1 and at most the number of available cores. Default=1')
    parser.add_argument('--create_kmer_database', dest='create_kmer_database', type=lambda x:bool(distutils.util.strtobool(x)))
    parser.add_argument('--number-hybrids', dest='number_hybrids', type=int, default=5, help='The number of hybrid alignments to keep per spectrum. Default=5')
    parser.add_argument('--number-natives', dest='number_natives', type=int, default=5, help='The number of native alignments to keep per spectrum. Default=5')
    parser.add_argument('--verbose', dest='verbose', type=lambda x:bool(distutils.util.strtobool(x)))
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--no-config', dest='config', action='store_false')
    parser.set_defaults(config=True)
    args = parser.parse_args()
    main(args)