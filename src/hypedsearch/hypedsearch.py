import sys
import argparse
import dataclasses
import utils, runner
from dataclasses import dataclass

#There are two ways to run hypedsearch:
#1) command line
#2) file system Config file

#There are two ways to pass in dependencies: (spectra_file(s); database_file; output_location)
#1) in-memory
#2) file system       

@dataclass
class Hypedsearch_Arguments:
    min_peptide_len: int = 0
    max_peptide_len: int = 0
    tolerance: float = 0.0
    precursor_tolerance: float = 0.0
    verbose: str = ''
    peak_filter: float = 0.0
    relative_abundance_filter: float = 0.0
    digest: str = ''
    DEBUG: bool = True
    cores: int = 0
    n: int = 0
    truth_set: str = ''

@dataclass
class In_Memory_Arguments:
    spectra_file: list = None
    database_file: list = None
    hypedsearch_arguments: Hypedsearch_Arguments = None        

@dataclass
class IO_Arguments:    
    spectra_folder_path: str = ''
    database_file_path: str = ''
    output_dir_path: str = ''
    hypedsearch_arguments: Hypedsearch_Arguments = None        
        
@dataclass
class Config_File_Arguments:
    config_file_path: string = ''
        
@dataclass
class Main_Arguments:
    In_Memory_Arguments: In_Memory_Arguments = None
    IO_Arguments: IO_Arguments = None
    Config_File_Arguments: Config_File_Arguments = None    

def set_args(args) -> dict:
    final_args = {
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
    if typeof(args) == In_Memory_Arguments:
        return final_args
    else if typeof(args) == IO_Arguments:
        return final_args
    else if typeof(args) == Config_File_Arguments:
        return final_args
    else:
        raise Exception("Sorry, no numbers below zero")
        
def main(args: Main_Arguments) -> None:
    arguments = set_args(args)
    #runner.run(arguments)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')
    parser.add_argument('--spectra-folder-path', dest='spectra_folder_path', type=str, default='./', help='Path to folder containing spectra files.')
    parser.add_argument('--database-file-path', dest='database_file_path', type=str, default='./', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir-path', dest='output_dir_path', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--config-file-path', dest='config_file_path', type=bool, default=True, help='Use the config.yaml file adjacent to main.py instead of using command line arguments. Default=True')
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
