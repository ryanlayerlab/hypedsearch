import sys
import argparse
import dataclasses

from dataclasses import dataclass

@dataclass
class Main_Arguments:
    spectra_folder: str = ''
    #database_file: str
    #output_dir: str
    #min_peptide_len: int
    #max_peptide_len: int
    #tolerance: float
    #precursor_tolerance: float
    #verbose: str
    #peak_filter: float
    #relative_abundance_filter: float
    #digest: str
    #DEBUG: bool
    #cores: int
    #n: int
    #truth_set: str

def main(args: Main_Arguments) -> None:
    print("Hello World!")

if __name__ == "__main__":
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
