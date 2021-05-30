import sys
import dataclasses

from dataclasses import dataclass

@dataclass
class Main_Arguments:
    spectra_folder: str
    database_file: str
    output_dir: str
    min_peptide_len: int
    max_peptide_len: int
    tolerance: float
    precursor_tolerance: float
    verbose: str
    peak_filter: float
    relative_abundance_filter: float
    digest: str
    DEBUG: bool
    cores: int
    n: int
    truth_set: str

def main(args: Main_Arguments) -> None:
    print("Hello World!")

if __name__ == "__main__":
    main()
