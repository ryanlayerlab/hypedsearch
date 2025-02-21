import time
from typing import List

from src.constants import ION_CHARGES_TO_CONSIDER
from src.erik import parse_mzml
from src.fasta_utils import get_proteins_from_fasta

# from src.identification import get_matched_fragment
from src.lookups.protein_product_ion_db import create_protein_product_ion_db


def hypedsearch(
    ppm_tol: int,
    mzml_path: str,
    fasta_path: str,
    db_path: str,
    max_kmer_len: int,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
):
    spectra = list(parse_mzml(mzml_path=mzml_path))
    t0 = time.time()
    db = create_protein_product_ion_db(
        db_path=db_path,
        fasta_path=fasta_path,
        max_kmer_len=max_kmer_len,
        charges_to_consider=charges_to_consider,
    )
    t1 = time.time()
    duration = round(t1 - t0, 2)
    print(f"Preparing the database took {duration} secs = {round(duration/60, 2)} mins")
    # for spectrum in spectra:
    #     spectrum = process_spectrum()
    #     get_matched_fragment()

    spectrum = spectra[0]
    peak = spectrum.peaks[0]
    # get_matched_fragment
    pass
