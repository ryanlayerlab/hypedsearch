import sys
from pathlib import Path
from time import time

repo_dir = Path(__file__).parents[3]
print(f"Adding {repo_dir} to sys.path")
assert repo_dir.name == "hypedsearch"
sys.path.append(str(repo_dir))

import multiprocessing as mp

from src.constants import SPECTRA_DIR
from src.hybrids_via_all_hybrids import run_hypedsearch_on_one_spectrum
from src.mass_spectra import Spectrum
from src.utils import get_time_in_diff_units, setup_logger

# Set constants
TESTING = True
FASTA_PATH = Path("fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta").absolute()
OUTPUT_DIR = Path("tmp/slurm_results").absolute()
PROTEIN_NAMES = Path("results/new_fasta/top_10_proteins.txt").absolute()
PRECURSOR_MZ_PPM_TOL = 20
NUM_PSMS = 10

assert FASTA_PATH.exists()
assert OUTPUT_DIR.exists()
assert PROTEIN_NAMES.exists()


def process_spectrum(spectrum: Spectrum):
    run_hypedsearch_on_one_spectrum(
        fasta_path=FASTA_PATH,
        precursor_mz_ppm_tol=PRECURSOR_MZ_PPM_TOL,
        output_dir=OUTPUT_DIR,
        spectrum=spectrum,
        protein_names=PROTEIN_NAMES,
        other_prots_for_fasta=FASTA_PATH,
        num_psms=NUM_PSMS,
    )


if __name__ == "__main__":

    logger = setup_logger()
    logger.info(f"Running with TESTING = {TESTING}")

    mzmls = list(SPECTRA_DIR.glob("*.mzML"))
    mzmls = [mzml.absolute() for mzml in mzmls]

    # For testing
    if TESTING:
        mzmls = [mzmls[0]]

    for mzml in mzmls:
        logger.info(f"Processing mzML file: {mzml}")
        t0 = time()

        spectra = Spectrum.from_mzml(mzml_path=mzml)

        if TESTING:
            spectra = spectra[:10]

        # # Non-parallel way
        # logger.info(f"Processing spectra serially")
        # for spectrum in spectra:
        #     process_spectrum(spectrum=spectrum)

        # Paralllel way
        num_procs = mp.cpu_count()  # Will use all 192 if available
        logger.info(f"Processing spectra in parallel with {num_procs} processes")
        with mp.Pool(processes=num_procs) as pool:
            results = pool.map(process_spectrum, spectra)

        logger.info(
            f"Finished processing {len(spectra)} spectra from mzml {mzml.stem} in {get_time_in_diff_units(time() - t0)}"
        )
