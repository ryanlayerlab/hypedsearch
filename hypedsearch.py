import logging
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from time import time
from typing import Dict, List, Optional, Required

import pandas as pd

from src.comet_utils import CometPSM, get_comet_protein_counts, run_comet_on_one_mzml
from src.config import AppConfig, Required, generate_cli
from src.constants import (
    COMET_EXECUTABLE,
    COMET_PARAMS,
    DEFAULT_MAX_K,
    DEFAULT_MIN_K,
    DEFAULT_NUM_PEAKS,
    DEFAULT_PPM_TOLERANCE,
    LOGS_DIR,
    MOUSE_PROTEOME,
    SCAN,
)
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.peptide_spectrum_comparison import (
    create_hybrids_fasta,
    extend_clusters,
    get_clusters_from_ions,
    get_possible_hybrids,
)
from src.peptides_and_ions import (
    Peptide,
    get_proteins_by_name,
    get_proteins_from_fasta,
    get_uniq_kmer_to_protein_map,
)
from src.protein_product_ion_database import (
    ProteinProductIonDb,
    create_db,
    get_positions_in_proteins_of_peak_matching_ions,
    get_product_ions_matching_spectrum,
)
from src.utils import make_directory, setup_logger

setup_logger(log_file=LOGS_DIR / f"{__name__}.log")
logger = logging.getLogger(__name__)


@dataclass
class HypedsearchConfig(AppConfig):
    mzml_dir: str
    mzml_path: str
    output_dir: str
    db_path: str
    # spectrum_file: Required[str]
    scan_num: int = -1
    fasta_path: str = MOUSE_PROTEOME
    top_n_proteins: int = 100
    num_peaks: int = DEFAULT_NUM_PEAKS
    comet_exe_path: str = COMET_EXECUTABLE
    comet_params_path: str = COMET_PARAMS
    cleanup: bool = True
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE
    precursor_ppm_tol: float = DEFAULT_PPM_TOLERANCE

    # min_k: int = DEFAULT_MIN_K
    # max_k: int = DEFAULT_MAX_K
    # mzml_paths: List[str]

    def __post_init__(self):
        self.mzml_dir = Path(self.mzml_dir).absolute()
        self.mzml_path = Path(self.mzml_path).absolute()
        self.output_dir = Path(self.output_dir).absolute()
        self.fasta_path = Path(self.fasta_path).absolute()
        self.comet_exe_path = Path(self.comet_exe_path).absolute()
        self.comet_params_path = Path(self.comet_params_path).absolute()
        self.db_path = Path(self.db_path).absolute()

    @staticmethod
    def get_variable_help_strings() -> Dict[str, str]:
        help_strings = {}
        help_strings["db_path"] = "Path to database file."
        help_strings["fasta_path"] = "Path to FASTA file."
        help_strings["mzml_path"] = "Path to mzML file."
        help_strings["num_peaks"] = (
            "Filter the number of spectrum peaks considered to the top N=num_peaks intensity peaks. "
            "If 0, no peak filtering."
        )
        help_strings["scan_num"] = "Scan number to analyze. If -1, analyze all scans."
        help_strings["cleanup"] = (
            "If True, remove intermediate files leaving only the Hypedsearch .csv result file"
        )
        help_strings["peak_to_ion_ppm_tol"] = (
            "PPM tolerance when matching spectrum peaks to product ions"
        )
        help_strings["precursor_ppm_tol"] = (
            "Proposed peptides must be within this PPM tolerance of the spectrum precursor m/z"
        )
        return help_strings

    @staticmethod
    def get_variable_flags() -> Dict[str, str]:
        variable_flags = {}
        variable_flags["db_path"] = "d"
        variable_flags["mzml_path"] = "m"
        variable_flags["fasta_path"] = "f"
        variable_flags["num_peaks"] = "np"
        variable_flags["comet_params_path"] = "cpp"
        variable_flags["scan_num"] = "s"
        return variable_flags

    def validate() -> bool:
        # Add validation logic if needed
        return True


def hypedsearch(config: HypedsearchConfig):
    logger.info(f"Configuration: {config}")
    hs_start_time = time()

    # Run Comet on all the MZML files in the specified directory and its subdirectories
    logger.info("Starting Comet run 1...")
    start_time = time()
    mzml_paths = list(config.mzml_dir.rglob("*.mzML"))
    comet_run_1_txts = []
    for mzml in mzml_paths:
        logger.info(f"Running Comet on MZML file: {mzml}")
        # Create the output directory for this mzML file
        mzml_output_dir = config.output_dir / mzml.stem
        comet_run_1_txts.append(
            run_comet_on_one_mzml(
                template_comet_params_path=config.comet_params_path,
                fasta_path=config.fasta_path,
                output_dir=mzml_output_dir,
                comet_exe_path=config.comet_exe_path,
                mzml_path=mzml,
                overwrite=False,
                output_file_stem="run_1",
                keep_params=False,
            )
        )
    logger.info(f"Comet run 1 completed in {time() - start_time:.2f} seconds.")

    # Get the number of times each protein appears in the Comet run 1 output
    logger.info("Counting proteins from Comet run 1...")
    comet_run_1_rows = []
    for comet_run_1_txt in comet_run_1_txts:
        comet_run_1_rows += CometPSM.from_txt(file_path=comet_run_1_txt)
    protein_counts = get_comet_protein_counts(
        comet_rows=comet_run_1_rows, shorten_names=False
    )

    # Create database of top N most common proteins
    most_common_proteins = [
        protein_and_count[0]
        for protein_and_count in protein_counts.most_common(config.top_n_proteins)
    ]
    db_proteins = get_proteins_by_name(
        protein_names=most_common_proteins, fasta_path=config.fasta_path
    )
    # db_proteins = get_proteins_from_fasta(fasta_path=config.fasta_path)
    protein_id_to_name_map = {prot.id: prot.name for prot in db_proteins}
    logger.info("Getting unique kmer to protein map...")
    uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(proteins=db_proteins)
    db = create_db(
        db_path=config.db_path,
        db_proteins=db_proteins,
        uniq_kmer_to_protein_map=uniq_kmer_to_protein_map,
        overwrite=False,
    )

    # Get spectrum and filter peaks
    spectra = Spectrum.from_mzml(mzml_path=config.mzml_path)

    if config.scan_num == -1:
        for spectrum in spectra:
            logger.info(f"Processing spectrum {spectrum.scan_num}...")
            run_hypedsearch_on_one_spectrum(
                config=config,
                db=db,
                spectrum=spectrum,
                protein_id_to_name_map=protein_id_to_name_map,
                uniq_kmer_to_protein_map=uniq_kmer_to_protein_map,
                db_proteins=db_proteins,
            )
    else:
        logger.info(f"Processing spectrum {config.scan_num}...")
        spectrum = list(
            filter(lambda spectrum: spectrum.scan_num == config.scan_num, spectra)
        )
        assert (
            len(spectrum) == 1
        ), f"Scan number must be unique. There were {len(spectrum)} spectra with scan number {config.scan_num}."
        spectrum = spectrum[0]

        run_hypedsearch_on_one_spectrum(
            config=config,
            db=db,
            spectrum=spectrum,
            protein_id_to_name_map=protein_id_to_name_map,
            uniq_kmer_to_protein_map=uniq_kmer_to_protein_map,
            db_proteins=db_proteins,
        )
    logger.info(f"Hypedsearch completed in {time() - hs_start_time:.2f} seconds.")

    # Cleanup
    if config.cleanup:
        logger.info("Cleaning up intermediate files...")
        # Remove the Comet run 1 output files
        (config.output_dir / f"{config.mzml_path.stem}/run_1.txt").unlink()
        (config.output_dir / f"{config.mzml_path.stem}/run_1.pep.xml").unlink()
        config.db_path.unlink()


def run_hypedsearch_on_one_spectrum(
    config: HypedsearchConfig,
    db: ProteinProductIonDb,
    spectrum: Spectrum,
    protein_id_to_name_map: Dict[str, str],
    uniq_kmer_to_protein_map: Dict[str, List[str]],
    db_proteins: List[Peptide],
):
    # Perform peak filtering
    if config.num_peaks > 0:
        spectrum.filter_to_top_n_peaks(n=config.num_peaks)

    # Find ions that match spectrum peaks
    peaks_with_matches = get_product_ions_matching_spectrum(
        spectrum=spectrum,
        db=db,
        peak_product_ion_ppm_tolerance=config.peak_to_ion_ppm_tol,
    )
    positioned_ions = get_positions_in_proteins_of_peak_matching_ions(
        peaks_with_matches=peaks_with_matches,
        kmer_to_protein_map=uniq_kmer_to_protein_map,
        db=db,
    )

    # Get clusters
    logger.info("Getting clusters...")
    t0 = time()
    clusters = get_clusters_from_ions(ions=positioned_ions)
    logger.info(f"Getting clusters took {time() - t0:.2f} seconds")
    logger.info("Extending clusters...")
    t0 = time()
    extended_clusters = extend_clusters(
        spectrum_clusters=clusters, spectrum=spectrum, db=db
    )
    logger.info(f"Extending clusters took {time() - t0:.2f} seconds")

    # Get possible hybrids
    logger.info("Creating possible hybrids...")
    t0 = time()
    possible_hybrids = get_possible_hybrids(
        extended_clusters=extended_clusters,
        spectrum=spectrum,
        precursor_mz_ppm_tolerance=config.precursor_ppm_tol,
    )
    logger.info(f"Getting possible hybrids took {time() - t0:.2f} seconds")

    scan_dir = config.output_dir / f"{config.mzml_path.stem}/scan={spectrum.scan_num}"
    make_directory(scan_dir)
    # Create new FASTA with hybrids
    logger.info("Creating new FASTA with hybrids...")
    t0 = time()
    new_fasta_path = scan_dir / "hybrids.fasta"
    create_hybrids_fasta(
        db_proteins=db_proteins,
        hybrids=possible_hybrids,
        protein_id_to_name_map=protein_id_to_name_map,
        fasta_path=new_fasta_path,
    )
    logger.info(f"Creating new FASTA took {time() - t0:.2f} seconds")

    # Run Comet
    logger.info("Running Comet using new FASTA...")
    comet_run_2_txt = run_comet_on_one_mzml(
        template_comet_params_path=config.comet_params_path,
        fasta_path=new_fasta_path,
        output_dir=scan_dir,
        comet_exe_path=config.comet_exe_path,
        mzml_path=config.mzml_path,
        output_file_stem=f"run_2",
        keep_params=False,
    )

    # Get Comet run 1 vs 2 results
    logger.info("Getting Comet run 1 and 2 results...")
    t0 = time()
    run_1_path = config.output_dir / f"{config.mzml_path.stem}/run_1.txt"
    comet_run_1 = CometPSM.from_txt(file_path=run_1_path, as_df=True)
    comet_run_1 = comet_run_1[comet_run_1[SCAN] == spectrum.scan_num]
    comet_run_1.reset_index(inplace=True, drop=True)
    comet_run_1["run"] = 1

    comet_run_2 = CometPSM.from_txt(file_path=comet_run_2_txt, as_df=True)
    comet_run_2 = comet_run_2[comet_run_2[SCAN] == spectrum.scan_num]
    comet_run_2.reset_index(inplace=True, drop=True)
    comet_run_2["run"] = 2

    df = pd.concat([comet_run_1, comet_run_2], ignore_index=True)
    df.to_csv(
        scan_dir / "hs_results.csv",
        index=False,
    )
    logger.info(f"Getting Comet results took {time() - t0:.2f} seconds")


if __name__ == "__main__":
    cli = generate_cli(cls=HypedsearchConfig, main_fcn=hypedsearch)
    cli()
