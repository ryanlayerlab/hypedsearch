import logging
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from time import time
from typing import Dict, List, Optional, Required

import pandas as pd

from config import AppConfig, Required, generate_cli
from src.comet_utils import CometPSM, get_comet_protein_counts, run_comet
from src.constants import (
    COMET_EXECUTABLE,
    COMET_PARAMS,
    DEFAULT_MAX_K,
    DEFAULT_MIN_K,
    DEFAULT_NUM_PEAKS,
    LOGS_DIR,
    MOUSE_PROTEOME,
    SCAN,
)
from src.mass_spectra import Spectrum, get_specific_spectrum_by_sample_and_scan_num
from src.peptide_spectrum_comparison import (
    create_new_fasta_including_hybrids,
    get_clusters_from_ions,
    get_extended_clusters,
    get_possible_hybrids,
)
from src.peptides_and_ions import (
    get_proteins_by_name,
    get_proteins_from_fasta,
    get_uniq_kmer_to_protein_map,
)
from src.protein_product_ion_database import (
    create_db,
    get_product_ions_matching_spectrum,
    process_peak_matching_ions,
)
from src.utils import setup_logger

setup_logger(log_file=LOGS_DIR / f"{__name__}.log")
logger = logging.getLogger(__name__)


@dataclass
class HypedsearchConfig(AppConfig):
    mzml_dir: str
    # db_path: Required[str]
    # spectrum_file: Required[str]
    # scan_num: int = -1
    # fasta_path: str = MOUSE_PROTEOME
    # num_peaks: int = DEFAULT_NUM_PEAKS
    # comet_exe_path: str = COMET_EXECUTABLE
    # comet_params_path: str = COMET_PARAMS
    # min_k: int = DEFAULT_MIN_K
    # max_k: int = DEFAULT_MAX_K
    # mzml_paths: List[str]


@dataclass
class HypedsearchConfig(AppConfig):
    mzml_dir: Path
    mzml_path: Path
    output_dir: Path
    db_path: Path
    # spectrum_file: Required[str]
    scan_num: int = -1
    fasta_path: Path = MOUSE_PROTEOME
    top_n_proteins: int = 100
    num_peaks: int = DEFAULT_NUM_PEAKS
    comet_exe_path: Path = COMET_EXECUTABLE
    comet_params_path: Path = COMET_PARAMS
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

    # Run Comet on all the MZML files in the specified directory
    logger.info("Starting Comet run 1...")
    start_time = time()
    mzml_paths = config.mzml_dir.glob("*.mzML")
    comet_run_1_txts = []
    for mzml in mzml_paths:
        logger.info(f"Running Comet on MZML file: {mzml}")
        comet_run_1_txts.append(
            run_comet(
                template_comet_params_path=config.comet_params_path,
                fasta_path=config.fasta_path,
                parent_output_dir=config.output_dir,
                comet_exe_path=config.comet_exe_path,
                mzml_path=mzml,
                overwrite=False,
                output_file_stem="run_1",
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
        uniq_kmers=uniq_kmer_to_protein_map.keys(),
        overwrite=False,
    )

    # Get spectrum and filter peaks
    spectra = Spectrum.from_mzml(mzml_path=config.mzml_path)
    spectrum = list(
        filter(lambda spectrum: spectrum.scan_num == config.scan_num, spectra)
    )
    assert (
        len(spectrum) == 1
    ), f"Scan number must be unique. There were {len(spectrum)} spectra with scan number {config.scan_num}."
    spectrum = spectrum[0]
    if config.num_peaks > 0:
        spectrum.filter_to_top_n_peaks(n=config.num_peaks)

    # Find ions that match spectrum peaks
    peaks_with_matches = get_product_ions_matching_spectrum(
        spectrum=spectrum,
        db=db,
    )
    positioned_ions = process_peak_matching_ions(
        peaks_with_matches=peaks_with_matches,
        kmer_to_protein_map=uniq_kmer_to_protein_map,
        db=db,
    )

    # Get clusters
    clusters = get_clusters_from_ions(ions=positioned_ions)
    extended_clusters = get_extended_clusters(
        spectrum_clusters=clusters, spectrum=spectrum, db=db
    )

    # Get possible hybrids
    possible_hybrids = get_possible_hybrids(
        extended_clusters=extended_clusters, spectrum=spectrum
    )

    # Create new FASTA with hybrids
    sample = config.mzml_path.stem
    new_fasta_path = (
        config.output_dir / f"{sample}/scan={config.scan_num}_hybrids.fasta"
    )
    create_new_fasta_including_hybrids(
        db_proteins=db_proteins,
        hybrids=possible_hybrids,
        protein_id_to_name_map=protein_id_to_name_map,
        fasta_path=new_fasta_path,
    )

    # Run Comet
    comet_run_2_txt = run_comet(
        template_comet_params_path=config.comet_params_path,
        fasta_path=new_fasta_path,
        parent_output_dir=config.output_dir,
        comet_exe_path=config.comet_exe_path,
        mzml_path=config.mzml_path,
        output_file_stem="run_2",
    )

    # Get Comet run 1 vs 2 results
    run_1_path = config.output_dir / f"{config.mzml_path.stem}/run_1.txt"
    comet_run_1 = CometPSM.from_txt(file_path=run_1_path, as_df=True)
    comet_run_1 = comet_run_1[comet_run_1[SCAN] == config.scan_num]
    comet_run_1.reset_index(inplace=True, drop=True)
    comet_run_1["run"] = 1

    comet_run_2 = CometPSM.from_txt(file_path=comet_run_2_txt, as_df=True)
    comet_run_2 = comet_run_1[comet_run_2[SCAN] == config.scan_num]
    comet_run_2.reset_index(inplace=True, drop=True)
    comet_run_2["run"] = 2

    df = pd.concat([comet_run_1, comet_run_2], ignore_index=True)
    df.to_csv(
        config.output_dir
        / f"{config.mzml_path.stem}/scan={config.scan_num}_hs_results.csv",
        index=False,
    )
    logger.info("Hypedsearch completed!")


if __name__ == "__main__":
    cli = generate_cli(cls=HypedsearchConfig, main_fcn=hypedsearch)
    cli()
