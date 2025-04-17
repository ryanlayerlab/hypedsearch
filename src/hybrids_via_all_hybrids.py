import logging
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from time import time
from typing import List, Optional, Union
from venv import logger

import click

from src.click_utils import ClickOptions
from src.comet_utils import run_comet_on_one_mzml
from src.constants import (
    DEFAULT_MAX_K,
    DEFAULT_MIN_K,
    DEFAULT_NUM_COMET_PSMS,
    DEFAULT_NUM_PEAKS,
    DEFAULT_PPM_TOLERANCE,
    PROTON_MASS,
    WATER_MASS,
    HybridFormingMethods,
)
from src.hybrids_via_clusters import HybridPeptide, SeqWithMass, create_hybrids_fasta
from src.mass_spectra import Spectrum
from src.peptides_and_ions import (
    Peptide,
    compute_peptide_mz,
    get_proteins_by_name,
    get_uniq_kmer_to_protein_map,
)
from src.protein_abundance import get_most_common_proteins, get_protein_comet_counts
from src.protein_product_ion_database import (
    ProteinProductIonDb,
    get_positions_in_proteins_of_peak_matching_ions,
    get_product_ions_matching_spectrum,
)
from src.sql_database import Sqlite3Database
from src.utils import (
    get_time_in_diff_units,
    hypedsearch_output_stem,
    log_params,
    relative_ppm_tolerance_in_daltons,
    setup_logger,
)

logger = logging.getLogger(__name__)


@dataclass
class HSOutput:
    hybrids: List[HybridPeptide]
    comet_txt: Optional[Path]
    run_time: Optional[float] = None


def form_all_hybrids(
    precursor_charge: int,
    precursor_mz: float,
    proteins: List[Peptide],
    precursor_mz_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    min_k: int = DEFAULT_MIN_K,
    max_k: int = DEFAULT_MAX_K,
) -> List[HybridPeptide]:
    logger.info("Forming all hybrids...")
    start_time = time()
    uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(
        proteins=proteins, min_k=min_k, max_k=max_k
    )
    db_rows = [
        SeqWithMass.from_seq(seq=seq, charge=precursor_charge)
        for seq in uniq_kmer_to_protein_map
    ]
    db = Sqlite3Database()
    table_name = "kmers"
    db.create_table_from_dataclass(table_name=table_name, obj=SeqWithMass)
    db.insert_dataclasses(table_name=table_name, data_classes=db_rows)
    db.add_index(table_name=table_name, index_name="mass", colms_to_index=["mz"])

    # For each b-sequence search the y-sequence database to find the hybrids
    # that would produce a peptide within the given PPM tolerance of the precursor m/z
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
    )
    adjusted_precursor_mz = precursor_mz + (WATER_MASS / precursor_charge) + PROTON_MASS
    potench_hybrids = []
    for b_seq in uniq_kmer_to_protein_map:
        b_seq_mz = compute_peptide_mz(aa_seq=b_seq, charge=precursor_charge)
        lower_bdd = adjusted_precursor_mz - mz_tolerance - b_seq_mz
        upper_bdd = adjusted_precursor_mz + mz_tolerance - b_seq_mz
        query = f"""
            SELECT
                *
            FROM {table_name} as ion
            WHERE ion.mz BETWEEN {lower_bdd} AND {upper_bdd}
        """
        matches = [match["seq"] for match in db.read_query(query=query)]
        for y_seq in matches:
            potench_hybrids.append(
                HybridPeptide(
                    b_seq=b_seq,
                    y_seq=y_seq,
                    b_prot_ids=uniq_kmer_to_protein_map[b_seq],
                    y_prot_ids=uniq_kmer_to_protein_map[y_seq],
                )
            )

    logger.info(
        f"Forming all hybrids took {get_time_in_diff_units(time() - start_time)}"
    )
    return potench_hybrids


def determe_hybrid_forming_method(
    protein_names: Optional[Path],
) -> HybridFormingMethods:
    if protein_names is not None:
        return HybridFormingMethods.all
    else:
        raise RuntimeError("Unrecognized hybrid forming method!")


def run_hypedsearch_on_spectrum(
    mzml: Path,
    fasta_path: Path,
    precursor_mz_ppm_tol: float,
    scan: int,
    output_dir: Path,
    num_peaks: int = 0,
    min_k: int = DEFAULT_MIN_K,
    max_k: int = DEFAULT_MAX_K,
    protein_names: Optional[Path] = None,
    other_prots_for_fasta: Optional[Union[List[Peptide], Path]] = None,
    num_psms: int = DEFAULT_NUM_COMET_PSMS,
) -> HSOutput:
    assert output_dir.exists(), f"Output directory is expected to exist but doesn't:\n{output_dir}"  
    logger.info(f"Running HS on scan={scan} of MZML:\n{mzml}")
    start_time = time()

    # Load spectra
    spectra = Spectrum.from_mzml(mzml_path=mzml)
    spectra = list(filter(lambda spectrum: spectrum.scan_num == scan, spectra))
    assert (
        len(spectra) == 1
    ), f"There are {len(spectra)} spectra corresponding to scan {scan} in {mzml}"
    spectrum = spectra[0]

    # Preprocess spectrum
    if num_peaks > 0:
        spectrum.filter_to_top_n_peaks(n=num_peaks)

    # Form hybrids
    hybrid_method = determe_hybrid_forming_method(protein_names)
    if hybrid_method == HybridFormingMethods.all:
        proteins = get_proteins_by_name(
            protein_names=protein_names, fasta_path=fasta_path
        )
        potench_hybrids = form_all_hybrids(
            precursor_charge=spectrum.precursor_charge,
            precursor_mz=spectrum.precursor_mz,
            proteins=proteins,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            min_k=min_k,
            max_k=max_k,
        )

        # For setting the hybrids FASTA information below
        prot_id_to_name_map = {prot.id: prot.name for prot in proteins}

    if len(potench_hybrids) == 0:
        logger.info("No hybrids were formed! Exiting...")
        sys.exit(0)

    # Prepare the hybrids to be written to a FASTA file by setting their description and name
    logger.info
    for hybrid in potench_hybrids:
        hybrid.set_fasta_info(prot_id_to_name_map=prot_id_to_name_map)

    # Create FASTA file and run Comet
    logger.info("Creating FASTA file and running Comet...")
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir).absolute()

        # Create new FASTA with hybrids and other proteins the user wants included
        new_fasta_path = (
            temp_dir / f"{spectrum.mzml.stem}_{spectrum.scan_num}_hybrids.fasta"
        )
        create_hybrids_fasta(
            hybrids=potench_hybrids,
            fasta_path=new_fasta_path,
            other_prots=other_prots_for_fasta,
        )

        # Run Comet
        output_stem = hypedsearch_output_stem(mzml_stem=spectrum.mzml.stem)
        comet_txt = None
        try:
            comet_txt = run_comet_on_one_mzml(
                fasta=new_fasta_path,
                mzml=spectrum.mzml,
                output_dir=temp_dir,
                scan=spectrum.scan_num,
                stem=output_stem,
                num_psms=num_psms,
            )
            shutil.copy2(comet_txt, output_dir)
            comet_txt = output_dir / comet_txt.name
        except RuntimeError as err:
            logger.info(f"Comet failed! Here's the error message:\n{err}")

    elapsed_time = time() - start_time
    logger.info(f"Running HS took {get_time_in_diff_units(elapsed_time)}")
    return HSOutput(hybrids=potench_hybrids, comet_txt=comet_txt, run_time=elapsed_time)


options = ClickOptions()


def run_hypedsearch_on_scan_cli_options():
    def decorator(func):
        for opt in reversed(
            [
                options.scan(),
                options.num_psms(),
                options.max_k(),
                options.min_k(),
                options.precursor_mz_ppm_tol(),
                options.num_peaks(),
                options.output_dir(),
                options.mzml(),
                options.fasta_path(),
                options.protein_names(),
            ]
        ):
            func = opt(func)
        return func

    return decorator


@click.command(
    name="run-comet-on-all-hybrids",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@run_hypedsearch_on_scan_cli_options()
@log_params
def cli(
    protein_names: Path,
    fasta_path: Path,
    mzml: Path,
    output_dir: Path,
    num_peaks: int,
    precursor_mz_ppm_tol: float,
    min_k: int,
    max_k: int,
    num_psms: int,
    scan: Optional[int] = None,
):
    run_hypedsearch_on_spectrum(
        protein_names=protein_names,
        fasta_path=fasta_path,
        mzml=mzml,
        output_dir=output_dir,
        num_peaks=num_peaks,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        min_k=min_k,
        max_k=max_k,
        num_psms=num_psms,
        scan=scan,
    )


if __name__ == "__main__":
    setup_logger()
    cli()
