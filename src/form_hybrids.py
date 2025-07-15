import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from time import time
from typing import Callable, Dict, List, Literal, Optional

import click
import numpy as np

from src.constants import (
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    MZML,
    SCAN_HYBRIDS,
)
from src.hybrids_via_clusters import get_hybrids_via_clusters
from src.hypedsearch_utils import (
    HybridPeptide,
    create_hybrids_fasta,
    postprocess_hybrids,
)
from src.mass_spectra import Mzml, Spectrum, get_spectrum_from_mzml
from src.protein_product_ion_database import ProteinProductIonDb
from src.utils import (
    PathType,
    flatten_list_of_lists,
    get_time_in_diff_units,
    load_json,
    log_params,
    log_time,
    setup_logger,
    to_json,
)


@dataclass
class HybridFormingMethod:
    """
    API class for forming hybrids
    """

    method: Literal["all", "clustering"]
    form_hybrids: Callable
    info: Dict
    prot_id_to_name_map: Dict[int, str]

    @classmethod
    def from_hybrid_formation_relevant_params(
        cls,
        db_path: Optional[Path] = None,
        fasta: Optional[Path] = None,
        protein_names: Optional[Path] = None,
        precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
        min_k: int = DEFAULT_MIN_KMER_LEN,
        max_k: int = DEFAULT_MAX_KMER_LEN,
    ):

        if db_path is not None:
            hybrid_forming_fcn = lambda spectrum: get_hybrids_via_clusters(
                spectrum=spectrum,
                db_path=db_path,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            )

            info = {
                "method": "clustering",
                "db_path": db_path,
                "precursor_mz_ppm_tol": precursor_mz_ppm_tol,
                "peak_to_ion_ppm_tol": peak_to_ion_ppm_tol,
            }
            return HybridFormingMethod(
                method="clustering",
                info=info,
                form_hybrids=hybrid_forming_fcn,
                prot_id_to_name_map=ProteinProductIonDb(
                    db_path=db_path
                ).get_protein_id_to_name_map(),
            )


@log_time(level=logging.INFO)
def form_hybrids_for_spectrum(
    db_path: Path,
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
    fasta: Path,
    mzml: Optional[Path] = None,
    scan: Optional[int] = None,
    spectrum: Optional[Spectrum] = None,
    out_path: Optional[Path] = None,
) -> Dict[str, List[HybridPeptide]]:
    """
    Form hybrids for the spectrum.
    Also postprocess the hybrids to remove native sequences and add protein names.
    """
    if spectrum is None:
        spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
    hybrid_former = HybridFormingMethod.from_hybrid_formation_relevant_params(
        db_path=db_path,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
    )
    hybrids = hybrid_former.form_hybrids(spectrum=spectrum)
    seq_to_hybrid_peptides = postprocess_hybrids(
        hybrids=hybrids,
        fasta_path=fasta,
        remove_native=True,
        prot_id_to_name_map=hybrid_former.prot_id_to_name_map,
    )
    if out_path is not None:
        if out_path.is_dir():
            out_path = out_path / f"{spectrum.scan}.json"
        serialized_hy_peps = {
            key: [hy_pep.to_dict() for hy_pep in hy_peps]
            for key, hy_peps in seq_to_hybrid_peptides.items()
        }
        to_json(data=serialized_hy_peps, out_path=out_path)
    return seq_to_hybrid_peptides


@click.command(
    name="form-hybrids",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help="Form hybrids for the given spectra",
)
@click.option(
    "--mzml",
    "-m",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the MZML file.",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    required=False,
    help="Scan number of the spectrum to form hybrids for. If not provided, all scans will be processed.",
)
@click.option(
    "--db_path",
    "-d",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the product-ion database file.",
)
@click.option(
    "--precursor_mz_ppm_tol",
    "-pmpt",
    type=float,
    default=DEFAULT_PRECURSOR_MZ_PPM_TOL,
    show_default=True,
    help="Precursor m/z PPM tolerance. Hybrids will be within this PPM of the precursor m/z.",
)
@click.option(
    "--peak_to_ion_ppm_tol",
    "-pipt",
    type=float,
    default=DEFAULT_PEAK_TO_ION_PPM_TOL,
    show_default=True,
    help="The PPM tolerance within which a spectrum peak will match a fragment ion",
)
@click.option(
    "--fasta",
    "-f",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the FASTA file",
)
@click.option(
    "--out_dir",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Where hybrid .json files will be saved as '<scan>.json'",
)
@log_params
@log_time(level=logging.INFO)
def cli_form_hybrids(
    mzml: Path,
    db_path: Path,
    fasta: Path,
    out_dir: Path,
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
    scan: Optional[int],
):
    if scan is None:
        logger.info(f"Processing all scans in {mzml}")
        for spectrum in Spectrum.parse_ms2_from_mzml(mzml):
            logger.info(f"Processing scan {spectrum.scan} in {mzml}")
            form_hybrids_for_spectrum(
                db_path=db_path,
                fasta=fasta,
                out_path=out_dir,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
                spectrum=spectrum,
            )
    else:
        logger.info(f"Processing scan {scan} in {mzml}")
        form_hybrids_for_spectrum(
            db_path=db_path,
            fasta=fasta,
            out_path=out_dir,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            scan=scan,
            mzml=mzml,
        )


def load_hybrid_seqs_from_json(hybrid_json: Path) -> List[str]:
    seq_with_hybrids_dict = load_json(hybrid_json)
    return list(seq_with_hybrids_dict.keys())


def load_hybrid_peptides_from_json(hybrid_json: Path) -> List[HybridPeptide]:
    seq_with_hybrids_dict = load_json(hybrid_json)
    scan = int(hybrid_json.stem)
    hybrids = []
    for _, hybrid_list in seq_with_hybrids_dict.items():
        for hybrid_dict in hybrid_list:
            hybrid_dict["scan"] = scan
            hybrid_dict["sample"] = hybrid_json.parent.stem
            hybrids.append(HybridPeptide(**hybrid_dict))
    return hybrids


def combine_hybrids_from_scans(
    hybrids_dir: Path,
    sample: Optional[str] = None,
) -> List[HybridPeptide]:
    hybrids = []
    hybrid_jsons = list(hybrids_dir.glob("*.json"))

    for idx, hybrid_json in enumerate(hybrid_jsons):
        logger.debug(
            f"Loading hybrids from {hybrid_json} ({idx + 1} of {len(hybrid_jsons)})"
        )
        seq_with_hybrids_dict = load_json(hybrid_json)
        scan = int(hybrid_json.stem)
        for seq, hybrid_list in seq_with_hybrids_dict.items():
            for hybrid_dict in hybrid_list:
                hybrid_dict["scan"] = scan
                if sample is not None:
                    hybrid_dict["sample"] = sample
                else:
                    hybrid_dict["sample"] = hybrids_dir.stem
                hybrids.append(HybridPeptide(**hybrid_dict))

    return hybrids


def hybrids_to_fasta(
    hybrids_path: Path,
    new_fasta_path: Path,
    old_fasta: Optional[Path] = None,
):
    if hybrids_path.is_dir():
        logger.info(f"Loading hybrids from directory {hybrids_path}")
        hybrid_seqs = flatten_list_of_lists(
            [
                load_hybrid_seqs_from_json(hybrid_json=hybrid_json)
                for hybrid_json in hybrids_path.glob("*.json")
            ]
        )
    else:
        hybrid_seqs = load_hybrid_seqs_from_json(hybrid_json=hybrids_path)
    create_hybrids_fasta(
        hybrid_seqs=hybrid_seqs,
        new_fasta_path=new_fasta_path,
        old_fasta=old_fasta,
    )


@click.command(
    name="hybrids-to-fasta",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help="Convert hybrids from a directory to a FASTA file.",
)
@click.option(
    "--hybrids_path",
    "-hp",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--new_fasta_path",
    "-nfp",
    type=PathType(),
    required=True,
    help="Path to save the new FASTA file with hybrid sequences.",
)
@click.option(
    "--old_fasta",
    "-of",
    type=PathType(),
    required=False,
    help="Path to an old FASTA file to include existing proteins.",
)
def cli_hybrids_to_fasta(
    hybrids_path: Path,
    new_fasta_path: Path,
    old_fasta: Optional[Path] = None,
):
    hybrids_to_fasta(
        hybrids_path=hybrids_path,
        new_fasta_path=new_fasta_path,
        old_fasta=old_fasta,
    )


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


if __name__ == "__main__":
    logger = setup_logger()
    cli.add_command(cli_form_hybrids)
    cli.add_command(cli_hybrids_to_fasta)
    cli()
