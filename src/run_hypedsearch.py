import logging
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from time import time
from typing import List, Optional

import click

from src.comet_utils import run_comet_on_one_mzml
from src.constants import DEFAULT_NUM_PEAKS, DEFAULT_PPM_TOLERANCE, MOUSE_PROTEOME
from src.mass_spectra import Spectrum, get_spectrum_from_mzml
from src.peptide_spectrum_comparison import (  # get_possible_hybrids,
    ExtendedCluster,
    HybridPeptide,
    SpectrumExtendedClusters,
    create_hybrids_fasta,
    extend_clusters,
    get_clusters_from_ions,
    get_hybrids,
)
from src.protein_product_ion_database import (
    PositionedIon,
    ProteinProductIonDb,
    get_positions_in_proteins_of_peak_matching_ions,
    get_product_ions_matching_spectrum,
)
from src.utils import (
    PathType,
    get_time_in_diff_units,
    log_params,
    make_directory,
    setup_logger,
)

logger = logging.getLogger(__name__)


# @dataclass
# class HSIntermediates:
#     positioned_ions: List[PositionedIon]
#     b_ext_clusters: List[ExtendedCluster]
#     y_ext_clusters: List[ExtendedCluster]
#     hybrids: List[HybridPeptide]


# def get_potential_hybrids(
#     spectrum: Spectrum,
#     db_path: Path,
#     peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
#     precursor_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
# ):
#     # Connect to the database
#     db = ProteinProductIonDb(db_path=db_path, overwrite=False)

#     # Find ions that match spectrum peaks
#     peaks_with_matches = get_product_ions_matching_spectrum(
#         spectrum=spectrum,
#         db=db,
#         peak_product_ion_ppm_tolerance=peak_to_ion_ppm_tol,
#     )
#     positioned_ions = get_positions_in_proteins_of_peak_matching_ions(
#         peaks_with_matches=peaks_with_matches,
#         db=db,
#     )

#     # Get clusters
#     logger.info("Getting clusters...")
#     t0 = time()
#     clusters = get_clusters_from_ions(ions=positioned_ions)
#     logger.info(f"Getting clusters took {time() - t0:.2f} seconds")
#     logger.info("Extending clusters...")
#     t0 = time()
#     extended_clusters = extend_clusters(
#         spectrum_clusters=clusters, spectrum=spectrum, db=db
#     )
#     logger.info(f"Extending clusters took {time() - t0:.2f} seconds")

#     # Get possible hybrids
#     logger.info("Creating possible hybrids...")
#     t0 = time()
#     possible_hybrids = get_possible_hybrids(
#         extended_clusters=extended_clusters,
#         spectrum=spectrum,
#         precursor_mz_ppm_tolerance=precursor_ppm_tol,
#     )
#     logger.info(
#         f"Getting possible hybrids took {get_time_in_diff_units(time() - t0)}. There are {len(possible_hybrids)} possible hybrids."
#     )
#     return HSIntermediates(
#         positioned_ions=positioned_ions,
#         hybrids=possible_hybrids,
#         b_ext_clusters=extended_clusters.b,
#         y_ext_clusters=extended_clusters.y,
#     )


def run_hs_on_one_spectrum(
    db_path: Path,
    spectrum: Spectrum,
    output_dir: Path,
    num_peaks: int = DEFAULT_NUM_PEAKS,
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    precursor_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    fasta_path: Optional[Path] = None,
    keep_fasta: bool = False,
    overwrite: bool = False,
):
    # Check if run already exists. If it does, skip it.
    output_file = (
        output_dir
        / f"hs_{spectrum.mzml.stem}.{spectrum.scan_num}-{spectrum.scan_num}.txt"
    )
    if output_file.exists() and not overwrite:
        logger.info(
            f"Output file {output_file} already exists. Skipping spectrum {spectrum.scan_num} from MZML {spectrum.mzml.name}.\n\n"
        )
        return

    # Run HS on the spectrum
    fcn_start_time = time()
    logger.info(
        f"Running HS on spectrum {spectrum.scan_num} from MZML {spectrum.mzml.name}..."
    )

    # Perform peak filtering
    if num_peaks > 0:
        spectrum.filter_to_top_n_peaks(n=num_peaks)

    # Get hybrids
    results = get_hybrids(
        db_path=db_path,
        spectrum=spectrum,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        precursor_mz_ppm_tol=precursor_ppm_tol,
    )

    if len(results.hybrids) == 0:
        logger.info("No possible hybrids found. Exiting...\n\n")
        return

    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_path = Path(temp_dir).absolute()

        # Create new FASTA with hybrids
        logger.info("Creating new FASTA with hybrids...")
        if keep_fasta:
            new_fasta_path = (
                output_dir / f"{spectrum.mzml.stem}_{spectrum.scan_num}_hybrids.fasta"
            )
        else:
            new_fasta_path = (
                tmp_path / f"{spectrum.mzml.stem}_{spectrum.scan_num}_hybrids.fasta"
            )
        db = ProteinProductIonDb(db_path=db_path, overwrite=False)
        create_hybrids_fasta(
            hybrids=results.hybrids,
            fasta_path=new_fasta_path,
            db=db,
            other_prots=fasta_path,
        )

        # Run Comet
        logger.info("Running Comet using new FASTA...")
        output_stem = f"hs_{spectrum.mzml.stem}"
        comet_txt = None
        try:
            comet_txt = run_comet_on_one_mzml(
                fasta=new_fasta_path,
                mzml=spectrum.mzml,
                output_dir=tmp_path,
                scan=spectrum.scan_num,
                stem=output_stem,
            )
            shutil.copy2(comet_txt, output_dir)
        except RuntimeError as err:
            logger.info(f"Comet failed! Here's the error message:\n{err}")

    logger.info(
        f"Running HS on spectrum {spectrum.scan_num} from MZML {spectrum.mzml.name} took {get_time_in_diff_units(time() - fcn_start_time)}\n\n"
    )
    # return (results, comet_txt)


def run_hs_on_mzml(
    db_path: Path,
    mzml: Path,
    output_dir: Path,
    num_peaks: int = DEFAULT_NUM_PEAKS,
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    precursor_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    scan: Optional[int] = None,
    keep_fasta: bool = False,
    fasta_path: Optional[Path] = None,
    overwrite: bool = False,
):
    # Create the output directory if it doesn't exist
    make_directory(output_dir)

    # If scan number is passed, then run HS on that scan only
    if scan is not None:
        spectrum = get_spectrum_from_mzml(mzml_path=mzml, scan_num=scan)
        run_hs_on_one_spectrum(
            db_path=db_path,
            spectrum=spectrum,
            output_dir=output_dir,
            num_peaks=num_peaks,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            precursor_ppm_tol=precursor_ppm_tol,
            keep_fasta=keep_fasta,
            fasta_path=fasta_path,
            overwrite=overwrite,
        )

    # Otherwise, run HS on all scans
    else:
        # Get all the spectra from the mzml file
        spectra = Spectrum.from_mzml(mzml_path=mzml)

        # Run HS on each spectrum
        for spectrum in spectra:
            run_hs_on_one_spectrum(
                db_path=db_path,
                spectrum=spectrum,
                output_dir=output_dir,
                num_peaks=num_peaks,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
                precursor_ppm_tol=precursor_ppm_tol,
                keep_fasta=keep_fasta,
                fasta_path=fasta_path,
                overwrite=overwrite,
            )


@click.command(
    name="run-hs",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--db_path",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the already-created protein and product ion database.",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="Path to the mzML file.",
)
@click.option(
    "--output_dir",
    "-o",
    type=PathType(),
    required=True,
    help="Path to the output directory.",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    default=None,
    help="Scan number to run HS on. If not provided, HS will be run on all scans.",
)
@click.option(
    "--num_peaks",
    "-n",
    type=int,
    default=DEFAULT_NUM_PEAKS,
    show_default=True,
    help="Number of peaks to keep in the spectrum.",
)
@click.option(
    "--peak_to_ion_ppm_tol",
    "-p",
    type=float,
    default=DEFAULT_PPM_TOLERANCE,
    show_default=True,
    help="PPM tolerance for peak to ion matching.",
)
@click.option(
    "--precursor_ppm_tol",
    "-P",
    type=float,
    default=DEFAULT_PPM_TOLERANCE,
    show_default=True,
    help="PPM tolerance for precursor m/z matching.",
)
@click.option(
    "--keep_fasta",
    "-kf",
    type=bool,
    default=False,
    show_default=True,
    help="Whether or not to keep the FASTA file with the hybrids",
)
@click.option(
    "--fasta_path",
    "-fp",
    type=PathType(),
    help="If provided, the proteins in the given FASTA file will be included with the potential hybrid sequences in the FASTA file that Comet searches",
)
@click.option(
    "--overwrite",
    "-ow",
    type=bool,
    default=False,
    show_default=True,
    help="If True, checks if Hypedsearch output already exists and, if it does, does not run Hypedsearch",
)
@log_params
def run_hs_on_mzml_cli(
    db_path: Path,
    mzml: Path,
    output_dir: Path,
    num_peaks: int,
    peak_to_ion_ppm_tol: float,
    precursor_ppm_tol: float,
    keep_fasta: bool,
    scan: Optional[int],
    fasta_path: Optional[Path],
    overwrite: bool,
):
    t0 = time()
    logger.info("Starting to run Hypedsearch")
    run_hs_on_mzml(
        db_path=db_path,
        mzml=mzml,
        output_dir=output_dir,
        num_peaks=num_peaks,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        precursor_ppm_tol=precursor_ppm_tol,
        keep_fasta=keep_fasta,
        scan=scan,
        fasta_path=fasta_path,
        overwrite=overwrite,
    )
    logger.info(
        f"Running Hypedsearch in total took {get_time_in_diff_units(time() - t0)}"
    )


if __name__ == "__main__":
    setup_logger()
    run_hs_on_mzml_cli()
