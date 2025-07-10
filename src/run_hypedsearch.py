import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from src.comet_utils import CometPSM
from src.constants import (
    DEFAULT_CRUX_PARAMS,
    DEFAULT_CRUX_PATH,
    DEFAULT_NUM_PSMS,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
)
from src.form_hybrids import form_hybrids_for_spectrum
from src.hypedsearch_utils import create_hybrids_fasta
from src.mass_spectra import Mzml, Spectrum
from src.run_comet import run_comet_via_crux
from src.utils import setup_logger, to_json


@dataclass
class HSOutput:
    hybrids_json: Path
    hybrids_native_fasta: Path
    hybrids_fasta: Path
    native_target: Path
    native_decoy: Path
    hybrids_native_target: Path
    hybrids_native_decoy: Path
    hybrids_target: Path

    @classmethod
    def from_out_dir(cls, out_dir: Path) -> "HSOutput":
        return cls(
            hybrids_json=out_dir / "hybrids.json",
            hybrids_native_fasta=out_dir / "hybrids_native.fasta",
            hybrids_fasta=out_dir / "hybrids.fasta",
            native_target=out_dir / "native_target.txt",
            native_decoy=out_dir / "native_decoy.txt",
            hybrids_native_target=out_dir / "hybrids_native_target.txt",
            hybrids_native_decoy=out_dir / "hybrids_native_decoy.txt",
            hybrids_target=out_dir / "hybrids_target.txt",
        )


def process_spectrum(
    spectrum: Spectrum,
    db_path: Path,
    fasta: Path,
    out_dir: Path,
    crux_path: Path = DEFAULT_CRUX_PATH,
    crux_params: Path = DEFAULT_CRUX_PARAMS,
    num_psms: int = DEFAULT_NUM_PSMS,
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
):
    # Constants
    hs_output = HSOutput.from_out_dir(out_dir)

    # Run Comet
    try:
        comet_output = run_comet_via_crux(
            mzml=spectrum.mzml,
            fasta=fasta,
            out_dir=out_dir,
            crux_path=crux_path,
            crux_params=crux_params,
            decoy_search=2,
            min_scan=spectrum.scan,
            max_scan=spectrum.scan,
            num_psms=num_psms,
        )
        # Rename output files
        comet_output.target_output.rename(hs_output.native_target)
        comet_output.decoy_output.rename(hs_output.native_decoy)

    # If the native search returns no PSMs, create empty files
    except AssertionError:
        return
    #     hs_output.native_target.touch()
    #     hs_output.native_decoy.touch()

    # Form hybrids
    seq_to_hybrids = form_hybrids_for_spectrum(
        spectrum=spectrum,
        db_path=db_path,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        fasta=fasta,
    )
    serialized_hy_peps = {
        key: [hy_pep.to_dict() for hy_pep in hy_peps]
        for key, hy_peps in seq_to_hybrids.items()
    }
    to_json(data=serialized_hy_peps, out_path=hs_output.hybrids_json)

    # Create FASTA containing hybrids. Add the hybrids to the proteome (as opposed to
    # running Comet just on the hybrids) so that the decoy search makes sense
    _ = create_hybrids_fasta(
        hybrid_seqs=list(seq_to_hybrids.keys()),
        new_fasta_path=hs_output.hybrids_native_fasta,
        old_fasta=fasta,
    )

    # Run Comet on hybrids+native FASTA
    comet_output = run_comet_via_crux(
        mzml=spectrum.mzml,
        fasta=hs_output.hybrids_native_fasta,
        out_dir=out_dir,
        crux_path=crux_path,
        crux_params=crux_params,
        decoy_search=2,
        min_scan=spectrum.scan,
        max_scan=spectrum.scan,
        num_psms=num_psms,
    )
    # Rename output files
    comet_output.target_output.rename(hs_output.hybrids_native_target)
    comet_output.decoy_output.rename(hs_output.hybrids_native_decoy)

    # Create FASTA containing just hybrids to look at xcorr of hybrids
    _ = create_hybrids_fasta(
        hybrid_seqs=list(seq_to_hybrids.keys()),
        new_fasta_path=hs_output.hybrids_fasta,
    )
    comet_output = run_comet_via_crux(
        mzml=spectrum.mzml,
        fasta=hs_output.hybrids_fasta,
        out_dir=out_dir,
        crux_path=crux_path,
        crux_params=crux_params,
        decoy_search=0,
        min_scan=spectrum.scan,
        max_scan=spectrum.scan,
        num_psms=num_psms,
    )
    # Rename output files
    comet_output.target_output.rename(hs_output.hybrids_target)

    # Delete new FASTA to save space
    hs_output.hybrids_native_fasta.unlink()
    hs_output.hybrids_fasta.unlink()


def process_mzml(
    mzml: Path,
    db_path: Path,
    fasta: Path,
    out_dir: Path,
    num_psms: int,
    testing: bool = False,
):
    # Make directory for MZML
    mzml_dir = out_dir / mzml.stem
    mzml_dir.mkdir(parents=True, exist_ok=True)
    for spectrum_idx, spectrum in enumerate(Mzml(path=mzml).ms2_spectra):
        scan_dir = mzml_dir / f"scan_{spectrum.scan}"
        scan_dir.mkdir(parents=True, exist_ok=True)
        process_spectrum(
            spectrum=spectrum,
            db_path=db_path,
            fasta=fasta,
            out_dir=scan_dir,
            num_psms=num_psms,
        )
        if testing and (spectrum_idx == 1):
            break
    # Output for snakemake
    (mzml_dir / "done.txt").touch()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process the spectra in an MZML file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mzml", "-m", type=Path, required=True, help="Path to the MZML file."
    )
    parser.add_argument(
        "--db_path",
        "-d",
        type=Path,
        required=True,
        help="Path to the product-ion database file.",
    )
    parser.add_argument(
        "--fasta", "-f", type=Path, required=True, help="Path to the FASTA file."
    )
    parser.add_argument(
        "--out_dir", "-o", type=Path, required=True, help="Where files will be saved."
    )
    parser.add_argument(
        "--num_psms",
        "-n",
        type=int,
        required=True,
        help="Number of Comet PSMs to return.",
    )
    parser.add_argument(
        "--testing",
        "-t",
        action="store_true",
        help="For testing, only process two spectra",
    )

    return parser.parse_args()


if __name__ == "__main__":
    setup_logger()
    args = parse_args()
    process_mzml(
        mzml=args.mzml,
        db_path=args.db_path,
        fasta=args.fasta,
        out_dir=args.out_dir,
        num_psms=args.num_psms,
        testing=args.testing,
    )
