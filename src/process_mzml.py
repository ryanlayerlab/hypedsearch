import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import click

from src.click_utils import PathType
from src.comet_utils import CometPSM, CometTxt
from src.constants import DEFAULT_PEAK_TO_ION_PPM_TOL, DEFAULT_PRECURSOR_MZ_PPM_TOL
from src.form_hybrids import form_hybrids_for_spectrum
from src.hypedsearch_utils import HybridPeptide, create_hybrids_fasta
from src.mass_spectra import Spectrum
from src.run_comet import CruxCometOutput, run_comet_on_one_mzml, run_comet_via_crux
from src.utils import get_default_comet_executable_path, setup_logger, to_json


def hybrids_json_name(scan: int) -> str:
    return f"hybrids_{scan}.json"


@dataclass
class Hybrids:
    seq_to_hybrids: Dict[str, List[HybridPeptide]]

    def to_json(self, out_dir: Path, scan: int):
        """
        Serialize the hybrids to a JSON file.
        """
        serialized_hy_peps = {
            key: [hy_pep.to_dict() for hy_pep in hy_peps]
            for key, hy_peps in self.seq_to_hybrids.items()
        }
        to_json(
            data=serialized_hy_peps,
            out_path=self.json_path(out_dir=out_dir, scan=scan),
        )

    def to_fasta(self, out_dir: Path, scan: int, old_fasta: Optional[Path] = None):
        create_hybrids_fasta(
            hybrid_seqs=list(self.seq_to_hybrids.keys()),
            new_fasta_path=self.fasta_path(
                out_dir=out_dir,
                scan=scan,
            ),
            old_fasta=old_fasta,
        )

    @staticmethod
    def json_path(out_dir: Path, scan: int) -> Path:
        return out_dir / f"hybrids_{scan}.json"

    @staticmethod
    def fasta_path(
        out_dir: Path,
        scan: int,
        #    old_fasta: Optional[Path] = None
    ) -> Path:
        # if old_fasta is not None:
        #     return out_dir / f"native_hybrids_{scan}.fasta"
        # else:
        return out_dir / f"hybrids_{scan}.fasta"


def process_spectrum(
    mzml: Path,
    template_comet_params: Path,
    spectrum: Spectrum,
    db_path: Path,
    fasta: Path,
    out_dir: Path,
    num_psms: int,
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    comet_exe: Path = get_default_comet_executable_path(),
) -> CruxCometOutput:
    # Form hybrids and save them
    hybrids = Hybrids(
        seq_to_hybrids=form_hybrids_for_spectrum(
            spectrum=spectrum,
            db_path=db_path,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            fasta=fasta,
        )
    )
    hybrids.to_json(out_dir=out_dir, scan=spectrum.scan)

    # Run Comet separate decoy search on native + hybrid FASTA
    # I think we want to include the original FASTA in the hybrids FASTA so that the decoy database
    # is larger
    hybrids.to_fasta(
        out_dir=out_dir,
        scan=spectrum.scan,
        old_fasta=fasta,  # Uncomment this line if you want to include the original FASTA in the hybrids FASTA
    )
    hybrid_comet_output = run_comet_on_one_mzml(
        template_comet_params=template_comet_params,
        num_psms=num_psms,
        precursor_mz_ppm_tolerance=precursor_mz_ppm_tol,
        out_dir=out_dir,
        fasta=hybrids.fasta_path(out_dir=out_dir, scan=spectrum.scan),
        mzml=mzml,
        comet_exe=comet_exe,
        min_scan=spectrum.scan,
        max_scan=spectrum.scan,
        decoy_search=2,
        stem="hybrids",
    )

    return hybrid_comet_output


def process_mzml(
    mzml: Path,
    template_comet_params: Path,
    num_psms: int,
    fasta: Path,
    out_dir: Path,
    db_path: Path,
    comet_exe: Path = get_default_comet_executable_path(),
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    min_scan: int = 0,
    max_scan: int = 0,
):
    # TODO: handle when len(native_psms) == 0,
    # Run Comet on the FASTA to get native PSMs

    # native_comet_output = run_comet_via_crux(
    #     mzml=mzml,
    #     fasta=fasta,
    #     out_dir=out_dir,
    #     decoy_search=0,
    #     num_psms=num_psms,
    #     min_scan=min_scan,
    #     max_scan=max_scan,
    # )

    native_comet_output = run_comet_on_one_mzml(
        template_comet_params=template_comet_params,
        num_psms=num_psms,
        precursor_mz_ppm_tolerance=precursor_mz_ppm_tol,
        out_dir=out_dir,
        fasta=fasta,
        mzml=mzml,
        comet_exe=comet_exe,
        min_scan=min_scan,
        max_scan=max_scan,
        stem="native",
        decoy_search=0,
    )
    scans_with_psms = list(
        set(psm.scan for psm in CometPSM.from_txt(txt_path=native_comet_output))
    )
    to_json(data=scans_with_psms, out_path=out_dir / "scans_with_psms.json")
    spectra = Spectrum.parse_ms2_from_mzml(spectra_file=mzml)

    for spectrum in spectra:
        if spectrum.scan not in scans_with_psms:
            continue

        _ = process_spectrum(
            template_comet_params=template_comet_params,
            mzml=mzml,
            comet_exe=comet_exe,
            spectrum=spectrum,
            db_path=db_path,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            fasta=fasta,
            out_dir=out_dir,
            num_psms=num_psms,
        )


def aggregate_best_target_psms_and_best_decoy_psms(mzml_dir: Path):
    target_psm_pattern = re.compile(r"^hybrids\.comet\.\d+-\d+\.txt$")
    target_psm_files = [
        f
        for f in mzml_dir.iterdir()
        if f.is_file() and target_psm_pattern.match(f.name)
    ]
    target_psms, decoy_psms = [], []
    for target_psm_file in target_psm_files:
        scan_num = int(target_psm_file.name.split(".")[2].split("-")[0])
        decoy_file = list(
            mzml_dir.glob(f"hybrids.comet.{scan_num}-{scan_num}.decoy.txt")
        )
        assert (
            len(decoy_file) == 1
        ), f"Expected exactly one decoy file for scan {scan_num}, found {len(decoy_file)}"

        # The PSMs in the file should be ordered by xcorr so just grab the first one
        target_psms.append(
            (scan_num, CometTxt(path=target_psm_file).get_first_psm_line())
        )
        decoy_psms.append((scan_num, CometTxt(path=decoy_file[0]).get_first_psm_line()))

    # Save the best PSMs to a TXT files
    # Sort PSMs by scan
    target_psms.sort(key=lambda x: x[0])
    decoy_psms.sort(key=lambda x: x[0])
    header = CometTxt(path=target_psm_file).get_header()
    # Save to TXT files
    with open(mzml_dir / "target.txt", "w") as f:
        f.write(f"{header}\n")
        for _, psm in target_psms:
            f.write(f"{psm}\n")
    with open(mzml_dir / "decoy.txt", "w") as f:
        f.write(f"{header}\n")
        for _, psm in decoy_psms:
            f.write(f"{psm}\n")


@click.command(
    name="process-mzml",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=(
        "Process the given MZML where 'process' means run Comet on the MZML and FASTA, "
        "form hybrids, and run Comet with a separate decoy run on the hybrids+FASTA."
    ),
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="Path to the MZML file.",
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--out_dir",
    "-o",
    type=PathType(),
    required=True,
)
@click.option(
    "--db_path",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the product-ion database file.",
)
@click.option(
    "--num_psms",
    "-n",
    default=1,
    show_default=True,
    help="Number of PSMs to keep per scan.",
)
def cli_process_mzml(
    mzml: Path,
    num_psms: int,
    fasta: Path,
    out_dir: Path,
    db_path: Path,
):
    """
    CLI entry point for processing MZML files.
    """
    process_mzml(
        mzml=mzml,
        num_psms=num_psms,
        fasta=fasta,
        out_dir=out_dir,
        db_path=db_path,
    )


if __name__ == "__main__":
    setup_logger()
    cli_process_mzml()
