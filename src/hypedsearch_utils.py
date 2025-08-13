import logging
import os
import shlex
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import click
import yaml
from pydantic import BaseModel, field_validator

from src.constants import (
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    HS_PREFIX,
)
from src.hybrids_via_clusters import (
    form_spectrum_hybrids_via_clustering,
    serialize_hybrids,
)
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.utils import PathType, log_time, run_command_line_cmd, setup_logger, to_json

CRUX_BIN_PATH = Path("comet/crux-4.3.Darwin.x86_64/bin").absolute()
logger = logging.getLogger(__name__)


@dataclass
class FormHybridsConfig:
    # mzml_dir: Path
    # mzmls: List[Path]
    mzml_to_scans: Dict[Path, Union[str, List[int]]]
    database: Path
    results_dir: Path
    fasta: Path
    precursor_mz_ppm_tol: float
    peak_to_ion_ppm_tol: float

    def __post_init__(self):
        for mzml, scans in self.mzml_to_scans.items():
            assert Path(mzml).exists(), f"mzML file {mzml} does not exist."
            if isinstance(scans, str):
                assert scans == "all", "Only 'all' is allowed as a string for scans."
                self.mzml_to_scans[mzml] = Mzml(path=mzml).scans
        self.database = Path(self.database)
        self.results_dir = Path(self.results_dir)
        self.fasta = Path(self.fasta)

    def form_hybrids(self):
        self.results_dir.mkdir(parents=True, exist_ok=True)
        for mzml_path, scans in self.mzml_to_scans.items():
            mzml = Mzml(path=mzml_path)
            for scan in scans:
                spectrum = mzml.get_spectrum(scan=scan)
                seq_to_hybrids = form_spectrum_hybrids_via_clustering(
                    database=self.database,
                    fasta=self.fasta,
                    precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
                    peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
                    spectrum=spectrum,
                )
                # Save the hybrids to a JSON file
                to_json(
                    data=serialize_hybrids(seq_to_hybrids=seq_to_hybrids),
                    out_path=self.results_dir
                    / f"{Path(mzml_path).stem}_scan={spectrum.scan}.json",
                )

    @staticmethod
    def persistent_hypedsearch_outputs_for_spectrum(
        out_dir: Union[str, Path], sample: str, scan: int
    ) -> List[str]:
        out_dir = Path(out_dir)
        return [out_dir / f"{sample}_scan={scan}_native_psm.txt"]

    @staticmethod
    def temporary_hypedsearch_outputs_for_spectrum(
        sample: str, scan: int
    ) -> List[Path]:
        pass

    # @field_validator("mzmls", mode="before")
    # @classmethod
    # def set_mzml_paths(cls, value, info):
    #     mzml_dir = Path(info.data.get("mzml_dir"))
    #     return [mzml_dir / mzml for mzml in value]


@dataclass
class CometResultFiles:
    params: Path
    log: Path
    target: Path
    decoy: Optional[Path] = None


@dataclass
class CometRunner:
    env: Dict = field(init=False)

    def __post_init__(self):
        # Check that crux is available on the system and set self.env
        env = os.environ.copy()
        env["PATH"] = os.pathsep.join([str(CRUX_BIN_PATH), env["PATH"]])
        self.env = env
        cmd = "crux version"
        result = subprocess.run(
            cmd, capture_output=True, text=True, shell=True, env=self.env
        )
        assert (
            result.returncode == 0
        ), "It appears `crux` is NOT available on the system."

    @staticmethod
    def comet_outputs(
        out_dir: Path,
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        decoy_search: int = 0,
    ) -> CometResultFiles:
        if file_root != "":
            file_root = f"{file_root}."
        scan_range = ""
        if (scan_min != 0) and (scan_max != 0):
            scan_range = f"{scan_min}-{scan_max}."
        if (decoy_search == 0) or (decoy_search == 1):
            return CometResultFiles(
                target=out_dir / "".join([file_root, "comet.", scan_range, "txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )
        elif decoy_search == 2:
            return CometResultFiles(
                target=out_dir
                / "".join([file_root, "comet.", scan_range, "target.txt"]),
                decoy=out_dir / "".join([file_root, "comet.", scan_range, "decoy.txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )

    def run_comet(
        self,
        mzml: Path,
        fasta: Path,
        crux_comet_params: Path,
        decoy_search: int = 0,
        out_dir: Path = "./",
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        spectrum_batch_size: int = 500,
    ):
        cmd_parts = [
            "crux comet",
            "--num_threads 1",
            "--verbosity 60",
            "--overwrite T",
            f"--spectrum_batch_size {spectrum_batch_size}",
            f"--parameter-file {crux_comet_params}",
            f"--decoy_search {decoy_search}",
            f"--fileroot '{file_root}'",
            f"--scan_range '{scan_min} {scan_max}'",
            f"--output-dir {out_dir}",
            f"{mzml}",
            f"{fasta}",
        ]
        cmd = " ".join(cmd_parts)
        run_args = shlex.split(cmd)
        result = subprocess.run(
            cmd, capture_output=True, text=True, shell=True, env=self.env
        )
        # Validate the response from Comet
        if result.returncode != 0:
            if (result.returncode == 1) and "no spectra searched" in result.stderr:
                logger.info(
                    f"Warning: `crux comet` finished with return code 1 and 'no spectra searched' in stderr. "
                    "This can happen and generally is not an error even though the return code is 1"
                )
            else:
                raise RuntimeError(
                    f"`crux comet` failed with return code {result.returncode}.\n"
                    f"STDOUT: {result.stdout}\n"
                    f"STDERR: {result.stderr}"
                )


# class HypedSearchRunner:
#     spectrum: Spectrum


def hybrid_fasta_name(hybrid_seq: str) -> str:
    return f"{HS_PREFIX}{hybrid_seq}"


def create_hybrids_fasta(
    hybrid_seqs: List[str],
    new_fasta_path: Path,
    old_fasta: Optional[Path] = None,
    protein_names: Optional[Union[List[str], Path]] = None,
) -> List[Peptide]:
    """
    Writes a list of hybrids (hybrids) to a FASTA file (new_fasta_path).
    If old_fasta is provided, it will also include the proteins from the old FASTA.
    If protein_names is provided, it will only include those proteins from the old FASTA.
    """
    prots = []
    if old_fasta is not None:
        if protein_names is not None:
            # Get specific proteins from FASTA by name
            prots = get_proteins_by_name(
                protein_names=protein_names, fasta_path=old_fasta
            )
        else:
            # Get all the proteins in the FASTA
            prots = Peptide.from_fasta(fasta_path=old_fasta)

    for _, hybrid_seq in enumerate(hybrid_seqs):
        new_peptide = Peptide(
            seq=hybrid_seq,
            name=hybrid_fasta_name(hybrid_seq=hybrid_seq),
        )

        prots.append(new_peptide)

    Fasta.write_fasta(peptides=prots, out_path=new_fasta_path)
    return prots


@dataclass
class HypedsearchOutput:
    native_target: Path
    native_decoy: Path
    hybrid_target: Path
    hybrid_decoy: Path
    hybrids: Path

    @classmethod
    def get_hypedsearch_outputs(
        cls, out_dir: Path, mzml: Path, scan: int
    ) -> "HypedsearchOutput":
        mzml_name = f"{mzml.name[:-5]}"  # remove .mzML extension
        native_file_root = f"native_{mzml_name}"
        hybrid_file_root = f"hybrid_{mzml_name}"
        native_outputs = CometRunner.comet_outputs(
            out_dir=out_dir,
            file_root=native_file_root,
            scan_min=scan,
            scan_max=scan,
            decoy_search=2,
        )
        hybrid_outputs = CometRunner.comet_outputs(
            out_dir=out_dir,
            file_root=hybrid_file_root,
            scan_min=scan,
            scan_max=scan,
            decoy_search=2,
        )
        hybrids_json = out_dir / f"{mzml_name}_scan={scan}.json"
        return cls(
            native_target=native_outputs.target,
            native_decoy=native_outputs.decoy,
            hybrid_target=hybrid_outputs.target,
            hybrid_decoy=hybrid_outputs.decoy,
            hybrids=hybrids_json,
        )


@dataclass
class HypedsearchRunner:

    @staticmethod
    def run_hypedsearch(
        spectrum: Spectrum,
        database: Path,
        fasta: Path,
        crux_comet_params: Path,
        out_dir: Path,
        precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> HypedsearchOutput:
        comet = CometRunner()
        hs_outputs = HypedsearchOutput.get_hypedsearch_outputs(
            out_dir=out_dir, mzml=spectrum.mzml, scan=spectrum.scan
        )
        # Native Comet run
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            comet.run_comet(
                mzml=spectrum.mzml,
                fasta=fasta,
                crux_comet_params=crux_comet_params,
                decoy_search=2,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
                out_dir=tmp_path,
            )
            comet_output_files = CometRunner.comet_outputs(
                out_dir=tmp_path,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
                decoy_search=2,
            )

            # If the PSM files weren't created, create empty files. This is for snakemake
            if comet_output_files.target.exists() is False:
                comet_output_files.target.touch()
                comet_output_files.decoy.touch()

            # Move target and decoy files from temporary directory to output directory
            comet_output_files.target.rename(hs_outputs.native_target)
            comet_output_files.decoy.rename(hs_outputs.native_decoy)

        # Form hybrids
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            seq_to_hybrids = form_spectrum_hybrids_via_clustering(
                database=database,
                fasta=fasta,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
                spectrum=spectrum,
            )
            to_json(
                data=serialize_hybrids(seq_to_hybrids=seq_to_hybrids),
                out_path=hs_outputs.hybrids,
            )

            # Create hybrids+native FASTA in temporary directory
            hybrid_fasta = tmp_path / f"{spectrum.mzml.stem}_scan={spectrum.scan}.fasta"
            create_hybrids_fasta(
                hybrid_seqs=list(seq_to_hybrids.keys()),
                new_fasta_path=hybrid_fasta,
                old_fasta=fasta,
            )

            # Run Comet on hybrids FASTA
            comet.run_comet(
                mzml=spectrum.mzml,
                fasta=hybrid_fasta,
                crux_comet_params=crux_comet_params,
                decoy_search=2,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
                out_dir=tmp_path,
            )
            comet_output_files = CometRunner.comet_outputs(
                out_dir=tmp_path,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
                decoy_search=2,
            )

            # If the PSM files weren't created, create empty files. This is for snakemake
            if comet_output_files.target.exists() is False:
                comet_output_files.target.touch()
                comet_output_files.decoy.touch()

            # Move target and decoy files from temporary directory to output directory
            comet_output_files.target.rename(hs_outputs.hybrid_target)
            comet_output_files.decoy.rename(hs_outputs.hybrid_decoy)

        return hs_outputs


@click.command(
    name="run-hypedsearch",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=("Run Hypedsearch on a spectrum"),
)
@click.option(
    "--database",
    "-db",
    type=PathType(),
    required=True,
    help="Path to the database.",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="Path to the mzML file.",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    required=True,
    help="Scan number in mzML of spectrum to run Hypedsearch on",
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--crux_comet_params",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the crux.comet.params file to use with `crux comet`",
)
@click.option(
    "--out_dir",
    "-o",
    type=PathType(),
    required=True,
    help="Output directory for Hypedsearch results.",
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
@log_time(level=logging.INFO)
def cli_run_hypedsearch(
    mzml: Path,
    scan: int,
    database: Path,
    fasta: Path,
    crux_comet_params: Path,
    out_dir: Path,
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
):
    spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
    HypedsearchRunner.run_hypedsearch(
        spectrum=spectrum,
        database=database,
        fasta=fasta,
        crux_comet_params=crux_comet_params,
        out_dir=out_dir,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
    )


if __name__ == "__main__":
    setup_logger()
    cli_run_hypedsearch()
