import logging
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Optional

import click
from pydantic import BaseModel

from src.click_utils import PathType
from src.comet_utils import CometPSM
from src.constants import (
    COMET_PARAMS,
    CRUX,
    DEFAULT_COMET_PARAMS_FILE,
    DEFAULT_COMET_PRECURSOR_MZ_PPM_TOL,
    DEFAULT_COMET_SCAN_RANGE,
    DEFAULT_CRUX_PARAMS,
    DEFAULT_CRUX_PATH,
    DEFAULT_NUM_PSMS,
    MOUSE_PROTEOME,
)
from src.mass_spectra import Spectrum
from src.utils import (
    ExistingPath,
    get_default_comet_executable_path,
    log_params,
    log_time,
    setup_logger,
)

logger = logging.getLogger(__name__)


@dataclass
class CometParamsFile:
    """ """

    template: Path = DEFAULT_COMET_PARAMS_FILE
    num_psms: int = DEFAULT_NUM_PSMS
    precursor_mz_ppm_tol: float = DEFAULT_COMET_PRECURSOR_MZ_PPM_TOL
    decoy_search: Literal[0, 1, 2] = (
        0,
    )  # 0=no (default), 1=internal decoy concatenated, 2=internal decoy separate
    file_lines: List[str] = field(init=False)

    def __post_init__(self):
        # Check that template comet.params actually exists
        self.template = Path(self.template).absolute()
        assert (
            self.template.exists()
        ), "The path provided to comet.params file does not exist!"

        # Read the file and set self.file_lines
        with open(self.template, "r") as f:
            self.file_lines = f.readlines()

        self._update_params_obj()

    def _update_params_obj(self):
        self._update_num_output_psms_per_scan(num_output_psms=self.num_psms)
        self._update_precursor_mz_ppm_tol(tol=self.precursor_mz_ppm_tol)
        self._update_decoy_search(decoy_search=self.decoy_search)

    def _update_num_output_psms_per_scan(self, num_output_psms: int):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("num_output_lines ="):
                self.file_lines[line_idx] = f"num_output_lines = {num_output_psms}\n"
                break

    def _update_precursor_mz_ppm_tol(self, tol: float):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("peptide_mass_tolerance_upper ="):
                self.file_lines[line_idx] = f"peptide_mass_tolerance_upper = {tol}\n"
                break
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("peptide_mass_tolerance_lower ="):
                self.file_lines[line_idx] = f"peptide_mass_tolerance_lower = -{tol}\n"
                break

    def _update_decoy_search(self, decoy_search: Literal[0, 1, 2]):
        for line_idx, line in enumerate(self.file_lines):
            if line.strip().startswith("decoy_search ="):
                self.file_lines[line_idx] = f"decoy_search = {decoy_search}\n"
                break

    def write(self, output_path: Path):
        with open(output_path, "w") as f:
            f.writelines(self.file_lines)


class CometRunner(BaseModel):
    """
    Comet class to serve as API to Comet
    """

    comet_exe: ExistingPath  # path to comet executable
    comet_params: CometParamsFile  # path to base comet.params file
    fasta: ExistingPath
    out_dir: ExistingPath
    min_scan: int = DEFAULT_COMET_SCAN_RANGE[0]
    max_scan: int = DEFAULT_COMET_SCAN_RANGE[1]
    stem: str = ""

    def _comet_output_stem_without_scan_range(self, mzml: Path) -> str:
        if len(self.stem) == 0:
            return "comet"
        else:
            return f"{self.stem}.comet"

    def _comet_output_path(self, mzml: Path) -> Path:
        """
        Returns the path to the Comet output file for the given mzml file.
        """
        output_stem = self._comet_output_stem_without_scan_range(mzml=mzml)
        if (self.min_scan != 0) or (self.max_scan != 0):
            output_stem += f".{self.min_scan}-{self.max_scan}"
        return self.out_dir / f"{output_stem}.txt"

    def run_comet_on_mzml(self, mzml: Path, overwrite: bool = False) -> Path:
        """
        Returns path to the Comet txt output.
        """
        # If the Comet output already exists, only run Comet again if overwrite=True
        comet_output_path = self._comet_output_path(mzml=mzml)
        if comet_output_path.exists() and not overwrite:
            logger.info(
                f"Comet output files for {mzml} already exist in {self.out_dir} and "
                "overwrite is set to {overwrite} so skipping Comet run."
            )
        else:
            # Create the Comet params file in a temporary directory so that it doesn't persist
            # and then run Comet
            with tempfile.TemporaryDirectory() as temp_dir:
                comet_params_path = (Path(temp_dir) / "comet.params").absolute()
                self.comet_params.write(output_path=comet_params_path)

                # Run Comet
                cmd_parts = [
                    str(self.comet_exe),
                    f"-P{comet_params_path}",  # comet.params path
                    f"-F{self.min_scan} -L{self.max_scan}",  # scan range
                    f"-D{self.fasta}",  # fasta path
                    f"-N{self.out_dir / self._comet_output_stem_without_scan_range(mzml=mzml)}",
                    str(mzml),
                ]
                cmd = " ".join(cmd_parts)
                logger.info(f"Running Comet with command:\n{cmd}")
                cmd_result = subprocess.run(
                    cmd,
                    shell=True,
                    capture_output=True,
                    text=True,
                )

        return comet_output_path


@dataclass
class CruxCometOutput:
    min_scan: int
    max_scan: int
    file_root: str
    out_dir: Path
    decoy_search: Literal[0, 1, 2] = 0
    log_output: Path = field(init=False)
    params_output: Path = field(init=False)
    target_output: Path = field(init=False)
    decoy_output: Optional[Path] = field(init=False, default=None)

    def __post_init__(self):
        first_section = (
            f"{self.file_root}.comet." if len(self.file_root) > 0 else "comet."
        )
        self.log_output = self.out_dir / f"{first_section}log.txt"
        self.params_output = self.out_dir / f"{first_section}params.txt"

        # If min_scan and max_scan are both 0, we don't include the scan range in the output file name
        scan_section = f"{self.min_scan}-{self.max_scan}." if self.min_scan > 0 else ""

        # If decoy_search is 2, we add target and decoy to the output name
        if self.decoy_search == 2:
            self.target_output = (
                self.out_dir / f"{first_section}{scan_section}target.txt"
            )
            self.decoy_output = self.out_dir / f"{first_section}{scan_section}decoy.txt"
        else:
            self.target_output = self.out_dir / f"{first_section}{scan_section}txt"


def run_comet_via_crux(
    mzml: Path,
    fasta: Path,
    out_dir: Path,
    crux_path: Path = DEFAULT_CRUX_PATH,
    crux_params: Path = DEFAULT_CRUX_PARAMS,
    decoy_search: Literal[0, 1, 2] = 0,
    min_scan: int = 0,
    max_scan: int = 0,
    num_psms: int = DEFAULT_NUM_PSMS,
    file_root: str = "",
    overwrite: bool = False,
    keep_log_and_params: bool = False,
) -> CruxCometOutput:
    """
    Returns path to the Comet txt output.
    We want to run Comet in a temporary directory because it creates comet.log.txt and comet.params.txt
    and I don't want those to conflict with one another if Comet is being run in the same
    directory many times
    """
    assert ((max_scan > 0) and (min_scan > 0)) or ((max_scan == 0) and (min_scan == 0))

    # Run Comet
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir).absolute()
        file_root_str = f"--fileroot {file_root}" if file_root else ""
        cmd_parts = [
            str(crux_path),
            "comet",
            f"--parameter-file {crux_params}",
            f"--decoy_search {decoy_search}",
            f"--scan_range '{min_scan} {max_scan}'",
            f"--output-dir {temp_dir}",
            f"--num_output_lines {num_psms}",
            f"--overwrite {'T' if overwrite else 'F'}",
            file_root_str,
            f"{mzml} {fasta}",
        ]
        cmd = " ".join(cmd_parts)
        logger.info(f"Running Comet with command:\n{cmd}")
        cmd_result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        assert (
            cmd_result.returncode == 0
        ), f"Comet seems to have failed with this stderr:\n{cmd_result.stderr}"

        # Move the output file from the temporary directory to
        comet_output = CruxCometOutput(
            min_scan=min_scan,
            max_scan=max_scan,
            file_root=file_root,
            out_dir=temp_dir,
            decoy_search=decoy_search,
        )
        my_output = CruxCometOutput(
            min_scan=min_scan,
            max_scan=max_scan,
            file_root=file_root,
            out_dir=out_dir,
            decoy_search=decoy_search,
        )

        comet_output.target_output.rename(my_output.target_output)
        if decoy_search == 2:
            comet_output.decoy_output.rename(my_output.decoy_output)

    return my_output


def run_comet_on_one_spectrum(
    template_comet_params: Path,
    num_psms: int,
    precursor_mz_ppm_tolerance: float,
    decoy_search: Literal[0, 1, 2],
    out_dir: Path,
    fasta: Path,
    spectrum: Spectrum,
    stem: str,
    comet_exe: Path,
    mzml: Path,
    overwrite: bool,
) -> Path:
    assert spectrum.mzml.exists(), "The spectrum's mzML file does not exist!"
    return run_comet_on_one_mzml(
        template_comet_params=template_comet_params,
        num_psms=num_psms,
        precursor_mz_ppm_tolerance=precursor_mz_ppm_tolerance,
        decoy_search=decoy_search,
        out_dir=out_dir,
        fasta=fasta,
        min_scan=spectrum.scan,
        max_scan=spectrum.scan,
        stem=stem,
        comet_exe=comet_exe,
        mzml=mzml,
        overwrite=overwrite,
    )


def run_comet_on_one_mzml(
    mzml: Path,
    template_comet_params: Path,
    num_psms: int,
    precursor_mz_ppm_tolerance: float,
    out_dir: Path,
    fasta: Path,
    min_scan: int = 0,
    max_scan: int = 0,
    decoy_search: Literal[0, 1, 2] = 0,
    comet_exe: Path = get_default_comet_executable_path(),
    overwrite: bool = True,
    stem: str = "",
) -> Path:
    comet_params = CometParamsFile(
        template=template_comet_params,
        num_psms=num_psms,
        precursor_mz_ppm_tol=precursor_mz_ppm_tolerance,
        decoy_search=decoy_search,
    )
    comet_runner = CometRunner(
        comet_exe=comet_exe,
        comet_params=comet_params,
        fasta=fasta,
        out_dir=out_dir,
        stem=stem,
        min_scan=min_scan,
        max_scan=max_scan,
    )
    comet_output_path = comet_runner.run_comet_on_mzml(
        mzml=mzml,
        overwrite=overwrite,
    )
    return comet_output_path


def get_scans_comet_runs_on(
    out_dir: Path,
    fasta: Path,
    mzml: Path,
    min_scan: int = 0,
    max_scan: int = 0,
) -> List[int]:
    out_txt = run_comet_via_crux(
        mzml=mzml,
        fasta=fasta,
        out_dir=out_dir,
        min_scan=min_scan,
        max_scan=max_scan,
        num_psms=1,
        overwrite=True,
    )
    return list(set(psm.scan for psm in CometPSM.from_txt(txt_path=out_txt)))


@click.command(
    name="run-comet-on-mzml",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help=(
        "Path to .mzML file on which to run Comet. "
        "Or, if the path provided is to a directory, runs Comet on all the .mzML files found there"
    ),
)
@click.option(
    "--output_dir",
    "-o",
    type=PathType(),
    required=True,
    help="Path to directory where Comet outputs will be saved.",
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    default=MOUSE_PROTEOME,
    show_default=True,
    help="Path to FASTA file with which to run Comet",
)
@click.option(
    "--comet_exe",
    "-ce",
    type=PathType(),
    show_default=True,
    default=get_default_comet_executable_path(),
    help="Path to Comet executable",
)
@click.option(
    "--comet_params",
    "-cp",
    type=PathType(),
    show_default=True,
    default=COMET_PARAMS,
    help="Path to Comet parameters file, typically named comet.params",
)
@click.option(
    "--stem",
    "-n",
    type=str,
    help="The Comet output files will have this file stem. Defaults to the stem of the .mzML file",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    help=(
        "The scan number to run Comet on. "
        "Defaults to the scan_range in the provided Comet params file which is usually all scans in the .mzML file"
    ),
)
@click.option(
    "--num_psms",
    "-np",
    type=int,
    help=(
        "The number of PSMs per scan that Comet reports. "
        "Defaults to the num_output_lines in the provided Comet params file."
    ),
)
@click.option(
    "--precursor_mz_ppm_tol",
    "-pt",
    type=float,
    default=20,
    show_default=True,
    help=("The precursor mass tolerance in PPM."),
)
@click.option(
    "--overwrite",
    "-ow",
    is_flag=True,
    help="If Comet outputs already exist, this controls whether or not to overwrite them.",
)
@log_params
def cli_run_comet(
    output_dir: Path,
    mzml: Path,
    precursor_mz_ppm_tol: float,
    fasta: Path = MOUSE_PROTEOME,
    stem: Optional[str] = None,
    comet_exe: Path = get_default_comet_executable_path(),
    comet_params: Path = COMET_PARAMS,
    overwrite: bool = False,
    scan: Optional[int] = None,
    num_psms: Optional[int] = None,
):
    if mzml.is_dir():
        mzmls = list(mzml.glob("*.mzML"))
        logger.info(
            "You set the mzml argument to a directory. So running Comet on all .mzML files found there:\n"
            f"{mzmls}"
        )
        for mzml in mzmls:
            run_comet_on_one_mzml(
                mzml=mzml.absolute(),
                fasta=fasta.absolute(),
                comet_exe=comet_exe.absolute(),
                template_comet_params=comet_params.absolute(),
                out_dir=output_dir.absolute(),
                stem=stem,
                overwrite=overwrite,
                scan=scan,
                num_psms=num_psms,
                precursor_mz_ppm_tolerance=precursor_mz_ppm_tol,
            )

    else:
        run_comet_on_one_mzml(
            mzml=mzml.absolute(),
            fasta=fasta.absolute(),
            comet_exe=comet_exe.absolute(),
            template_comet_params=comet_params.absolute(),
            out_dir=output_dir.absolute(),
            stem=stem,
            overwrite=overwrite,
            scan=scan,
            num_psms=num_psms,
            precursor_mz_ppm_tolerance=precursor_mz_ppm_tol,
        )


if __name__ == "__main__":
    setup_logger()
    cli_run_comet()
