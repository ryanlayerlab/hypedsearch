import logging
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import yaml
from pydantic import BaseModel

from src.constants import (
    DECOY,
    DEFAULT_NUM_COMET_RETRIES,
    DEFAULT_NUM_COMET_THREADS,
    LINUX_CRUX_EXECUTABLE,
    MAC_CRUX_EXECUTABLE,
    SINGULARITY_IMAGE,
    TARGET,
)
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide
from src.psm import CometPSM
from src.utils import (
    CmdLineResult,
    CmdLineRunner,
    PathType,
    load_json,
    move_file,
    read_new_line_separated_file,
    save_dict,
    setup_logger,
)

ASSIGN_CONFIDENCE_NAME = "assign-confidence.target.txt"
logger = logging.getLogger(__name__)


class CometOutputs(BaseModel):
    params: Path
    log: Path
    target: Path
    decoy: Optional[Path] = None

    @classmethod
    def crux_comet_outputs(
        cls,
        out_dir: Union[str, Path],
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        decoy_search: int = 0,
    ) -> "CometOutputs":
        """
        No matter the parameters, running Comet via crux always produces a params file
        and a log file, comet.params.txt and comet.log.txt, respectively.
        When decoy == 0 (no decoy search) or 1 (concatenated decoy search), there
        will only be a single output file: <file_root>.comet.[<scan_min>-<scan_max>.]txt.
        When decoy == 2, there will be a target output file and a decoy output file:
        <file_root>.comet.[<scan_min>-<scan_max>.]target.txt and
        <file_root>.comet.[<scan_min>-<scan_max>.]decoy.txt.
        The scan range part is not included if scan_min=scan_max=0.
        """
        out_dir = Path(out_dir)
        if file_root != "":
            file_root = f"{file_root}."
        scan_range = ""
        if (scan_min != 0) and (scan_max != 0):
            scan_range = f"{scan_min}-{scan_max}."
        if (decoy_search == 0) or (decoy_search == 1):
            return cls(
                target=out_dir / "".join([file_root, "comet.", scan_range, "txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )
        elif decoy_search == 2:
            return cls(
                target=out_dir
                / "".join([file_root, "comet.", scan_range, "target.txt"]),
                decoy=out_dir / "".join([file_root, "comet.", scan_range, "decoy.txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )

    @staticmethod
    def parse_standardized_comet_txt(comet_txt: Union[str, Path]) -> Tuple[str, str]:
        ouput_regex = r"^(?P<mzml>.+?)\.comet\.(?P<start>\d+)-(?P<end>\d+)\.(?P<psm_type>[^.]+)\.txt$"
        name = Path(comet_txt).name
        match = re.match(ouput_regex, name)
        if not match:
            raise ValueError(
                f"Trying to parse Comet TXT's name {name} using regex {ouput_regex}"
            )
        mzml, start, end, psm_type = match.groups()
        assert start == end, f"Expected start and end scan to be the same in {name}"
        assert psm_type in [
            TARGET,
            DECOY,
        ], f"Expected psm_type to be {TARGET} or {DECOY} in {name}"
        return mzml, start, psm_type

    @classmethod
    def standardized_comet_outputs(
        cls,
        out_dir: Union[str, Path],
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        decoy_search: int = 0,
    ) -> "CometOutputs":
        """
        For output file name consistency which will help snakemake check output files reliably,
        standardize the Comet output file names to always include the scan range and to
        include `.target.txt` for the target output files when decoy_search is 0 or 2
        """
        out_dir = Path(out_dir)
        if file_root != "":
            file_root = f"{file_root}."
        scan_range = f"{scan_min}-{scan_max}."
        if decoy_search == 0:
            return cls(
                target=out_dir
                / "".join([file_root, "comet.", scan_range, "target.txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )
        elif decoy_search == 1:
            return cls(
                target=out_dir / "".join([file_root, "comet.", scan_range, "txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )
        elif decoy_search == 2:
            return cls(
                target=out_dir
                / "".join([file_root, "comet.", scan_range, "target.txt"]),
                decoy=out_dir / "".join([file_root, "comet.", scan_range, "decoy.txt"]),
                log=out_dir / "".join([file_root, "comet.log.txt"]),
                params=out_dir / "".join([file_root, "comet.params.txt"]),
            )


class CometConfig(BaseModel):
    mzml_to_scans: Dict[Path, List[int]]
    crux_comet_params: Path
    decoy_search: int
    fasta: Path
    out_dir: Path
    num_threads: int = DEFAULT_NUM_COMET_THREADS

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "CometConfig":
        """
        Load a HypedsearchConfig from a YAML file.
        """
        with open(path, "r") as file:
            return cls(**yaml.safe_load(file))

    @classmethod
    def from_json(cls, path: Union[str, Path]) -> "CometConfig":
        return cls(**load_json(path=path))

    def to_dict(self) -> Dict:
        return self.model_dump(mode="json")

    def save(self, path: Union[str, Path]) -> None:
        data = self.to_dict()
        save_dict(data=data, path=path)

    def expected_outputs(
        self,
        psm_type: Literal[TARGET, DECOY, "both"] = "both",
    ) -> List[Path]:
        return get_expected_comet_outputs_for_mzml_to_scans(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=self.out_dir,
            decoy_search=self.decoy_search,
            psm_type=psm_type,
        )

    def missing_outputs(
        self, psm_type: Literal[TARGET, DECOY, "both"] = "both"
    ) -> List[Path]:
        expected_outputs = self.expected_outputs(psm_type=psm_type)
        missing_outputs = [Path(p) for p in expected_outputs if not Path(p).exists()]
        return missing_outputs

    def run_comet_on_mzml(
        self, mzml: Union[str, Path], scan_min: int = 0, scan_max: int = 0
    ) -> CometOutputs:
        return Crux().run_comet(
            mzml=mzml,
            fasta=self.fasta,
            crux_comet_params=self.crux_comet_params,
            decoy_search=self.decoy_search,
            num_threads=self.num_threads,
            out_dir=self.out_dir,
            file_root=Mzml.get_mzml_name(mzml=mzml),
            scan_min=scan_min,
            scan_max=scan_max,
        )


class CometRun(BaseModel):
    fasta: Path
    mzml: Path
    crux_comet_params: Path
    out_dir: Path
    decoy_search: Literal[0, 1, 2] = 0
    scan_min: int = 0
    scan_max: int = 0
    num_threads: Optional[int] = None

    @property
    def file_root(self) -> str:
        return Mzml(path=self.mzml).name

    @staticmethod
    def get_run_comet_command(
        crux_path: Union[str, Path],
        comet_run: "CometRun",
    ) -> str:
        cmd_parts = [
            f"{crux_path} comet",
            "--verbosity 60",
            f'--parameter-file "{comet_run.crux_comet_params}"',
            f"--decoy_search {comet_run.decoy_search}",
            f"--fileroot '{comet_run.file_root}'",
            f"--scan_range '{comet_run.scan_min} {comet_run.scan_max}'",
            f'--output-dir "{comet_run.out_dir}"',
            (
                f"--num_threads {comet_run.num_threads}"
                if comet_run.num_threads is not None
                else ""
            ),
            f'"{comet_run.mzml}"',
            f'"{comet_run.fasta}"',
        ]
        return " ".join(cmd_parts)

    @property
    def nonstandardized_comet_outputs(self) -> CometOutputs:
        return CometOutputs.crux_comet_outputs(
            out_dir=self.out_dir,
            file_root=self.file_root,
            scan_min=self.scan_min,
            scan_max=self.scan_max,
            decoy_search=self.decoy_search,
        )

    def run_comet_locally(self, crux_path: Union[str, Path]) -> CmdLineResult:
        logger.debug(f"Running Comet via crux located at {crux_path}")
        cmd_result = CmdLineRunner.run_cmd(
            cmd=self.get_run_comet_command(crux_path=crux_path, comet_run=self)
        )
        return cmd_result

    def run_comet_in_singularity(
        self,
        singularity_image: Union[str, Path],
        singularity_crux_path: str = "crux",
        singularity_num_threads: int = 1,
    ):
        logger.debug("Running Comet via Singularity")
        singularity_run = self.__class__(
            fasta=f"/data/{self.fasta.name}",
            mzml=f"/data/{self.mzml.name}",
            crux_comet_params="/data/crux.comet.params",
            out_dir="/outdir",
            decoy_search=self.decoy_search,
            scan_min=self.scan_min,
            scan_max=self.scan_max,
            num_threads=singularity_num_threads,
        )
        singularity_cmd_parts = [
            "singularity exec",
            f'--bind "{self.mzml}":"{singularity_run.mzml}"',
            f'--bind "{self.fasta}":"{singularity_run.fasta}"',
            f'--bind "{self.crux_comet_params}":"{singularity_run.crux_comet_params}"',
            f'--bind "{self.out_dir}":"{singularity_run.out_dir}"',
            f"{singularity_image}",
            self.get_run_comet_command(
                crux_path=singularity_crux_path, comet_run=singularity_run
            ),
        ]
        singularity_cmd = " ".join(singularity_cmd_parts)
        return CmdLineRunner.run_cmd(cmd=singularity_cmd)

    @property
    def standardized_comet_outputs(self):
        return CometOutputs.standardized_comet_outputs(
            out_dir=Path(self.out_dir),
            file_root=self.file_root,
            scan_min=self.scan_min,
            scan_max=self.scan_max,
            decoy_search=self.decoy_search,
        )

    def run_comet_and_keep_only_results(
        self,
        crux_path: Optional[str | Path] = None,
    ):
        on_singularity = False
        if crux_path is None:
            on_singularity = True
        # Run Comet. Run in a temporary directory so comet.log.txt and comet.params.txt are not kept
        if on_singularity:
            num_retries = DEFAULT_NUM_COMET_RETRIES
        else:
            num_retries = 0
        num_calls = 1 + num_retries
        for run_idx in range(num_calls):
            try:
                with tempfile.TemporaryDirectory() as tmp_dir:
                    tmp_path = Path(tmp_dir)
                    tmp_run = deepcopy(self)
                    tmp_run.out_dir = tmp_path
                    if on_singularity:
                        process = tmp_run.run_comet_in_singularity(
                            singularity_image=SINGULARITY_IMAGE,
                        )
                    else:
                        process = tmp_run.run_comet_locally(crux_path=crux_path)

                    # Sometimes Comet won't produce an output file when the spectrum
                    # that's being run on produces this Comet output:
                    # "- Load spectra: Warning - no spectra searched" with return code 1
                    if (
                        process.returncode == 1
                        and "no spectra searched" in process.stderr
                    ):
                        logger.warning(
                            "Running Comet finished with return code 1 and 'no spectra searched' in stderr. "
                            "This can happen when the spectrum you're trying to run Comet on doesn't "
                            "meet some criteria Comet is following. "
                            f"Creating the expected Comet outputs as empty files and continuing. Here's the STDERR and STDOUT:\nSTDERR:\n{process.stderr}\nSTDOUT:\n{process.stdout}"
                        )
                        self.standardized_comet_outputs.target.touch()  # create empty file to indicate that this was run and produced no PSMs
                        if self.standardized_comet_outputs.decoy:
                            self.standardized_comet_outputs.decoy.touch()
                        return process

                    # Move files from temp directory to final resting place
                    move_file(
                        src=tmp_run.nonstandardized_comet_outputs.target,
                        dest=self.standardized_comet_outputs.target,
                    )
                    if tmp_run.nonstandardized_comet_outputs.decoy:
                        move_file(
                            src=tmp_run.nonstandardized_comet_outputs.decoy,
                            dest=self.standardized_comet_outputs.decoy,
                        )
            except OSError as e:
                logger.warning(
                    f"Failed running Comet. This is try {run_idx + 1} of {num_calls}. Here's the run info:\n{self}.\nAnd here's the error:\n{e}"
                )
                continue
            return process


def run_comet(
    fasta: Path,
    mzml: Path,
    crux_comet_params: Path,
    out_dir: Path,
    decoy_search: Literal[0, 1, 2] = 0,
    scan_min: int = 0,
    scan_max: int = 0,
    num_threads: Optional[int] = None,
    crux_path: Optional[str | Path] = None,
) -> CometRun:
    run = CometRun(
        fasta=fasta,
        mzml=mzml,
        crux_comet_params=crux_comet_params,
        out_dir=out_dir,
        decoy_search=decoy_search,
        scan_min=scan_min,
        scan_max=scan_max,
        num_threads=num_threads,
    )
    run.run_comet_and_keep_only_results(crux_path=crux_path)
    return run


@dataclass
class Crux:
    # crux_path: Path = MAC_CRUX_EXECUTABLE
    # env: Dict = field(init=False)

    def check_if_crux_is_available(self):
        # Check that crux is available on the system and set self.env
        result = subprocess.run(
            f"{self.crux_path} version",
            capture_output=True,
            text=True,
            shell=True,
        )
        assert (
            result.returncode == 0
        ), f"It appears `crux` is NOT available on the system {platform.system()} at path {self.crux_path}."

    @cached_property
    def crux_path(self) -> str:
        if sys.platform == "darwin":
            return MAC_CRUX_EXECUTABLE
        elif sys.platform == "linux":
            return LINUX_CRUX_EXECUTABLE
        else:
            return None

    @staticmethod
    def combine_crux_comet_files(
        files: List[Union[str, Path]], out_path: Path
    ) -> Optional[List[str]]:
        # Get header
        header = None
        for f in files:
            if f.stat().st_size != 0:
                header = read_new_line_separated_file(path=f)[0]
                break
        if header is None:
            logger.info("Couldn't find a header in any of the files.")
            return None
        # Get body lines
        body_lines = []
        for f in files:
            lines = read_new_line_separated_file(path=f)
            body_lines.extend(lines[1:])
        # Write combined file
        file_lines = [header] + body_lines
        out_path.write_text("\n".join(file_lines))
        return file_lines

    def run_assign_confidence(
        self, target_txts: List[Path], out_path: Path, overwrite: bool = False
    ) -> CmdLineResult:
        if out_path.exists() and not overwrite:
            logger.info(
                f"{out_path} already exists and overwrite is set to False. Skipping running `crux assign-confidence`."
            )
            return
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            cmd_parts = [
                f"{self.crux_path} assign-confidence",
                # "--overwrite T",
                f"--output-dir {tmp_path}",
                "--list-of-files T",
            ]
            for txt in target_txts:
                cmd_parts.append(str(txt))
            cmd = " ".join(cmd_parts)
            result = CmdLineRunner().run_cmd(cmd=cmd)
            if result.returncode != 0:
                raise RuntimeError(
                    f"Running `crux assign-confidence` failed with return code {result.returncode}. Here's the STDERR and STDOUT:\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
                )
            # Move `crux assign-conidence` to out_path
            tmp_file = tmp_path / ASSIGN_CONFIDENCE_NAME
            shutil.copy(tmp_file, out_path)
        return result


def run_comet_on_custom_seqs(
    seqs: Union[Set[str], List[str]],
    spectra: List[Spectrum],
    comet_params: Union[str, Path],
    crux_path: Union[str, Path] = MAC_CRUX_EXECUTABLE,
) -> Dict[str, Optional[CometPSM]]:
    seqs = set(seqs)
    spectrum_to_psms = {}
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        fasta_path = tmp_path / "fasta.fasta"
        Fasta.write_fasta(
            peptides=[
                Peptide(seq=seq, name=f"seq{idx}") for idx, seq in enumerate(seqs)
            ],
            path=fasta_path,
        )
        for spectrum in spectra:
            run = CometRun(
                fasta=fasta_path,
                mzml=spectrum.mzml,
                crux_comet_params=comet_params,
                out_dir=tmp_path,
                decoy_search=0,
                scan_min=spectrum.scan,
                scan_max=spectrum.scan,
            )
            process = run.run_comet_and_keep_only_results(crux_path=crux_path)
            psms = CometPSM.from_txt(txt=run.standardized_comet_outputs.target)
            spectrum_to_psms[spectrum.uid] = psms[0] if len(psms) > 0 else None

    return spectrum_to_psms


def get_expected_comet_outputs_for_mzml_to_scans(
    mzml_to_scans: Dict[str, Set[int]],
    decoy_search: Literal[0, 1, 2],
    psm_type: Literal["both", TARGET, DECOY],
    out_dir: Union[str, Path],
) -> List[Path]:

    expected_outputs = {TARGET: [], DECOY: []}
    for mzml, scans in mzml_to_scans.items():
        for scan in scans:
            comet_outputs = CometOutputs.standardized_comet_outputs(
                out_dir=out_dir,
                decoy_search=decoy_search,
                file_root=Mzml.get_mzml_name(mzml=mzml),
                scan_min=scan,
                scan_max=scan,
            )
            expected_outputs[TARGET].append(comet_outputs.target)
            if comet_outputs.decoy is not None:
                expected_outputs[DECOY].append(comet_outputs.decoy)

    if psm_type == "both":
        expected_outputs = expected_outputs[TARGET] + expected_outputs[DECOY]
    else:
        expected_outputs = expected_outputs[psm_type]

    return expected_outputs


@click.command(
    "run-comet",
    help=("Run `crux comet`"),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--params",
    "-p",
    type=PathType(),
    required=True,
    help="Path to comet.params file",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=False,
    help="Path to MZML",
)
@click.option(
    "--out_dir",
    "-o",
    type=PathType(),
    required=False,
    help="Path to output directory.",
)
def cli_run_comet(mzml: Path, fasta: Path, params: Path, out_dir: Path):
    Crux().run_comet(
        mzml=mzml,
        fasta=fasta,
        crux_comet_params=params,
        out_dir=out_dir,
        file_root=Mzml.get_mzml_name(mzml=mzml),
    )


@click.command(
    "run-assign-confidence",
    help=("Run `crux assign-confidene` on all '*.target.txt' files in a directory."),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--in_dir",
    "-id",
    type=PathType(),
    required=True,
    help="Path to folder containing '*.target.txt' files.",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=False,
    help=(
        f"Path to folder where `crux assign-confidence` output will be saved as {ASSIGN_CONFIDENCE_NAME}. "
        + "Defaults to <in_dir>"
    ),
)
def cli_run_assign_confidence(in_dir: Path, out_dir: Optional[Path]):
    if out_dir is None:
        out_dir = in_dir
    Crux().run_assign_confidence(
        target_txts=list(in_dir.glob("*.target.txt")),
        out_path=out_dir / ASSIGN_CONFIDENCE_NAME,
    )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_run_comet)
    cli.add_command(cli_run_assign_confidence)
    cli()
