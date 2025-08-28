import os
import platform
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union
from venv import logger

import click
import yaml
from pydantic import BaseModel

from src.comet_utils import CometPSM
from src.constants import COMET_DIR, DECOY, RUN_COMET_SMK, TARGET
from src.kmer_database import create_kmer_database
from src.mass_spectra import Mzml
from src.protein_abundance import (
    get_most_common_proteins,
    get_protein_counts_from_comet_results,
)
from src.utils import (
    PathType,
    read_new_line_separated_file,
    save_dict,
    setup_logger,
    to_json,
    to_yaml,
    write_new_line_separated_file,
)

COMET_SMK = Path("snakefiles/run_comet.smk")


class CometConfig(BaseModel):
    mzml_to_scans: Dict[Path, List[int]]
    crux_comet_params: Path
    decoy_search: int
    fasta: Path
    out_dir: Path
    crux_path: Path

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "CometConfig":
        """
        Load a HypedsearchConfig from a YAML file.
        """
        with open(path, "r") as file:
            return cls(**yaml.safe_load(file))

    def to_dict(self) -> Dict:
        return self.model_dump(mode="json")

    def save(self, path: Union[str, Path]) -> None:
        data = self.to_dict()
        save_dict(data=data, path=path)

    def expected_outputs(
        self,
        out_dir: Union[str, Path],
        psm_type: Literal[TARGET, DECOY, "both"] = "both",
    ) -> List[str]:
        return get_expected_comet_outputs(
            mzml_to_scans=self.mzml_to_scans,
            out_dir=out_dir,
            decoy_search=self.decoy_search,
            psm_type=psm_type,
        )

    def run_comet_snakemake_command(self, config_path: Union[str, Path]) -> str:
        return f"snakemake -s {RUN_COMET_SMK} --configfile {config_path} ..."


def get_expected_comet_outputs(
    mzml_to_scans: Dict[str, List[int]],
    out_dir: Path,
    decoy_search: Literal[0, 1, 2],
    psm_type: Literal["both", TARGET, DECOY],
    remove_existing_files: bool = True,
) -> List[str]:
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

    if remove_existing_files:
        expected_outputs = [out for out in expected_outputs if not out.exists()]

    return [str(out) for out in expected_outputs]


class CometOutputs(BaseModel):
    params: Path
    log: Path
    target: Path
    decoy: Optional[Path] = None

    @classmethod
    def crux_comet_outputs(
        cls,
        out_dir: Path,
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        decoy_search: int = 0,
    ) -> "CometOutputs":
        """ """
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

    @classmethod
    def standardized_comet_outputs(
        cls,
        out_dir: Union[str, Path],
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        decoy_search: int = 0,
    ) -> "CometOutputs":
        """ """
        out_dir = Path(out_dir)
        if file_root != "":
            file_root = f"{file_root}."
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


@dataclass
class Crux:
    path: Path
    # env: Dict = field(init=False)

    def __post_init__(self):
        # Check that crux is available on the system and set self.env
        result = subprocess.run(
            f"{self.path} version",
            capture_output=True,
            text=True,
            shell=True,
        )
        assert (
            result.returncode == 0
        ), f"It appears `crux` is NOT available on the system {platform.system()} at path {self.path}."

    @staticmethod
    def validate_comet_output(result: subprocess.CompletedProcess):
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

    def run_comet(
        self,
        mzml: Union[str, Path],
        fasta: Union[str, Path],
        crux_comet_params: Union[str, Path],
        decoy_search: Literal[0, 1, 2],
        out_dir: Union[str, Path],
        file_root: str = "",
        scan_min: int = 0,
        scan_max: int = 0,
        num_threads: Optional[int] = None,
        # run_method: Literal["background", "foreground"],
    ) -> CometOutputs:
        # Check if expected outputs already exist and skip Comet run if they do
        expected_outputs = CometOutputs.standardized_comet_outputs(
            out_dir=Path(out_dir),
            file_root=file_root,
            scan_min=scan_min,
            scan_max=scan_max,
            decoy_search=decoy_search,
        )
        if expected_outputs.target.exists():
            logger.info(
                f"File {expected_outputs.target} already exists. Skipping `crux comet`."
            )
            return expected_outputs

        # Run Comet
        mzml = Path(mzml)
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            cmd_parts = [
                f"{self.path} comet",
                "--verbosity 60",
                f"--parameter-file {crux_comet_params}",
                f"--decoy_search {decoy_search}",
                f"--fileroot '{file_root}'",
                f"--scan_range '{scan_min} {scan_max}'",
                f"--output-dir {tmp_path}",
                f"--num_threads {num_threads}" if num_threads is not None else "",
                f"{mzml}",
                f"{fasta}",
            ]
            cmd = " ".join(cmd_parts)
            logger.info(f"Running Comet with command:\n{cmd}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                shell=True,
            )
            self.validate_comet_output(result=result)
            # Move the Comet outputs to out_dir
            tmp_outputs = CometOutputs.crux_comet_outputs(
                out_dir=tmp_path,
                file_root=file_root,
                scan_min=scan_min,
                scan_max=scan_max,
                decoy_search=decoy_search,
            )
            try:
                logger.info(f"Moving {tmp_outputs.target} to {expected_outputs.target}")
                shutil.move(tmp_outputs.target, expected_outputs.target)
                if tmp_outputs.decoy is not None:
                    shutil.move(tmp_outputs.decoy, expected_outputs.decoy)
            except FileNotFoundError:
                expected_outputs.target.touch()
                if tmp_outputs.decoy is not None:
                    expected_outputs.decoy.touch()
        return expected_outputs

    @staticmethod
    def combine_crux_comet_files(
        files: List[Union[str, Path]], out_path: Path
    ) -> List[str]:
        # Get header
        header = None
        for f in files:
            if f.stat().st_size != 0:
                header = read_new_line_separated_file(path=f)[0]
                break
        assert header is not None
        # Get body lines
        body_lines = []
        for f in files:
            lines = read_new_line_separated_file(path=f)
            body_lines.extend(lines[1:])
        # Write combined file
        file_lines = [header] + body_lines
        out_path.write_text("\n".join(file_lines))
        return file_lines

    def run_assign_confidence(self, target_txts: List[Path], out_path: Path) -> Path:
        logger.info("Running 'crux assign-confidence'...")
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            cmd_parts = [
                f"{self.path} assign-confidence",
                "--overwrite T",
                f"--output-dir {tmp_path}",
                "--list-of-files T",
            ]
            for txt in target_txts:
                cmd_parts.append(str(txt))
            cmd = " ".join(cmd_parts)
            # result = subprocess.run(
            #     cmd, capture_output=True, text=True, shell=True, env=self.env
            # )
            # assert result.returncode == 0
            logger.info(f"Running command:\n{cmd}")
            _ = subprocess.run(
                cmd,
                stdout=sys.stdout,
                stderr=sys.stderr,
                text=True,
                shell=True,
            )
            # Move `crux assign-conidence` to out_path
            tmp_file = tmp_path / "assign-confidence.target.txt"
            shutil.copy(tmp_file, out_path)
        logger.info("Finished running 'crux assign-confidence'")

    @staticmethod
    def native_comet_snakemake_config_and_cmd(
        mzml_to_scans: Dict[Path, List[int]],
        crux_comet_params: Path,
        decoy_search: int,
        fasta: Path,
        native_run_dir: Path,
        native_run_smk_config: Path,
        crux_path: Path,
    ) -> Tuple[str, CometConfig]:
        # Create config for snakemake
        config = CometConfig(
            mzml_to_scans=self.mzml_to_scans,
            crux_comet_params=self.crux_comet_params,
            decoy_search=2,
            fasta=self.fasta,
            out_dir=self.native_run_dir,
            crux_path=self.crux_path,
        )
        config.save(path=self.native_run_smk_config)

        # Print command to run snakemake
        cmd = f"snakemake -s {RUN_COMET_SMK} --configfile {self.native_run_smk_config} ..."
        logger.info(
            f"Saved native Comet run snakemake config to {self.native_run_smk_config}. To run snakemake, use this command:\n{cmd}"
        )
        return (cmd, config)


@dataclass
class HSConfig:
    name: str
    mzml_to_scans: Dict[Path, List[int]]
    crux_comet_params: Dict
    fasta: Path
    decoy_search: int
    parent_out_dir: Path
    top_n_proteins: int
    config_dir: Path = field(init=False)
    native_run_dir: Path = field(init=False)
    hybrid_run_dir: Path = field(init=False)
    hybrid_scan_result_dir: Path = field(init=False)
    native_run_config: Path = field(init=False)
    native_assign_confidence_txt: Path = field(init=False)
    top_prots_txt: Path = field(init=False)
    db_dir: Path = field(init=False)
    db_path: Path = field(init=False)
    kmer_to_prot_json: Path = field(init=False)

    def __post_init__(self):
        # Make directories
        self.name_dir = self.parent_out_dir / f"{self.name}"
        self.name_dir.mkdir(parents=True, exist_ok=True)
        self.config_dir = self.name_dir / "configs"
        self.config_dir.mkdir(parents=True, exist_ok=True)
        self.native_run_dir = self.name_dir / "native_run"
        self.native_run_dir.mkdir(parents=True, exist_ok=True)
        self.db_dir = self.name_dir / f"top{self.top_n_proteins}prots"
        self.hybrid_run_dir = self.db_dir / "hybrid_run"
        self.hybrid_run_dir.mkdir(parents=True, exist_ok=True)
        self.hybrid_scan_result_dir = self.hybrid_run_dir / f"scan_results"
        self.hybrid_scan_result_dir.mkdir(parents=True, exist_ok=True)
        self.native_run_config = self.config_dir / "native_run.yaml"
        self.native_assign_confidence_txt = (
            self.native_run_dir / "assign-confidence.txt"
        )
        self.top_prots_txt = self.db_dir / f"proteins.txt"
        self.db_path = self.db_dir / "kmers.db"
        self.kmer_to_prot_json = self.db_dir / "kmer_to_proteins.json"

    def create_native_run_comet_config(self) -> CometConfig:
        config = CometConfig(
            mzml_to_scans=self.mzml_to_scans,
            crux_comet_params=self.crux_comet_params,
            fasta=self.fasta,
            decoy_search=self.decoy_search,
            out_dir=self.native_run_dir,
        )
        config.to_yaml(self.native_run_config)

    def create_hybrid_run_config(self) -> CometConfig:
        pass

    def create_kmer_database(self):
        if self.db_path.exists():
            logger.info(f"Database {self.db_path} already exists. Skipping creation.")
        else:
            psms = CometPSM.from_txt(txt=self.native_assign_confidence_txt)
            prot_counts = get_protein_counts_from_comet_results(psms=psms)
            most_common_proteins = get_most_common_proteins(
                protein_counts=prot_counts, top_n=self.top_n_proteins
            )
            logger.info(f"Top {self.top_n_proteins} proteins:\n{most_common_proteins}")
            _ = write_new_line_separated_file(
                lines=most_common_proteins, path=self.top_prots_txt
            )
            _ = create_kmer_database(
                kmer_to_proteins_path=self.kmer_to_prot_json,
                fasta=self.fasta,
                proteins=most_common_proteins,
                db_path=self.db_path,
            )
