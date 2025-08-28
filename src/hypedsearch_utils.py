import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union

import click
import yaml
from pydantic import BaseModel, field_validator

from src.comet_utils import CometPSM, CometPSMs
from src.constants import (
    BACKGROUND,
    DECOY,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    FOREGROUND,
    HS_PREFIX,
    HYBRID,
    MZML_EXT,
    NATIVE,
    TARGET,
)
from src.hybrids_via_clusters import (
    form_spectrum_hybrids_via_clustering,
    serialize_hybrids,
)
from src.kmer_database import create_kmer_database
from src.mass_spectra import Mzml, Spectrum
from src.peptides_and_ions import Fasta, Peptide, get_proteins_by_name
from src.protein_abundance import (
    get_most_common_proteins,
    get_protein_counts_from_comet_results,
)
from src.utils import (
    PathType,
    load_json,
    log_time,
    read_new_line_separated_file,
    run_command_line_cmd,
    setup_logger,
    to_json,
    write_new_line_separated_file,
)

CRUX_BIN_PATH = Path("comet/crux-4.3.Darwin.x86_64/bin").absolute()
logger = logging.getLogger(__name__)


@dataclass
class Hypedsearch:
    name: str
    mzmls: List[Path]
    parent_out_dir: Path
    fasta: Path
    crux_comet_params: Path
    precursor_mz_ppm_tol: float
    peak_to_ion_ppm_tol: float

    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> "HypedsearchConfig":
        """
        Load a HypedsearchConfig from a YAML file.
        """
        with open(yaml_path, "r") as file:
            config_dict = yaml.safe_load(file)

        # For those mzML files that have "all" scans, we will replace it with the actual scans
        mzmls_to_scans = {}
        key_name = "mzml_to_scans"
        for mzml, scans in config_dict[key_name].items():
            assert Path(mzml).exists(), f"mzML file {mzml} does not exist."
            if isinstance(scans, str):
                assert scans == "all", "Only 'all' is allowed as a string for scans."
                scans = Mzml(mzml=mzml).scans
            mzmls_to_scans[mzml] = scans
        config_dict[key_name] = mzmls_to_scans
        return cls(**config_dict)

    def to_yaml(self, path: Path) -> None:
        data = self.model_dump(mode="json")  # Converts Paths to strings
        with open(path, "w") as f:
            yaml.safe_dump(
                data, f, sort_keys=False
            )  # preserve field order so mzml_to_scans is last

    # def form_hybrids(self):
    #     self.results_dir.mkdir(parents=True, exist_ok=True)
    #     for mzml_path, scans in self.mzml_to_scans.items():
    #         mzml = Mzml(path=mzml_path)
    #         for scan in scans:
    #             spectrum = mzml.get_spectrum(scan=scan)
    #             seq_to_hybrids = form_spectrum_hybrids_via_clustering(
    #                 database=self.database,
    #                 fasta=self.fasta,
    #                 precursor_mz_ppm_tol=self.precursor_mz_ppm_tol,
    #                 peak_to_ion_ppm_tol=self.peak_to_ion_ppm_tol,
    #                 spectrum=spectrum,
    #             )
    #             # Save the hybrids to a JSON file
    #             to_json(
    #                 data=serialize_hybrids(seq_to_hybrids=seq_to_hybrids),
    #                 out_path=self.results_dir
    #                 / f"{Path(mzml_path).stem}_scan={spectrum.scan}.json",
    #             )


@click.command(
    name="create-hs-config",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=("Run Hypedsearch on a spectrum"),
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="Path to the directory containing mzML files",
)
@click.option(
    "--database",
    "-db",
    type=PathType(),
    required=True,
    help="Path to the database.",
)
@click.option(
    "--mzml_dir",
    "-m",
    type=PathType(),
    required=True,
    help="Path to the mzML file.",
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
    "-od",
    type=PathType(),
    required=True,
    help="Output directory for Hypedsearch results.",
)
@click.option(
    "--out_config",
    "-oc",
    type=PathType(),
    required=True,
    help="Where to save the Hypedsearch config file",
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
def cli_create_hypdedsearch_config(
    mzml_dir: Path,
    database: Path,
    fasta: Path,
    crux_comet_params: Path,
    out_dir: Path,
    out_config: Path,
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
):
    mzml_to_scans = {}
    for mzml in mzml_dir.glob("*.mzML"):
        mzml_to_scans[str(mzml)] = "all"
    config = HypedsearchConfig(
        mzml_to_scans=mzml_to_scans,
        database=database,
        out_dir=out_dir,
        fasta=fasta,
        crux_comet_params=crux_comet_params,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
    )
    config.to_yaml(path=out_config)


@dataclass
class CometResultFiles:
    params: Path
    log: Path
    target: Path
    decoy: Optional[Path] = None


def hybrid_fasta_name(hybrid_seq: str) -> str:
    return f"{HS_PREFIX}{hybrid_seq}"


def create_hybrids_fasta(
    hybrid_seqs: List[str],
    output_fasta_path: Path,
    fasta_to_include: Optional[Path] = None,
    protein_names: Optional[Union[List[str], Path]] = None,
) -> List[Peptide]:
    """
    Writes a list of hybrids (hybrids) to a FASTA file (new_fasta_path).
    If old_fasta is provided, it will also include the proteins from the old FASTA.
    If protein_names is provided, it will only include those proteins from the old FASTA.
    """
    prots = []
    if fasta_to_include is not None:
        if protein_names is not None:
            # Get specific proteins from FASTA by name
            prots = get_proteins_by_name(
                protein_names=protein_names, fasta_path=fasta_to_include
            )
        else:
            # Get all the proteins in the FASTA
            prots = Peptide.from_fasta(fasta_path=fasta_to_include)

    for _, hybrid_seq in enumerate(hybrid_seqs):
        new_peptide = Peptide(
            seq=hybrid_seq,
            name=hybrid_fasta_name(hybrid_seq=hybrid_seq),
        )

        prots.append(new_peptide)

    Fasta.write_fasta(peptides=prots, out_path=output_fasta_path)
    return prots


@dataclass
class SnakemakeRunner:

    @staticmethod
    def run_snakemake(
        smk_file: Path,
        run_config: Path,
        num_cores: int,
        expected_outputs: List[str] = [],
        keep_going: bool = True,
        run_method: Literal["foreground", "background"] = "background",
        log_file: Optional[str] = None,
    ):
        cmd_parts = [
            "snakemake",
            "-s",
            f"{smk_file}",
            "--configfile",
            f"{run_config}",
            "--cores",
            f"{num_cores}",
        ]
        if keep_going:
            cmd_parts.append("--keep-going")
        cmd_parts.extend(expected_outputs)
        if log_file is not None:
            cmd_parts.append(f"&> {log_file}")
        cmd = " ".join(cmd_parts)
        logger.info(f"Running Snakemake with this command:\n{cmd}")
        if run_method == "foreground":
            logger.info("Running in foreground")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                shell=True,
                env=Crux().env,
            )
            assert result.returncode == 0, "Appears snakemake run failed"

        elif run_method == "background":
            logger.info("Running in background")
            subprocess.Popen(
                cmd,
                # cmd_parts,
                stdout=sys.stdout,
                stderr=sys.stderr,
                start_new_session=True,  # Detach from parent session
                shell=True,
            )
        else:
            raise ValueError(f"Invalid run_method: {run_method}")


@dataclass
class MzmlSet:
    name: str
    mzmls: List[Path]

    @cached_property
    def num_scans(self) -> int:
        num_scans = 0
        for mzml in self.mzmls:
            num_scans += len(Mzml(mzml=mzml).scans)
        return num_scans

    def get_psms_by_type(
        self,
        out_dir: Path,
        native_or_hybrid: Literal[NATIVE, HYBRID],
        comet_or_assign_confidence: Literal["comet", "assign-confidence"],
        psm_type: Literal[TARGET, DECOY] = TARGET,
    ) -> CometPSMs:
        psms = []
        for mzml in self.mzmls:
            mzml_name = Mzml.get_mzml_name(mzml)
            if comet_or_assign_confidence == "comet":
                txt_path = out_dir / HypedsearchOutput.comet_txt_name_for_mzml(
                    mzml_name=mzml_name,
                    native_or_hybrid=native_or_hybrid,
                    psm_type=psm_type,
                )
            elif comet_or_assign_confidence == "assign-confidence":
                txt_path = out_dir / HypedsearchOutput.assign_confidence_txt_name(
                    mzml_set_name=self.name, native_or_hybrid=native_or_hybrid
                )
            else:
                raise ValueError(
                    f"Invalid 'comet_or_assign_confidence' value: {comet_or_assign_confidence}. Must be 'comet' or 'assign-confidence'"
                )
            psms.extend(CometPSM.from_txt(txt=txt_path))
        return CometPSMs(psms=psms)

    # def hybrid_target_psms(self, out_dir: Path) -> List[CometPSM]:

    #         assert len(matches) == 1, f"{mzml_name}: {len(matches)}"
    #         yield CometPSM.parse_from_file(matches[0])
    # def hybrid_target_psms(self, out_dir: Path) -> List[CometPSM]:


@dataclass
class HypedsearchConfig:
    name: str
    mzmls: List[Path]
    parent_out_dir: Path
    fasta: Path
    crux_comet_params: Path
    precursor_mz_ppm_tol: float
    peak_to_ion_ppm_tol: float

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "HypedsearchConfig":
        """
        Load a HypedsearchConfig from a YAML file.
        """
        with open(path, "r") as file:
            return cls(**yaml.safe_load(file))


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
    "--kmer_to_protein_map",
    "-kpm",
    type=PathType(),
    required=True,
    help="kmer-to-protein map.",
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
    kmer_to_protein_map: Path,
):
    setup_logger()
    spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
    HypedsearchRunner.run_hypedsearch(
        spectrum=spectrum,
        database=database,
        fasta=fasta,
        crux_comet_params=crux_comet_params,
        out_dir=out_dir,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        kmer_to_protein_map=kmer_to_protein_map,
    )


def get_missing_hypedsearch_outputs(
    results_dir: Path,
    config: Path,
) -> List[Path]:
    hs_config = HypedsearchConfig.from_yaml(yaml_path=config)
    # Get expected outputs
    expected_outputs = []
    for mzml, scans in hs_config.mzml_to_scans.items():
        for scan in scans:
            outputs = HypedsearchOutput.get_hypedsearch_outputs_for_scan(
                out_dir=results_dir, mzml=mzml, scan=scan
            )
            expected_outputs.extend(
                [
                    outputs.native_target,
                    outputs.native_decoy,
                    outputs.hybrid_target,
                    outputs.hybrid_decoy,
                ]
            )

    # Check which expected outputs are extant and missing
    missing_outputs = []
    for expected_output in expected_outputs:
        if not expected_output.exists():
            missing_outputs.append(expected_output)
    logger.info(
        f"Config: {config.name}\n \t - Expected outputs: {len(expected_outputs)}\n \t - Missing outputs: {len(missing_outputs)}"
    )
    return missing_outputs


def create_config_for_missing_scans(
    results_dir: Path, config: Path, out_path: Optional[Path] = None
) -> HypedsearchConfig:

    # Get missing output file paths
    missing_outputs = get_missing_hypedsearch_outputs(
        results_dir=results_dir, config=config
    )
    # Group missing outputs by MZML: {<mzml name>: [<mising scans>]}
    mzml_names_to_scans = defaultdict(set)
    for out in missing_outputs:
        mzml_stem = HypedsearchOutput.get_info_from_hypedsearch_comet_output_file(
            out.name, info="mzml"
        )
        scan = int(
            HypedsearchOutput.get_info_from_hypedsearch_comet_output_file(
                out.name, info="scan"
            )
        )
        mzml_names_to_scans[mzml_stem].add(scan)

    # Get paths of MZMLs from config to get {<mzml path>: [<mising scans>]}
    hs_config = HypedsearchConfig.from_yaml(config)
    mzml_path_to_scans = defaultdict(list)
    for mzml, missing_scans in mzml_names_to_scans.items():
        # Find path of MZML
        mzml_path = None
        for path in hs_config.mzml_to_scans.keys():
            if Mzml.get_mzml_name(path) == mzml:
                mzml_path = path
                break
        assert mzml_path is not None, f"Could not find path for {mzml}"
        mzml_path_to_scans[mzml_path] = list(missing_scans)

    new_config = HypedsearchConfig(
        mzml_to_scans=dict(mzml_path_to_scans),
        database=hs_config.database,
        out_dir=results_dir,
        crux_comet_params=hs_config.crux_comet_params,
        fasta=hs_config.fasta,
        precursor_mz_ppm_tol=hs_config.precursor_mz_ppm_tol,
        peak_to_ion_ppm_tol=hs_config.peak_to_ion_ppm_tol,
        kmer_to_protein_map=hs_config.kmer_to_protein_map,
    )
    if out_path is not None:
        new_config.to_yaml(path=out_path)
    return new_config


# def find_seq_in_psms(psms: List[CometPSM], seq: )


@dataclass
class NativeHybridPSMComparison:
    native_psms: CometPSMs
    hybrid_psms: CometPSMs

    @property
    def scans_with_native_and_hybrid_psm(self) -> Set[Tuple[str, int]]:
        return self.native_psms.scans.intersection(self.hybrid_psms.scans)

    def filter_hybrids(self, kmers: Dict, min_ion_len: int):
        filtered_hybrids = []
        for psm in self.hybrid_psms.psms:
            possible_hybrids = find_possible_hybrids(
                seq=psm.seq, kmers=kmers, min_side_len=min_ion_len
            )
            if len(possible_hybrids) > 0:
                filtered_hybrids.append(psm)
        self.hybrid_psms.psms = filtered_hybrids


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    cli.add_command(cli_run_hypedsearch)
    cli.add_command(cli_create_hypdedsearch_config)
    # cli.add_command(cli_run_comet)
    # cli.add_command(cli_run_assign_confidence)
    cli()
