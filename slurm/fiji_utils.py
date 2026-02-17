import datetime
import logging
import re
import shutil
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Union

import click
import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import LINUX_CRUX_EXECUTABLE, MAC_CRUX_EXECUTABLE
from src.hypedsearch import HybridPSMScorer, HypedsearchRunConfig
from src.utils import PathType, copy_file, log_params, setup_logger

logger = logging.getLogger(__name__)

LOCALSCRATCH = "/localscratch"
DEFAULT_FIJI_CONFIG_NAME = "hs.fiji.config.json"
node_to_data_dir_map = {"": "/localscratch"}


def collect_benchmarking_data(dir: Path) -> pd.DataFrame:
    benchmark_file_regex = r"^(?P<mzml>(.+?))\.(?P<scan>\d+)\.log$"
    dfs = []
    mzmls = []
    scans = []
    dates = []
    for f in Path(dir).glob("*.log"):
        match = re.match(benchmark_file_regex, f.name).groupdict()
        mzmls.append(match["mzml"])
        scans.append(match["scan"])
        dates.append(datetime.datetime.fromtimestamp(f.stat().st_ctime))
        dfs.append(pd.read_csv(f, sep="\t"))
    df = pd.concat(dfs)
    df["mzml"] = mzmls
    df["scan"] = scans
    df["date"] = dates
    return df


def create_config_for_fiji_run(
    config: Union[str, Path],
    node_data_dir: Union[str, Path] = LOCALSCRATCH,
    dry_run: bool = False,
) -> HypedsearchRunConfig:
    hs_config = HypedsearchRunConfig.from_json(path=config)

    # Create directories that need to exist
    name_dir = Path(node_data_dir) / hs_config.name
    if not dry_run:
        name_dir.mkdir(parents=True, exist_ok=True)

    # Copy MZMLs
    new_mzml_to_scans = {}
    for mzml, scans in hs_config.mzml_to_scans.items():
        new_mzml_path = name_dir / Path(mzml).name
        new_mzml_to_scans[new_mzml_path] = scans
        if not dry_run:
            copy_file(src=mzml, dest=new_mzml_path)

    # Copy other files
    # comet.params
    new_params_path = name_dir / hs_config.crux_comet_params.name
    if not dry_run:
        copy_file(
            src=hs_config.crux_comet_params,
            dest=new_params_path,
        )

    # k-mer database
    new_kmer_db_path = name_dir / hs_config.kmer_db.name
    if not dry_run:
        copy_file(src=hs_config.kmer_db, dest=new_kmer_db_path)

    # FASTA
    new_fasta_path = name_dir / hs_config.fasta.name
    if not dry_run:
        copy_file(src=hs_config.fasta, dest=new_fasta_path)

    # FM-index
    new_fm_index_path = name_dir / hs_config.fasta_fm_index.name
    if not dry_run:
        copy_file(src=hs_config.fasta_fm_index, dest=new_fm_index_path)

    # Create new config
    new_config = deepcopy(hs_config)
    new_config.mzml_to_scans = new_mzml_to_scans
    new_config.fasta = new_fasta_path
    new_config.kmer_db = new_kmer_db_path
    new_config.crux_comet_params = new_params_path
    new_config.parent_output_dir = name_dir
    new_config.fasta_fm_index = new_fm_index_path

    return new_config


def move_scan_results_from_node_to_persistent_storage(
    config: Path,
    node_data_dir: Path,
):
    non_fiji_config = HypedsearchRunConfig.from_json(path=config)
    fiji_config = create_config_for_fiji_run(
        config=config, node_data_dir=node_data_dir, dry_run=True
    )
    # Move hybrid scan results
    src = fiji_config.hybrid_run_scan_results_dir
    dest = non_fiji_config.hybrid_run_scan_results_dir
    logger.info(f"Copying hybrid scan results from {src} to {dest}...")
    shutil.copytree(src, dest, dirs_exist_ok=True)


@click.command(
    name="prep-files-on-fiji",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=("Prepare files on Fiji node for Hypedsearch run"),
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch JSON config",
)
@click.option(
    "--data_dir",
    "-d",
    type=PathType(),
    default=Path("/localscratch"),
    show_default=True,
    help=(
        "Path to the data directory on the Fiji node. "
        "A Hypedsearch config file will be created here: <data_dir>/<config.name>"
    ),
)
@log_params
def cli_prep_files_on_fiji(
    config: Path,
    data_dir: Path,
):
    logger = setup_logger()
    logger.info("Setting up files on Fiji node...")
    fiji_config = create_config_for_fiji_run(config=config, node_data_dir=data_dir)
    fiji_config.save(path=fiji_config.parent_output_dir / DEFAULT_FIJI_CONFIG_NAME)
    logger.info("Finished preparing files on Fiji node")


@click.command(
    name="move-scan-results",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=(
        "Move scan results files from a Fiji node data directory back to the main output directory"
    ),
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch JSON config",
)
@click.option(
    "--data_dir",
    "-d",
    type=PathType(),
    default=Path("/localscratch"),
    show_default=True,
    help=(
        "Path to the data directory on the Fiji node. "
        "A Hypedsearch config file will be created here: <data_dir>/<config.name>"
    ),
)
@log_params
def cli_move_scan_results(
    config: Path,
    data_dir: Path,
):
    move_scan_results_from_node_to_persistent_storage(
        config=config, node_data_dir=data_dir
    )


@click.command(
    name="collect-benchmark-data",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=("Collect snakemake benchmarking data from a directory into a file"),
)
@click.option(
    "--dir",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the directory with the benchmarking files",
)
@click.option(
    "--out_path",
    "-o",
    type=PathType(),
    required=True,
    help="Path to the output file where the collected benchmarking data will be saved",
)
def cli_collect_benchmark_data(
    dir: Path,
    out_path: Path,
):
    df = collect_benchmarking_data(dir=Path(dir))
    df.to_csv(out_path, index=False)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_prep_files_on_fiji)
    cli.add_command(cli_collect_benchmark_data)
    cli.add_command(cli_move_scan_results)
    cli()
