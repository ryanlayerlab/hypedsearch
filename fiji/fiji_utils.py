import datetime
import logging
import re
import shutil
from pathlib import Path
from typing import Optional, Union

import click
import pandas as pd

from src.constants import CRUX_PATH_IN_SINGULARITY
from src.hypedsearch import HybridRunConfig, HypedsearchConfig
from src.utils import PathType, copy_file, log_params, setup_logger

logger = logging.getLogger(__name__)


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


def prep_hybrid_run_on_fiji(
    hybrid_run_config_path: Union[str, Path],
    node_data_dir: Path,
    crux_path: Path = CRUX_PATH_IN_SINGULARITY,
) -> HybridRunConfig:
    hybrid_run_config = HybridRunConfig.from_json(path=hybrid_run_config_path)

    # The out directory is path/to/<name>/hybrid_run/scan_results. We want to preserve
    # this except replace "path/to" with node_data_dir
    scan_results_dir = node_data_dir / hybrid_run_config.out_dir.relative_to(
        hybrid_run_config.out_dir.parents[2]
    )
    scan_results_dir.mkdir(parents=True, exist_ok=True)
    hybrid_run_dir = scan_results_dir.parents[0]
    name_dir = scan_results_dir.parents[1]

    # Handle MZMLs
    logger.info(f"Copying MZMLs to {name_dir}...")
    new_mzml_to_scans = {}
    for mzml, scans in hybrid_run_config.mzml_to_scans.items():
        new_mzml_path = str(name_dir / Path(mzml).name)
        new_mzml_to_scans[new_mzml_path] = scans
        copy_file(src=mzml, dest=new_mzml_path)

    # Handle files that need to be copied
    logger.info(f"Copying other needed files to {name_dir}...")
    fiji_config = {}
    for attr in ["kmer_db", "fasta", "crux_comet_params", "kmer_to_proteins_map"]:
        old_path = Path(getattr(hybrid_run_config, attr))
        new_out_path = str(name_dir / old_path.name)
        fiji_config[attr] = new_out_path
        copy_file(src=old_path, dest=new_out_path)

    # Create the new config
    fiji_node_config = hybrid_run_config.__class__(
        mzml_to_scans=new_mzml_to_scans,
        out_dir=str(scan_results_dir),
        fasta=fiji_config["fasta"],
        kmer_db=fiji_config["kmer_db"],
        kmer_to_proteins_map=fiji_config["kmer_to_proteins_map"],
        crux_path=crux_path,
        peak_to_ion_ppm_tol=hybrid_run_config.peak_to_ion_ppm_tol,
        precursor_mz_ppm_tol=hybrid_run_config.precursor_mz_ppm_tol,
        crux_comet_params=fiji_config["crux_comet_params"],
        min_cluster_len=hybrid_run_config.min_cluster_len,
        min_cluster_support=hybrid_run_config.min_cluster_support,
        max_allowed_ion_charge=hybrid_run_config.max_allowed_ion_charge,
        log_dir=hybrid_run_config.log_dir,
    )

    # Save new config
    out_path = hybrid_run_dir / hybrid_run_config_path.name
    logger.info(f"Saving config to {out_path}...")
    fiji_node_config.save(path=out_path)
    return fiji_node_config


@click.command(
    name="prep-hybrid-run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=(
        "Run Hypedsearch via snakemake on a Fiji node using the given snakemake config"
    ),
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the snakemake config file",
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
def cli_prep_hybrid_run_on_fiji(
    config: Path,
    data_dir: Path,
):
    logger = setup_logger()
    logger.info("Setting up files on Fiji node...")
    prep_hybrid_run_on_fiji(hybrid_run_config_path=config, node_data_dir=data_dir)
    logger.info("Finished preparing files on Fiji node")


def prep_native_run_on_fiji(
    hypedsearch_config: Union[HypedsearchConfig, str, Path],
    data_dir: Path,
):
    if isinstance(hypedsearch_config, (str, Path)):
        hypedsearch_config = HypedsearchConfig.from_json(path=hypedsearch_config)

    # Create output directory
    out_dir = data_dir / hypedsearch_config.name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Move files that need to be moved to the output directory
    fiji_config = {}
    for attr in ["fasta", "crux_comet_params"]:
        old_path = Path(getattr(hypedsearch_config, attr))
        new_out_path = str(out_dir / old_path.name)
        fiji_config[attr] = new_out_path
        copy_file(src=old_path, dest=new_out_path)

    # Handle MZMLs
    new_mzml_to_scans = {}
    for mzml, scans in old_config.mzml_to_scans.items():
        new_mzml_path = str(node_data_dir / Path(mzml).name)
        new_mzml_to_scans[new_mzml_path] = scans
        copy_file(src=mzml, dest=new_mzml_path)


@click.command(
    name="prep-native-run-on-fiji",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=("Native run via snakemake on a Fiji node"),
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the HypedSearch JSON config file",
)
@click.option(
    "--data_dir",
    "-d",
    type=PathType(),
    default=Path("/localscratch"),
    show_default=True,
    help=("Path to the data directory on the Fiji node. "),
)
@log_params
def cli_prep_native_run_on_fiji(
    hypedsearch_config: Path,
    data_dir: Path,
):
    logger = setup_logger()
    logger.info("Setting up files on Fiji node")
    original_config = HypedsearchConfig.from_json(path=hypedsearch_config)
    logger.info("Finished preparing files on Fiji node")


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


# @click.command(
#     name="cleanup",
#     context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
#     help=(
#         "Cleanup the files on the Fiji node (moving outputs to /scratch and deleting temporary others)"
#     ),
# )
# @click.option(
#     "--config",
#     "-c",
#     type=PathType(),
#     required=True,
#     help="Path to the snakemake config file",
# )
# def cli_cleanup(
#     config: Path,
# ):
#     setup_logger()
#     hs_on_fiji_config = HypedsearchOnFijiConfig.prepare_files_on_fiji_node(
#         run_config=config,
#         copy_files=False,
#     )
#     # Copy results directory to original output directory
#     shutil.copytree(
#         hs_on_fiji_config.fiji_config.out_dir,
#         hs_on_fiji_config.original_config.out_dir,
#         dirs_exist_ok=True,
#     )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    cli.add_command(cli_prep_hybrid_run_on_fiji)
    cli.add_command(cli_collect_benchmark_data)
    # cli.add_command(cli_cleanup)
    cli()
