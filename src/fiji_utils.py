import logging
import re
import shutil
from pathlib import Path

import click
import pandas as pd

from src.hypedsearch import HybridRunConfig
from src.utils import PathType, log_params, setup_logger

logger = logging.getLogger(__name__)


def collect_benchmarking_data(dir: Path, out_path: Path):
    benchmark_file_regex = r"^(?P<mzml>(.+?))\.(?P<scan>\d+)\.log$"
    dfs = []
    mzmls = []
    scans = []
    for f in dir.glob("*.log"):
        match = re.match(benchmark_file_regex, f.name).groupdict()
        mzmls.append(match["mzml"])
        scans.append(match["scan"])
        dfs.append(pd.read_csv(f, sep="\t"))
    df = pd.concat(dfs)
    df["mzml"] = mzmls
    df["scan"] = scans
    df.to_csv(out_path, index=False)


@click.command(
    name="prep-files-on-fiji",
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
def cli_prep_files_on_fiji(
    config: Path,
    data_dir: Path,
):
    logger = setup_logger()
    logger.info("Setting up files on Fiji node")
    original_config = HybridRunConfig.from_yaml(path=config)
    fiji_node_config = original_config.create_config_for_fiji_node(
        node_data_dir=data_dir
    )
    fiji_node_config.save(path=data_dir / config.name)
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
    collect_benchmarking_data(dir=dir, out_path=out_path)


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
    cli.add_command(cli_prep_files_on_fiji)
    cli.add_command(cli_collect_benchmark_data)
    # cli.add_command(cli_cleanup)
    cli()
