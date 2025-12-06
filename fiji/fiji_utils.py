import datetime
import logging
import re
import shutil
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import click
import pandas as pd
from pydantic import BaseModel, field_validator

from src.constants import LINUX_CRUX_EXECUTABLE, MAC_CRUX_EXECUTABLE
from src.hypedsearch import HybridPSMScorer, HypedsearchRunConfig
from src.utils import PathType, copy_file, log_params, setup_logger

logger = logging.getLogger(__name__)

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


@dataclass
class HypedsearchOnFijiConfig:
    hs_config: HypedsearchRunConfig

    @classmethod
    def from_json(cls, path: Union[str, Path]):
        hs_config = HypedsearchRunConfig.from_json(path=path)
        return cls(hs_config=hs_config)

    def prepare_files_for_hybrid_run(
        self,
        node: Optional[str] = None,
        node_data_dir: Optional[Path] = None,
        out_path: Optional[Union[str, Path]] = None,
    ) -> HypedsearchRunConfig:
        if node_data_dir is None:
            node_data_dir = node_to_data_dir_map[node]

        # Create directories that need to exist
        name_dir = node_data_dir / self.hs_config.name
        name_dir.mkdir(parents=True, exist_ok=True)

        # Copy MZMLs
        new_mzml_to_scans = {}
        for mzml, scans in self.hs_config.mzml_to_scans.items():
            new_mzml_path = str(name_dir / Path(mzml).name)
            new_mzml_to_scans[new_mzml_path] = scans
            copy_file(src=mzml, dest=new_mzml_path)

        # Copy other files
        # comet.params
        old_comet_params = self.hs_config.psm_scorer.comet_params
        new_comet_params = name_dir / old_comet_params.name
        copy_file(src=old_comet_params, dest=new_comet_params)

        # k-mer database
        old_kmer_db = self.hs_config.hybrid_former.kmer_db
        new_kmer_db = name_dir / old_kmer_db.name
        copy_file(src=old_kmer_db, dest=new_kmer_db)

        # FASTA
        old_fasta = self.hs_config.hybrid_former.fasta
        new_fasta = name_dir / old_fasta.name
        copy_file(src=old_fasta, dest=new_fasta)

        # Create new config
        hybrid_former = deepcopy(self.hs_config.hybrid_former)
        hybrid_former.kmer_db = new_kmer_db
        hybrid_former.fasta = new_fasta

        psm_scorer = deepcopy(self.hs_config.psm_scorer)
        if psm_scorer.fasta is not None:
            psm_scorer.fasta = new_fasta
        psm_scorer.comet_params = new_comet_params

        fiji_config = HypedsearchRunConfig(
            mzml_to_scans=new_mzml_to_scans,
            parent_out_dir=name_dir,
            name=self.hs_config.name,
            spectrum_selector=self.hs_config.spectrum_selector,
            spectrum_preprocessor=self.hs_config.spectrum_preprocessor,
            psm_scorer=psm_scorer,
            hybrid_former=hybrid_former,
        )
        if out_path is None:
            out_path = name_dir / f"hs.fiji.config.json"
        fiji_config.to_json(path=out_path)
        return fiji_config

    def move_scan_results(
        self,
        node_data_dir: Path,
    ):
        fiji_config = self.prepare_files_for_hybrid_run(node_data_dir=node_data_dir)
        # Move native scan results
        src = fiji_config.native_run_scan_results_dir
        dest = self.hs_config.native_run_scan_results_dir
        logger.info(f"Copying native scan results from {src} to {dest}...")
        shutil.copytree(src, dest, dirs_exist_ok=True)

        # Move hybrid scan results
        src = fiji_config.hybrid_run_scan_results_dir
        dest = self.hs_config.hybrid_run_scan_results_dir
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
    fiji_runner = HypedsearchOnFijiConfig.from_json(path=config)
    _ = fiji_runner.prepare_files_for_hybrid_run(node_data_dir=data_dir)
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
    logger = setup_logger()
    fiji_runner = HypedsearchOnFijiConfig.from_json(path=config)
    fiji_runner.move_scan_results(node_data_dir=data_dir)


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
    cli.add_command(cli_prep_files_on_fiji)
    cli.add_command(cli_collect_benchmark_data)
    cli.add_command(cli_move_scan_results)
    cli()
