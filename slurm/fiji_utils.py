import datetime
import logging
import re
import shutil
import textwrap
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Union

import click
import pandas as pd

from src.constants import DEFAULT_JCT_LEN, DEFAULT_Q_THRESHOLD
from src.hypedsearch import HypedsearchRunConfig, run_hypedsearch
from src.utils import (
    PathType,
    check_if_file_is_empty,
    copy_file,
    log_params,
    setup_logger,
    write_new_line_separated_file,
)

logger = logging.getLogger(__name__)

RUN_CMD_SH = "./slurm/run_command_via_slurm.sh"
SHEBANG = "#!/bin/bash"
LOCALSCRATCH = "/localscratch"
DEFAULT_FIJI_CONFIG_NAME = "hs.fiji.config.json"
DEFAULT_MEM = "500GB"
DEFAULT_N_CORES = 180
DEFAULT_TIME = "24:00:00"
DEFAULT_LOG_DIR = "logs/hypedsearch"
DEFAULT_PARTITION = "highmem"
node_to_data_dir_map = {"": "/localscratch"}
HIGH_MEM_NODES = [67, 68, 69, 70, 71, 73]
SHORT_PARTITION_NUM_CPUS = 64
SHORT_PARTITION_RAM = "450GB"


@dataclass
class SbatchConfig:
    name: str
    time: str = DEFAULT_TIME
    mem: str = DEFAULT_MEM
    n_cores: int = DEFAULT_N_CORES
    partition: str = DEFAULT_PARTITION
    log_dir: Union[str, Path] = DEFAULT_LOG_DIR
    n_nodes: int = 1
    nodelist: Optional[str] = None

    def get_sbatch_directives_lines(self) -> List[str]:
        header = [
            f"#SBATCH --job-name={self.name}",
            f"#SBATCH --mem={self.mem}",
            f"#SBATCH --ntasks={self.n_cores}",
            f'#SBATCH --partition="{self.partition}"',
            f"#SBATCH --nodes={self.n_nodes}",
            f"#SBATCH --time={self.time}",
            f"#SBATCH --output={self.log_dir}/{self.name}.out",
            f"#SBATCH --error={self.log_dir}/{self.name}.err",
            "#SBATCH --mail-type=BEGIN,FAIL,END",
            "#SBATCH --mail-user=erjo3868@colorado.edu",
        ]
        if self.nodelist is not None:
            header.append(f"#SBATCH --nodelist={self.nodelist}")
        return header


def create_sbatch_script_to_run_hypedsearch(
    sbatch_config: SbatchConfig,
    config: Union[str, Path],
    out_path: Optional[Union[str, Path]] = None,
    n_cores: int = 80,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    fiji_config = create_config_for_fiji_run(config=config, dry_run=True)
    fiji_config_path = fiji_config.parent_output_dir / DEFAULT_FIJI_CONFIG_NAME
    cmds = [
        # f"python -m slurm.fiji_utils prep-files-on-fiji -c {config}",
        # f"python -m src.hypedsearch run-in-parallel -c {fiji_config_path} -n {n_cores} -os",
        f"python -m slurm.hypedsearch run-in-parallel -c {config} -n {n_cores} -os",
    ]
    lines = ["#!/bin/bash"] + sbatch_config.get_sbatch_directives_lines() + [""] + cmds
    if out_path is None:
        out_path = hs_config.parent_output_dir / f"{hs_config.name}.hybrid_run.sbatch"
    write_new_line_separated_file(lines=lines, path=out_path)


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
    logger.info(
        f"Moving and preparing files on Fiji node in directory {node_data_dir}..."
    )
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
    new_kmer_db_path = name_dir / hs_config.kmer_db_path.name
    if not dry_run:
        copy_file(src=hs_config.kmer_db_path, dest=new_kmer_db_path)

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
    new_config.kmer_db_path = new_kmer_db_path
    new_config.crux_comet_params = new_params_path
    new_config.parent_output_dir = name_dir
    new_config.fasta_fm_index = new_fm_index_path

    return new_config


@dataclass
class ConfigAndPath:
    config_path: Path
    hs_config: HypedsearchRunConfig = field(init=False)

    def __post_init__(self):
        self.hs_config = HypedsearchRunConfig.from_json(path=self.config_path)


@dataclass
class HypedsearchRunScripts:
    config_paths: List[str | Path]
    script_dir: Path
    configs: List[ConfigAndPath] = field(init=False)

    def __post_init__(self):
        self.script_dir = Path(self.script_dir)
        self.script_dir.mkdir(parents=True, exist_ok=True)
        self.configs = [
            ConfigAndPath(config_path=config_path) for config_path in self.config_paths
        ]

    def create_script_to_create_custom_psm_dfs(self, time: str = "08:00:00"):
        lines = [SHEBANG]
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name psmDFs_{config.hs_config.name}",
                '--mem "2GB"',
                '--part "short"',
                f'--time "{time}"',
                "--cores 2",
                f'--cmd "python cli.py create-psm-dfs -c {config.config_path}"',
            ]
            lines.append(" ".join(cmd_parts))
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "create_custom_psm_dfs_via_slurm.sh"
        )

    def create_script_to_run_native_run_via_slurm(
        self, n_cores: int = 60, time: str = "10:00:00"
    ):
        lines = [SHEBANG]
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name {config.hs_config.name}_nativeRun",
                '--mem "5GB"',
                '--part "short"',
                f'--time "{time}"',
                f"--cores {n_cores}",
                f'--cmd "python -m cli native-run -os -c {config.config_path}"',
            ]
            lines.append(" ".join(cmd_parts))
            write_new_line_separated_file(
                lines=lines, path=self.script_dir / "native_run_via_slurm.sh"
            )

    def create_script_to_run_native_run_locally(self):
        lines = [SHEBANG]
        for config in self.configs:
            lines.append(f"python -m cli native-run -c {config.config_path}")
            write_new_line_separated_file(
                lines=lines, path=self.script_dir / "native_run_locally.sh"
            )

    def create_script_for_hybrid_run_via_slurm(
        self, n_cores: int = 10, nodes: List[int] = HIGH_MEM_NODES
    ):
        fiji_nodes_to_use = [f"fijinode-{num}" for num in nodes]
        lines = [SHEBANG]
        idx = 0
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name {config.hs_config.name}_hybridRun",
                '--mem "500GB"',
                '--part "highmem"',
                '--time "24:00:00"',
                "--cores 180",
                f'--nodelist "{fiji_nodes_to_use[idx % len(nodes)]}"',
                f'--cmd "python -m cli run-hypedsearch -c {config.config_path} -n {n_cores} -os -p"',
            ]
            lines.append(" ".join(cmd_parts))
            idx += 1
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "hybrid_run_via_slurm.sh"
        )

    def slurm_script_for_hybrid_run_on_short_nodes(
        self, n_cores: int = 10, script_name: str = "short_partition_hybrid_run.sh"
    ):
        lines = [SHEBANG]
        idx = 0
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name {config.hs_config.name}_hybridRun",
                f'--mem "{SHORT_PARTITION_RAM}"',
                '--part "short"',
                '--time "23:00:00"',
                f"--cores {SHORT_PARTITION_NUM_CPUS}",
                f'--cmd "python -m cli run-hypedsearch -c {config.config_path} -n {n_cores} -os -p"',
            ]
            lines.append(" ".join(cmd_parts))
            idx += 1
        write_new_line_separated_file(lines=lines, path=self.script_dir / script_name)

    def create_script_for_native_plots_via_slurm(self):
        lines = [SHEBANG]
        idx = 0
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name natPlots_{config.hs_config.name}",
                '--mem "1GB"',
                '--part "short"',
                '--time "01:00:00"',
                "--cores 1",
                f'--cmd "python cli.py create-native-plots -c {config.config_path}"',
            ]
            lines.append(" ".join(cmd_parts))
            idx += 1
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "create_native_plots_via_slurm.sh"
        )

    def create_script_for_hybrid_plots_via_slurm(self):
        lines = [SHEBANG]
        idx = 0
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name hyPlots_{config.hs_config.name}",
                '--mem "1GB"',
                '--part "short"',
                '--time "01:00:00"',
                "--cores 1",
                f'--cmd "python cli.py create-hybrid-plots -c {config.config_path}"',
            ]
            lines.append(" ".join(cmd_parts))
            idx += 1
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "create_hybrid_plots_via_slurm.sh"
        )

    def create_script_for_local_hybrid_run(self):
        lines = [SHEBANG]
        for config in self.configs:
            cmd_parts = [
                "python -m cli run-in-parallel",
                f"-c {config.config_path} -n 8",
            ]
            lines.append(" ".join(cmd_parts))
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "hybrid_run_locally.sh"
        )

    def create_script_to_combine_hybrid_txts_via_slurm(self):
        lines = [SHEBANG]
        for config in self.configs:
            cmd_parts = [
                RUN_CMD_SH,
                f"--name combineHyTxts_{config.hs_config.name}",
                '--mem "2GB"',
                '--part "short"',
                '--time "01:00:00"',
                "--cores 10",
                f'--cmd "python -m cli combine-comet-txts -c {config.config_path}"',
            ]
            lines.append(" ".join(cmd_parts))
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "combine_hybrid_txts_via_slurm.sh"
        )

    def create_python_script_to_see_if_hybrid_run_txts_were_combined(self):
        config_paths_str = ", ".join(
            [f'Path("{str(config.config_path)}")' for config in self.configs]
        )
        script = f"""
        import sys
        from pathlib import Path
        repo_dir = Path(__file__).parents[3]
        assert (
            repo_dir.name == "hypedsearch"
        ), f"Expected repo_dir to be 'hypedsearch'"
        sys.path.append(str(repo_dir))
        from src.hypedsearch import HypedsearchRunConfig
        finished_configs, not_finished_configs = [], []
        for config_path in [{config_paths_str}]:
            hs_config = HypedsearchRunConfig.from_json(path=config_path)
            if len(list(hs_config.hybrid_run_dir.glob("*.txt"))) != 1:
                not_finished_configs.append(hs_config.name)
            else:
                finished_configs.append(hs_config.name)
        print(f"Finished configs: {{finished_configs}}")
        print(f"Not finished configs: {{not_finished_configs}}")
        """
        script = textwrap.dedent(script).lstrip(
            "\n"
        )  # remove common indent and optional leading blank line
        out_path = Path(self.script_dir / "check_if_hybrid_run_txts_were_combined.py")
        out_path.write_text(script, encoding="utf-8")

    def create_script_to_check_for_missing_comet_txts_locally(self):
        lines = [SHEBANG]
        for config in self.configs:
            cmd_parts = [
                "python -m cli get-missing-hs-outputs"
                f"--config {config.config_path}"
                "--verbose",
                f"> tmp/{config.hs_config.name}_missing_outputs.txt",
            ]
            lines.append(" ".join(cmd_parts))
        write_new_line_separated_file(
            lines=lines,
            path=self.script_dir / "check_for_missing_hs_outputs_locally.sh",
        )

    def create_script_to_process_native_and_hybrid_results(
        self, jct_len: int, q_threshold: float
    ):
        lines = [SHEBANG]
        for config in self.configs:
            lines.append(
                " ".join(
                    [
                        "python -m cli process-hs",
                        f"--config {config.config_path}",
                        f"--jct_len {jct_len}",
                        f"--q_threshold {q_threshold}",
                    ]
                )
            )
        write_new_line_separated_file(
            lines=lines,
            path=self.script_dir / "process_native_and_hybrid_results_locally.sh",
        )

    def create_script_to_run_param_medic_locally(self):
        lines = [SHEBANG]
        for config in self.configs:
            lines.append(
                " ".join(
                    [
                        "python -m cli run-param-medic",
                        f"--config {config.config_path}",
                    ]
                )
            )
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "run_param_medic.sh"
        )

    def create_script_for_hybrid_run_on_layer_lab_computer(self, cores: int):
        lines = [SHEBANG]
        for config in self.configs:
            lines.append(
                " ".join(
                    [
                        "uv run python cli.py run-hypedsearch ",
                        f"-n {cores}",
                        "-cp ../crux-5.0.0.Linux.x86_64/bin/crux",
                        "-p",
                        f"-c {config.config_path}",
                    ]
                )
            )
        write_new_line_separated_file(
            lines=lines, path=self.script_dir / "hybrid_run_on_layer_lab_comp.sh"
        )

    def create_all_scripts(
        self, jct_len: int = DEFAULT_JCT_LEN, q_threshold: float = DEFAULT_Q_THRESHOLD
    ):
        self.create_script_to_run_native_run_via_slurm()
        self.create_script_to_run_native_run_locally()
        self.create_script_for_local_hybrid_run()
        self.create_script_for_hybrid_run_via_slurm()
        self.create_script_to_combine_hybrid_txts_via_slurm()
        self.create_script_to_check_for_missing_comet_txts_locally()
        self.create_script_to_run_param_medic_locally()
        self.create_script_to_process_native_and_hybrid_results(
            jct_len=jct_len, q_threshold=q_threshold
        )


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


@click.command(
    name="run-hs-on-slurm",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help=(""),
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    required=True,
    help="",
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
def cli_run_hs_on_slurm(
    config: Path,
    n_cores: int,
    data_dir: Path,
):
    # Move files to local directory
    fiji_config = create_config_for_fiji_run(config=config, node_data_dir=data_dir)
    fiji_config_path = fiji_config.parent_output_dir / DEFAULT_FIJI_CONFIG_NAME
    fiji_config.save(path=fiji_config_path)

    # Run Hypedsearch
    run_hypedsearch(config=fiji_config_path, n_cores=n_cores, on_singularity=True)

    # Move files back to persistent storage
    move_scan_results_from_node_to_persistent_storage(
        config=config, node_data_dir=data_dir
    )


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
    cli.add_command(cli_run_hs_on_slurm)
    cli()
