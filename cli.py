from pathlib import Path
from typing import Optional

import click

from src.constants import DEFAULT_Q_THRESHOLD, MAC_CRUX_EXECUTABLE
from src.hypedsearch import HypedsearchRunConfig, run_hypedsearch_in_parallel
from src.utils import PathType, log_params, setup_logger


@click.command(
    name="native-run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--on_singularity",
    "-os",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@log_params
def cli_native_comet_run(
    config: Path,
    on_singularity: bool,
):
    if on_singularity:
        crux_path = None
    else:
        crux_path = MAC_CRUX_EXECUTABLE
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.run_native_comet(crux_path=crux_path, on_singularity=on_singularity)


@click.command(
    name="analyze-native-run",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--q_threshold",
    "-q",
    type=float,
    required=True,
    default=DEFAULT_Q_THRESHOLD,
    show_default=True,
    help="q-value threshold",
)
@log_params
def cli_analyze_native_results(config: Path, q_threshold: float):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    # try:
    #     hs_config.run_native_assign_confidence()
    # except:
    #     pass
    hs_config.analyze_native_results()


@click.command(
    name="run-in-parallel",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    required=True,
    help="",
)
@click.option(
    "--on_singularity",
    "-os",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@log_params
def cli_run_in_parallel(config: Path, n_cores: int, on_singularity: bool):
    if on_singularity:
        crux_path = None
    else:
        crux_path = MAC_CRUX_EXECUTABLE
    run_hypedsearch_in_parallel(
        config=config,
        n_cores=n_cores,
        crux_path=crux_path,
        on_singularity=on_singularity,
    )


@click.command(
    name="",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    required=True,
    help="",
)
@click.option(
    "--on_singularity",
    "-os",
    is_flag=True,
    help="If outputs already exist, this controls whether or not to overwrite them.",
)
@log_params
def cli_run_in_parallel(config: Path, n_cores: int, on_singularity: bool):
    if on_singularity:
        crux_path = None
    else:
        crux_path = MAC_CRUX_EXECUTABLE
    run_hypedsearch_in_parallel(
        config=config,
        n_cores=n_cores,
        crux_path=crux_path,
        on_singularity=on_singularity,
    )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_native_comet_run)
    cli.add_command(cli_analyze_native_results)
    cli.add_command(cli_run_in_parallel)
    cli()
