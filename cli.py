from pathlib import Path
from typing import Optional

import click

from src.constants import DEFAULT_JCT_LEN, DEFAULT_Q_THRESHOLD, MAC_CRUX_EXECUTABLE
from src.hypedsearch import HypedsearchOutputs, HypedsearchRunConfig, run_hypedsearch
from src.kmer_database import KmerDatabase
from src.postprocess_hs_results import process_native_and_hybrid_runs_via_config
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
    name="process-native-run",
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
def cli_process_native_run(
    config: Path,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.create_native_run_plots()


@click.command(
    name="run-hypedsearch",
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
    help="",
)
@click.option(
    "--parallel",
    "-p",
    is_flag=True,
    help="",
)
def cli_run_hypedsearch(
    config: Path, n_cores: int, on_singularity: bool, parallel: bool
):
    if on_singularity:
        crux_path = None
    else:
        crux_path = MAC_CRUX_EXECUTABLE
    run_hypedsearch(
        config=config,
        n_cores=n_cores,
        crux_path=crux_path,
        on_singularity=on_singularity,
        run_in_parallel=parallel,
    )


@click.command(
    name="process-hs",
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
    default=DEFAULT_Q_THRESHOLD,
    show_default=True,
    help="",
)
@click.option(
    "--jct_len",
    "-j",
    type=int,
    default=DEFAULT_JCT_LEN,
    show_default=True,
    help="",
)
@log_params
def cli_process_native_and_hybrid_runs_via_config(
    config: Path, q_threshold: float, jct_len: int
):
    process_native_and_hybrid_runs_via_config(
        q_threshold=q_threshold,
        jct_len=jct_len,
        hs_config=config,
    )


@click.command(
    name="combine-comet-txts",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--results_dir",
    "-rd",
    type=PathType(),
    required=False,
    help="Path to the directory containing scan results",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=False,
    help="Path to the output directory where the combined results for each mzML will be saved",
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=False,
    help="Path to the Hypedsearch config JSON",
)
def cli_combine_comet_txts(
    config: Optional[Path], results_dir: Optional[Path], out_dir: Optional[Path]
):
    if config is None:
        HypedsearchOutputs.combine_comet_results_by_mzml_and_psm_type(
            folder=results_dir, out_dir=out_dir
        )
    else:
        hs_config = HypedsearchRunConfig.from_json(path=config)
        hs_config.combine_hybrid_run_scan_outputs(overwrite=True)


@click.command(
    name="get-missing-hs-outputs",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=False,
    help="Path to the Hypedsearch config JSON",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="",
)
def cli_get_spectra_with_missing_hs_outputs(config: Path, verbose: bool):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.check_for_missing_scans(print_missing=verbose)


@click.command(
    name="run-param-medic",
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
def cli_run_param_medic(config: Path):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.run_param_medic()


@click.command(
    name="kmer-db-proteins",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--kmer_db",
    "-kdb",
    type=PathType(),
    required=True,
    help="Path to the kmer database",
)
def cli_kmer_db_proteins(kmer_db: Path):
    db = KmerDatabase(db_path=kmer_db)
    msg = f"Proteins in kmer database {kmer_db}:\n" + "\n".join(db.proteins)
    print(msg)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    logger = setup_logger()
    # Miscellaneous stuff
    cli.add_command(cli_run_param_medic)
    cli.add_command(cli_kmer_db_proteins)

    # Native run stuff
    cli.add_command(cli_native_comet_run)
    cli.add_command(cli_process_native_run)

    # Hybrid run stuff
    cli.add_command(cli_run_hypedsearch)
    cli.add_command(cli_combine_comet_txts)
    # cli.add_command(cli_process_native_and_hybrid_runs_via_config)
    cli.add_command(cli_get_spectra_with_missing_hs_outputs)

    cli()
