import os

from src.mass_spectra import Spectrum

os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
from typing import Optional

import click
import numpy as np  # Must import after setting OPENBLAS_NUM_THREADS

from src.constants import (
    DEFAULT_JCT_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_THRESHOLD,
    HS_CONFIG_SUFFIX,
    MAC_CRUX_EXECUTABLE,
)
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
    run_hypedsearch,
)
from src.kmer_database import KmerDatabase
from src.postprocess_hs_results import NativeVsHybridComparison
from src.utils import (
    PathType,
    check_if_file_is_empty,
    log_params,
    setup_logger,
    to_json,
)


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
    hs_config.native_comet_run_on_all_spectra(crux_path=crux_path)


@click.command(
    name="check-native-run-complete",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config_dir",
    "-cd",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--config_suffix",
    "-cs",
    type=str,
    default=HS_CONFIG_SUFFIX,
    show_default=True,
    required=True,
    help="",
)
def cli_check_if_native_run_is_done(config_dir: Path, config_suffix: str):
    configs = list(config_dir.glob(f"*{config_suffix}"))
    logger.info(
        f"Found {len(configs)} configs in {config_dir} with suffix {config_suffix}. Checking for missing or empty comet output files..."
    )
    configs_with_missing_results = []
    for config in configs:
        hs_config = HypedsearchRunConfig.from_json(path=config)
        comet_runs = hs_config.native_comet_run_on_all_spectra(dry_run=True)
        for comet_run in comet_runs:
            if check_if_file_is_empty(path=comet_run.target) or check_if_file_is_empty(
                path=comet_run.decoy
            ):
                configs_with_missing_results.append(config)
    if len(configs_with_missing_results) == 0:
        logger.info(
            "All comet runs have non-empty output files. Native run appears to be complete!"
        )
    else:
        logger.warning(
            f"{len(configs_with_missing_results)} configs have missing or empty comet output files. These configs are:\n"
            + "\n".join([str(c) for c in configs_with_missing_results])
        )


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
    name="create-native-plots",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="",
)
def cli_create_native_plots(
    config: Path,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.create_spectra_plots()
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
    "--crux_path",
    "-cp",
    type=PathType(),
    required=False,
    help="",
)
@click.option(
    "--n_cores",
    "-n",
    type=int,
    default=8,
    show_default=True,
    required=False,
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
def cli_run_hypedsearch_via_config(
    config: Path,
    n_cores: int,
    on_singularity: bool,
    parallel: bool,
    crux_path: Path | None,
):
    if on_singularity:
        crux_path = None
    else:
        if crux_path is None:
            crux_path = MAC_CRUX_EXECUTABLE
    # if parallel:
    #     import multiprocessing as mp
    #     mp.set_start_method("spawn", force=True)
    run_hypedsearch(
        config=config,
        n_cores=n_cores,
        crux_path=crux_path,
        run_in_parallel=parallel,
    )


@click.command(
    name="run-hypedsearch-on-spectrum",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--hypedsearch_config",
    "-hsc",
    type=PathType(),
    required=False,
    help="",
)
@click.option(
    "--hybrid_run_config",
    "-hrc",
    type=PathType(),
    required=False,
    help="",
)
@click.option(
    "--crux_path",
    "-cp",
    type=PathType(),
    required=False,
    help="",
)
@click.option(
    "--out_dir",
    "-od",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--mzml",
    "-m",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--scan",
    "-s",
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
def cli_run_hypedsearch_on_spectrum(
    out_dir: Path,
    hypedsearch_config: Path | None,
    hybrid_run_config: Path | None,
    on_singularity: bool,
    crux_path: Path | None,
    mzml: Path,
    scan: int,
):
    if on_singularity:
        crux_path = None
    else:
        if crux_path is None:
            crux_path = MAC_CRUX_EXECUTABLE
    # if parallel:
    #     import multiprocessing as mp
    #     mp.set_start_method("spawn", force=True)
    if hypedsearch_config is not None:
        hybrid_run_params = HypedsearchRunConfig.from_json(
            path=hypedsearch_config
        ).hybrid_run_params
    elif hybrid_run_config is not None:
        hybrid_run_params = HybridRunParams.load(path=hybrid_run_config)
    else:
        raise ValueError(
            "Must provide either a Hypedsearch config or a hybrid run config."
        )
    hybrid_run_on_spectrum(
        spectrum=Spectrum.get_spectrum(mzml=mzml, scan=scan),
        params=hybrid_run_params,
        fasta_dir=out_dir,
        crux_path=crux_path,
        overwrite=True,
        out_dir=out_dir,
    )


@click.command(
    name="combine-comet-txts",
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
def cli_combine_comet_txts(config: Path):
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
    name="kmer-db-info",
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
def cli_kmer_db_info(kmer_db: Path):
    db = KmerDatabase(db_path=kmer_db)
    msg = (
        f"k-mer database info:\n"
        f"- min_k = {db.min_k}\n"
        f"- max_k = {db.max_k}\n"
        f"- proteins in database: {db.proteins}"
    )
    print(msg)


@click.command(
    name="neofusion",
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
    "--ppm_tol",
    "-p",
    type=float,
    default=DEFAULT_PEAK_TO_ION_PPM_TOL,
    show_default=True,
    required=False,
    help="",
)
def cli_get_spectra_with_missing_hs_outputs(config: Path, verbose: bool):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_config.check_for_missing_scans(print_missing=verbose)


@click.command(
    name="create-psm-dfs",
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
# @click.option(
#     "--out_dir",
#     "-o",
#     type=PathType(),
#     required=True,
#     help="",
# )
@click.option(
    "--ppm_tol",
    "-p",
    type=float,
    default=DEFAULT_PEAK_TO_ION_PPM_TOL,
    show_default=True,
    required=False,
    help="",
)
def cli_create_psm_dfs(
    config: Path,
    ppm_tol: float,
    # out_dir: Path
):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    comp = NativeVsHybridComparison.from_config(config=config)
    comp.create_psm_dataframes(ppm_tol=ppm_tol, out_dir=hs_config.name_dir)


@click.command(
    name="create-hybrid-plots",
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
    "--jct_len",
    "-jl",
    type=int,
    required=False,
    default=DEFAULT_JCT_LEN,
    show_default=True,
    help="",
)
def cli_create_hybrid_run_plots(
    config: Path,
    jct_len: int,
):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    comp = NativeVsHybridComparison.from_config(config=config)
    comp.create_hybrid_run_plots(
        psm_df_dir=hs_config.name_dir,
        out_dir=hs_config.hybrid_run_dir,
        title=hs_config.name,
    )
    to_json(
        data=[psm.uid for psm in comp.neofusion_output.accepted_hybrids],
        path=hs_config.hybrid_run_dir / "neofusion_accepted_hybrid_uids.json",
    )
    comp.junction_analysis(
        jct_len=jct_len,
        out_dir=hs_config.hybrid_run_dir,
        title=hs_config.name,
    )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    logger = setup_logger()
    # Miscellaneous stuff
    cli.add_command(cli_run_param_medic)
    cli.add_command(cli_kmer_db_info)
    cli.add_command(cli_check_if_native_run_is_done)

    # Native run stuff
    cli.add_command(cli_native_comet_run)
    cli.add_command(cli_process_native_run)
    cli.add_command(cli_create_native_plots)

    # Hybrid run stuff
    cli.add_command(cli_run_hypedsearch_via_config)
    cli.add_command(cli_combine_comet_txts)
    cli.add_command(cli_get_spectra_with_missing_hs_outputs)
    cli.add_command(cli_run_hypedsearch_on_spectrum)

    # Postprocessing
    cli.add_command(cli_create_psm_dfs)
    cli.add_command(cli_create_hybrid_run_plots)

    cli()
