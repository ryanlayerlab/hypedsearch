import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import click
import pandas as pd

from src.constants import (
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_THRESHOLD,
    PSMS_DF_NAME,
)
from src.hypedsearch import HypedsearchRunConfig
from src.hypedsearch_run_analysis import ResultsAnalysis, create_junction_support_df
from src.mass_spectra import Spectrum, create_spectra_plots
from src.plot_utils import fig_setup, save_fig
from src.psm import CometPSM, ProteinAbundance, convert_comet_psms_to_psms
from src.utils import PathType, flatten_list_of_lists, setup_logger

logger = logging.getLogger(__name__)


@click.command(
    "spectra-plots",
    help=("Create spectra plots for spectra defined in run config"),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to run config",
)
def cli_spectra_plots(config: Path):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    create_spectra_plots(
        spectra=list(hs_config.spectrum_uid_to_spectrum.values()),
        sample=hs_config.name,
        out_dir=hs_config._results_dir,
    )


@click.command(
    "create-psms",
    help=("Create and save spectra PSM objects"),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to run config",
)
@click.option(
    "--hy_side_len",
    "-s",
    type=int,
    default=DEFAULT_MIN_SIDE_LEN,
    show_default=True,
    required=True,
    help="Path to run config",
)
def cli_create_psms(config: Union[Path, str], hy_side_len: int):
    hs_config = HypedsearchRunConfig.from_json(path=config)

    # Spectrum plots
    create_spectra_plots(
        spectra=list(hs_config.spectrum_uid_to_spectrum.values()),
        sample=hs_config.name,
        out_dir=hs_config._results_dir,
    )

    # Create and save SpectrumPSMs objects
    # if not (hs_config._results_dir / SPECTRUM_PSMS_NAME).exists():
    hs_out = ResultsAnalysis(
        hs_config=hs_config,
        min_side_len=hy_side_len,
        remove_carbamidomethylation=True,
    )
    hs_out.collect_outputs()
    hs_out.set_q_values()
    psms = hs_out.save_spectrum_psms()
    hs_out.save_hybrid_peptides()

    # Save dataframes
    psm_df = pd.DataFrame(flatten_list_of_lists(psm.to_rows() for psm in psms))
    psm_df.to_csv(
        hs_config._results_dir / PSMS_DF_NAME,
        index=False,
    )


@click.command(
    "process",
    help=("Load all spectra and PSMs, set q-values, and save to json files"),
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to run config",
)
@click.option(
    "--hy_side_len",
    "-s",
    type=int,
    default=DEFAULT_MIN_SIDE_LEN,
    show_default=True,
    required=True,
    help="Path to run config",
)
@click.option(
    "--q_threshold",
    "-q",
    type=float,
    default=DEFAULT_Q_THRESHOLD,
    show_default=True,
    required=True,
    help="q-value threshold",
)
def cli_process(config: Path, hy_side_len: int, q_threshold: float):
    hs_config = HypedsearchRunConfig.from_json(path=config)
    hs_out = ResultsAnalysis(
        hs_config=hs_config, min_side_len=hy_side_len, remove_carbamidomethylation=True
    )
    hs_out.collect_outputs()
    hs_out.set_q_values()

    # Spectra plots
    spectrum_df = create_spectra_plots(
        spectra=list(hs_out.spectrum_uid_to_spectrum.values()),
        sample=hs_config.name,
        out_dir=hs_config._results_dir,
    )

    # Native targets
    # Protein abundance plot
    logger.info("Computing protein abundances...")
    prot_ab = ProteinAbundance.from_comet_psms(
        quality_psms=list(hs_out.top_native_targets.values()),
        q_threshold=q_threshold,
    )
    fig, axs = fig_setup(h=8, w=8)
    prot_ab.plot_sorted_prot_cnts(top_n_prots=20, ax=axs[0])
    save_fig(
        path=hs_config._results_dir / "protein_abundance.png",
        fig=fig,
        title=f"{hs_config.name} (approx) protein abundance (q <= {q_threshold})",
    )

    logger.info("Processing native targets...")
    n_t_df = process_psms(
        psms=list(hs_out.top_native_targets.values()),
        uid_to_spectrum_map=hs_out.spectrum_uid_to_spectrum,
        psm_type="native_target",
        prot_ab=prot_ab,
    )
    n_t_df.to_csv(hs_config._results_dir / "native_target_psms.csv", index=False)

    # Native decoys
    logger.info("Processing native decoys...")
    n_d_df = process_psms(
        psms=list(hs_out.top_native_decoys.values()),
        uid_to_spectrum_map=hs_out.spectrum_uid_to_spectrum,
        psm_type="native_decoy",
    )
    n_d_df.to_csv(hs_config._results_dir / "native_decoy_psms.csv", index=False)

    # Hybrid targets
    logger.info("Processing hybrid targets...")
    h_t_df = process_psms(
        psms=list(hs_out.top_hybrid_targets.values()),
        uid_to_spectrum_map=hs_out.spectrum_uid_to_spectrum,
        psm_type="hybrid_target",
    )
    h_t_df.to_csv(hs_config._results_dir / "hybrid_target_psms.csv", index=False)

    # Combine all PSMs
    df = pd.concat(
        [
            n_t_df.merge(spectrum_df, on="uid"),
            n_d_df.merge(spectrum_df, on="uid"),
            h_t_df.merge(spectrum_df, on="uid"),
        ]
    )
    df["seq_len"] = df["seq"].apply(len)
    df.to_csv(hs_config._results_dir / "all_psms.csv", index=False)

    # Junction analysis
    jct_df = create_junction_support_df(
        hybrid_psms=list(hs_out.top_hybrid_targets.values()),
        protein_name_to_seq_map=hs_out.protein_name_to_seq_map,
        min_aa_jct_len=3,
        seq_to_hybrids_map=hs_out.seq_to_hybrids_map,
        spectrum_uid_to_spectrum_map=hs_out.spectrum_uid_to_spectrum,
    )
    jct_df.to_csv(hs_config._results_dir / "hy_jct.csv", index=False)

    # Save spectra and ComePSMs
    Spectrum.save_to_json(
        spectra=list(hs_out.spectrum_uid_to_spectrum.values()),
        path=hs_config._results_dir / "spectra.json",
    )
    CometPSM.save_to_json(
        psms=list(hs_out.top_native_targets.values()),
        path=hs_config._results_dir / "top_native_targets.json",
    )
    CometPSM.save_to_json(
        psms=list(hs_out.top_native_decoys.values()),
        path=hs_config._results_dir / "top_native_decoys.json",
    )
    CometPSM.save_to_json(
        psms=list(hs_out.top_hybrid_targets.values()),
        path=hs_config._results_dir / "top_hybrid_targets.json",
    )


# def create_spectrum_psm_df(
#     spectrum_psms: List[SpectrumPSMs],
#     prot_ab: Optional[ProteinAbundance] = None,
# ):
#     # Create dataframe
#     rows = []
#     for spectrum_psm in spectrum_psms:
#         if spectrum_psm.native_target:

#         for psm in spectrum_psm.psms:
#             row = psm.to_row()
#             if prot_ab is not None:
#                 # Get protein abundance
#                 row["prot_ab"] = max(
#                     prot_ab.get_rel_ab(protein=prot) for prot in psm.positions
#                 )
#             rows.append(row)
#     df = pd.DataFrame(rows)
#     return df


def process_psms(
    psms: List[CometPSM],
    uid_to_spectrum_map: Dict[str, Spectrum],
    psm_type: str,
    prot_ab: Optional[ProteinAbundance] = None,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
) -> pd.DataFrame:
    psms = convert_comet_psms_to_psms(
        comet_psms=psms,
        uid_to_spectrum_map=uid_to_spectrum_map,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
    )
    # Create dataframe
    rows = []
    for psm in psms:
        row = psm.to_row()
        if prot_ab is not None:
            # Get protein abundance
            row["prot_ab"] = max(
                prot_ab.get_rel_ab(protein=prot) for prot in psm.positions
            )
        rows.append(row)
    df = pd.DataFrame(rows)
    df["psm_type"] = psm_type
    return df


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_spectra_plots)
    cli.add_command(cli_process)
    cli.add_command(cli_create_psms)
    cli()
