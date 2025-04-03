import logging
from collections import Counter
from pathlib import Path
from typing import List, Optional

import click
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from src.comet_utils import CometPSM, load_comet_data
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
from src.utils import PathType, flatten_list_of_lists

logger = logging.getLogger(__name__)


def common_arguments(func):
    """Decorator for the common arguments to the click commands"""
    func = click.option(
        "--comet_results_dir",
        "-d",
        type=PathType(),
        required=True,
        help="Path to the directory containing Comet .txt result files",
    )(func)
    func = click.option(
        "--top_n_psms",
        "-n",
        type=int,
        help="If provided, only consider the top N PSMs per spectrum",
    )(func)
    return func


def get_protein_comet_counts(
    comet_results_dir: Path,
    top_n_psms: Optional[int] = None,
) -> Counter:
    # Get all Comet PSMs
    comet_txts = list(comet_results_dir.glob("*.txt"))
    psms = flatten_list_of_lists([CometPSM.from_txt(txt) for txt in comet_txts])

    # Filter PSMs to just those in the top N
    if top_n_psms is not None:
        psms = filter(lambda psm: psm.num <= top_n_psms, psms)

    # Get all the proteins in the resulting PSMs and then count the number of occurences for each protein
    all_comet_proteins = flatten_list_of_lists([psm.proteins for psm in psms])
    comet_protein_counts = Counter(all_comet_proteins)

    return comet_protein_counts


@click.command(
    name="plot-counts",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@common_arguments
def plot_protein_counts(
    comet_results_dir: Path,
    top_n_psms: Optional[int] = None,
):

    prot_counts = get_protein_comet_counts(
        comet_results_dir=comet_results_dir, top_n_psms=top_n_psms
    )

    df = pd.DataFrame(prot_counts.items(), columns=["prot", "count"])
    df.sort_values("count", inplace=True, ignore_index=True, ascending=False)

    fig, axs = fig_setup(2, 1)
    ax = axs[0]
    _ = ax.plot(
        df["count"],
        "bo",
        ms=1,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="protein index",
        ylabel="count",
        # title="Number of times protein appears in Comet run 1 data"
    )
    # finalize(axs)

    # fig, axs = fig_setup(1, 1)
    ax = axs[1]
    _ = ax.plot(
        df["count"],
        "bo",
        ms=1,
    )

    # Label the left-most N points
    num_pts_to_label = 8
    fs = 7
    for i in range(num_pts_to_label):
        _ = ax.text(
            i,
            df["count"].iloc[i],
            df["prot"].iloc[i],
            fontsize=fs,
            verticalalignment="bottom",
            horizontalalignment="left",
        )

    set_title_axes_labels(
        ax=ax,
        xlabel="protein index",
        ylabel="count",
        # title="Number of times protein appears in Comet run 1 data"
    )
    ax.set_xlim(left=-1, right=100)
    finalize(axs)
    params = f"topNpsms={top_n_psms}.cometDir={comet_results_dir.stem}"
    fig.suptitle(params)
    plt.tight_layout()
    fig.savefig(f"plots/protein_counts.{params}.png", dpi=200)


def get_most_common_proteins(protein_counts: Counter, top_n: int) -> List[str]:
    most_common_proteins = [
        protein_and_count[0] for protein_and_count in protein_counts.most_common(top_n)
    ]
    return most_common_proteins


@click.command(
    name="get-common-proteins",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
)
@common_arguments
@click.option(
    "--top_n_proteins",
    "-t",
    required=True,
    type=int,
    help="Get the top_n most common proteins",
)
def get_most_common_proteins_cli(
    comet_results_dir: Path,
    top_n_proteins: int,
    top_n_psms: Optional[int] = None,
):
    prot_counts = get_protein_comet_counts(
        comet_results_dir=comet_results_dir, top_n_psms=top_n_psms
    )
    most_common_proteins = get_most_common_proteins(
        protein_counts=prot_counts, top_n=top_n_proteins
    )
    for prot in most_common_proteins:
        print(prot)


@click.group()
def cli():
    pass


if __name__ == "__main__":
    cli.add_command(get_most_common_proteins_cli)
    cli.add_command(plot_protein_counts)
    cli()
