import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import click
import pandas as pd
import seaborn as sns

from src.comet_utils import CometPSM
from src.plot_utils import fig_setup, finalize, save_fig, set_title_axes_labels
from src.utils import PathType, flatten_list_of_lists, setup_logger, to_json

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


def load_comet_psms(
    comet_results_dir: Path,
    q_value_threshold: Optional[float] = None,
    top_n_psms: Optional[int] = None,
):
    if q_value_threshold is not None:
        # Load PSMs from the assign-confidence.target.txt file
        assign_confidence_file = comet_results_dir / "assign-confidence.target.txt"
        psms = CometPSM.from_txt(txt_path=assign_confidence_file)
        # Filter PSMs by q_value threshold
        psms = filter(lambda psm: psm.q_value <= q_value_threshold, psms)
        return psms
    elif top_n_psms is not None:
        # Load PSMs from all *comet.target.txt files
        comet_target_files = comet_results_dir.glob("*.comet.target.txt")
        psms = []
        for comet_target_file in comet_target_files:
            psms.extend(CometPSM.from_txt(txt_path=comet_target_file))
        return psms
    else:
        raise RuntimeError("Either q_value_threshold or top_n_psms must be provided.")


def get_protein_counts_from_comet_results(
    psms: List[CometPSM],
) -> Counter[str, int]:
    """ """
    all_comet_proteins = flatten_list_of_lists([psm.proteins for psm in psms])
    comet_protein_counts = Counter(all_comet_proteins)
    return comet_protein_counts


def get_prefix_counts_by_length(
    seqs: List[str],
) -> Dict[int, Dict[str, int]]:
    """
    Get the number of times each unique k-mer in the given sequences appears across all sequences
    and organize the counts by k-mer length.
    """
    all_prefixes = flatten_list_of_lists(
        [[seq[:i] for i in range(1, len(seq) + 1)] for seq in seqs]
    )
    prefix_counts = Counter(all_prefixes)
    prefix_counts_by_length = defaultdict(dict)
    for prefix, count in prefix_counts.items():
        prefix_counts_by_length[len(prefix)][prefix] = count
    return prefix_counts_by_length


def plot_protein_counts(
    prot_counts: Dict[str, int],
    ms: int = 10,
    top_n_prots: int = 10,
):
    df = pd.DataFrame(prot_counts.items(), columns=["prot", "count"])
    df.sort_values("count", inplace=True, ignore_index=True, ascending=False)

    fig, axs = fig_setup(1, 2)

    # 1st plot
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
    )

    # 2nd plot - same as first but only showing the most 100 abundant proteins
    ax = axs[1]
    _ = ax.plot(
        df["count"],
        "bo",
        ms=1,
    )

    # Label the left-most N points
    colors = sns.color_palette("hsv", top_n_prots)
    for idx in range(top_n_prots):
        # _ = ax.text(
        _ = ax.scatter(
            idx,
            df["count"].iloc[idx],
            label=df["prot"].iloc[idx],
            color=colors[idx],
            s=ms,
        )

    set_title_axes_labels(
        ax=ax,
        xlabel="protein index",
    )
    _ = ax.set_xlim(left=-1, right=100)
    finalize(axs)
    return fig, axs


def get_and_plot_most_common_proteins(
    comet_results_dir: Path,
    out_path: Path,
    top_n_proteins: int = 10,
    q_value_threshold: Optional[float] = None,
    top_n_psms: Optional[int] = None,
):
    psms = load_comet_psms(
        comet_results_dir=comet_results_dir,
        q_value_threshold=q_value_threshold,
        top_n_psms=top_n_psms,
    )

    prot_counts = get_protein_counts_from_comet_results(psms=psms)
    most_common_proteins = get_most_common_proteins(
        protein_counts=prot_counts, top_n=top_n_proteins
    )

    # Save top proteins
    with open(out_path, "w") as f:
        for prot in most_common_proteins:
            print(prot)
            f.write(f"{prot}\n")

    # Plot protein abundances
    fig, _ = plot_protein_counts(
        prot_counts=prot_counts,
    )
    save_fig(
        path=out_path.parent / "protein_abundances.png",
        fig=fig,
    )


def get_most_common_proteins(protein_counts: Counter, top_n: int) -> List[str]:
    most_common_proteins = [
        protein_and_count[0] for protein_and_count in protein_counts.most_common(top_n)
    ]
    return most_common_proteins


@click.command(
    name="protein-abundances",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help=(
        "Get the most common proteins from Comet results.\n\n"
        "\t - If q_value_threshold is provided, then it's assumed there's an assign-confidence.target.txt "
        "file in the comet_results_dir and the PSMs will be loaded from that file and "
        "filtered to those with q_value <= q_value_threshold.\n\n"
        "\t - If top_n_psms is provided, then all PSMs will be loaded from all *.comet.target.txt "
        "files in comet_results_dir and filtered to just the top_n_psms per spectrum.\n\n"
        "The top_n_proteins will be saved in a TXT as a new-line separated list in the given out_path.\n\n"
        "A plot of the protein counts will be saved in the same directory as the TXT."
    ),
)
@click.option(
    "--comet_results_dir",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the directory containing Comet .txt result files",
)
@click.option(
    "--q_value_threshold",
    "-q",
    type=float,
    help=("If provided, only consider PSMs with q_value <= q_value_threshold. "),
)
@click.option(
    "--top_n_psms",
    "-n",
    type=int,
    help="If provided, only consider the top N PSMs per spectrum",
)
@click.option(
    "--top_n_proteins",
    "-t",
    required=True,
    type=int,
    help="Get the top_n most common proteins",
)
@click.option(
    "--out_path",
    "-o",
    type=PathType(),
    required=True,
    help=("The top_n_proteins will be saved here as a new-line separated list."),
)
def cli_get_and_plot_most_common_proteins(
    comet_results_dir: Path,
    top_n_proteins: int,
    out_path: Path,
    q_value_threshold: Optional[float] = None,
    top_n_psms: Optional[int] = None,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    get_and_plot_most_common_proteins(
        comet_results_dir=comet_results_dir,
        top_n_proteins=top_n_proteins,
        q_value_threshold=q_value_threshold,
        top_n_psms=top_n_psms,
        out_path=out_path,
    )


def get_and_save_prefix_counts_by_length(
    comet_results_dir: Path,
    out_path: Path,
    q_value_threshold: Optional[float] = None,
    top_n_psms: Optional[int] = None,
):
    psms = load_comet_psms(
        comet_results_dir=comet_results_dir,
        q_value_threshold=q_value_threshold,
        top_n_psms=top_n_psms,
    )
    prefix_counts_by_length = get_prefix_counts_by_length(
        seqs=[psm.seq for psm in psms]
    )
    to_json(data=prefix_counts_by_length, out_path=out_path)


@click.command(
    name="prefix-abundances",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help=(
        "Get the abundances/counts of all the prefixes of the PSM sequences.\n\n"
        "\t - If q_value_threshold is provided, then it's assumed there's an assign-confidence.target.txt "
        "file in the comet_results_dir and the PSMs will be loaded from that file and "
        "filtered to those with q_value <= q_value_threshold.\n\n"
        "\t - If top_n_psms is provided, then all PSMs will be loaded from all *.comet.target.txt "
        "files in comet_results_dir and filtered to just the top_n_psms per spectrum.\n\n"
        "The prefix counts will be saved in a JSON file in the given out_path.\n\n"
    ),
)
@click.option(
    "--comet_results_dir",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the directory containing Comet .txt result files",
)
@click.option(
    "--q_value_threshold",
    "-q",
    type=float,
    help=("If provided, only consider PSMs with q_value <= q_value_threshold. "),
)
@click.option(
    "--top_n_psms",
    "-n",
    type=int,
    help="If provided, only consider the top N PSMs per spectrum",
)
@click.option(
    "--out_path",
    "-o",
    type=PathType(),
    required=True,
    help=("The prefix counts will be saved here as a JSON."),
)
def cli_get_prefix_counts_by_length(
    comet_results_dir: Path,
    out_path: Path,
    q_value_threshold: Optional[float] = None,
    top_n_psms: Optional[int] = None,
):
    get_and_save_prefix_counts_by_length(
        comet_results_dir=comet_results_dir,
        out_path=out_path,
        q_value_threshold=q_value_threshold,
        top_n_psms=top_n_psms,
    )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_get_and_plot_most_common_proteins)
    cli.add_command(cli_get_prefix_counts_by_length)
    cli()
    # cli_get_and_plot_most_common_proteins()
