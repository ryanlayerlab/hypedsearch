import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import click
from matplotlib.axes import Axes
from pydantic import BaseModel

from src.constants import DEFAULT_Q_THRESHOLD, Q_VAL
from src.plot_utils import fig_setup, finalize, save_fig, set_title_axes_labels
from src.psm import CometPSM, get_high_confidence_psms
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
        psms = CometPSM.from_txt(txt=assign_confidence_file)
        # Filter PSMs by q_value threshold
        psms = filter(lambda psm: psm.q_value <= q_value_threshold, psms)
        return psms
    elif top_n_psms is not None:
        # Load PSMs from all *comet.target.txt files
        comet_target_files = comet_results_dir.glob("*.comet.target.txt")
        psms = []
        for comet_target_file in comet_target_files:
            psms.extend(CometPSM.from_txt(txt=comet_target_file))
        return psms
    else:
        raise RuntimeError("Either q_value_threshold or top_n_psms must be provided.")


def get_protein_counts_from_comet_psms(
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


class ProteinAbundance(BaseModel):
    protein_counts: Counter

    @classmethod
    def from_comet_txt(
        cls, txt: Union[str, Path], q_val_thresh: float = DEFAULT_Q_THRESHOLD
    ):
        psms = CometPSM.from_txt(txt=txt)
        return cls.from_comet_psms(psms=psms, q_threshold=q_val_thresh)

    @classmethod
    def from_comet_psms(
        cls, psms: List[CometPSM], q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> "ProteinAbundance":
        psms = get_high_confidence_psms(psms=psms, score=Q_VAL, threshold=q_threshold)
        all_comet_proteins = flatten_list_of_lists([psm.proteins for psm in psms])
        protein_counts = Counter(all_comet_proteins)
        return cls(protein_counts=protein_counts)

    def top_n_prots(
        self, n: int, with_cnts: bool = False
    ) -> Union[Set[str], Dict[str, int]]:
        most_common_proteins = {
            prot: cnt for prot, cnt in self.protein_counts.most_common(n)
        }
        if with_cnts:
            return most_common_proteins
        else:
            return set(most_common_proteins.keys())

    def plot(
        self, top_n_prots: Optional[int] = None, ax: Optional[Axes] = None
    ) -> Axes:
        # Define data
        items = sorted(self.protein_counts.items(), key=lambda x: x[1], reverse=True)
        if top_n_prots is not None:
            items = items[:top_n_prots]
        keys, values = zip(*items)

        # Plot
        if ax is None:
            fig, axs = fig_setup(h=8, w=10)
            ax = axs[0]
        ax.scatter(range(len(keys)), values)
        ax.set_xticks(range(len(keys)), keys, rotation=90, fontsize=8)
        set_title_axes_labels(
            ax=ax,
            # title="Protein counts",
            xlabel="Protein",
            ylabel="PSM counts",
        )
        finalize(ax)
        return ax

    def get_ab(self, protein: str) -> int:
        return self.protein_counts[protein]

    def get_rel_ab(self, protein: str) -> float:
        max_count = max(self.protein_counts.values())
        return self.protein_counts[protein] / max_count


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

    prot_counts = get_protein_counts_from_comet_psms(psms=psms)
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


def get_most_common_proteins(protein_counts: Counter, top_n: int) -> Set[str]:
    most_common_proteins = [
        protein_and_count[0] for protein_and_count in protein_counts.most_common(top_n)
    ]
    return set(most_common_proteins)


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
    to_json(data=prefix_counts_by_length, path=out_path)


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
