import logging

import seaborn as sns

from src.plot_utils import fig_setup, finalize

logger = logging.getLogger(__name__)


def plot_fpr_and_hybrid_score(df, min_hybrid_score):
    df.reset_index(drop=True, inplace=True)
    _, axs = fig_setup()
    _ = sns.scatterplot(
        x=(df.index.to_series() + 1),
        y=df["fpr"],
        ax=axs[0],
        # label="FPR",
        s=5,
        color="tab:blue",
    )
    axs[0].set_xlabel("Index in order of decreasing hybrid score")
    axs[0].set_ylabel(r"$\hat{FPR}$", color="tab:blue")
    axs[0].tick_params(axis="y", labelcolor="tab:blue")

    ax_copy = axs[0].twinx()
    # Plot horizontal line at minimum allowed hybrid score
    _ = ax_copy.axhline(
        y=min_hybrid_score,
        color="tab:green",
        linestyle="--",
        label="Min allowed hybrid score",
    )
    _ = sns.scatterplot(
        x=(df.index.to_series() + 1),
        y=df["hybrid_score"],
        ax=ax_copy,
        # label="hybrid score",
        color="tab:red",
        s=7,
    )
    ax_copy.set_ylabel("Hybrid score", color="tab:red")
    ax_copy.tick_params(axis="y", labelcolor="tab:red")
    finalize(axs)
