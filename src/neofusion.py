import logging
from typing import Optional

import pandas as pd
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


def neofusion_iteration(
    df: pd.DataFrame,
    max_fpr: float,
    min_score_delta: float,
    max_native_q_value: float,
    plot: bool = False,
) -> Optional[int]:
    # Filter to only those spectra with <hybrid score> - <native_score> >= min_score_delta
    tmp = df[(df["hybrid_score"] - df["native_score"]) > min_score_delta]
    tmp.reset_index(drop=True, inplace=True)
    logger.info(
        f"Of {df.shape[0]} spectra, {tmp.shape[0]} have hybrid score - native score >= {min_score_delta}"
    )

    # Classify native PSMs as gold-standards or not
    tmp["gold_standard"] = tmp["native_q"] < max_native_q_value
    logger.info(
        f"Of {tmp.shape[0]} spectra, {tmp['gold_standard'].sum()} are gold-standards with native q-value < {max_native_q_value}"
    )

    # Set minimimum allowed hybrid score
    min_hybrid_score = min(tmp[tmp["gold_standard"]]["native_score"])
    logger.info(
        f"Minimum allowed hybrid score (which is the minimum gold-standard native PSM score) is {min_hybrid_score}"
    )

    # Sort in decreasing hybrid score order
    tmp.sort_values(by="hybrid_score", ascending=False, inplace=True, ignore_index=True)
    # tmp.head()

    # For each row, compute the false positive rate and whether the hybrid score is below the threshold
    tmp["hybrid_score_too_low"] = tmp["hybrid_score"] < min_hybrid_score
    tmp["fpr"] = tmp["gold_standard"].cumsum() / (tmp.index.to_series() + 1)
    # tmp.head()

    # Plot
    if plot:
        plot_fpr_and_hybrid_score(df=tmp, min_hybrid_score=min_hybrid_score)

    # Find the point when the hybrid score dips below the minimum hybrid score
    idx = tmp.index[
        (tmp["hybrid_score"] > min_hybrid_score)
        & (tmp["hybrid_score"].shift(-1) <= min_hybrid_score)
    ].to_list()
    if len(idx) == 1:
        if tmp.iloc[idx[0]].fpr < max_fpr:
            true_positives = (~tmp.iloc[0 : idx[0] + 1]["gold_standard"]).sum()
            return true_positives
        else:
            return None
    else:
        logger.info(
            f"Didn't find a point when the hybrid score dipped below the minimum hybrid score. {idx}"
        )
        return None
