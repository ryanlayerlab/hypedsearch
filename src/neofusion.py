import logging
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from typing import Dict, List, Optional, Tuple

import pandas as pd
from matplotlib.axes import Axes

from src.constants import DEFAULT_FPR, DEFAULT_Q_RANGE, DEFAULT_SCORE_CHANGE_RANGE
from src.plot_utils import finalize, plot_sorted_1d_data, set_title_axes_labels
from src.psm import CometPSM

logger = logging.getLogger(__name__)

HYBRID_SCORE = "hybrid_score"
NATIVE_SCORE = "native_score"


@dataclass
class NeoFusionIteration:
    q: float
    min_score_delta: float
    fpr: float
    tp: int
    min_hybrid_score: float
    accepted_hybrid_psm_spectrum_uids: List[str]

    @property
    def num_accepted(self) -> int:
        return len(self.accepted_hybrid_psm_spectrum_uids)

    @property
    def fp(self) -> int:
        return self.num_accepted - self.tp

    @property
    def info_dict(self) -> dict:
        return {
            "q": round(self.q, 3),
            "min_score_delta": round(self.min_score_delta, 3),
            "min_hybrid_score": round(self.min_hybrid_score, 3),
            "fpr": round(self.fpr, 3),
            "tp": self.tp,
            "fp": self.fp,
            "num_accepted_psms": self.num_accepted,
        }

    @property
    def param_str(self) -> str:
        return f"q{self.q}_delta{self.min_score_delta}_fpr{self.fpr}_minHybridScore{self.min_hybrid_score}"

    # def xcorr_plot(self, )


@dataclass
class NeoFusionOutput:
    iterations: List[NeoFusionIteration]
    best_iteration: NeoFusionIteration
    accepted_hybrids: List[CometPSM]


@dataclass
class NeoFusionRunner:
    native_targets: List[CometPSM]
    hybrid_targets: List[CometPSM]
    q_vals: List[float] = field(default_factory=lambda: DEFAULT_Q_RANGE.copy())
    score_deltas: List[float] = field(
        default_factory=lambda: DEFAULT_SCORE_CHANGE_RANGE.copy()
    )
    fpr_threshold: float = DEFAULT_FPR

    @cached_property
    def native_df(self) -> pd.DataFrame:
        return pd.DataFrame(psm.to_dict() for psm in self.native_targets)

    @cached_property
    def hybrid_df(self) -> pd.DataFrame:
        return pd.DataFrame(psm.to_dict() for psm in self.hybrid_targets)

    @cached_property
    def df(self) -> pd.DataFrame:
        return self.native_df.merge(
            right=self.hybrid_df, on="uid", suffixes=("_native", "_hybrid")
        )

    @cached_property
    def uid_to_native_target(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.native_targets}

    @cached_property
    def uid_to_hybrid_target(self) -> Dict[str, CometPSM]:
        return {psm.uid: psm for psm in self.hybrid_targets}

    @staticmethod
    def create_neofusion_df(
        native_assign_conf: Dict[str, CometPSM],
        top_hybrid_targets: Dict[str, CometPSM],
    ) -> pd.DataFrame:
        # Create dataframe for NeoFusion analysis
        df = pd.DataFrame(
            [
                [
                    spectrum_uid,
                    native_assign_conf[spectrum_uid].xcorr,
                    top_hybrid_targets[spectrum_uid].xcorr,
                    native_assign_conf[spectrum_uid].q_value,
                ]
                for spectrum_uid in set(top_hybrid_targets.keys()).intersection(
                    native_assign_conf.keys()
                )
            ],
            columns=[
                "uid",
                "n_score",
                "h_score",
                "n_q",
            ],
        )
        df["delta"] = df.h_score - df.n_score
        return df

    @staticmethod
    def neofusion_iteration(
        df: pd.DataFrame,
        q_val: float,
        score_delta: float,
        fpr_thresh: float,
    ) -> Optional[NeoFusionIteration]:
        for colm in ["n_score", "h_score", "n_q", "delta"]:
            assert colm in df.columns
        neo_df = deepcopy(df)

        # Remove rows where hybrid score isn't high enough compared to native score
        neo_df = neo_df[neo_df.delta >= score_delta].copy()

        # Set which native PSMs are "gold-standard"
        neo_df["gold"] = neo_df.n_q <= q_val
        neo_df["fp"] = neo_df["gold"].copy()

        # Min hybrid score is lowest gold-standard native score
        min_hybrid_score = neo_df[neo_df.gold].n_score.min()
        if pd.isna(min_hybrid_score):
            return None

        # Sort in descending order
        neo_df.sort_values(by="h_score", ascending=False, inplace=True)
        neo_df.reset_index(drop=True, inplace=True)

        # Iterate through rows and, for each row i, find number number of false and true positives
        # in rows 1, 2, ..., i.
        fpr_colm = []
        tp_colm = []
        for row_idx, row in neo_df.iterrows():
            if row.h_score < min_hybrid_score:
                fpr_colm.append(None)
                tp_colm.append(None)
                continue

            tmp = neo_df.iloc[: row_idx + 1, :]
            fpr_colm.append(tmp.fp.sum() / tmp.shape[0])
            tp_colm.append((~tmp.fp).sum())

        neo_df["fpr"] = fpr_colm
        neo_df["tp"] = tp_colm
        if neo_df[neo_df["fpr"] < fpr_thresh].shape[0] > 0:
            # Get row that maximizes the number of true positives
            try:
                tmp = neo_df[neo_df["fpr"] < fpr_thresh]
                tp_maximizing_idx = tmp.tp.idxmax()
                tp_maximizing_row = tmp.loc[tp_maximizing_idx]
                accepted_psm = neo_df.uid.iloc[: tp_maximizing_idx + 1].tolist()

                return NeoFusionIteration(
                    q=float(q_val),
                    min_score_delta=float(score_delta),
                    fpr=float(tp_maximizing_row.fpr),
                    tp=int(tp_maximizing_row.tp),
                    min_hybrid_score=float(min_hybrid_score),
                    accepted_hybrid_psm_spectrum_uids=accepted_psm,
                )
            except:
                logger.debug(f"Issue with q_val={q_val}, score_delta={score_delta}")
        return None

    def xcorr_plot(self, neofusion_iteration: NeoFusionIteration):
        df["hybrid_accepted"] = df.uid.isin(
            set(psm.uid for psm in neofusion_out.accepted_hybrids)
        )

        fig, axs = fig_setup(w=8)
        ax = axs[0]
        tmp = df
        _ = sns.scatterplot(
            x=tmp.xcorr_native,
            y=tmp.xcorr_hybrid,
            hue=df.native_accepted,
            s=5,
            # color="black",
            ax=ax,
        )
        tmp = df[df.hybrid_accepted]
        _ = sns.scatterplot(
            x=tmp.xcorr_native,
            y=tmp.xcorr_hybrid,
            s=8,
            marker=(8, 1, 0),
            color="red",
            ax=ax,
        )
        plot_line(
            ax=ax,
            b=neofusion_iteration.min_score_delta,
            label=f"min score delta ({round(neofusion_iteration.best_iteration.min_score_delta, 3)})",
        )
        plot_line(
            ax=ax,
            m=0,
            b=neofusion_iteration.min_hybrid_score,
            label=f"min hybrid score ({round(neofusion_iteration.best_iteration.min_hybrid_score, 3)})",
        )
        finalize(axs)
        ax.legend(
            loc="center left",
            bbox_to_anchor=(0.95, 0.8),
            borderaxespad=0,
        )

    def run_neofusion(
        self,
    ) -> NeoFusionOutput:
        logger.debug("Running NeoFusion")
        df = self.create_neofusion_df(
            native_assign_conf=self.uid_to_native_target,
            top_hybrid_targets=self.uid_to_hybrid_target,
        )

        neofusion_iterations = []
        for q_val in self.q_vals:
            for min_score_delta in self.score_deltas:
                result = self.neofusion_iteration(
                    df=df,
                    q_val=q_val,
                    score_delta=min_score_delta,
                    fpr_thresh=self.fpr_threshold,
                )
                if result is not None:
                    neofusion_iterations.append(result)
        best_iteration, accepted_hybrids = self.select_hybrid_psms_from_best_iteration(
            neofusion_iterations=neofusion_iterations
        )

        return NeoFusionOutput(
            iterations=neofusion_iterations,
            best_iteration=best_iteration,
            accepted_hybrids=accepted_hybrids,
        )

    @staticmethod
    def plot_neofusion_true_positive_data(
        neofusion_iterations: List[NeoFusionIteration],
        title: str = "",
    ) -> Axes:
        data = {res.param_str: res.tp for res in neofusion_iterations}
        ax = plot_sorted_1d_data(data=data)
        set_title_axes_labels(
            ax=ax,
            title=title,
            xlabel="NeoFusion parameters\n(sorted in decreasing TP order)",
            ylabel='"True positives (TPs)"',
        )
        finalize(ax)
        return ax

    def select_hybrid_psms_from_best_iteration(
        self,
        neofusion_iterations: List[NeoFusionIteration],
    ) -> Tuple[NeoFusionIteration, List[CometPSM]]:
        best_iteration = max(neofusion_iterations, key=lambda x: x.tp)
        accepted_hybrids = [
            self.uid_to_hybrid_target[spectrum_uid]
            for spectrum_uid in best_iteration.accepted_hybrid_psm_spectrum_uids
        ]
        return best_iteration, accepted_hybrids
