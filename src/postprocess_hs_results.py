import logging
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union
from venv import logger

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pydantic import BaseModel
from scipy.interpolate import PchipInterpolator

from src.constants import (
    DEFAULT_FPR,
    DEFAULT_JCT_LEN,
    DEFAULT_Q_RANGE,
    DEFAULT_Q_THRESHOLD,
    DEFAULT_SCORE_CHANGE_RANGE,
    NAT_DECOY,
    NAT_TARGET,
    NEOFUSION,
    Q_VAL,
    XCORR,
)
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import HypedsearchRunConfig
from src.mass_spectra import Spectrum
from src.neofusion import NeoFusionOutput, NeoFusionRunner
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.psm import CometPSM, CometRunAnalysis, create_xcorr_dists_plot
from src.utils import to_json

logger = logging.getLogger(__name__)
DEFAULT_XCORR_PLOT_NAME = "native_xcorr.png"
Q_ACCEPTANCE_METHOD = "q<="
NEOFUSION_METHOD = "NeoFusion"


def fit_xcorr_to_qval_interpolator(
    psms: Union[pd.DataFrame, List[CometPSM]],
    # ax: Optional[Axes] = None,
) -> PchipInterpolator:
    # Get data
    if isinstance(psms, pd.DataFrame):
        df = psms.copy()
    else:
        df = pd.DataFrame(
            [(psm.xcorr, psm.q_value) for psm in psms],
            columns=[XCORR, Q_VAL],
        )
    xy = df.drop_duplicates(subset=XCORR)
    xy.sort_values(by=XCORR, inplace=True)
    x = xy[XCORR].to_numpy()
    y = xy[Q_VAL].to_numpy()
    interpolator = PchipInterpolator(x, y)
    return interpolator


def add_qvalue_interpolator_to_xcorr_plot(
    q_value_psms: Union[List[CometPSM], pd.DataFrame],
    ax: Axes,
    q_threshold: Optional[float] = None,
    ylabel: str = "Native log10(q-value)",
):
    q_interpolator = fit_xcorr_to_qval_interpolator(
        psms=q_value_psms,
    )
    ax_copy = ax.twinx()
    xmin, xmax = ax.get_xlim()
    x_new = np.linspace(xmin, xmax, 500)
    _ = ax_copy.plot(x_new, q_interpolator(x_new), "r--", label="Native q-value")
    if q_threshold is not None:
        _ = ax_copy.axhline(
            y=q_threshold, color="red", linestyle="--", label=f"q={q_threshold}"
        )
    ax_copy.set_yscale("log")  # set y-axis to log10 scale
    ax_copy.set_ylabel(ylabel, color="tab:red")
    ax_copy.tick_params(axis="y", labelcolor="tab:red")


def create_native_xcorr_plot(
    native_target_psms: List[CometPSM],
    native_decoy_psms: List[CometPSM],
    assign_conf_psms: Optional[List[CometPSM]] = None,
):
    fig, axs = fig_setup()
    ax = axs[0]
    psms_by_type = {
        NAT_TARGET: [psm for psm in native_target_psms if psm.num == 1],
        NAT_DECOY: [psm for psm in native_decoy_psms if psm.num == 1],
    }
    _ = score_histogram(psms_by_type=psms_by_type, score=XCORR, ax=ax)
    if assign_conf_psms:
        add_qvalue_interpolator_to_xcorr_plot(
            q_value_psms=assign_conf_psms,
            ax=ax,
        )
    set_title_axes_labels(ax=ax, xlabel="xcorr", ylabel="Density")
    finalize(axs)
    return fig, axs


def score_histogram(
    psms_by_type: Dict[str, List[Any]],
    score: str,
    ax: Optional[Axes] = None,
) -> Axes:
    if ax is None:
        _, axs = fig_setup()
        ax = axs[0]
    for key, psms in psms_by_type.items():
        if isinstance(psms, pd.DataFrame):
            data = psms[score]
        else:
            try:
                data = [getattr(psm, score) for psm in psms]
            except:
                data = psms
        _ = sns.kdeplot(
            data,
            ax=ax,
            label=f"{key} (n = {len(data)})",
        )
    return ax


def group_hybrid_psms_by_junction(
    hybrid_psms: List[CometPSM],
    jct_len: int,
):
    hybrid_psms_with_multiple_explanations = [
        hybrid_psm for hybrid_psm in hybrid_psms if len(hybrid_psm.proteins) > 1
    ]
    uids_to_skip = set(psm.uid for psm in hybrid_psms_with_multiple_explanations)
    logger.info(
        f"Found {len(hybrid_psms_with_multiple_explanations)} ({len(uids_to_skip)} UIDs) hybrid PSMs with > 1 explanations. Ignoring them..."
    )
    jct_to_psms = defaultdict(list)
    for hybrid_psm in hybrid_psms:
        if hybrid_psm.uid in uids_to_skip:
            continue
        for prot in hybrid_psm.proteins:
            hy_pep = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
            jct_to_psms[hy_pep.get_junction_str(jct_len=jct_len)].append(hybrid_psm)
    return dict(jct_to_psms), hybrid_psms_with_multiple_explanations


@dataclass
class JunctionAnalysis:
    hybrid_psms: List[CometPSM]
    jct_len: int
    spectra: List[Spectrum]

    @cached_property
    def uid_to_spectrum(self) -> Dict[str, Spectrum]:
        return {spectrum.uid: spectrum for spectrum in self.spectra}

    @property
    def hybrid_psms_with_multiple_explanations(self) -> List[CometPSM]:
        return [
            hybrid_psm
            for hybrid_psm in self.hybrid_psms
            if len(hybrid_psm.proteins) > 1
        ]

    @property
    def jct_to_psms_map(self) -> Dict[str, List[CometPSM]]:
        uids_to_skip = set(
            psm.uid for psm in self.hybrid_psms_with_multiple_explanations
        )
        logger.info(
            f"Found {len(self.hybrid_psms_with_multiple_explanations)} ({len(uids_to_skip)} UIDs) hybrid PSMs with > 1 explanations. Ignoring them..."
        )
        jct_to_psms = defaultdict(list)
        for hybrid_psm in self.hybrid_psms:
            if hybrid_psm.uid in uids_to_skip:
                continue
            for prot in hybrid_psm.proteins:
                hy_pep = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
                jct_to_psms[hy_pep.get_junction_str(jct_len=self.jct_len)].append(
                    hybrid_psm
                )
        return dict(jct_to_psms)

    def create_jct_df(
        self,
    ) -> pd.DataFrame:
        rows = []
        for jct, psms in self.jct_to_psms_map.items():
            hybrid_seqs = []
            spectra_info = []
            for psm in psms:
                assert (
                    len(psm.proteins) == 1
                ), "Hybrid PSM has > 1 possible hybrid explanation. These should be filtered out by now so something is up."
                hy_pep = HybridPeptide.parse_hybrid_peptide_str(
                    hybrid_str=psm.proteins[0]
                )
                hybrid_seqs.append(hy_pep.hyphen_seq)
                spectrum = self.uid_to_spectrum[psm.uid]
                spectra_info.append(
                    f"{spectrum.uid}, seq={psm.seq}, m/z={round(spectrum.precursor_mz, 3)}, rt={round(spectrum.retention_time, 3)}, charge={spectrum.precursor_charge}, intensity={round(spectrum.precursor_intensity, 3)}"
                )
            # jct, hybrid_seqs, spectra_info
            seq_cntr = dict(Counter(hybrid_seqs))
            rows.append(
                (
                    jct,
                    len(psms),
                    len(seq_cntr),
                    seq_cntr,
                    spectra_info,
                )
            )
        df = pd.DataFrame(
            rows,
            columns=[
                "jct",
                "num_spectra_supporting",
                "num_uniq_seqs",
                "num_spectra_by_seq",
                "spectra_info",
            ],
        )
        df.sort_values(
            by="num_uniq_seqs", ascending=False, ignore_index=True, inplace=True
        )
        return df

    @cached_property
    def junction_df(
        self,
    ) -> pd.DataFrame:
        rows = []
        for jct, psms in self.jct_to_psms_map.items():
            hybrid_seqs = []
            for psm in psms:
                assert (
                    len(psm.proteins) == 1
                ), "Hybrid PSM has > 1 possible hybrid explanation. These should be filtered out by now so something is up."
                hy_pep = HybridPeptide.parse_hybrid_peptide_str(
                    hybrid_str=psm.proteins[0]
                )
                hybrid_seqs.append(hy_pep.hyphen_seq)
            rows.append(
                (
                    jct,
                    len(psms),
                    list(set(psm.uid for psm in psms)),
                    np.mean([psm.xcorr for psm in psms]),
                    max([psm.xcorr for psm in psms]),
                    min([psm.q_value for psm in psms]),
                    dict(Counter(hybrid_seqs)),
                )
            )
        df = pd.DataFrame(
            rows,
            columns=[
                "jct",
                "num_psm_supporting",
                "spectra_supporting",
                "mean_xcorr",
                "max_xcorr",
                "min_q_value",
                "psm_seq_cnter",
            ],
        )
        df["num_uniq_seqs"] = df.psm_seq_cnter.apply(lambda cnter: len(cnter.keys()))
        df.sort_values(
            by="num_uniq_seqs", ascending=False, ignore_index=True, inplace=True
        )
        return df

    def create_jct_num_supporting_spectra_vs_num_uniq_seqs(
        self,
        out_path: Optional[str | Path] = None,
        title: str = "",
    ):
        fig, axs = fig_setup()
        sns.scatterplot(
            data=self.junction_df,
            x="num_psm_supporting",
            y="num_uniq_seqs",
            ax=axs[0],
            s=7,
            # label=f"n = {len(self.junction_df)} junctions",
        )
        counts = (
            self.junction_df.groupby(["num_psm_supporting", "num_uniq_seqs"])
            .size()
            .reset_index(name="n")  # n = number of rows with that (x, y)
        )
        for _, row in counts.iterrows():
            axs[0].text(
                row["num_psm_supporting"],  # small x offset
                row["num_uniq_seqs"],  # small y offset
                str(row["n"]),
                fontsize=9,
                ha="left",
                va="bottom",
            )
        set_title_axes_labels(
            ax=axs[0],
            title=title,
            xlabel="Number of PSMs/spectra supporting junction",
            ylabel="Number of unique peptide\nsequences supporting junction",
        )

        finalize(axs)
        if out_path:
            save_fig(path=out_path)


def aggregate_junction_dfs_over_samples(
    top_dir: Union[str, Path],
):
    top_dir = Path(top_dir)
    dfs = []
    for sample_dir in top_dir.glob("*_hybrid_results"):
        mzml_name = sample_dir.name[: -len("_hybrid_results")]
        csv = sample_dir / "NeoFusion_accepted_hybrid_junctions.csv"
        if not csv.exists():
            logger.info(
                f"Junction CSV for MZML {mzml_name} not found at {csv}. Skipping..."
            )
            continue
        df = pd.read_csv(sample_dir / "NeoFusion_accepted_hybrid_junctions.csv")
        df["mzml"] = mzml_name
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    df.sort_values(by="num_uniq_seqs", ascending=False, ignore_index=True, inplace=True)
    return df


@dataclass
class SpectrumPSMs:
    native_targets: List[CometPSM]
    native_decoys: List[CometPSM]
    hybrid_targets: Optional[List[CometPSM]] = None

    def __post_init__(self):
        uid = self.native_targets[0].uid
        for psm in self.native_targets:
            assert (
                psm.uid == uid
            ), "All native target PSMs for a spectrum must have the same UID"
        for psm in self.native_decoys:
            assert (
                psm.uid == uid
            ), "All native decoy PSMs for a spectrum must have the same UID"
        for psm in self.hybrid_targets:
            assert (
                psm.uid == uid
            ), "All hybrid target PSMs for a spectrum must have the same UID"

    @property
    def top_native_target(self) -> CometPSM:
        psms = self.get_top_psms(psms=self.native_targets)
        assert (
            len(psms) == 1
        ), f"Expected exactly 1 top native target PSM for spectrum {self.native_targets[0].uid}, but found {len(psms)}"
        return psms[0]

    @property
    def top_native_decoy(self) -> CometPSM:
        psms = self.get_top_psms(psms=self.native_decoys)
        assert (
            len(psms) == 1
        ), f"Expected exactly 1 top native decoy PSM for spectrum {self.native_decoys[0].uid}, but found {len(psms)}"
        return psms[0]

    @property
    def top_hybrid_target(self) -> Optional[CometPSM]:
        if self.hybrid_targets is None:
            return None
        psms = self.get_top_psms(psms=self.hybrid_targets)
        assert (
            len(psms) <= 1
        ), f"Expected at most 1 top hybrid target PSM for spectrum {self.hybrid_targets[0].uid}, but found {len(psms)}"
        return psms[0]

    @staticmethod
    def get_top_psms(psms: List[CometPSM]) -> List[CometPSM]:
        return [psm for psm in psms if psm.num == 1]

    def to_dict(self):
        data = {
            "uid": self.top_native_target.uid,
            "nt_xcorr": self.top_native_target.xcorr,
            "nt_q": self.top_native_target.q_value,
            "nd_xcorr": self.top_native_decoy.xcorr,
            "nd_q": self.top_native_decoy.q_value,
            "ht_xcorr": (
                self.top_hybrid_target.xcorr if self.top_hybrid_target else None
            ),
            "ht_q": self.top_hybrid_target.q_value if self.top_hybrid_target else None,
        }
        return data


def plot_xcorr_of_accepted_vs_not_accepted_hybrids(
    native_targets: List[CometPSM],
    accepted_hybrids: List[CometPSM],
    hybrid_targets: List[CometPSM],
    ax: Optional[Axes] = None,
) -> Tuple[Figure, List[Axes]]:
    if ax is None:
        fig, axs = fig_setup()
        ax = axs[0]
    accepted_uids = set(psm.uid for psm in accepted_hybrids)
    native_targets = {psm.uid: psm for psm in native_targets}
    hybrid_targets = {psm.uid: psm for psm in hybrid_targets}
    df = pd.DataFrame(
        {
            "uid": uid,
            "native_xcorr": native_targets[uid].xcorr if uid in native_targets else 0,
            "hybrid_xcorr": hybrid_targets[uid].xcorr if uid in hybrid_targets else 0,
        }
        for uid in set(native_targets.keys()) & set(hybrid_targets.keys())
    )
    df["accepted"] = df["uid"].apply(lambda uid: uid in accepted_uids)
    s = 7
    sns.scatterplot(
        data=df[~df.accepted],
        x="native_xcorr",
        y="hybrid_xcorr",
        s=s,
        marker="o",
        color="blue",
        # label=f"All spectra with a native and hybrid PSM (n={len(native_and_hybrid_psms)})",
        ax=ax,
    )
    sns.scatterplot(
        data=df[df.accepted],
        x="native_xcorr",
        y="hybrid_xcorr",
        s=s,
        marker="X",
        color="red",
        label=f"Accepted hybrid PSMs (n={len(accepted_hybrids)})",
        ax=ax,
    )
    plot_line(ax=ax, label="y=x")
    set_title_axes_labels(
        ax=ax,
        xlabel="Native target xcorr",
        ylabel="Hybrid target xcorr",
    )
    finalize(ax)


class HybridJunction(BaseModel):
    jct: str
    supporting_psms: List[CometPSM]
    supporting_spectra: List[Spectrum]

    @property
    def max_precursor_intensity(self) -> float:
        return max(spectrum.precursor_intensity for spectrum in self.supporting_spectra)


class NativeToHybridComparison(BaseModel):
    native_run: CometRunAnalysis
    hybrid_run: CometRunAnalysis
    spectra: List[Spectrum]
    name: Optional[str] = None

    @property
    def uid_to_spectrum_psms(self) -> Dict[str, SpectrumPSMs]:
        uid_to_spectrum_psms = defaultdict(
            lambda: {"native_targets": [], "native_decoys": [], "hybrid_targets": []}
        )

        for psm in self.native_run.top_targets:
            uid_to_spectrum_psms[psm.uid]["native_targets"].append(psm)
        for psm in self.native_run.top_decoys:
            uid_to_spectrum_psms[psm.uid]["native_decoys"].append(psm)
        for psm in self.hybrid_run.top_targets:
            uid_to_spectrum_psms[psm.uid]["hybrid_targets"].append(psm)

        for uid, psms_as_dict in uid_to_spectrum_psms.items():
            uid_to_spectrum_psms[uid] = SpectrumPSMs(**psms_as_dict)

        return dict(uid_to_spectrum_psms)

    @cached_property
    def spectra_df(self) -> pd.DataFrame:
        return Spectrum.to_df(spectra=self.spectra)

    @classmethod
    def from_config(cls, config: str | Path | HypedsearchRunConfig):
        if isinstance(config, (str, Path)):
            config = HypedsearchRunConfig.from_json(path=config)
        native_run = CometRunAnalysis(
            targets=config.native_target_psms,
            decoys=config.native_decoy_psms,
            assign_conf=config.native_assign_confidence_psms,
            interpolate=True,
        )
        hybrid_run = CometRunAnalysis(
            targets=config.hybrid_target_psms,
            assign_conf=native_run.assign_conf,
            interpolate=True,
        )
        return cls(
            name=config.name,
            native_run=native_run,
            hybrid_run=hybrid_run,
        )

    def get_protein_abundance_df(
        self, q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> pd.DataFrame:
        return self.native_run.get_protein_abundance(q_threshold=q_threshold).df

    def xcorr_target_and_decoy_distributions(
        self, title: Optional[str] = None, ax: Optional[Axes] = None
    ) -> Axes:
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        create_xcorr_dists_plot(
            psms_by_type={
                "Top native targets": self.native_run.top_targets,
                "Top native decoys": self.native_run.top_decoys,
                "Top hybrid targets": self.hybrid_run.top_targets,
            },
            q_interpolating_psms=(
                self.native_run.assign_conf
                if len(self.native_run.assign_conf) > 0
                else None
            ),
            title=title,
            ax=ax,
        )

    def xcorr_accepted_vs_not_accepted_plot(
        self, accepted_hybrids: Optional[List[CometPSM]] = None
    ):
        if accepted_hybrids is None:
            logger.info(
                "No accepted hybrids provided so accepting hybrids via NeoFusion using default settings"
            )
            accepted_hybrids = self.run_neofusion().accepted_hybrids
        plot_xcorr_of_accepted_vs_not_accepted_hybrids(
            native_targets=self.native_run.top_targets,
            accepted_hybrids=accepted_hybrids,
            hybrid_targets=self.hybrid_run.top_targets,
        )

    def run_neofusion(self) -> NeoFusionOutput:
        neofusion = NeoFusionRunner(
            native_targets=self.native_run.top_targets,
            hybrid_targets=self.hybrid_run.top_targets,
        )
        neofusion_output = neofusion.run_neofusion()
        neofusion.plot_neofusion_true_positive_data(
            neofusion_iterations=neofusion_output.iterations
        )
        return neofusion_output

    def get_spectrum_psms(self, uid: str) -> SpectrumPSMs:
        return self.uid_to_spectrum_psms[uid]

    def to_df(self):
        return pd.DataFrame(
            spectrum_psms.to_dict()
            for spectrum_psms in self.uid_to_spectrum_psms.values()
        )


def interpolate_psm_q_values_from_assign_conf(
    interpolating_psms: List[CometPSM], q_valueless_psms: List[CometPSM]
) -> List[CometPSM]:
    interpolator = fit_xcorr_to_qval_interpolator(psms=interpolating_psms)
    for psm in q_valueless_psms:
        psm.q_value = float(interpolator(psm.xcorr))
    return q_valueless_psms


def get_native_and_hybrid_runs_via_config(
    hs_config: Union[str, Path, HypedsearchRunConfig],
):
    if isinstance(hs_config, (str, Path)):
        hs_config = HypedsearchRunConfig.from_json(path=hs_config)


def process_native_and_hybrid_runs_via_config(
    hs_config: str | Path | HypedsearchRunConfig,
    q_threshold: float = DEFAULT_Q_THRESHOLD,
    perform_native_analysis: bool = True,
    acceptance_method: str = NEOFUSION,
    jct_len: int = DEFAULT_JCT_LEN,
):
    # Load config and make directories
    if isinstance(hs_config, (str, Path)):
        hs_config = HypedsearchRunConfig.from_json(path=hs_config)
    native_results_dir = (
        hs_config.parent_output_dir / f"results/{hs_config.name}_native_results"
    )
    native_results_dir.mkdir(parents=True, exist_ok=True)
    hybrid_results_dir = (
        hs_config.parent_output_dir / f"results/{hs_config.name}_hybrid_results"
    )
    hybrid_results_dir.mkdir(parents=True, exist_ok=True)

    # Run analysis
    native_run = process_native_run_from_hs_config(
        config=hs_config,
        q_threshold=q_threshold,
        analyze=perform_native_analysis,
        out_dir=native_results_dir,
    )
    try:
        hybrid_run = CometRunAnalysis(
            targets=hs_config.hybrid_target_psms, assign_conf=native_run.assign_conf
        )
        comp = NativeToHybridComparison(
            name=hs_config.name,
            q_threshold=q_threshold,
            native_run=native_run,
            hybrid_run=hybrid_run,
        )
        comp.create_xcorr_plot(
            out_path=hybrid_results_dir / "xcorr.png",
        )
        _, accepted_hybrids = comp.accept_hybrids_via_neofusion(
            best_iteration_out_path=hybrid_results_dir / "best_neofusion_iteration.json"
        )
        NativeToHybridComparison.accepted_vs_not_accepted_plot(
            uid_to_top_native_target=native_run.uid_to_top_target,
            uid_to_top_hybrid_target=hybrid_run.uid_to_top_target,
            accepted_hybrid_uids=set(psm.uid for psm in accepted_hybrids),
            title=f"{hs_config.name} spectra",
            out_path=hybrid_results_dir
            / f"{acceptance_method}_accepted_vs_not_accepted.png",
        )
        jct_analyzer = JunctionAnalysis(hybrid_psms=accepted_hybrids, jct_len=jct_len)
        CometPSM.save_psms_to_json(
            psms=jct_analyzer.hybrid_psms_with_multiple_explanations,
            path=hybrid_results_dir
            / f"{acceptance_method}_accepted_hybrid_psms_with_multiple_explanations.json",
        )
        jct_analyzer.junction_df.to_csv(
            hybrid_results_dir / f"{acceptance_method}_accepted_hybrid_junctions.csv",
            index=False,
        )
        jct_analyzer.create_jct_num_supporting_spectra_vs_num_uniq_seqs(
            out_path=hybrid_results_dir
            / f"{acceptance_method}_accepted_hybrid_jct_support.png",
            title=f"{hs_config.name}\n{acceptance_method} accepted hybrids",
        )
        return accepted_hybrids, jct_analyzer
    except:
        logger.info("Problem processing the hybrid run. Continuing")
