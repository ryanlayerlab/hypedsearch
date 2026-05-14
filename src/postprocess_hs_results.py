import logging
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple, Union
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
    DEFAULT_PEAK_TO_ION_PPM_TOL,
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
    create_joint_plot,
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.psm import (
    CometPSM,
    CometRunAnalysis,
    PeptideSeqSpectrumComparer,
    convert_comet_psms_to_custom_psms,
    create_xcorr_dists_plot,
    hybrid_psm_plot,
    spectrum_peptide_plot,
)
from src.utils import to_json

logger = logging.getLogger(__name__)
DEFAULT_XCORR_PLOT_NAME = "native_xcorr.png"
Q_ACCEPTANCE_METHOD = "q<="
NEOFUSION_METHOD = "NeoFusion"
NTS, NDS, HTS = "native_targets", "native_decoys", "hybrid_targets"
PSM_TYPES = [NTS, NDS, HTS]
PSM_TYPE_TO_CSV_NAME = {psm_type: f"{psm_type}_psms.csv" for psm_type in PSM_TYPES}


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


def get_jct_len_from_jct(jct: str):
    left, right = jct.split("-")
    left = left.split(")")[-1]
    right = right.split("(")[0]
    assert len(left) == len(
        right
    ), f"Expected left and right sides of junction to be same length, but got {left} and {right}"
    return len(left)


def make_sure_hybrid_psm_corresponds_to_jct(psm: CometPSM, jct: str) -> bool:
    for prot in psm.proteins:
        hy_pep = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
        if hy_pep.get_junction_str(jct_len=get_jct_len_from_jct(jct=jct)) == jct:
            return True
    return False


def get_hybrid_seq_supporting_jct(psm: CometPSM, jct: str) -> Optional[str]:
    for prot in psm.proteins:
        hy_pep = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
        if hy_pep.get_junction_str(jct_len=get_jct_len_from_jct(jct=jct)) == jct:
            return hy_pep.hyphen_seq
    return None


@dataclass
class SupportedHybridJunction:
    psms: List[CometPSM]
    jct: str

    def __post_init__(self):
        # Make sure each PSM has a protein that corresponds to the junction
        for psm in self.psms:
            assert make_sure_hybrid_psm_corresponds_to_jct(
                psm=psm, jct=self.jct
            ), f"PSM with uid {psm.uid} does not have a protein that corresponds to junction {self.jct}"
        # if self.spectra is not None:
        #     assert len(self.psms) == len(
        #         self.spectra
        #     ), "If spectra are provided, there must be the same number of spectra and PSMs"

    @property
    def hybrid_seqs(self) -> List[str]:
        hybrid_seqs = [
            get_hybrid_seq_supporting_jct(psm=psm, jct=self.jct) for psm in self.psms
        ]
        assert len(hybrid_seqs) == len(
            self.psms
        ), "Expected to get a hybrid seq for each PSM"
        return hybrid_seqs

    @property
    def uids(self):
        return [psm.uid for psm in self.psms]

    def get_spectra_supporting_jct(
        self, uid_to_spectrum: Dict[str, Spectrum]
    ) -> List[Spectrum]:
        return [uid_to_spectrum[uid] for uid in self.uids]


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
        jct_to_psms = defaultdict(list)
        for hybrid_psm in self.hybrid_psms:
            for prot in hybrid_psm.proteins:
                hy_pep = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=prot)
                jct_to_psms[hy_pep.get_junction_str(jct_len=self.jct_len)].append(
                    hybrid_psm
                )
        return dict(jct_to_psms)

    @property
    def junctions(self) -> List[SupportedHybridJunction]:
        return [
            SupportedHybridJunction(jct=jct, psms=psms)
            for jct, psms in self.jct_to_psms_map.items()
        ]

    def create_jct_df(
        self,
    ) -> pd.DataFrame:
        rows = []
        for junction in self.junctions:
            spectra_info = defaultdict(list)
            # Get information from spectra that support the junction
            for psm in junction.psms:
                spectrum = self.uid_to_spectrum[psm.uid]
                spectra_info["rts"].append(spectrum.retention_time)
                spectra_info["mzs"].append(spectrum.precursor_mz)
                spectra_info["zs"].append(spectrum.precursor_charge)
                spectra_info["intensities"].append(spectrum.precursor_intensity)
            rows.append(
                (
                    junction.jct,
                    junction.hybrid_seqs,
                    junction.uids,
                    spectra_info["mzs"],
                    spectra_info["rts"],
                    spectra_info["zs"],
                    spectra_info["intensities"],
                )
            )
        df = pd.DataFrame(
            rows,
            columns=[
                "jct",
                "hybrid_seqs",
                "uids",
                "spectra_mzs",
                "spectra_rts",
                "spectra_zs",
                "spectra_intensities",
            ],
        )
        df["num_uniq_seqs"] = df.hybrid_seqs.apply(lambda seqs: len(set(seqs)))
        df["num_spectra_supporting"] = df.spectra_mzs.apply(lambda mzs: len(mzs))
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


class PSMs(BaseModel):
    psms: List[CometPSM]

    @cached_property
    def top_psm_by_xcorr(self) -> CometPSM | None:
        """Gets PSM for which xcorr is highest"""
        top_psms = [psm for psm in self.psms if psm.xcorr == self.max_xcorr]
        assert (
            len(top_psms) <= 1
        ), f"Expected at most 1 top PSM by xcorr, but found {len(top_psms)}"
        return top_psms[0] if len(top_psms) == 1 else None

    @property
    def min_xcorr(self) -> float:
        return min(psm.xcorr for psm in self.psms)

    @property
    def max_xcorr(self) -> float:
        return max(psm.xcorr for psm in self.psms)

    @property
    def num_psms(self) -> int:
        return len(self.psms)

    @property
    def psms_sorted_by_xcorr(self) -> List[CometPSM]:
        return sorted(self.psms, key=lambda psm: psm.xcorr, reverse=True)


@dataclass
class HypedsearchSpectrumResults:
    spectrum: Spectrum
    native_targets: PSMs
    native_decoys: PSMs
    hybrid_targets: PSMs

    def __post_init__(self):
        # If the PSMs are passed as List[CometPSM] convert them to PSMs
        if isinstance(self.native_targets, list):
            self.native_targets = PSMs(psms=self.native_targets)
        if isinstance(self.native_decoys, list):
            self.native_decoys = PSMs(psms=self.native_decoys)
        if isinstance(self.hybrid_targets, list):
            self.hybrid_targets = PSMs(psms=self.hybrid_targets)

    @property
    def native_targets_combined_with_hybrid_targets(self) -> PSMs:
        """
        Combine the native and hybrid targets because HypedSearch doesn't include natives
        in the same competition as hybrids. Return a list of PSMs that is not longer
        than the number of native targets
        """
        psms = self.native_targets.psms + self.hybrid_targets.psms
        psms = sorted(psms, key=lambda psm: psm.xcorr, reverse=True)
        return PSMs(psms=psms[: self.native_targets.num_psms])

    @property
    def xcorr_ranges(self) -> Dict[str, float]:
        return {
            "uid": self.spectrum.uid,
            "hybrid_min_xcorr": self.native_targets_combined_with_hybrid_targets.min_xcorr,
            "hybrid_max_xcorr": self.native_targets_combined_with_hybrid_targets.max_xcorr,
            "native_min_xcorr": self.native_targets.min_xcorr,
            "native_max_xcorr": self.native_targets.max_xcorr,
            "native_decoy_min_xcorr": self.native_decoys.min_xcorr,
            "native_decoy_max_xcorr": self.native_decoys.max_xcorr,
        }

    @property
    def top_native_beats_top_hybrid(self) -> bool:
        if self.native_targets.max_xcorr >= self.hybrid_targets.max_xcorr:
            return True
        else:
            return False

    @property
    def top_decoy_beats_top_native(self) -> bool:
        if self.native_decoys.max_xcorr >= self.native_targets.max_xcorr:
            return True
        else:
            return False

    @property
    def uid(self):
        return self.spectrum.uid

    def create_top_hybrid_psm_plot(
        self, ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL, ax: Axes | None = None
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        hy_psm = self.hybrid_targets.top_psm_by_xcorr
        hybrid_psm_plot(
            spectrum=self.spectrum, hy_comet_psm=hy_psm, ppm_tol=ppm_tol, ax=ax
        )


class NativeVsHybridComparison(BaseModel):
    native_run: CometRunAnalysis
    hybrid_run: CometRunAnalysis
    spectra: List[Spectrum]
    name: Optional[str] = None

    @cached_property
    def uid_to_spectrum(self):
        return {spectrum.uid: spectrum for spectrum in self.spectra}

    @cached_property
    def uid_to_spectrum_results(self) -> Dict[str, HypedsearchSpectrumResults]:
        uid_to_spectrum_psms = defaultdict(lambda: {NTS: [], NDS: [], HTS: []})

        for psm in self.native_run.targets:
            uid_to_spectrum_psms[psm.uid][NTS].append(psm)
        for psm in self.native_run.decoys:
            uid_to_spectrum_psms[psm.uid][NDS].append(psm)
        for psm in self.hybrid_run.targets:
            uid_to_spectrum_psms[psm.uid][HTS].append(psm)

        for uid, data in uid_to_spectrum_psms.items():
            data["spectrum"] = self.uid_to_spectrum[uid]
            for key in [NTS, NDS, HTS]:
                data[key] = PSMs(psms=data[key])
            uid_to_spectrum_psms[uid] = HypedsearchSpectrumResults(**data)

        return dict(uid_to_spectrum_psms)

    @property
    def xcorr_ranges_df(self):
        return pd.DataFrame(
            [res.xcorr_ranges for res in self.uid_to_spectrum_results.values()]
        )

    @property
    def num_psms_per_spectrum(self):
        return {
            NTS: Counter(
                [
                    results.native_targets.num_psms
                    for results in self.uid_to_spectrum_results.values()
                ]
            ),
            NDS: Counter(
                [
                    results.native_decoys.num_psms
                    for results in self.uid_to_spectrum_results.values()
                ]
            ),
            HTS: Counter(
                [
                    results.hybrid_targets.num_psms
                    for results in self.uid_to_spectrum_results.values()
                ]
            ),
        }

    def create_native_vs_hybrid_xcorr_min_and_max_change_plot(
        self,
        ax: Axes | None = None,
        title: str = "",
        accepted_uids: Set[str] | List[str] = [],
    ):
        df = self.xcorr_ranges_df
        if ax is None:
            fig, axs = fig_setup()
            ax = axs[0]
        s = 7
        _ = sns.scatterplot(
            y=df.hybrid_max_xcorr - df.native_max_xcorr,
            x=df.hybrid_min_xcorr - df.native_min_xcorr,
            s=s,
            ax=ax,
        )
        if len(accepted_uids) > 0:
            df["accepted_hybrid"] = df.uid.isin(accepted_uids)
            tmp = df[df.accepted_hybrid]
            _ = sns.scatterplot(
                y=tmp.hybrid_max_xcorr - tmp.native_max_xcorr,
                x=tmp.hybrid_min_xcorr - tmp.native_min_xcorr,
                s=3 * s,
                ax=ax,
                marker="X",
                color="red",
                label="NeoFusion accepted hybrids",
            )
        set_title_axes_labels(
            ax=ax,
            title=title,
            ylabel="hybrid max xcorr - native max xcorr",
            xlabel="hybrid min xcorr - native min xcorr",
        )
        finalize(ax)
        # _ = ax.legend(
        #     loc="upper right",
        #     frameon=False,
        #     bbox_to_anchor=(1.75, 0.5),
        # )

    def create_native_vs_hybrid_xcorr_range_change_plot(
        self,
        # ax: Axes | None = None,
        title: str = "",
        top_n_psms: Optional[int] = None,
    ):
        # Create dataframe
        df = []
        for uid, results in self.uid_to_spectrum_results.items():
            if top_n_psms is not None:
                nts = PSMs(
                    psms=results.native_targets.psms_sorted_by_xcorr[:top_n_psms]
                )
                hts = PSMs(
                    psms=results.native_targets_combined_with_hybrid_targets.psms_sorted_by_xcorr[
                        :top_n_psms
                    ]
                )
            else:
                nts = results.native_targets
                hts = results.native_targets_combined_with_hybrid_targets
            df.append(
                {
                    "uid": uid,
                    "hybrid_min_xcorr": hts.min_xcorr,
                    "hybrid_max_xcorr": hts.max_xcorr,
                    "native_min_xcorr": nts.min_xcorr,
                    "native_max_xcorr": nts.max_xcorr,
                }
            )
        df = pd.DataFrame(df)

        # Make plot
        # if ax is None:
        #     _, axs = fig_setup()
        #     ax = axs[0]
        s = 7
        # _ = sns.scatterplot(
        #     y=df.hybrid_max_xcorr - df.native_max_xcorr,
        #     x=df.hybrid_min_xcorr - df.native_max_xcorr,
        #     s=s,
        #     ax=ax,
        # )
        if top_n_psms is not None:
            title += f" (top {top_n_psms} PSMs per spectrum)"
        p = create_joint_plot(
            x=df.hybrid_min_xcorr - df.native_max_xcorr,
            y=df.hybrid_max_xcorr - df.native_max_xcorr,
            s=s,
            xlabel="hybrid min xcorr - native max xcorr",
            ylabel="hybrid max xcorr - native max xcorr",
            title=title,
        )

    @cached_property
    def spectra_df(self) -> pd.DataFrame:
        return Spectrum.to_df(spectra=self.spectra)

    @cached_property
    def neofusion_output(self) -> NeoFusionOutput:
        neofusion = NeoFusionRunner(
            native_targets=self.native_run.top_targets,
            hybrid_targets=self.hybrid_run.top_targets,
        )
        neofusion_output = neofusion.run_neofusion()
        return neofusion_output

    @classmethod
    def from_config(cls, config: str | Path | HypedsearchRunConfig):
        if isinstance(config, (str, Path)):
            config = HypedsearchRunConfig.from_json(path=config)
        return cls(
            name=config.name,
            native_run=config.native_comet_run,
            hybrid_run=config.hybrid_comet_run,
            spectra=config.spectra,
        )

    @classmethod
    def from_txts_and_mzml(
        cls,
        native_targets: str | Path,
        native_decoys: str | Path,
        native_assign_conf: str | Path,
        hybrid_targets: str | Path,
        mzml: str | Path,
    ):
        native_targets = CometPSM.from_txt(txt=native_targets)
        native_decoys = CometPSM.from_txt(txt=native_decoys)
        native_assign_conf = CometPSM.from_txt(txt=native_assign_conf)
        hybrid_targets = CometPSM.from_txt(txt=hybrid_targets)
        return cls.from_comet_psms(
            native_targets=native_targets,
            native_decoys=native_decoys,
            native_assign_conf=native_assign_conf,
            hybrid_targets=hybrid_targets,
            mzml=mzml,
        )

    @classmethod
    def from_comet_psms(
        cls,
        native_targets: List[CometPSM],
        native_assign_conf: List[CometPSM],
        native_decoys: List[CometPSM],
        hybrid_targets: List[CometPSM],
        mzml: str | Path,
    ):
        native_run = CometRunAnalysis(
            targets=native_targets,
            decoys=native_decoys,
            assign_conf=native_assign_conf,  # this will extrapolate the native q(xcorr) fcn to hybrids
        )
        hybrid_run = CometRunAnalysis(
            targets=hybrid_targets,
            decoys=[],
            assign_conf=native_assign_conf,  # this will extrapolate the native q(xcorr) fcn to hybrids
        )
        return cls(
            native_run=native_run,
            hybrid_run=hybrid_run,
            spectra=Spectrum.parse_ms2_from_mzml(mzml=mzml),
        )

    def create_psm_dataframes(
        self,
        psm_types: List[str] = PSM_TYPES,
        ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
        out_dir: Optional[Path] = None,
    ) -> Dict[str, pd.DataFrame]:
        psm_type_to_psms = {}
        for psm_type in psm_types:
            if psm_type == NTS:
                psm_type_to_psms[psm_type] = self.native_run.targets
            elif psm_type == NDS:
                psm_type_to_psms[psm_type] = self.native_run.decoys
            elif psm_type == HTS:
                psm_type_to_psms[psm_type] = self.hybrid_run.targets
            else:
                raise ValueError(
                    f"Invalid psm_type {psm_type}. Must be one of 'native_targets', 'native_decoys', or 'hybrid_targets'."
                )
        psm_type_to_df = {}
        for psm_type, psms in psm_type_to_psms.items():
            if out_dir is not None:
                out_path = out_dir / PSM_TYPE_TO_CSV_NAME[psm_type]
                if out_path.exists():
                    logger.info(
                        f"PSM dataframe for {psm_type} already exists at {out_path}. Skipping creation..."
                    )
                    continue
            logger.info(f"Creating custom PSM dataframe for {psm_type}")
            df = pd.DataFrame(
                convert_comet_psms_to_custom_psms(
                    spectra=self.spectra,
                    comet_psms=psms,
                    ppm_tol=ppm_tol,
                )
            )
            if out_dir is not None:
                df.to_csv(out_path, index=False)
            psm_type_to_df[psm_type] = df
        return psm_type_to_df

    def get_protein_abundance_df(
        self, q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> pd.DataFrame:
        return self.native_run.get_protein_abundance(q_threshold=q_threshold).df

    def create_hybrid_run_plots(
        self, psm_df_dir: str | Path, title: str = "", out_dir: str | Path | None = None
    ):
        # Load data
        psm_df_dir = Path(psm_df_dir)
        # native_targets = pd.read_csv(psm_df_dir / PSM_TYPE_TO_CSV_NAME[NTS])
        # native_decoys = pd.read_csv(psm_df_dir / PSM_TYPE_TO_CSV_NAME[NDS])
        hybrid_targets = pd.read_csv(psm_df_dir / PSM_TYPE_TO_CSV_NAME[HTS])
        accepted_uids = [psm.uid for psm in self.neofusion_output.accepted_hybrids]
        hybrid_targets["accepted_hybrid"] = hybrid_targets.uid.isin(accepted_uids)
        title_y = 1.05

        # Plots
        fig, axs = fig_setup()
        self.create_native_vs_hybrid_xcorr_min_and_max_change_plot(
            accepted_uids=accepted_uids,
            ax=axs[0],
            # title=title,
        )
        if out_dir is not None:
            save_fig(
                fig=fig,
                title=title,
                path=Path(out_dir) / f"native_vs_hybrid_run_xcorr_change.png",
            )

        # Plot top hybrid targets
        tmp = hybrid_targets[(hybrid_targets.num == 1)]
        g = sns.JointGrid(
            data=tmp,
            x="hybrid_left_support",
            y="hybrid_right_support",
            marginal_ticks=True,
        )

        # Main scatter, size mapped by
        _ = g.plot_joint(
            sns.scatterplot,
            size=np.log10(tmp.precursor_intensity),
            sizes=(20, 200),
            legend="brief",
        )
        _ = g.plot_marginals(
            sns.histplot,
            kde=True,
        )
        _ = sns.move_legend(
            g.ax_joint,
            "upper left",
            bbox_to_anchor=(1.2, 0.5),
            frameon=True,
        )
        _ = g.fig.suptitle(f"{title}\nTop hybrid PSMs (num=1)", y=title_y)
        if out_dir is not None:
            g.savefig(Path(out_dir) / f"top_hybrid_psms_support_jointplot.png")

        # Plot accepted hybrid targets
        tmp = hybrid_targets[(hybrid_targets.num == 1) & hybrid_targets.accepted_hybrid]
        g = sns.JointGrid(
            data=tmp,
            x="hybrid_left_support",
            y="hybrid_right_support",
            marginal_ticks=True,
        )
        _ = g.plot_joint(
            sns.scatterplot,
            size=np.log10(tmp.precursor_intensity),
            sizes=(20, 200),
            legend="brief",
        )
        _ = g.plot_marginals(
            sns.histplot,
            kde=True,
        )
        _ = sns.move_legend(
            g.ax_joint,
            "upper left",
            bbox_to_anchor=(1.2, 0.5),
            frameon=True,
        )
        _ = g.fig.suptitle(f"{title}\nAccepted hybrids via Neo-Fusion", y=title_y)
        if out_dir is not None:
            g.savefig(
                Path(out_dir) / f"neofusion_accepted_hybrids_support_jointplot.png"
            )

    def create_xcorr_target_and_decoy_distributions(
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

    def create_xcorr_accepted_vs_not_accepted_plot(
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
        # neofusion.plot_neofusion_true_positive_data(
        #     neofusion_iterations=neofusion_output.iterations
        # )
        return neofusion_output

    def create_xcorr_score_range_jointplots(
        self, top_n_psms: int | None = None, title: str = ""
    ):
        # Create dataframe from which to plot
        df = []
        for uid, results in self.uid_to_spectrum_results.items():
            if top_n_psms is not None:
                nts = PSMs(
                    psms=results.native_targets.psms_sorted_by_xcorr[:top_n_psms]
                )
                hts = PSMs(
                    psms=results.native_targets_combined_with_hybrid_targets.psms_sorted_by_xcorr[
                        :top_n_psms
                    ]
                )
            else:
                nts = results.native_targets
                hts = results.native_targets_combined_with_hybrid_targets
            df.append(
                {
                    "uid": uid,
                    "hybrid_min_xcorr": hts.min_xcorr,
                    "hybrid_max_xcorr": hts.max_xcorr,
                    "native_min_xcorr": nts.min_xcorr,
                    "native_max_xcorr": nts.max_xcorr,
                }
            )
        df = pd.DataFrame(df)

        # Plots
        s = 4
        if top_n_psms is not None:
            title += f"\nTop {top_n_psms} PSMs per spectrum"
        p = create_joint_plot(
            x=df.hybrid_min_xcorr - df.native_max_xcorr,
            y=df.hybrid_max_xcorr - df.native_max_xcorr,
            s=s,
            xlabel="hybrid min xcorr - native max xcorr",
            ylabel="hybrid max xcorr - native max xcorr",
            title=title,
        )
        p = create_joint_plot(
            x=df.native_max_xcorr - df.native_min_xcorr,
            y=df.hybrid_max_xcorr - df.hybrid_min_xcorr,
            s=s,
            xlabel="native max xcorr - native min xcorr",
            ylabel="hybrid max xcorr - hybrid min xcorr",
            title=title,
        )
        p = create_joint_plot(
            x=df.native_max_xcorr,
            y=df.native_max_xcorr - df.native_min_xcorr,
            s=s,
            xlabel="native max xcorr",
            ylabel="native max xcorr - native min xcorr",
            title=title,
        )
        p = create_joint_plot(
            x=df.hybrid_max_xcorr,
            y=df.hybrid_max_xcorr - df.hybrid_min_xcorr,
            s=s,
            xlabel="hybrid max xcorr",
            ylabel="hybrid max xcorr - hybrid min xcorr",
            title=title,
        )

    def junction_analysis(
        self,
        jct_len: int = DEFAULT_JCT_LEN,
        out_dir: str | Path | None = None,
        title: str = "",
    ):
        jct_analyzer = JunctionAnalysis(
            hybrid_psms=self.neofusion_output.accepted_hybrids,
            jct_len=jct_len,
            spectra=self.spectra,
        )
        jct_df = jct_analyzer.create_jct_df()
        if out_dir is not None:
            # Save dataframe
            jct_df.to_csv(Path(out_dir) / "neofusion_jct_df.csv", index=False)

            # Create and save plot
            fig, axs = fig_setup()
            _ = sns.scatterplot(
                data=jct_df,
                x="num_spectra_supporting",
                y="num_uniq_seqs",
                # s=7,
                ax=axs[0],
            )
            set_title_axes_labels(
                ax=axs[0],
                xlabel="Number spectra\nsupporting junction",
                ylabel="Number unique peptide seqs\nsupporting junction",
            )
            finalize(axs)
            save_fig(
                path=Path(out_dir) / "junction_support_plot.png",
                fig=fig,
                title=f"{title}\nNeoFusion accepted hybrids",
            )
        return jct_df


def neofusion_hybrid_acceptance_postprocessing(
    hs_config: Path,
    out_dir: str | Path,
    native_q_threshold: float = DEFAULT_Q_THRESHOLD,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    jct_len: int = DEFAULT_JCT_LEN,
):
    comp = NativeVsHybridComparison.from_config(config=hs_config)
    neofusion_output = comp.run_neofusion()
    if not (out_dir / "neofusion_psm_df.csv").exists():
        logger.info("Creating PSM dataframe with NeoFusion acceptance info")
        accepted_uids = [psm.uid for psm in neofusion_output.accepted_hybrids]
        df = comp.get_custom_psm_df(peak_to_ion_ppm_tol=peak_to_ion_ppm_tol)
        df["accepted_hybrids"] = df["uid"].isin(accepted_uids)
        df[f"accepted_natives{native_q_threshold}"] = (
            df.q_value_nt <= native_q_threshold
        )

        # Add spectrum info
        df["retention_time"] = df.uid.apply(
            lambda uid: comp.uid_to_spectrum[uid].retention_time
        )
        df["precursor_mz"] = df.uid.apply(
            lambda uid: comp.uid_to_spectrum[uid].precursor_mz
        )
        df["precursor_charge"] = df.uid.apply(
            lambda uid: comp.uid_to_spectrum[uid].precursor_charge
        )
        df["precursor_intensity"] = df.uid.apply(
            lambda uid: comp.uid_to_spectrum[uid].precursor_intensity
        )
        df.to_csv(out_dir / "neofusion_psm_df.csv", index=False)

    logger.info("Creating junction dataframe")
    jct_runner = JunctionAnalysis(
        hybrid_psms=neofusion_output.accepted_hybrids,
        jct_len=jct_len,
        spectra=comp.spectra,
    )
    jct_df = jct_runner.create_jct_df()
    jct_df.to_csv(out_dir / f"neofusion_jct_df.csv", index=False)

    return df


def interpolate_psm_q_values_from_assign_conf(
    interpolating_psms: List[CometPSM], q_valueless_psms: List[CometPSM]
) -> List[CometPSM]:
    interpolator = fit_xcorr_to_qval_interpolator(psms=interpolating_psms)
    for psm in q_valueless_psms:
        psm.q_value = float(interpolator(psm.xcorr))
    return q_valueless_psms
