import itertools
import os
import time
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pydantic import BaseModel, Field, model_validator
from scipy.interpolate import PchipInterpolator
from scipy.stats import ecdf
from statsmodels.distributions.empirical_distribution import ECDF

from src.comet_utils import CometPSM
from src.constants import (
    ASSIGN_CONFIDENCE,
    DATA_DIR,
    DECOY,
    DEFAULT_FPR,
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_RANGE,
    DEFAULT_Q_VAL_THRESH,
    DEFAULT_SCORE_CHANGE_RANGE,
    HUMAN_PROTEOME,
    HY_TARGET,
    HYBRID,
    NAT_DECOY,
    NAT_TARGET,
    NATIVE,
    NEOFUSION,
    Q_VAL,
    SPECTRA_PSMS_FILE_NAME,
    TARGET,
    TRUE_HYBRIDS_PATH,
    XCORR,
)
from src.hybrids_via_clusters import HybridJunction, HybridPeptide, HybridPosition
from src.hypedsearch import (
    HypedsearchRunConfig,
    SpectrumCometRunResults,
    TrueHybrid,
    find_possible_hybrids_for_seq,
    get_seq_to_hybrids_map,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml, Spectrum, organize_by_spectrum_uid
from src.peptide_spectrum_comparison import PSM
from src.peptides_and_ions import Fasta, ProteinRange
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.protein_abundance import ProteinAbundance
from src.utils import (
    PathType,
    Position,
    decompress_and_depickle,
    flatten_list_of_lists,
    get_positions_of_subseq_in_seq,
    pickle_and_compress,
    setup_logger,
    to_json,
    write_new_line_separated_file,
)

NONXCORR_PSM_SCORES = [
    "prop_intensity_supported",
    "prop_ions_matched",
    "prop_prefixes_supported",
    "prop_suffixes_supported",
]

PSM_SCORES = NONXCORR_PSM_SCORES + ["xcorr"]


@dataclass
class HypedsearchOutputs:
    hs_config: HypedsearchRunConfig
    min_side_len: int = DEFAULT_MIN_SIDE_LEN
    remove_carbamidomethylation: bool = True
    native_targets: Dict[str, List[CometPSM]] = field(init=False)
    native_assign_conf: Dict[str, CometPSM] = field(init=False)
    native_decoys: Dict[str, List[CometPSM]] = field(init=False)
    hybrid_targets: Dict[str, List[CometPSM]] = field(init=False)
    hybrid_decoys: Dict[str, List[CometPSM]] = field(init=False)

    def __post_init__(self):
        self.native_targets = self._get_native_targets()
        self.native_assign_conf = self._get_native_assign_conf()
        self.native_decoys = self._get_native_decoys()
        self.hybrid_targets = self._get_hybrid_targets()
        self.hybrid_decoys = self._get_hybrid_decoys()

    @cached_property
    def hybrid_seqs(self) -> Set[str]:
        return set(
            psm.seq for psm in flatten_list_of_lists(self.hybrid_targets.values())
        )

    @cached_property
    def seq_to_hybrids_map(self) -> Dict[str, List[HybridPeptide]]:
        return get_seq_to_hybrids_map(
            seqs=self.hybrid_seqs,
            db_path=self.hs_config.kmer_db_path,
            min_side_len=self.min_side_len,
        )

    @cached_property
    def q_value_interpolator(self) -> PchipInterpolator:
        return fit_xcorr_to_qval_interpolator(psms=self.native_assign_conf.values())

    @cached_property
    def spectrum_uid_to_spectrum(self) -> Dict[str, Spectrum]:
        spectrum_uid_to_spectrum = {}
        for mzml in self.hs_config.mzml_to_scans.keys():
            mzml = Mzml(mzml=mzml)
            spectrum_uid_to_spectrum.update(mzml.id_to_spectrum)
        return spectrum_uid_to_spectrum

    @cached_property
    def all_spectrum_uids(self) -> List[str]:
        return list(
            set(self.native_targets.keys()).union(set(self.hybrid_targets.keys()))
        )

    @property
    def top_hybrid_targets(self):
        top_hybrid_targets = organize_by_spectrum_uid(
            [
                psm
                for psm in flatten_list_of_lists(self.hybrid_targets.values())
                if psm.num == 1
            ]
        )
        for spectrum_uid, psms in top_hybrid_targets.items():
            assert len(psms) == 1, f"More than one top hybrid PSM for {spectrum_uid}"
            top_hybrid_targets[spectrum_uid] = psms[0]
        return top_hybrid_targets

    def get_native_beating_low_q_value_hybrids(
        self, q_value_threshold: float
    ) -> List[CometPSM]:
        psms = []
        for spectrum_uid, hybrid_target in self.top_hybrid_targets.items():
            if spectrum_uid in self.native_assign_conf:
                native_target = self.native_assign_conf[spectrum_uid]
                if hybrid_target.xcorr < native_target.xcorr:
                    continue
            if hybrid_target.q_value <= q_value_threshold:
                psms.append(hybrid_target)
        return psms

    def _get_native_targets(self) -> Dict[str, List[CometPSM]]:
        return CometPSM.from_txts(
            txts=self.hs_config.get_output_txts(run_type=NATIVE, psm_type=TARGET),
            by_spectrum=True,
        )

    def _get_native_assign_conf(self) -> Dict[str, List[CometPSM]]:
        spectrum_to_psm = CometPSM.from_txts(
            txts=self.hs_config.get_output_txts(
                run_type=NATIVE, psm_type=ASSIGN_CONFIDENCE
            ),
            by_spectrum=True,
        )
        for key, psms in spectrum_to_psm.items():
            assert len(psms) == 1, f"More than one assign confidence PSM for {key}"
            spectrum_to_psm[key] = psms[0]
        return spectrum_to_psm

    def _get_native_decoys(self) -> Dict[str, List[CometPSM]]:
        return CometPSM.from_txts(
            txts=self.hs_config.get_output_txts(run_type=NATIVE, psm_type=DECOY),
            by_spectrum=True,
        )

    def _get_hybrid_targets(self) -> Dict[str, List[CometPSM]]:
        psms = [
            psm
            for psm in CometPSM.from_txts(
                txts=self.hs_config.get_output_txts(run_type=HYBRID, psm_type=TARGET),
                by_spectrum=False,
            )
            if psm.is_hybrid
        ]
        seq_to_hybrids = get_seq_to_hybrids_map(
            seqs=set(psm.seq for psm in psms),
            db_path=self.hs_config.kmer_db_path,
            min_side_len=self.min_side_len,
            remove_carbamidomethylation=self.remove_carbamidomethylation,
        )
        hybrid_supported_psms = [psm for psm in psms if psm.seq in seq_to_hybrids]
        return organize_by_spectrum_uid(data=hybrid_supported_psms)

    def _get_hybrid_decoys(self) -> Dict[str, List[CometPSM]]:
        return CometPSM.from_txts(
            txts=self.hs_config.get_output_txts(run_type=HYBRID, psm_type=DECOY),
            by_spectrum=True,
        )

    def get_native_results_for_spectrum(
        self, spectrum_uid: str
    ) -> SpectrumCometRunResults:
        return SpectrumCometRunResults(
            targets=self.native_targets.get(spectrum_uid, []),
            decoys=self.native_decoys.get(spectrum_uid, []),
            assign_conf=self.native_assign_conf.get(spectrum_uid, [None])[0],
        )

    def get_hybrid_results_for_spectrum(
        self, spectrum_uid: str
    ) -> SpectrumCometRunResults:
        return SpectrumCometRunResults(
            targets=self._get_hybrid_targets.get(spectrum_uid, []),
            decoys=self._get_hybrid_decoys.get(spectrum_uid, []),
            assign_conf=None,
        )

    def set_hybrid_target_q_values(self):
        for psms in self.hybrid_targets.values():
            for psm in psms:
                psm.q_value = float(self.q_value_interpolator(psm.xcorr))

    def get_protein_abundances(
        self, q_threshold: float = DEFAULT_Q_VAL_THRESH
    ) -> ProteinAbundance:
        return ProteinAbundance.from_comet_psms(
            psms=self.native_assign_conf.values(),
            q_val_thresh=q_threshold,
        )

    def get_native_target_psm(
        self, spectrum_uid: str, peak_to_ion_ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        if spectrum_uid in self.native_assign_conf:
            return PSM.from_spectrum_and_comet_psm(
                comet_psm=self.native_assign_conf[spectrum_uid],
                spectrum=self.spectrum_uid_to_spectrum[spectrum_uid],
                peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            )
        else:
            return None

    def get_native_decoy_psm(
        self, spectrum_uid: str, peak_to_ion_ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        if spectrum_uid in self.native_decoys:
            native_decoy = [
                psm for psm in self.native_decoys[spectrum_uid] if psm.num == 1
            ][0]
            return PSM.from_spectrum_and_comet_psm(
                comet_psm=native_decoy,
                spectrum=self.spectrum_uid_to_spectrum[spectrum_uid],
                peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            )
        else:
            return None

    def get_hybrid_target_psm(
        self, spectrum_uid: str, peak_to_ion_ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        if spectrum_uid in self.top_hybrid_targets:
            return PSM.from_spectrum_and_comet_psm(
                comet_psm=self.top_hybrid_targets[spectrum_uid],
                spectrum=self.spectrum_uid_to_spectrum[spectrum_uid],
                peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            )
        else:
            return None

    def get_psm_results(
        self, peak_to_ion_ppm_tol: int = DEFAULT_PEAK_TO_ION_PPM_TOL
    ) -> List[Dict]:
        """


        For ~15,000 spectra, this took ~5m to run locally in a Jupyter notebook
        """
        results = []
        for idx, spectrum_uid in enumerate(self.all_spectrum_uids):
            print(f"Processing {idx+1}/{len(self.all_spectrum_uids)}", end="\r")
            data = self.get_native_target_psm(
                spectrum_uid=spectrum_uid, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
            )
            if data is not None:
                data = data.to_dict()
                data["type"] = NAT_TARGET
                results.append(data)
            data = self.get_hybrid_target_psm(
                spectrum_uid=spectrum_uid, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
            )
            if data is not None:
                data = data.to_dict()
                data["type"] = HY_TARGET
                results.append(data)
            data = self.get_native_decoy_psm(
                spectrum_uid=spectrum_uid, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
            )
            if data is not None:
                data = data.to_dict()
                data["type"] = NAT_DECOY
                results.append(data)
        return results


@dataclass
class SpectrumPSMs:
    native_target: Optional[PSM] = None
    native_decoy: Optional[PSM] = None
    hybrid_target: Optional[PSM] = None
    hybrid_decoy: Optional[PSM] = None

    def __post_init__(self):
        spectrum_uids = []
        if self.native_target is not None:
            spectrum_uids.append(self.native_target.spectrum_uid)
        if self.native_decoy is not None:
            spectrum_uids.append(self.native_decoy.spectrum_uid)
        if self.hybrid_target is not None:
            spectrum_uids.append(self.hybrid_target.spectrum_uid)
        if self.hybrid_decoy is not None:
            spectrum_uids.append(self.hybrid_decoy.spectrum_uid)
        assert len(set(spectrum_uids)) == 1, "All PSMs must have the same spectrum UID"

    @property
    def native_seq(self):
        return self.native_target.seq

    @property
    def native_q(self):
        if self.native_target is None:
            return None
        return self.native_target.q_value

    @property
    def native_xcorr(self):
        return self.native_target.xcorr

    @property
    def hybrid_seq(self):
        if self.hybrid_target is not None:
            return self.hybrid_target.seq
        else:
            return None

    @property
    def hybrid_xcorr(self):
        if self.hybrid_target is not None:
            return self.hybrid_target.xcorr
        else:
            return None

    @property
    def hybrid_hyphen_seqs(self) -> List[str]:
        if self.hybrids is not None:
            return list(hy.hyphen_seq for hy in self.hybrids)
        else:
            return None

    @property
    def mzml(self):
        return self.spectrum.mzml

    @classmethod
    def from_spectrum_and_comet_psms(
        cls,
        spectrum: Spectrum,
        native_target: Optional[CometPSM],
        native_decoy: Optional[CometPSM],
        hybrid_target: Optional[CometPSM],
        peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ) -> "SpectrumPSMs":
        return cls(
            native_target=(
                PSM.from_spectrum_and_comet_psm(
                    spectrum=spectrum,
                    comet_psm=native_target,
                    peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
                )
                if native_target is not None
                else None
            ),
            native_decoy=(
                PSM.from_spectrum_and_comet_psm(
                    spectrum=spectrum,
                    comet_psm=native_decoy,
                    peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
                )
                if native_decoy is not None
                else None
            ),
            hybrid_target=(
                PSM.from_spectrum_and_comet_psm(
                    spectrum=spectrum,
                    comet_psm=hybrid_target,
                    peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
                )
                if hybrid_target is not None
                else None
            ),
        )

    @classmethod
    def load(cls, path: Union[str, Path]) -> Dict[str, "SpectrumPSMs"]:
        data = decompress_and_depickle(path=path)
        data = [cls(**d) for d in data]
        # return {d.spectrum_uid: d for d in data}
        return data

    @staticmethod
    def default_save_path(out_dir: Path, min_side_len: int):
        return out_dir / f"spectra_psms_minSideLen{min_side_len}.pklz"

    @staticmethod
    def to_df(spectra_psms: List["SpectrumPSMs"]) -> pd.DataFrame:
        data = flatten_list_of_lists(
            [spectrum_psms.to_dicts() for spectrum_psms in spectra_psms]
        )
        return pd.DataFrame(data)

    @staticmethod
    def save(
        spectra_psms: List["SpectrumPSMs"],
        out_path: Union[str, Path],
        overwrite: bool = False,
    ):
        out_path = Path(out_path)
        if out_path.exists() and not overwrite:
            logger.info(
                f"Spectra PSMs file already exists at {out_path}, not overwriting. "
                + "Skipping saving..."
            )
        else:
            data = [spectrum_psms.model_dump() for spectrum_psms in spectra_psms]
            out_path.parent.mkdir(parents=True, exist_ok=True)
            pickle_and_compress(obj=data, path=out_path)


@dataclass
class HybridSupport:
    left_native_support: float
    left_hybrid_support: float
    right_native_support: float
    right_hybrid_support: float

    @property
    def native_mean(self):
        return np.mean([self.left_native_support, self.right_native_support])

    @property
    def native_min(self):
        return min(self.left_native_support, self.right_native_support)

    @property
    def native_max(self):
        return max(self.left_native_support, self.right_native_support)

    @property
    def hybrid_mean(self):
        return np.mean([self.left_hybrid_support, self.right_hybrid_support])

    @property
    def hybrid_min(self):
        return min(self.left_hybrid_support, self.right_hybrid_support)

    @property
    def hybrid_max(self):
        return max(self.left_hybrid_support, self.right_hybrid_support)

    @property
    def mean_support(self):
        return np.mean(
            [
                self.left_native_support,
                self.left_hybrid_support,
                self.right_native_support,
                self.right_hybrid_support,
            ]
        )


def get_support_in_pileup(
    pos: ProteinRange,
    pileup: Dict[str, Dict[int, int]],
) -> List[int]:
    if pos.protein in pileup:
        return [
            pileup[pos.protein][idx]
            for idx in range(pos.inclusive_start, pos.exclusive_end)
        ]
    else:
        return [0 for _ in range(pos.inclusive_start, pos.exclusive_end)]


def align_psms_to_proteome(psms: List[CometPSM], fasta: Path) -> List[ProteinRange]:
    fasta = Fasta(path=fasta)
    positions = []
    for psm in psms:
        for prot_name in psm.proteins:
            psm_seq_positions = get_positions_of_subseq_in_seq(
                subseq=psm.seq, seq=fasta.protein_name_to_seq_map[prot_name]
            )
            for pos in psm_seq_positions:
                positions.append(ProteinRange.from_pos(protein=prot_name, pos=pos))
    return positions


def get_pileup_from_positions(positions: List[ProteinRange]):
    pileup = defaultdict(lambda: defaultdict(int))
    for pos in positions:
        for idx in range(pos.inclusive_start, pos.exclusive_end):
            pileup[pos.protein][idx] += 1
    return dict(pileup)


def get_hybrid_psm_pileups(psms: List[SpectrumPSMs], fasta: Path):
    hybrid_positions = SpectrumPSMs.align_hybrid_psms_to_proteome(
        psms=psms, fasta=fasta
    )
    left_pileup = get_pileup_from_positions(
        positions=[pos.left for pos in hybrid_positions]
    )
    right_pileup = get_pileup_from_positions(
        positions=[pos.right for pos in hybrid_positions]
    )
    return left_pileup, right_pileup


def get_native_pileup(psms: List[CometPSM], fasta: Path):
    return get_pileup_from_positions(
        positions=align_psms_to_proteome(
            psms=psms,
            fasta=fasta,
        )
    )


def plot_psm_pileup(
    left_hybrid_pileup: Dict[str, Dict[int, int]],
    right_hybrid_pileup: Dict[str, Dict[int, int]],
    fasta: Path,
    native_pileup: Optional[Dict[str, Dict[int, int]]] = None,
):
    prot_names = set(left_hybrid_pileup.keys()).union(set(right_hybrid_pileup.keys()))
    fasta = Fasta(path=fasta)
    fig, axs = fig_setup(nrows=len(prot_names), ncols=1, w=8)
    max_hybrid_cnt = max(
        max(
            flatten_list_of_lists(
                pileup.values() for pileup in left_hybrid_pileup.values()
            )
        ),
        max(
            flatten_list_of_lists(
                pileup.values() for pileup in right_hybrid_pileup.values()
            )
        ),
    )
    max_native_cnt = max(
        flatten_list_of_lists(pileup.values() for pileup in native_pileup.values())
    )
    for idx, prot_name in enumerate(prot_names):
        ax = axs[idx]
        prot_seq = fasta.protein_name_to_seq_map[prot_name]
        # Left-side hybrid pileup
        if prot_name in left_hybrid_pileup:
            data = np.array(
                [(ii, left_hybrid_pileup[prot_name][ii]) for ii in range(len(prot_seq))]
            )
            data = data[data[:, 1] != 0]  # remove zeros
            _ = ax.plot(
                data[:, 0],
                data[:, 1],
                "o",
                color="red",
                label=f"left-side support",
                ms=2,
            )
        # Right-side hybrid pileup
        if prot_name in right_hybrid_pileup:
            data = np.array(
                [
                    (ii, right_hybrid_pileup[prot_name][ii])
                    for ii in range(len(prot_seq))
                ]
            )
            data = data[data[:, 1] != 0]  # remove zeros
            _ = ax.plot(
                data[:, 0],
                data[:, 1],
                "o",
                color="blue",
                label=f"right-side support",
                ms=2,
            )
        set_title_axes_labels(ax=ax, title=prot_name)
        _ = ax.set_ylim(bottom=0, top=max_hybrid_cnt + 1)
        _ = ax.set_xlim(left=0, right=len(prot_seq) + 1)
        # Plot native pileup on its own axis because it may have a different scale
        if (native_pileup is not None) and (prot_name in native_pileup):
            ax_copy = ax.twinx()
            data = np.array(
                [(ii, native_pileup[prot_name][ii]) for ii in range(len(prot_seq))]
            )
            data = data[data[:, 1] != 0]  # remove zeros
            _ = ax_copy.plot(
                data[:, 0],
                data[:, 1],
                "o",
                color="green",
                label=f"native support",
                ms=2,
            )
            ax_copy.set_ylabel("Native support", color="tab:green")
            ax_copy.set_ylim(bottom=0, top=max_native_cnt + 1)
            ax_copy.tick_params(axis="y", labelcolor="tab:green")
    finalize(axs)
    return fig, axs


@dataclass
class NeoFusionIteration:
    q: float
    min_score_delta: float
    fpr: float
    tp: int
    min_hybrid_score: float
    accepted_hybrid_psm_spectrum_uids: List[str]

    @property
    def info(self) -> str:
        return (
            f"q={self.q}, min_score_delta={self.min_score_delta}, fpr={self.fpr}, "
            f"tp={self.tp}, min_hybrid_score={self.min_hybrid_score}, "
            f"num_accepted_psm={len(self.accepted_hybrid_psm_spectrum_uids)}"
        )

    @property
    def param_str(self) -> str:
        return f"q{self.q}_delta{self.min_score_delta}_fpr{self.fpr}_minHybridScore{self.min_hybrid_score}"


@dataclass
class NeoFusionRunner:
    native_assign_conf: Dict[str, CometPSM]
    top_hybrid_targets: Dict[str, CometPSM]
    q_vals: List[float] = field(default_factory=lambda: DEFAULT_Q_RANGE.copy())
    score_deltas: List[float] = field(
        default_factory=lambda: DEFAULT_SCORE_CHANGE_RANGE.copy()
    )
    fpr_threshold: float = DEFAULT_FPR

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
                    q=q_val,
                    min_score_delta=score_delta,
                    fpr=tp_maximizing_row.fpr,
                    tp=tp_maximizing_row.tp,
                    min_hybrid_score=min_hybrid_score,
                    accepted_hybrid_psm_spectrum_uids=accepted_psm,
                )
            except:
                logger.debug(f"Issue with q_val={q_val}, score_delta={score_delta}")
        return None

    def run_neofusion(
        self,
    ) -> List[NeoFusionIteration]:
        df = self.create_neofusion_df(
            native_assign_conf=self.native_assign_conf,
            top_hybrid_targets=self.top_hybrid_targets,
        )

        neofusion_results = []
        for q_val in self.q_vals:
            for min_score_delta in self.score_deltas:
                result = self.neofusion_iteration(
                    df=df,
                    q_val=q_val,
                    score_delta=min_score_delta,
                    fpr_thresh=self.fpr_threshold,
                )
                if result is not None:
                    neofusion_results.append(result)

        return neofusion_results

    @staticmethod
    def plot_neofusion_true_positive_data(
        neofusion_results: List[NeoFusionIteration],
        title: str = "",
    ) -> Axes:
        data = {res.param_str: res.tp for res in neofusion_results}
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
        neofusion_results: List[NeoFusionIteration],
    ) -> Tuple[NeoFusionIteration, List[CometPSM]]:
        best_iteration = max(neofusion_results, key=lambda x: x.tp)
        accepted_hybrids = [
            self.top_hybrid_targets[spectrum_uid]
            for spectrum_uid in best_iteration.accepted_hybrid_psm_spectrum_uids
        ]
        return best_iteration, accepted_hybrids


@dataclass
class AcceptedHybridPSMs:
    accepted_psms: List[SpectrumPSMs]

    @property
    def seq_to_accepted_psms(self) -> Dict[str, List[SpectrumPSMs]]:
        seq_to_psms = defaultdict(list)
        for psm in self.accepted_psms:
            seq_to_psms[psm.hybrid_seq].append(psm)
        return dict(seq_to_psms)

    @property
    def summary_dict(self):
        return {
            "num_accepted_hybrid_psms": len(self.accepted_psms),
            "num_accepted_unique_hybrid_seqs": len(self.seq_to_accepted_psms),
            "accepted_hybrid_psms": [psm.spectrum_uid for psm in self.accepted_psms],
        }


@dataclass
class TrueHybridSpectrumPSMsComparison:
    true_hybrid_seqs: Set[str]
    supported_true_seqs: List[Dict]

    @property
    def num_spectra_supporting_true_hybrids(self) -> int:
        return sum(
            seq_data["num_supporting_hybrid_psms"]
            for seq_data in self.supported_true_seqs
        )

    @classmethod
    def from_accepted_psms(
        cls, accepted_psms: AcceptedHybridPSMs, true_hybrid_seqs: Set[str]
    ) -> Dict:
        supported_true_seqs = []
        for seq in set(accepted_psms.seq_to_accepted_psms.keys()).intersection(
            true_hybrid_seqs
        ):
            seq_psms = accepted_psms.seq_to_accepted_psms[seq]
            data = {
                "seq": seq,
                "num_supporting_hybrid_psms": len(seq_psms),
                "spectra_uids": [psm.spectrum_uid for psm in seq_psms],
            }
            supported_true_seqs.append(data)
        return cls(
            true_hybrid_seqs=true_hybrid_seqs,
            supported_true_seqs=supported_true_seqs,
        )

    @property
    def summary_dict(self):
        return {
            "num_true_hybrid_seqs": len(self.true_hybrid_seqs),
            "num_supported_true_hybrid_seqs": len(self.supported_true_seqs),
            "num_spectra_supporting_true_hybrids": self.num_spectra_supporting_true_hybrids,
            "supported_true_hybrid_seqs": self.supported_true_seqs,
        }


@dataclass
class PositionedHybrid:
    pos: HybridPosition
    hybrid: HybridPeptide
    spectrum_uid: str

    @property
    def junction(self):
        return HybridJunction.from_hybrid_position(pos=self.pos)

    @property
    def seq(self):
        return self.hybrid.seq

    @property
    def hyphen_seq(self):
        return self.hybrid.hyphen_seq

    def to_dict(self):
        return {
            "spectrum_uid": self.spectrum_uid,
            "hybrid_seq": self.seq,
            "hybrid_hyphen_seq": self.hyphen_seq,
            "position": asdict(self.pos),
        }


def get_junction_matching_hybrid(
    spectrum_uid: str,
    spectrum_psms: SpectrumPSMs,
    target_junction: str,
    protein_name_to_seq_map: Dict[str, str],
):
    found_hybrid = None
    for hy in spectrum_psms.hybrids:
        if found_hybrid is not None:
            break
        for pos in get_positions_of_hybrid(
            hybrid=hy,
            protein_name_to_seq_map=protein_name_to_seq_map,
        ):
            if str(HybridJunction.from_hybrid_position(pos=pos)) == target_junction:
                found_hybrid = PositionedHybrid(
                    pos=pos, hybrid=hy, spectrum_uid=spectrum_uid
                )
                break
    return found_hybrid


@dataclass
class Experiment:
    name: str
    results_dir: Path
    hs_config: Path
    psms_path: Path

    def __post_init__(self):
        assert self.hs_config.exists(), f"HS config does not exist at {self.hs_config}"
        assert self.psms_path.exists(), f"PSMs path does not exist at {self.psms_path}"
        self.results_dir.mkdir(parents=True, exist_ok=True)

    @cached_property
    def psms(self):
        return SpectrumPSMs.load(path=self.psms_path)

    @property
    def _hs_config(self):
        return HypedsearchRunConfig.from_json(path=self.hs_config)

    @property
    def native_assign_confidence_txt(self) -> Path:
        return self._hs_config.native_assign_confidence_path

    @property
    def fasta_path(self) -> Path:
        return self._hs_config.hybrid_former.fasta

    def get_pileups(
        self,
        accepted_hybrid_psms: List[SpectrumPSMs],
        q_thresh: float = DEFAULT_Q_VAL_THRESH,
    ) -> Tuple[Dict, Dict, Dict]:
        left_pileup, right_pileup = get_hybrid_psm_pileups(
            psms=accepted_hybrid_psms, fasta=self.fasta_path
        )
        native_psms = [
            psm
            for psm in CometPSM.from_txt(txt=self.native_assign_confidence_txt)
            if psm.q_value <= q_thresh
        ]
        native_pileup = get_native_pileup(psms=native_psms, fasta=self.fasta_path)
        return left_pileup, right_pileup, native_pileup

    def xcorr_plot(self) -> Axes:
        _, axs = fig_setup(w=8)
        ax = axs[0]
        _ = xcorr_plot(psms=self.psms, ax=ax)
        # Add true-hybrid-containing hPSMs rugplot
        true_hybrid_containing_psms = compare_psms_to_true_hybrids(
            psms=self.psms, results_dir=self.results_dir
        )
        ymin, ymax = ax.get_ylim()
        rug_height = 0.1 * (ymax - ymin)  # small line height
        lw = 0.2
        data = [psm.hybrid_xcorr for psm in true_hybrid_containing_psms]
        _ = ax.vlines(
            data[0],
            ymax,
            ymax + rug_height,
            color="black",
            linewidth=lw,
            label=f"true-hybrid-containing hPSMs (n={len(data)})",
        )
        for datum in data[1:]:
            _ = ax.vlines(datum, ymax - rug_height, ymax, color="black", linewidth=lw)
        _ = ax.legend(loc="center left", bbox_to_anchor=(1.1, 0.9), frameon=False)
        _ = ax.set_ylim(0, ymax)
        _ = ax.set_title(self.name)
        return ax

    def nonxcorr_score_plots(self) -> List[Axes]:
        _, axs = fig_setup(nrows=len(NONXCORR_PSM_SCORES))
        for idx, score in enumerate(NONXCORR_PSM_SCORES):
            logger.info(f"Plotting score: {score}")
            ax = score_plot(psms=self.psms, score=score, ax=axs[idx])
            _ = ax.set_title(self.name)
        finalize(axs)
        return axs

    def process_experiment(
        self,
        acceptance_method: str = Literal[NEOFUSION, Q_VAL],
        q_thresh: float = DEFAULT_Q_VAL_THRESH,
    ):
        # Score plots
        ax = self.xcorr_plot()
        save_fig(self.xcorr_plot_path)
        _ = self.nonxcorr_score_plots()
        save_fig(self.nonxcorr_score_plot_path)

        # Protein abundance
        prot_ab = ProteinAbundance.from_comet_psms(
            psms=CometPSM.from_txt(txt=self.native_assign_confidence_txt),
            q_val_thresh=q_thresh,
        )
        ax = prot_ab.plot()
        _ = ax.set_title(f"{self.name}\nq <= {q_thresh}")
        save_fig(self.prot_ab_plot_path)

        # Accept hybrids
        acceptance_method_to_fcn = {
            NEOFUSION: lambda psms: accept_hybrid_psms_via_neo_fusion(psms=psms),
            Q_VAL: lambda psms: accept_hybrid_psms_via_interpolated_q_value(
                psms=psms, q_val_threshold=q_thresh
            ),
        }
        acceptance_fcn = acceptance_method_to_fcn[acceptance_method]
        accepted_hybrid_psms = acceptance_fcn(psms=self.psms)

        # PSM pileups
        left_hybrid_pileup, right_hybrid_pileup, native_pileup = self.get_pileups(
            accepted_hybrid_psms=accepted_hybrid_psms, q_thresh=q_thresh
        )
        fig, axs = plot_psm_pileup(
            left_hybrid_pileup=left_hybrid_pileup,
            right_hybrid_pileup=right_hybrid_pileup,
            fasta=self.fasta_path,
            native_pileup=native_pileup,
        )
        save_fig(
            self.psm_plot_path(acceptance_method=acceptance_method),
        )
        psm_df = create_psm_pileup_evidene_df(
            psms=self.psms,
            left_hybrid_pileup=left_hybrid_pileup,
            right_hybrid_pileup=right_hybrid_pileup,
            native_pileup=native_pileup,
            protein_name_to_seq_map=Fasta(path=self.fasta_path).protein_name_to_seq_map,
        )
        accepted_hy_psm_spectra = set(psm.spectrum_uid for psm in accepted_hybrid_psms)
        psm_df["accepted"] = psm_df["spectrum_uid"].apply(
            lambda uid: uid in accepted_hy_psm_spectra
        )
        psm_df["left_prot_ab"] = psm_df["left_prot"].apply(
            lambda prot: prot_ab.get_ab(protein=prot)
        )
        psm_df["right_prot_ab"] = psm_df["right_prot"].apply(
            lambda prot: prot_ab.get_ab(protein=prot)
        )
        psm_df.to_csv(
            self.pileup_df_path(acceptance_method=acceptance_method),
            index=False,
        )

        # Junction analysis
        # Get true hybrid junctions
        true_hybrid_jcts = get_true_hybrid_junctions(fasta=self.fasta_path)
        true_hybrid_jct_strs = set(str(jct) for jct in true_hybrid_jcts)
        jct_to_pos_hybrids = get_junction_to_positioned_hybrids_map(
            accepted_hybrid_psms=accepted_hybrid_psms, fasta=self.fasta_path
        )
        df = self.create_junction_df(
            jct_to_positioned_hybrids=jct_to_pos_hybrids,
            true_hybrid_jct_strs=true_hybrid_jct_strs,
        )
        df.to_csv(
            self.junction_df_path(acceptance_method=acceptance_method),
            index=False,
        )
        ax = self.junction_plot(
            df=df,
            acceptance_method=acceptance_method,
        )
        save_fig(self.junction_plot_path(acceptance_method=acceptance_method))

    def junction_df_path(self, acceptance_method: str) -> Path:
        return self.results_dir / f"{acceptance_method}_hybrid_junctions.csv"

    def load_junction_df(self, acceptance_method: str) -> pd.DataFrame:
        return pd.read_csv(self.junction_df_path(acceptance_method=acceptance_method))

    @staticmethod
    def create_junction_df(
        jct_to_positioned_hybrids: Dict[str, List[PositionedHybrid]],
        true_hybrid_jct_strs: Set[str],
    ):
        df = []
        for jct, hybrids in jct_to_positioned_hybrids.items():
            spectra = [hy.spectrum_uid for hy in hybrids]
            df.append(
                [
                    jct,
                    len(hybrids),
                    jct in true_hybrid_jct_strs,
                    list(set([hy.hyphen_seq for hy in hybrids])),
                    spectra,
                    len(spectra),
                ]
            )
        df = pd.DataFrame(
            df,
            columns=[
                "junction",
                "support_count",
                "true_jct",
                "hybrid_seqs",
                "spectra",
                "num_spectra",
            ],
        )
        df.sort_values(by="support_count", ascending=False, inplace=True)
        return df

    def junction_plot(
        self,
        df: pd.DataFrame,
        acceptance_method: str,
    ) -> Axes:
        # Create dataframe
        _, axs = fig_setup()
        ax = axs[0]
        _ = sns.histplot(
            data=df,
            x="support_count",
            kde=True,
            ax=ax,
        )
        _ = sns.rugplot(
            df[df.true_jct].support_count,
            ax=ax,
            height=0.05,
            color="red",
            label="Found true hybrid junctions",
        )

        finalize(ax)
        set_title_axes_labels(
            ax=ax,
            xlabel="Hybrid junctions",
            ylabel="Number accepted hybrids\nsupporting junction",
            title=self.plot_title(acceptance_method=acceptance_method),
        )
        return ax

    def plot_title(self, acceptance_method: str) -> Path:
        return f"{self.name}\nAccept={acceptance_method}"

    def pileup_df_path(self, acceptance_method: str) -> Path:
        return (
            self.results_dir
            / f"{acceptance_method}_accepted_hybrids_pileup_evidence.csv"
        )

    def psm_plot_path(self, acceptance_method: str) -> Path:
        return self.results_dir / f"{acceptance_method}_accepted_hybrids_psm_plot.png"

    @property
    def prot_ab_plot_path(self):
        return self.results_dir / "protein_abundance.png"

    @property
    def nonxcorr_score_plot_path(self):
        return self.results_dir / "non_xcorr_score_plots.png"

    @property
    def xcorr_plot_path(self):
        return self.results_dir / "xcorr_plot.png"

    def junction_plot_path(self, acceptance_method: str):
        return (
            self.results_dir / f"{acceptance_method}_hybrid_junction_support_plot.png"
        )


def psm_score_plots(
    df: pd.DataFrame, title: str, ppm_tol: Optional[float] = None
) -> Tuple[Figure, List[Axes]]:
    scores = [
        "xcorr",
        "prop_intensity_supported",
        "prop_prefixes_supported",
        "prop_suffixes_supported",
        "mz_ppm_diff",
    ]
    if ppm_tol is not None:
        df = df[np.abs(df.mz_ppm_diff) <= ppm_tol].copy()
    fig, axs = fig_setup(nrows=len(scores))
    for idx, score in enumerate(scores):
        # _, axs = fig_setup()
        # ax = axs[0]
        ax = axs[idx]
        for name, group in df.groupby("type"):
            _ = score_histogram(psms_by_type={name: group}, score=score, ax=ax)
            if score == "xcorr" and name == NAT_TARGET:
                add_qvalue_interpolator_to_xcorr_plot(native_psms=group, ax=ax)
    if ppm_tol is not None:
        title = f"{title}\n(restricted to <=20 PPM PSMs)"
    _ = fig.suptitle(title)
    finalize(axs)
    return fig, axs


def compare_psms_to_true_hybrids(
    psms: List[SpectrumPSMs], results_dir: Optional[Path] = None
) -> List[SpectrumPSMs]:
    """Returns true hybrid-containing PSMs"""
    true_hybrid_containing_psms = get_true_hybrid_containing_psms(psms=psms)
    logger.info(
        f"There are {len(true_hybrid_containing_psms)} hybrid PSMs that exactly match a true hybrid sequence"
    )
    # Save true hybrid containing PSMs
    if results_dir is not None:
        data_to_save = defaultdict(list)
        for psm in true_hybrid_containing_psms:
            data_to_save[psm.hybrid_seq].append(
                {psm.spectrum_uid: list(psm.hybrid_hyphen_seqs)}
            )
        to_json(
            data=dict(data_to_save),
            path=results_dir / "true_hybrid_exact_matches_found_with_spectra.json",
        )
        to_json(
            data=list(data_to_save.keys()),
            path=results_dir / "true_hybrid_exact_matches_found.json",
        )
    return true_hybrid_containing_psms


def get_true_hybrid_containing_psms(
    psms: List[SpectrumPSMs], true_hybrids: Path = TRUE_HYBRIDS_PATH
) -> List[SpectrumPSMs]:
    true_hybrids = TrueHybrid.load(path=true_hybrids)
    true_hybrid_seqs = set([hy.seq for hy in true_hybrids])
    return [
        psm
        for psm in [x for x in psms if x.has_hybrid]
        if psm.hybrid_seq in true_hybrid_seqs
    ]


def add_qvalue_interpolator_to_xcorr_plot(
    native_psms: Union[List[CometPSM], pd.DataFrame],
    ax: Axes,
):
    native_q_interpolator = fit_xcorr_to_qval_interpolator(
        psms=native_psms,
    )
    ax_copy = ax.twinx()
    xmin, xmax = ax.get_xlim()
    x_new = np.linspace(xmin, xmax, 500)
    _ = ax_copy.plot(x_new, native_q_interpolator(x_new), "r--", label="Native q-value")
    ax_copy.set_yscale("log")  # set y-axis to log10 scale
    ax_copy.set_ylabel("Native log10(q-value)", color="tab:red")
    ax_copy.tick_params(axis="y", labelcolor="tab:red")


def xcorr_plot(xcorr_by_type: Dict[str, List[float]], ax: Axes) -> Axes:
    # _, axs = fig_setup()
    # ax = axs[0]

    # Hybrid targets
    _ = score_histogram(psms_by_type=xcorr_by_type, score="xcorr", ax=ax)

    # Add q-value
    if NAT_TARGET in xcorr_by_type:
        native_q_interpolator = fit_xcorr_to_qval_interpolator(
            psms=xcorr_by_type[NAT_TARGET]
        )
        ax_copy = ax.twinx()
        xmin, xmax = ax.get_xlim()
        x_new = np.linspace(xmin, xmax, 500)
        _ = ax_copy.plot(
            x_new, native_q_interpolator(x_new), "r--", label="Native q-value"
        )
        ax_copy.set_yscale("log")  # set y-axis to log10 scale
        ax_copy.set_ylabel("Native log10(q-value)", color="tab:red")
        ax_copy.tick_params(axis="y", labelcolor="tab:red")
    set_title_axes_labels(ax=ax, xlabel="xcorr", ylabel="Density")
    finalize(ax)
    return ax


def accept_hybrids_that_beat_natives(psms: List[SpectrumPSMs]) -> List[SpectrumPSMs]:
    accepted_psms = SpectrumPSMs.get_winners_by_xcorr(psms=psms, psm_type=HYBRID)
    logger.info(
        f"Number of accepted hybrid PSMs that beat natives by XCorr: {len(accepted_psms)}"
    )
    return accepted_psms


def accept_hybrid_psms_via_neo_fusion(
    psms: List[SpectrumPSMs],
    q_vals: List[float] = DEFAULT_Q_RANGE,
    score_deltas: List[float] = DEFAULT_SCORE_CHANGE_RANGE,
    fpr_threshold: float = DEFAULT_FPR,
) -> List[SpectrumPSMs]:
    best_iteration, accepted_psms = NeoFusionRunner(
        q_vals=q_vals,
        score_deltas=score_deltas,
        fpr_threshold=fpr_threshold,
    ).run_neofusion(psms=psms)
    logger.info(f"Best NeoFusion iteration had:\n{best_iteration.info}")
    return accepted_psms


def hybrid_support_plot(df: pd.DataFrame, title: str) -> Axes:
    assert "support" in df.columns
    assert "true" in df.columns
    fig, axs = fig_setup(w=8)
    ax = axs[0]
    s = 7
    tmp = df[~df.true]
    _ = sns.scatterplot(
        x=tmp.index.to_list(),
        y=tmp.support.to_list(),
        color="blue",
        s=s,
        ax=ax,
    )
    tmp = df[df.true]
    _ = sns.scatterplot(
        x=tmp.index.to_list(),
        y=tmp.support.to_list(),
        color="red",
        marker="X",
        s=2 * s,
        label="'True' hybrids",
        ax=ax,
    )
    set_title_axes_labels(
        ax=ax,
        title=title,
        xlabel="Accepted hybrids (sorted by support)",
        ylabel="Num spectra supporting hybrid",
    )
    found_trues_str = "\n".join(set(df[df.true].seq))
    ax.text(
        1.02,
        0.5,
        f"Found trues:\n{found_trues_str}",
        transform=ax.transAxes,
        va="center",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white"),
    )
    finalize(axs)
    return ax


def psm_score_histogram(ax, data, label):
    _ = sns.histplot(
        data,
        # element="step",
        kde=True,
        stat="density",
        common_norm=False,
        ax=ax,
        label=f"{label} (n={len(data)})",
        # fill=False,
        alpha=0.4,
    )


def plot_native_vs_hybrid_scores(
    psms: List[SpectrumPSMs],
    # true_hybrid_seqs: Set[str],
    score: str = XCORR,
    ax: Optional[Axes] = None,
) -> Axes:

    # Get SpectrumPSM objects with both a native and a hybrid target PSM
    native_and_hybrid_psms = [
        psm
        for psm in psms
        if (psm.native_target is not None) and (psm.hybrid_target is not None)
    ]

    if ax is None:
        _, axs = fig_setup()
        ax = axs[0]
    s = 7
    sns.scatterplot(
        x=[getattr(psm.native_target, score) for psm in native_and_hybrid_psms],
        y=[getattr(psm.hybrid_target, score) for psm in native_and_hybrid_psms],
        s=s,
        marker="o",
        color="blue",
        label=f"n={len(native_and_hybrid_psms)}",
        ax=ax,
    )
    # sns.scatterplot(
    #     data=data,
    #     x=f"native_target_{score}",
    #     y=f"hybrid_target_{score}",
    #     s=2 * s,
    #     marker="x",
    #     color="red",
    #     label=f"True hybrid seqs (n={data.shape[0]})",
    #     ax=ax,
    #     linewidth=1.5,
    # )
    plot_line(ax=ax, label="y=x")
    set_title_axes_labels(
        ax=ax,
        xlabel=f"Native target {score}",
        ylabel=f"Hybrid target\n{score}",
    )
    return ax


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


def filter_to_top_n_highest_precursor_intensity_psms_per_mz(
    psms: List[SpectrumPSMs], n: int
) -> List[SpectrumPSMs]:
    mz_to_psms = defaultdict(list)
    for psm in psms:
        mz_to_psms[psm.spectrum.precursor_mz].append(psm)
    filtered_psms = []
    for mz, psms in mz_to_psms.items():
        psms = sorted(
            psms, key=lambda psm: psm.spectrum.precursor_abundance, reverse=True
        )
        filtered_psms.extend(psms[:n])
    return filtered_psms


def precursor_mz_plots(
    psms: List[SpectrumPSMs],
    sample: str,
    # plot_dir: Path
):
    mz_to_psms = defaultdict(list)
    for psm in psms:
        mz_to_psms[psm.spectrum.precursor_mz].append(psm)
    mz_to_psms = dict(mz_to_psms)

    num_supporting_spectra = Counter(len(psms) for mz, psms in mz_to_psms.items())
    _, axs = fig_setup(
        # nrows=1, ncols=2,
        w=10
    )
    ax = axs[0]
    _ = plot_sorted_1d_data(
        data=num_supporting_spectra,
        ax=ax,
        sort_idx=0,
        pt_labels=num_supporting_spectra,
        ax_labels=True,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="x = num supporting spectra",
        ylabel="Number of m/z's with num\nsupporting spectra = x",
        title=(
            f"sample = {sample}\n" + f"Num unique m/z values = {len(mz_to_psms)}\n"
            f"Num spectra = {sum(len(psms) for psms in mz_to_psms.values())}"
        ),
    )
    finalize(axs)
    # save_fig(
    #     path=plot_dir / f"number_spectra_supporting_each_uniq_mz.png",
    # )

    _, axs = fig_setup()
    cdf = ECDF(list(num_supporting_spectra.values()))
    ax = axs[0]
    _ = ax.plot(
        cdf.x,
        1 - cdf.y,
        "o-",
        ms=2,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="x = num supporting spectra",
        ylabel="Proportion of m/z's with\nnum supporting spectra >= x",
        title=(
            f"sample = {sample}\n" + f"Num unique m/z values = {len(mz_to_psms)}\n"
            f"num spectra = {sum(len(psms) for psms in mz_to_psms.values())}"
        ),
    )
    finalize(axs)


def create_psm_pileup_evidene_df(
    psms: List[SpectrumPSMs],
    left_hybrid_pileup: Dict,
    right_hybrid_pileup: Dict,
    native_pileup: Dict,
    protein_name_to_seq_map: Dict[str, str],
):
    df = []
    for psm in psms:
        if not psm.has_hybrid:
            continue
        for hy in psm.hybrids:
            # A single hybrid A-B can appear in multiple places if either A or B appears in multiple places
            hy_positions = get_positions_of_hybrid(
                hybrid=hy,
                protein_name_to_seq_map=protein_name_to_seq_map,
            )
            for hy_pos in hy_positions:
                left_hybrid_support = get_support_in_pileup(
                    pos=hy_pos.left, pileup=left_hybrid_pileup
                )
                right_hybrid_support = get_support_in_pileup(
                    pos=hy_pos.right, pileup=right_hybrid_pileup
                )
                left_native_support = get_support_in_pileup(
                    pos=hy_pos.left, pileup=native_pileup
                )
                right_native_support = get_support_in_pileup(
                    pos=hy_pos.right, pileup=native_pileup
                )
                df.append(
                    [
                        psm.spectrum_uid,
                        hy.hyphen_seq,
                        hy_pos.left.protein,
                        hy_pos.right.protein,
                        np.mean(left_hybrid_support),
                        np.mean(right_hybrid_support),
                        np.mean(left_native_support),
                        np.mean(right_native_support),
                    ]
                )
    df = pd.DataFrame(
        df,
        columns=[
            "spectrum_uid",
            "hybrid_seq",
            "left_prot",
            "right_prot",
            "left_hybrid_mean_support",
            "right_hybrid_mean_support",
            "left_native_mean_support",
            "right_native_mean_support",
        ],
    )
    return df


def get_hybrids_for_hybrid_psms(
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    hybrid_psms: List[CometPSM],
) -> List[HybridPeptide]:
    hybrids = []
    for psm in hybrid_psms:
        hybrids.extend(seq_to_hybrids_map[psm.seq])
    return hybrids


def get_hybrid_psm_to_junctions_map(
    hybrid_psms: List[CometPSM], protein_name_to_seq_map: Dict[str, str]
) -> pd.DataFrame:
    # Get support by junction point
    hy_psm_uid_to_jcts = defaultdict(list)
    for psm in hybrid_psms:
        if not psm.has_hybrid:
            continue
        for hy in psm.hybrids:
            # A single hybrid A-B can appear in multiple places if either A or B appears in multiple places
            for pos in get_positions_of_hybrid(
                hybrid=hy,
                protein_name_to_seq_map=protein_name_to_seq_map,
            ):
                hy_psm_uid_to_jcts[psm.spectrum_uid].append(
                    HybridJunction.from_hybrid_position(pos=pos)
                )
    return hy_psm_uid_to_jcts


def get_hybrid_junctions(
    hybrids: List[HybridPeptide],
    protein_name_to_seq_map: Dict[str, str],
) -> List[HybridJunction]:
    jcts = flatten_list_of_lists(
        hybrid.get_junctions(protein_name_to_seq_map=protein_name_to_seq_map)
        for hybrid in hybrids
    )
    return jcts


def get_true_hybrid_junctions(
    protein_name_to_seq_map: Dict[str, str], true_hybrids: Path = TRUE_HYBRIDS_PATH
) -> List[HybridJunction]:
    true_hybrids = TrueHybrid.load(path=true_hybrids)
    return get_hybrid_junctions(
        hybrids=true_hybrids, protein_name_to_seq_map=protein_name_to_seq_map
    )


def hybrid_support_analysis(
    hybrid_psms: List[CometPSM],
    results_dir: Path,
    sample_name: str,
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    true_seqs: Set[str],
    acceptance_method: str,
    protein_name_to_seq_map: Dict[str, str],
):
    prefix = f"{sample_name}_{acceptance_method}_hpsms"
    CometPSM.save_psms(
        psms=hybrid_psms,
        path=results_dir / f"{prefix}_n{len(hybrid_psms)}.json",
    )
    df = get_hybrid_support_df(
        hybrid_psms=hybrid_psms,
        seq_to_hybrids_map=seq_to_hybrids_map,
        protein_name_to_seq_map=protein_name_to_seq_map,
    )
    df["true"] = df.seq.apply(lambda seq: seq in true_seqs)
    df.to_csv(
        results_dir / f"{prefix}_support_df.csv",
        index=False,
    )
    ax = hybrid_support_plot(df=df, title=f"{sample_name}\n{acceptance_method}")
    save_fig(path=results_dir / f"{prefix}_support_plot.png")


def get_hybrid_support_df(
    hybrid_psms: List[CometPSM],
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    protein_name_to_seq_map: Dict[str, str],
) -> pd.DataFrame:
    hybrid_to_psms = defaultdict(list)
    for psm in hybrid_psms:
        hybrids = seq_to_hybrids_map[psm.seq]
        for hybrid in hybrids:
            hybrid_to_psms[
                hybrid.to_str(protein_name_to_seq_map=protein_name_to_seq_map)
            ].append(psm)

    df = []
    for hybrid, psms in hybrid_to_psms.items():
        assert all(psm.seq == psms[0].seq for psm in psms)
        data = {
            "hybrid": hybrid,
            "supporting_spectra": [psm.spectrum_uid for psm in psms],
            "seq": psms[0].seq,
        }
        df.append(data)
    df = pd.DataFrame(df)
    df["support"] = df.supporting_spectra.apply(lambda x: len(x))
    df.sort_values("support", ascending=False, inplace=True, ignore_index=True)
    return df


def get_junction_to_positioned_hybrids_map(
    accepted_hybrid_psms: List[SpectrumPSMs], fasta: Path
) -> Dict[str, List[PositionedHybrid]]:
    spectrum_to_psms = {psm.spectrum_uid: psm for psm in accepted_hybrid_psms}
    prot_name_to_seq = Fasta(path=fasta).protein_name_to_seq_map
    accepted_hybrid_psms = accept_hybrid_psms_via_neo_fusion(psms=accepted_hybrid_psms)

    scan_to_hybrid_jcts = get_hybrid_psm_to_junctions_map(
        hybrid_psms=accepted_hybrid_psms, fasta=fasta
    )
    logger.info(
        f"<number of possible hybrid junctions explaining scan> : <num scans>\n{dict(Counter(len(jcts) for jcts in scan_to_hybrid_jcts.values()))}"
    )
    # Get <junction> : <scans supporting junction>
    jct_to_scans = defaultdict(list)
    for scan, jcts in scan_to_hybrid_jcts.items():
        for jct in jcts:
            jct_to_scans[str(jct)].append(scan)

    # Get <junction> : <positioned hybrid>
    jct_to_positioned_hybrids = defaultdict(list)
    for jct, spectrum_uids in jct_to_scans.items():
        for spectrum_uid in spectrum_uids:
            found_hybrid = get_junction_matching_hybrid(
                spectrum_uid=spectrum_uid,
                spectrum_psms=spectrum_to_psms[spectrum_uid],
                protein_name_to_seq_map=prot_name_to_seq,
                target_junction=jct,
            )
            assert found_hybrid is not None
            jct_to_positioned_hybrids[jct].append(found_hybrid)

    return jct_to_positioned_hybrids


def extract_spectrum_command(
    spectrum_scan: int,
    spectrum_idx: int,
    local_mzml_path: Path,
    seq: str,
    out_dir: Path = Path("./"),
    container_mzml_path: Optional[str] = None,
):
    """
    Create the command to extract a spectrum from an mzML file using `wine msconvert`.

    - container_mzml_path : should be the Path to the mzML file in the Docker container.
        Defaults to `local_mzml_path.relative_to(DATA_DIR.absolute())`
    """

    if container_mzml_path is None:
        container_mzml_path = local_mzml_path.relative_to(DATA_DIR.absolute())
    mzml = Mzml(mzml=local_mzml_path)
    file_name = f"mzml{mzml.name}_seq{seq}_idx{spectrum_idx}_scan{spectrum_scan}.mgf"
    return f'wine msconvert {local_mzml_path} --filter "index {spectrum_idx}" --outfile {file_name} -o {out_dir} --mgf'


def create_extract_spectra_from_mzml_bash_script(
    spectra_psms: List[SpectrumPSMs],
    local_mzml_path: Union[Path, str],
    container_out_dir: Union[Path, str],
    local_script_out_dir: Union[Path, str],
    container_mzml_path: Optional[str] = None,
):
    script_lines = [
        extract_spectrum_command(
            spectrum_scan=psm.spectrum.scan,
            spectrum_idx=psm.spectrum.mzml_index,
            local_mzml_path=local_mzml_path,
            container_mzml_path=container_mzml_path,
            out_dir=container_out_dir,
            seq=psm.hybrid_seq,
        )
        for psm in spectra_psms
    ]
    mzml = Mzml(mzml=local_mzml_path)
    write_new_line_separated_file(
        lines=script_lines,
        path=Path(local_script_out_dir)
        / f"{mzml.name}_get_true_hybrid_supporting_spectra.sh",
    )


@click.command(
    name="spectrum-psms",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    Create SpectrumPSMs objects for Hypedsearch experiment
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    multiple=True,
    help="Paths to one or more Hypedsearch JSON configs",
)
@click.option(
    "--min_side_len",
    "-m",
    type=int,
    required=False,
    show_default=True,
    default=DEFAULT_MIN_SIDE_LEN,
    help="Minimum hybrid side length",
)
def cli_create_spectrum_psms(
    config: Tuple[Path, ...],
    min_side_len: int,
):
    for hs_config in config:
        # Creating SpectrumPSMs
        spectra_psms = SpectrumPSMs.from_hs_config(
            hs_config=hs_config, min_side_len=min_side_len, save=True
        )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_create_spectrum_psms)
    cli()
