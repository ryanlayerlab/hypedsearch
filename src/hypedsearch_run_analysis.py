from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass, field, replace
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple, Union
from venv import logger

import click
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.interpolate import PchipInterpolator

from src.constants import (
    ASSIGN_CONFIDENCE,
    DATA_DIR,
    DECOY,
    DEFAULT_FPR,
    DEFAULT_JCT_LEN,
    DEFAULT_MIN_SIDE_LEN,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_Q_RANGE,
    DEFAULT_Q_THRESHOLD,
    DEFAULT_SCORE_CHANGE_RANGE,
    GOOD_HYBRID_BAD_NATIVE,
    HUMAN_PROTEOME,
    HY_TARGET,
    HYBRID_PEPTIDES_NAME,
    NAT_DECOY,
    NAT_TARGET,
    NEOFUSION,
    PSMS_DF_NAME,
    Q_VAL,
    SPECTRA_DF_NAME,
    SPECTRUM_PSMS_NAME,
    TARGET,
    TRUE_HYBRIDS_PATH,
    XCORR,
)
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import HypedsearchRunConfig, TrueHybrid, get_seq_to_hybrids_map
from src.mass_spectra import Mzml, Spectrum, organize_by_spectrum_uid
from src.peptides_and_ions import Fasta, ProteinRange
from src.plot_utils import (
    fig_setup,
    finalize,
    plot_line,
    plot_sorted_1d_data,
    save_fig,
    set_title_axes_labels,
)
from src.psm import PSM, CometPSM, ProteinAbundance, spectrum_peptide_plot
from src.utils import (
    PathType,
    flatten_list_of_lists,
    get_positions_of_subseq_in_seq,
    load_json,
    log_params,
    log_time,
    setup_logger,
    to_json,
    write_new_line_separated_file,
)

ALLOWED_ACCEPTANCE_METHODS = ["good_hybrid_bad_native"]


def num_hybrids_per_seq_histplot(seq_to_hybrids_map):
    df = pd.DataFrame(
        data=[
            (seq, len(hybrids), [hy.seq_with_hyphen for hy in hybrids])
            for seq, hybrids in seq_to_hybrids_map.items()
        ],
        columns=["seq", "num_hybrids", "hybrids"],
    )
    fig, axs = fig_setup()
    ax = axs[0]
    _ = sns.histplot(df.num_hybrids, ax=ax)
    set_title_axes_labels(
        ax=ax,
        xlabel="Number of hybrid explanations\nper hybrid sequence",
        ylabel="Count",
    )
    finalize(axs)
    return df, fig, axs


def group_hybrid_psms_by_junction(
    jct_len: int,
    hybrid_psms: List[CometPSM],
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    protein_name_to_seq_map: Dict[str, str],
) -> Dict[str, List[CometPSM]]:
    jct_to_psms = defaultdict(list)
    for psm in hybrid_psms:
        for hybrid in seq_to_hybrids_map[psm.seq]:
            for jct in hybrid.get_junctions(
                protein_name_to_seq_map=protein_name_to_seq_map,
                jct_len=jct_len,
            ):
                jct_to_psms[jct].append(psm)
    return dict(jct_to_psms)


def group_hybrid_psms_by_hybrid_str(
    hybrid_psms: List[CometPSM],
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
) -> Dict[str, List[CometPSM]]:
    hybrid_to_psms = defaultdict(list)
    for psm in hybrid_psms:
        for hybrid in seq_to_hybrids_map[psm.seq]:
            hybrid_to_psms[hybrid.to_str()].append(psm)
    return dict(hybrid_to_psms)


def get_top_psm_per_spectrum(
    psms: List[CometPSM],
) -> Dict[str, CometPSM]:
    top_psms = organize_by_spectrum_uid([psm for psm in psms if psm.num == 1])
    for spectrum_uid, psms in top_psms.items():
        assert len(psms) == 1, f"More than one top hybrid PSM for {spectrum_uid}"
        top_psms[spectrum_uid] = psms[0]
    return top_psms


@log_time()
def remove_native_psms(
    hybrid_psms: List[CometPSM],
    fasta: Union[str, Path, Fasta],
) -> List[CometPSM]:
    if isinstance(fasta, (str, Path)):
        fasta = Fasta(path=fasta)
    psm_seqs = set([psm.seq for psm in hybrid_psms])
    native_seqs = [seq for seq in psm_seqs if fasta.contains_seq(query_seq=seq)]
    return [psm for psm in hybrid_psms if psm.seq not in native_seqs]


@dataclass
class SpectrumPSMs:
    spectrum: Spectrum
    native_target: Optional[Union[CometPSM, PSM]] = None
    native_decoy: Optional[Union[CometPSM, PSM]] = None
    hybrid_target: Optional[Union[CometPSM, PSM]] = None

    @property
    def uid(self):
        return self.spectrum.uid

    def convert_comet_psms_to_psms(
        self, peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        return replace(
            self,
            native_target=(
                self.native_target.to_psm(
                    spectrum=self.spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
                if self.native_target is not None
                else None
            ),
            native_decoy=(
                self.native_decoy.to_psm(
                    spectrum=self.spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
                if self.native_decoy is not None
                else None
            ),
            hybrid_target=(
                self.hybrid_target.to_psm(
                    spectrum=self.spectrum, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
                )
                if self.hybrid_target is not None
                else None
            ),
        )

    def to_dict(self):
        return {
            "spectrum": self.spectrum.to_dict(),
            "native_target": (
                asdict(self.native_target) if self.native_target is not None else None
            ),
            "native_decoy": (
                asdict(self.native_decoy) if self.native_decoy is not None else None
            ),
            "hybrid_target": (
                asdict(self.hybrid_target) if self.hybrid_target is not None else None
            ),
        }

    @log_time()
    @staticmethod
    def save(psms: List["SpectrumPSMs"], path: Union[Path, str]):
        to_json(data=[psm.to_dict() for psm in psms], path=path)

    @classmethod
    def load(cls, path: Union[Path, str]):
        data = load_json(path=path)
        psms = []
        for item in data:
            psms.append(
                cls(
                    spectrum=Spectrum(**item["spectrum"]),
                    native_target=(
                        CometPSM(**item["native_target"])
                        if item["native_target"] is not None
                        else None
                    ),
                    native_decoy=(
                        CometPSM(**item["native_decoy"])
                        if item["native_decoy"] is not None
                        else None
                    ),
                    hybrid_target=(
                        CometPSM(**item["hybrid_target"])
                        if item["hybrid_target"] is not None
                        else None
                    ),
                )
            )
        return psms

    def to_rows(
        self, peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ) -> Dict[str, Any]:
        psm = self.convert_comet_psms_to_psms(peak_to_ion_ppm_tol=peak_to_ion_ppm_tol)
        rows = []
        if psm.native_target:
            row = psm.native_target.to_row()
            row["psm_type"] = "native_target"
            rows.append(row)
        if psm.native_decoy:
            row = psm.native_decoy.to_row()
            row["psm_type"] = "native_decoy"
            rows.append(row)
        if psm.hybrid_target:
            row = psm.hybrid_target.to_row()
            row["psm_type"] = "hybrid_target"
            rows.append(row)
        return rows

    def spectrum_peptide_plots(
        self, peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ):
        fig, axs = fig_setup(ncols=2)
        _ = spectrum_peptide_plot(
            spectrum=self.spectrum,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            seq=self.native_target.seq,
            ax=axs[0],
            title=f"Native target: {self.native_target.seq}\nxcorr={self.native_target.xcorr}, q={self.native_target.q_value}",
        )
        _ = spectrum_peptide_plot(
            spectrum=self.spectrum,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tol,
            seq=self.hybrid_target.seq,
            ax=axs[1],
            title=f"Hybrid target: {self.hybrid_target.seq}\nxcorr={self.hybrid_target.xcorr}, q={self.hybrid_target.q_value}",
        )
        fig.suptitle(
            f"Spectrum: {self.spectrum.uid}\nprecursor-intensity: {round(self.spectrum.precursor_intensity, 3)} (PPM tol={peak_to_ion_ppm_tol})"
        )


@dataclass
class ResultsAnalysis:
    hs_config: Union[str, Path, HypedsearchRunConfig]
    min_side_len: int = DEFAULT_MIN_SIDE_LEN
    remove_carbamidomethylation: bool = True
    native_targets: Dict[str, List[CometPSM]] = field(init=False)
    native_assign_conf: Dict[str, CometPSM] = field(init=False)
    native_decoys: Dict[str, List[CometPSM]] = field(init=False)
    hybrid_targets: Dict[str, List[CometPSM]] = field(init=False)

    def __post_init__(self):
        if isinstance(self.hs_config, (str, Path)):
            self.hs_config = HypedsearchRunConfig.from_json(path=self.hs_config)
        self.results_dir.mkdir(parents=True, exist_ok=True)

    @property
    def name(self) -> str:
        return self.hs_config.name

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
    def _native_output_txts(self) -> Dict[str, List[Path]]:
        return self.hs_config.get_expected_native_run_output_txts()

    @cached_property
    def q_value_interpolator(self) -> PchipInterpolator:
        return fit_xcorr_to_qval_interpolator(psms=self.native_assign_conf.values())

    @cached_property
    def spectrum_uid_to_spectrum(self) -> Dict[str, Spectrum]:
        return self.hs_config.spectrum_uid_to_spectrum

    @cached_property
    def all_spectrum_uids(self) -> List[str]:
        return list(
            set(self.native_targets.keys()).union(set(self.hybrid_targets.keys()))
        )

    @cached_property
    def top_native_targets(self) -> Dict[str, CometPSM]:
        return get_top_psm_per_spectrum(
            psms=flatten_list_of_lists(self.native_targets.values())
        )

    @cached_property
    def top_hybrid_targets(self) -> Dict[str, CometPSM]:
        return get_top_psm_per_spectrum(
            psms=flatten_list_of_lists(self.hybrid_targets.values())
        )

    @cached_property
    def top_native_decoys(self) -> Dict[str, CometPSM]:
        return get_top_psm_per_spectrum(
            psms=flatten_list_of_lists(self.native_decoys.values())
        )

    @cached_property
    def protein_name_to_seq_map(self) -> Dict[str, str]:
        fasta = Fasta(path=self.hs_config.fasta)
        return fasta.protein_name_to_seq_map

    @property
    def results_dir(self) -> Path:
        return self.hs_config.parent_output_dir / "results"

    def get_spectrum_psms_path(self, min_side_len: int) -> Path:
        return (
            self.results_dir
            / f"{self.name}.minSideLen={min_side_len}.spectrumPSMs.json"
        )

    def collect_outputs(self):
        self.native_targets = self._get_native_targets()
        self.native_assign_conf = self._get_native_assign_conf()
        self.native_decoys = self._get_native_decoys()
        self.hybrid_targets = self._get_hybrid_targets()
        logger.info("Loading spectra so they're quickly available...")
        _ = self.spectrum_uid_to_spectrum  # force loading spectra

    def get_spectrum_psms(self, spectrum_uid: str):
        return SpectrumPSMs(
            spectrum=self.spectrum_uid_to_spectrum[spectrum_uid],
            native_target=self.top_native_targets.get(spectrum_uid, None),
            native_decoy=self.top_native_decoys.get(spectrum_uid, None),
            hybrid_target=self.top_hybrid_targets.get(spectrum_uid, None),
        )

    def _get_native_targets(self) -> Dict[str, List[CometPSM]]:
        logger.info("Getting native target PSMs...")
        return CometPSM.from_txts(
            txts=self._native_output_txts[TARGET],
            by_spectrum=True,
        )

    def _get_native_decoys(self) -> Dict[str, List[CometPSM]]:
        logger.info("Getting native decoy PSMs...")
        return CometPSM.from_txts(
            txts=self._native_output_txts[DECOY],
            by_spectrum=True,
        )

    def _get_native_assign_conf(self) -> Dict[str, List[CometPSM]]:
        logger.info("Getting native assign confidence PSMs...")
        spectrum_to_psm = CometPSM.from_txt(
            txt=self._native_output_txts[ASSIGN_CONFIDENCE],
            by_spectrum=True,
        )
        for key, psms in spectrum_to_psm.items():
            assert len(psms) == 1, f"More than one assign confidence PSM for {key}"
            spectrum_to_psm[key] = psms[0]
        return spectrum_to_psm

    def _get_hybrid_targets(self) -> Dict[str, List[CometPSM]]:
        logger.info("Getting hybrid target PSMs...")
        psms = [
            psm
            for psm in CometPSM.from_txts(
                txts=self.hs_config.get_expected_hybrid_run_output_txts(),
                by_spectrum=False,
            )
            if psm.is_hybrid
        ]
        # Remove natives mistakenly labeled as hybrids. This can happen if the peptide
        # length is longer than the longest kmer in the kmer database.
        # Sadly, this takes a while: ~1.75 minutes for ~16,000 sequences
        logger.info("Removing native PSMs from hybrid PSMs...")
        psms = remove_native_psms(
            hybrid_psms=psms,
            fasta=self.hs_config.fasta,
        )

        seq_to_hybrids = get_seq_to_hybrids_map(
            seqs=set(psm.seq for psm in psms),
            db_path=self.hs_config.kmer_db_path,
            min_side_len=self.min_side_len,
            remove_carbamidomethylation=self.remove_carbamidomethylation,
        )
        hybrid_supported_psms = [psm for psm in psms if psm.seq in seq_to_hybrids]
        return organize_by_spectrum_uid(data=hybrid_supported_psms)

    def save_spectrum_psms(self) -> List[SpectrumPSMs]:
        psms = self.to_spectrum_psms()
        SpectrumPSMs.save(
            psms=psms, path=self.hs_config._results_dir / SPECTRUM_PSMS_NAME
        )
        return psms

    def load_spectrum_psms(self, by_uid: bool = False) -> List[SpectrumPSMs]:
        logger.info("Loading SpectrumPSMs objects...")
        psms = SpectrumPSMs.load(path=self.hs_config._results_dir / SPECTRUM_PSMS_NAME)
        if by_uid:
            uid_to_psms = {}
            for psm in psms:
                uid_to_psms[psm.spectrum.uid] = psm
            return uid_to_psms
        else:
            return psms

    def save_hybrid_peptides(self):
        data = [
            hy.to_dict()
            for hy in flatten_list_of_lists(self.seq_to_hybrids_map.values())
        ]
        to_json(data=data, path=self.hs_config._results_dir / HYBRID_PEPTIDES_NAME)

    def load_seq_to_hybrid_map(self):
        logger.info("Loading sequence-to-hybrid-peptides...")
        hy_peptides = [
            HybridPeptide(**hy)
            for hy in load_json(path=self.hs_config._results_dir / HYBRID_PEPTIDES_NAME)
        ]
        seq_to_hybrids_map = defaultdict(list)
        for hy in hy_peptides:
            seq_to_hybrids_map[hy.seq].append(hy)
        return dict(seq_to_hybrids_map)

    def load_psms_df(self) -> pd.DataFrame:
        logger.info("Loading PSM dataframe...")
        return pd.read_csv(self.hs_config._results_dir / PSMS_DF_NAME)

    def load_spectra_df(self) -> pd.DataFrame:
        return pd.read_csv(self.hs_config._results_dir / SPECTRA_DF_NAME)

    def load_stuff(self):
        psms_df = self.load_psms_df()
        spectra_df = self.load_spectra_df()
        seq_to_hybrids_map = self.load_seq_to_hybrid_map()
        psms = self.load_spectrum_psms()
        return psms, seq_to_hybrids_map, spectra_df, psms_df

    def set_q_values(self):
        logger.info(
            "Setting interpolated q-values for hybrid targets and native targets..."
        )
        for psms in self.hybrid_targets.values():
            for psm in psms:
                psm.q_value = float(self.q_value_interpolator(psm.xcorr))
        for psms in self.native_targets.values():
            for psm in psms:
                psm.q_value = float(self.q_value_interpolator(psm.xcorr))

    def get_protein_abundances(
        self, q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> ProteinAbundance:
        return ProteinAbundance.from_comet_psms(
            quality_psms=self.native_assign_conf.values(),
            q_threshold=q_threshold,
        )

    @log_time()
    def to_spectrum_psms(self) -> List[SpectrumPSMs]:
        logger.info("Converting all PSMs to SpectrumPSMs objects...")
        results = []
        for uid in self.spectrum_uid_to_spectrum.keys():
            results.append(self.get_spectrum_psms(spectrum_uid=uid))
        return results

    def xcorr_plot(self):
        psms = self.to_spectrum_psms()
        fig, axs = fig_setup(ncols=2)
        xcorr_plot(
            psms_by_type={
                NAT_TARGET: self.top_native_targets.values(),
                NAT_DECOY: self.top_native_decoys.values(),
                HY_TARGET: self.top_hybrid_targets.values(),
            },
            ax=axs[0],
        )
        plot_native_vs_hybrid_scores(psms=psms, ax=axs[1])
        finalize(axs)
        save_fig(
            fig=fig,
            title=f"{self.name}",
            path=self.hs_config.plots_dir / f"{self.name}.xcorr.png",
        )

    def create_plots(
        self,
    ):
        # Xcorr plot
        self.xcorr_plot()

        # Num of hybrid explanations per hybrid seq, e.g., hybrid seq = ABC with hybrid
        # explanations A-BC, AB-C
        df, fig, axs = num_hybrids_per_seq_histplot(
            seq_to_hybrids_map=self.seq_to_hybrids_map
        )
        save_fig(
            fig=fig,
            title=f"{self.name}",
            path=self.hs_config.plots_dir / f"{self.name}.numHybridsPerHybridSeq.png",
        )

        # Spectra plots
        fig, axs = Spectrum.plot_spectra_info(
            spectra=list(self.hs_config.spectrum_uid_to_spectrum.values()),
            add_counts=False,
        )
        save_fig(
            fig=fig,
            title=f"{self.name}",
            path=self.hs_config.plots_dir / f"{self.name}.spectraHistograms.png",
        )


@log_params
def process_hs_config(
    hs_config: Path,
    acceptance_method_name: Literal[
        GOOD_HYBRID_BAD_NATIVE, NEOFUSION
    ] = GOOD_HYBRID_BAD_NATIVE,
    q_threshold: Optional[float] = DEFAULT_Q_THRESHOLD,
    jct_len: int = DEFAULT_JCT_LEN,
    min_hybrid_side_len: int = DEFAULT_MIN_SIDE_LEN,
    remove_methylation: bool = True,
):
    results = ResultsAnalysis(
        hs_config=hs_config,
        min_side_len=min_hybrid_side_len,
        remove_carbamidomethylation=remove_methylation,
    )
    spectrum_psms_path = results.get_spectrum_psms_path(
        min_side_len=min_hybrid_side_len
    )
    if not spectrum_psms_path.exists():
        # Collect all the PSMs, and set q-values
        results.collect_outputs()
        results.set_q_values()

        # Create SpectrumPSMs
        psms = results.to_spectrum_psms()
        SpectrumPSMs.save(psms=psms, path=spectrum_psms_path)
        seq_to_hybrids_map = results.seq_to_hybrids_map

    else:
        logger.info("SpectrumPSMs already exist so loading them...")
        psms = SpectrumPSMs.load(path=spectrum_psms_path)
        seq_to_hybrids_map = get_seq_to_hybrids_map(
            seqs=set(psm.hybrid_target.seq for psm in psms if psm.hybrid_target),
            db_path=results.hs_config.kmer_db_path,
            min_side_len=min_hybrid_side_len,
            remove_carbamidomethylation=remove_methylation,
        )

    # Accept hybrids
    logger.info("Accepting hybrids...")
    identifier = AcceptanceMethod.get_acceptance_identifier(
        acceptance_method_name=acceptance_method_name,
        min_hybrid_side_len=min_hybrid_side_len,
        jct_len=jct_len,
    )
    if acceptance_method_name == GOOD_HYBRID_BAD_NATIVE:
        accepted_hybrids = (
            AcceptanceMethod.acccept_good_hybrids_with_no_good_native_explanation(
                psms=psms, q_threshold=q_threshold
            )
        )
    elif acceptance_method_name == NEOFUSION:
        accepted_hybrids = AcceptanceMethod.accept_hybrids_via_neofusion(
            psms=psms,
        )
    else:
        raise ValueError("Unrecognized acceptance method")

    # Plots
    logger.info("Plotting...")
    # Xcorr plot
    fig, axs = fig_setup(ncols=2)
    xcorr_plot(
        psms_by_type={
            NAT_TARGET: [
                psm.native_target for psm in psms if psm.native_target is not None
            ],
            NAT_DECOY: [
                psm.native_decoy for psm in psms if psm.native_decoy is not None
            ],
            HY_TARGET: [
                psm.hybrid_target for psm in psms if psm.hybrid_target is not None
            ],
        },
        ax=axs[0],
    )
    plot_xcorr_of_accepted_vs_not_accepted_hybrids(
        ax=axs[1],
        all_psms=psms,
        accepted_hybrid_psms=accepted_hybrids,
    )
    finalize(axs)
    save_fig(
        fig=fig,
        title=f"{results.name}",
        path=results.hs_config.plots_dir / f"{results.name}.{identifier}.xcorr.png",
    )

    # Num of hybrid explanations per hybrid seq, e.g., hybrid seq = ABC with hybrid
    # explanations A-BC, AB-C
    df, fig, axs = num_hybrids_per_seq_histplot(seq_to_hybrids_map=seq_to_hybrids_map)
    save_fig(
        fig=fig,
        title=f"{results.name}",
        path=results.hs_config.plots_dir / f"{results.name}.numHybridsPerHybridSeq.png",
    )

    # Spectra plots
    fig, axs = Spectrum.plot_spectra_info(
        spectra=list(results.hs_config.spectrum_uid_to_spectrum.values()),
        add_counts=False,
    )
    save_fig(
        fig=fig,
        title=f"{results.name}",
        path=results.hs_config.plots_dir / f"{results.name}.spectraHistograms.png",
    )

    # Junction analysis
    logger.info("Junction analysis...")
    df = get_junction_df(
        hybrid_psms=accepted_hybrids,
        protein_name_to_seq_map=results.protein_name_to_seq_map,
        seq_to_hybrids_map=seq_to_hybrids_map,
        spectrum_uid_to_spectrum_map=results.spectrum_uid_to_spectrum,
        jct_len=jct_len,
    )
    fig, axs = junction_plots(df=df)
    save_fig(
        fig=fig,
        title=f"{results.name}\n{identifier}",
        path=results.hs_config.plots_dir / f"{results.name}.{identifier}.png",
    )
    df.to_csv(results.results_dir / f"{results.name}.{identifier}.csv", index=False)


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

    # Sort each protein's pileup by location and turn the pileup into a dict from a default dict
    for prot, pile in pileup.items():
        pileup[prot] = dict(sorted(pile.items()))
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


# @dataclass
# class PSMPileup:
#     @staticmethod
def get_native_pileup(psms: List[CometPSM], fasta: Path = HUMAN_PROTEOME):
    return get_pileup_from_positions(
        positions=align_psms_to_proteome(
            psms=psms,
            fasta=fasta,
        )
    )


def plot_num_hybrid_explanations_per_psm(
    psms: List[CometPSM],
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    title: Optional[str] = None,
):
    df = pd.DataFrame(
        [
            [len(psm.seq), len(seq_to_hybrids_map[psm.seq]), psm.spectrum_uid]
            for psm in psms
        ],
        columns=["seq_len", "num_hybrid_explanations", "spectrum_uid"],
    )
    fig, axs = fig_setup(ncols=2)
    ax = axs[0]
    _ = sns.scatterplot(x=df.seq_len, y=df.num_hybrid_explanations, ax=ax, s=7)
    set_title_axes_labels(
        ax=ax,
        xlabel="PSM sequence length",
        ylabel="Number of possible\nhybrid explanations per PSM",
    )
    ax = axs[1]
    _ = sns.histplot(
        data=df.num_hybrid_explanations,
    )
    set_title_axes_labels(
        ax=ax,
        xlabel="Number of possible\nhybrid explanations per PSM",
        ylabel="Count",
    )
    if title is not None:
        fig.suptitle(title)
    finalize(axs)
    return df, fig, axs


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
                label="left-side support",
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
                label="right-side support",
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
                label="native support",
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
    q_threshold: Optional[float] = None,
):
    native_q_interpolator = fit_xcorr_to_qval_interpolator(
        psms=native_psms,
    )
    ax_copy = ax.twinx()
    xmin, xmax = ax.get_xlim()
    x_new = np.linspace(xmin, xmax, 500)
    _ = ax_copy.plot(x_new, native_q_interpolator(x_new), "r--", label="Native q-value")
    if q_threshold is not None:
        _ = ax_copy.axhline(
            y=q_threshold, color="red", linestyle="--", label=f"q={q_threshold}"
        )
    ax_copy.set_yscale("log")  # set y-axis to log10 scale
    ax_copy.set_ylabel("Native log10(q-value)", color="tab:red")
    ax_copy.tick_params(axis="y", labelcolor="tab:red")


def xcorr_plot(
    psms_by_type: Dict[str, List[CometPSM]], ax: Optional[Axes] = None
) -> Axes:
    if ax is None:
        fig, axs = fig_setup()
        ax = axs[0]
    _ = score_histogram(psms_by_type=psms_by_type, score=XCORR, ax=ax)
    if NAT_TARGET in psms_by_type:
        add_qvalue_interpolator_to_xcorr_plot(
            native_psms=psms_by_type[NAT_TARGET],
            ax=ax,
        )
    set_title_axes_labels(ax=ax, xlabel="xcorr", ylabel="Density")
    finalize(ax)
    return ax


def create_general_results_plots(psms: List[SpectrumPSMs]):
    fig, axs = fig_setup(ncols=2)
    xcorr_plot(
        psms_by_type={
            NAT_TARGET: [
                psm.native_target for psm in psms if psm.native_target is not None
            ],
            NAT_DECOY: [
                psm.native_decoy for psm in psms if psm.native_decoy is not None
            ],
            HY_TARGET: [
                psm.hybrid_target for psm in psms if psm.hybrid_target is not None
            ],
        },
        ax=axs[0],
    )
    plot_native_vs_hybrid_scores(psms=psms, ax=axs[1])
    finalize(axs)


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


def junction_support_plot(
    ax: Axes, df: pd.DataFrame, x_colm: str, y_colm: str, label_trues: bool = True
):
    if label_trues and "true" in df.columns:
        tmp = df[~df.true]
        _ = sns.scatterplot(
            data=tmp,
            x=x_colm,
            y=y_colm,
            color="blue",
            s=7,
            ax=ax,
        )
        tmp = df[df.true]
        _ = sns.scatterplot(
            x=tmp[x_colm].to_list(),
            y=tmp[y_colm].to_list(),
            color="red",
            marker="X",
            s=14,
            label="'Trues'",
            ax=ax,
        )
    else:
        sns.scatterplot(
            x=df[x_colm],
            y=df[y_colm],
            color="blue",
            s=7,
            ax=ax,
        )
    if y_colm == "min_q_value":
        ax.set_yscale("log")  # set y-axis to log10 scale
        set_title_axes_labels(ax=ax, xlabel=x_colm, ylabel=f"log10({y_colm})")
    else:
        set_title_axes_labels(
            ax=ax,
            xlabel=x_colm,
            ylabel=y_colm,
        )


def mean_precursor_abundance(
    psms: List[CometPSM], spectrum_uid_to_spectrum: Dict[str, Spectrum]
):
    abundances = []
    for psm in psms:
        spectrum = spectrum_uid_to_spectrum[psm.spectrum_uid]
        abundances.append(spectrum.precursor_intensity)
    return np.mean(abundances)


def max_precursor_abundance(
    psms: List[CometPSM], spectrum_uid_to_spectrum: Dict[str, Spectrum]
):
    abundances = []
    for psm in psms:
        spectrum = spectrum_uid_to_spectrum[psm.spectrum_uid]
        abundances.append(spectrum.precursor_intensity)
    return max(abundances)


def get_junction_df(
    hybrid_psms: List[CometPSM],
    protein_name_to_seq_map: Dict[str, str],
    jct_len: int,
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    spectrum_uid_to_spectrum_map: Dict[str, Spectrum],
):
    jct_to_psms = group_hybrid_psms_by_junction(
        hybrid_psms=hybrid_psms,
        protein_name_to_seq_map=protein_name_to_seq_map,
        jct_len=jct_len,
        seq_to_hybrids_map=seq_to_hybrids_map,
    )
    rows = []
    for jct, psms in jct_to_psms.items():
        rows.append(
            (
                jct,
                len(psms),
                list(set(psm.spectrum_uid for psm in psms)),
                mean_precursor_abundance(
                    psms=psms, spectrum_uid_to_spectrum=spectrum_uid_to_spectrum_map
                ),
                max_precursor_abundance(
                    psms=psms, spectrum_uid_to_spectrum=spectrum_uid_to_spectrum_map
                ),
                np.mean([psm.xcorr for psm in psms]),
                max([psm.xcorr for psm in psms]),
                min([psm.q_value for psm in psms]),
                dict(Counter(psm.seq for psm in psms)),
            )
        )
    df = pd.DataFrame(
        rows,
        columns=[
            "jct",
            "num_psm_supporting",
            "spectra_supporting",
            "mean_precursor_abundance",
            "max_precursor_abundance",
            "mean_xcorr",
            "max_xcorr",
            "min_q_value",
            "psm_seq_cnter",
        ],
    )
    df["num_uniq_seqs"] = df.psm_seq_cnter.apply(lambda cnter: len(cnter.keys()))
    df.sort_values(
        by="num_psm_supporting", ascending=False, ignore_index=True, inplace=True
    )
    return df


def junction_plots(
    df: pd.DataFrame,
):
    fig, axs = fig_setup(nrows=2, ncols=3)
    x_colm = "num_psm_supporting"
    for idx, y_colm in enumerate(
        [
            "max_precursor_abundance",
            "mean_precursor_abundance",
            "mean_xcorr",
            "max_xcorr",
            "min_q_value",
            "num_uniq_seqs",
        ]
    ):
        junction_support_plot(
            ax=axs[idx],
            df=df,
            x_colm=x_colm,
            y_colm=y_colm,
        )
    # if true_jcts is not None:
    #     found_trues_str = "\n".join(set(df[df.true].jct))
    #     _ = fig.text(
    #         0.1,
    #         -0.05,
    #         f"Found trues:\n{found_trues_str}",
    #         va="center",
    #         ha="left",
    #         bbox=dict(boxstyle="round", facecolor="white"),
    #     )
    finalize(axs=axs)
    return fig, axs


def create_junction_support_df(
    hybrid_psms: List[CometPSM],
    protein_name_to_seq_map: Dict[str, str],
    min_aa_jct_len: int,
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    spectrum_uid_to_spectrum_map: Dict[str, Spectrum],
    true_hybrids: Optional[List[TrueHybrid]] = None,
) -> pd.DataFrame:
    df = get_junction_df(
        hybrid_psms=hybrid_psms,
        protein_name_to_seq_map=protein_name_to_seq_map,
        jct_len=min_aa_jct_len,
        seq_to_hybrids_map=seq_to_hybrids_map,
        spectrum_uid_to_spectrum_map=spectrum_uid_to_spectrum_map,
    )
    if true_hybrids is not None:
        true_jcts = list(
            set(
                flatten_list_of_lists(
                    true.get_junctions(
                        protein_name_to_seq_map=protein_name_to_seq_map,
                        jct_len=min_aa_jct_len,
                    )
                    for true in true_hybrids
                )
            )
        )
        df["true"] = df.jct.apply(lambda x: x in true_jcts)
    return df


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


@dataclass
class HybridPSM:
    psm: CometPSM
    hybrids: List[HybridPeptide]

    @classmethod
    def create(cls, psm: CometPSM, seq_to_hybrids_map: Dict[str, List[HybridPeptide]]):
        hybrids = seq_to_hybrids_map[psm.seq]
        return cls(psm=psm, hybrids=hybrids)

    def get_junctions(self, jct_len: int = DEFAULT_JCT_LEN):
        return [hy.get_junction_str(jct_len=jct_len) for hy in self.hybrids]


def get_hybrids_for_hybrid_psms(
    seq_to_hybrids_map: Dict[str, List[HybridPeptide]],
    hybrid_psms: List[CometPSM],
) -> List[HybridPeptide]:
    hybrids = []
    for psm in hybrid_psms:
        hybrids.extend(seq_to_hybrids_map[psm.seq])
    return hybrids


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
    mzml = Mzml(path=local_mzml_path)
    file_name = f"mzml{mzml.name}_seq{seq}_idx{spectrum_idx}_scan{spectrum_scan}.mgf"
    return f'wine msconvert {local_mzml_path} --filter "index {spectrum_idx}" --outfile {file_name} -o {out_dir} --mgf'


@dataclass
class AcceptanceMethod:

    @staticmethod
    def get_acceptance_identifier(
        acceptance_method_name: str,
        min_hybrid_side_len: int,
        jct_len: int,
    ) -> str:
        return f"acceptanceMethod={acceptance_method_name}_minHySideLen={min_hybrid_side_len}_minJctLen{jct_len}"

    # def accept_all_hybrids(hs_out: ResultsAnalysis) -> List[CometPSM]:
    #     return list(hs_out.top_hybrid_targets.values())

    @staticmethod
    def accept_hybrids_that_beat_native(
        results: ResultsAnalysis,
    ) -> List[CometPSM]:
        accepted_hybrid_psms = []
        for spectrum_uid, hybrid_psm in results.top_hybrid_targets.items():
            native_psm = results.native_assign_conf.get(spectrum_uid, None)
            if native_psm is None:
                accepted_hybrid_psms.append(hybrid_psm)
                continue
            if hybrid_psm.xcorr > native_psm.xcorr:
                accepted_hybrid_psms.append(hybrid_psm)
        return accepted_hybrid_psms

    @staticmethod
    def acccept_good_hybrids_with_no_good_native_explanation(
        psms: List[SpectrumPSMs], q_threshold: float = DEFAULT_Q_THRESHOLD
    ) -> List[CometPSM]:
        accepted_hybrid_psms = []
        for psm in psms:
            if psm.hybrid_target is not None:
                if psm.hybrid_target.q_value > q_threshold:
                    # Skip hybrids with bad/high q-value
                    continue
                # Continue with hybrids with good/low q-value
                if psm.native_target is None:
                    accepted_hybrid_psms.append(psm.hybrid_target)
                elif psm.native_target.q_value > q_threshold:
                    accepted_hybrid_psms.append(psm.hybrid_target)
        return accepted_hybrid_psms

    @staticmethod
    def accept_hybrids_via_neofusion(
        psms: List[SpectrumPSMs],
        q_vals: List[float] = DEFAULT_Q_RANGE,
        score_deltas: List[float] = DEFAULT_SCORE_CHANGE_RANGE,
        fpr_threshold: float = DEFAULT_FPR,
    ) -> List[CometPSM]:
        neofusion_runner = NeoFusionRunner(
            native_assign_conf={
                psm.spectrum.uid: psm.native_target
                for psm in psms
                if psm.native_target is not None
            },
            top_hybrid_targets={
                psm.spectrum.uid: psm.hybrid_target
                for psm in psms
                if psm.hybrid_target is not None
            },
            q_vals=q_vals,
            score_deltas=score_deltas,
            fpr_threshold=fpr_threshold,
        )
        best_iteration, accepted_hybrids = (
            neofusion_runner.select_hybrid_psms_from_best_iteration(
                neofusion_results=neofusion_runner.run_neofusion()
            )
        )
        logger.info(f"Best NeoFusion iteration info:\n{best_iteration.info}")
        return accepted_hybrids


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
    mzml = Mzml(path=local_mzml_path)
    write_new_line_separated_file(
        lines=script_lines,
        path=Path(local_script_out_dir)
        / f"{mzml.name}_get_true_hybrid_supporting_spectra.sh",
    )


def plot_xcorr_of_accepted_vs_not_accepted_hybrids(
    ax: Axes, all_psms: List[SpectrumPSMs], accepted_hybrid_psms: List[CometPSM]
) -> Tuple[Figure, List[Axes]]:
    accepted_uids = set(psm.uid for psm in accepted_hybrid_psms)
    accepted_hybrid_psms = [
        psm for psm in all_psms if psm.spectrum.uid in accepted_uids
    ]

    native_and_hybrid_psms = [
        psm
        for psm in all_psms
        if (psm.native_target is not None) and (psm.hybrid_target is not None)
    ]
    score = "xcorr"
    s = 7
    sns.scatterplot(
        x=[getattr(psm.native_target, score) for psm in native_and_hybrid_psms],
        y=[getattr(psm.hybrid_target, score) for psm in native_and_hybrid_psms],
        s=s,
        marker="o",
        color="blue",
        label=f"All spectra with a native and hybrid PSM (n={len(native_and_hybrid_psms)})",
        ax=ax,
    )
    sns.scatterplot(
        x=[getattr(psm.native_target, score) for psm in accepted_hybrid_psms],
        y=[getattr(psm.hybrid_target, score) for psm in accepted_hybrid_psms],
        s=s,
        marker="X",
        color="red",
        label=f"Accepted hybrid spectra (n={len(accepted_hybrid_psms)})",
        ax=ax,
    )
    plot_line(ax=ax, label="y=x")
    set_title_axes_labels(
        ax=ax,
        xlabel="Native target xcorr",
        ylabel="Hybrid target xcorr",
    )


@click.command(
    name="process-config",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--config",
    "-c",
    type=PathType(),
    required=True,
    help="Path to config",
)
@click.option(
    "--acceptance_method",
    "-am",
    type=str,
    required=True,
    help=f"Name of the hybrid acceptance method to use. Accepted values: {GOOD_HYBRID_BAD_NATIVE} and {NEOFUSION}",
)
@click.option(
    "--min_hy_side_len",
    "-mhsl",
    type=int,
    required=False,
    show_default=True,
    default=DEFAULT_MIN_SIDE_LEN,
    help="Minimum hybrid side length",
)
@click.option(
    "--q_threshold",
    "-qt",
    type=float,
    required=False,
    show_default=True,
    default=DEFAULT_Q_THRESHOLD,
    help="",
)
@click.option(
    "--jct_len",
    "-jl",
    type=int,
    required=False,
    show_default=True,
    default=DEFAULT_JCT_LEN,
    help="",
)
def cli_process_config(
    config: Path,
    acceptance_method: str,
    q_threshold: float,
    min_hy_side_len: int,
    jct_len: int,
):
    process_hs_config(
        hs_config=config,
        acceptance_method_name=acceptance_method,
        q_threshold=q_threshold,
        min_hybrid_side_len=min_hy_side_len,
        jct_len=jct_len,
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
    cli.add_command(cli_process_config)
    cli()
