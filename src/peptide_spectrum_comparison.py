import logging
import math
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import List, Literal, Optional, Set, Union

import pandas as pd
import seaborn as sns
from matplotlib.pyplot import Axes

from src.comet_utils import CometPSM
from src.constants import (
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PPM_TOLERANCE,
    ION_INT_TO_TYPE,
    IONS_MATCHED,
    XCORR,
    IonTypes,
)

# from src.hypedsearch_utils import HybridPeptide
from src.kmer_database import DbKmer
from src.mass_spectra import Peak, Spectrum, plot_peaks
from src.peptides_and_ions import (
    Peptide,
    UnpositionedProductIon,
    compute_peptide_precursor_mz,
)
from src.plot_utils import fig_setup, finalize
from src.utils import flatten_list_of_lists, log_time, mass_difference_in_ppm

logger = logging.getLogger(__name__)


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: UnpositionedProductIon
    peaks: List[Peak]


@dataclass
class PeptideSpectrumComparison:
    spectrum: Spectrum
    peptide: Union[str, Peptide]
    hybrid_seq: Optional[str] = None
    peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL
    ion_types: List[IonTypes] = field(
        default_factory=lambda: [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
    )

    def __post_init__(self):
        if isinstance(self.peptide, str):
            self.peptide = Peptide(seq=self.peptide)

    @classmethod
    def from_str(
        cls,
        spectrum: Spectrum,
        peptide: str,
        peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    ):
        hybrid_seq = None
        if "-" in peptide:
            hybrid_seq = peptide
            peptide = peptide.replace("-", "")
        peptide = Peptide(seq=peptide)
        return cls(
            spectrum=spectrum,
            peptide=peptide,
            hybrid_seq=hybrid_seq,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
        )

    def get_ion_seqs_supported(self, ion_type: Literal["b", "y"]):
        tmp = self.product_ion_seqs_with_matching_peaks
        ion_seqs_supported = set(
            tmp[(tmp["ion_type"] == ion_type) & (tmp["num_matching_peaks"] > 0)]["seq"]
        )
        return ion_seqs_supported

    def _get_proportion_of_intensity_supported_by_seqs(
        self, ion_type: Literal["b", "y", "b+y", "hybrid"]
    ):
        tmp = self.product_ion_seqs_with_matching_peaks  # so that name is shorter
        if ion_type == "b":
            peaks = flatten_list_of_lists(
                tmp[tmp["ion_type"] == "b"]["matching_peaks"].to_list()
            )
        elif ion_type == "y":
            peaks = flatten_list_of_lists(
                tmp[tmp["ion_type"] == "y"]["matching_peaks"].to_list()
            )
        elif ion_type == "b+y":
            peaks = flatten_list_of_lists(tmp["matching_peaks"].to_list())
        elif ion_type == "hybrid":
            peaks = flatten_list_of_lists(
                tmp[tmp["hybrid_seq"] == 1]["matching_peaks"].to_list()
            )
        else:
            raise ValueError(
                f"Unknown ion type: {ion_type}. Should be 'b', 'y', 'b+y', or 'hybrid'."
            )
        uniq_peak_ids = set(peak.id for peak in peaks)
        uniq_matching_peaks = list(filter(lambda peak: peak.id in uniq_peak_ids, peaks))

        prop_intensity_supported_by_seqs = (
            sum(peak.intensity for peak in uniq_matching_peaks)
            / self.spectrum.total_intensity
        )
        return prop_intensity_supported_by_seqs

    @property
    def seq(self):
        return self.peptide.seq

    @property
    def sample(self):
        return self.spectrum.sample

    @property
    def scan(self):
        return self.spectrum.scan

    @cached_property
    def product_ions_with_matching_peaks(self):
        return get_peaks_that_match_peptide_product_ions(
            spectrum=self.spectrum,
            peptide=self.peptide,
            ion_types=self.ion_types,
            peak_ppm_tolerance=self.peak_to_ion_ppm_tolerance,
        )

    @cached_property
    def product_ion_seqs_with_matching_peaks(self):
        df = group_product_ions_and_matching_peaks_by_charge_and_ion_type(
            product_ions_with_matching_peaks=self.product_ions_with_matching_peaks
        )
        if self.hybrid_seq is not None:
            df["hybrid_seq"] = 0
            b_seq, y_seq = self.hybrid_seq.split("-")
            seqs_containing_hybrid_junction = [
                b_seq + y_seq[:idx] for idx in range(1, len(y_seq) + 1)
            ]
            seqs_containing_hybrid_junction += [
                b_seq[idx:] + y_seq for idx in range(1, len(b_seq))
            ]
            df.loc[(df.seq.isin(seqs_containing_hybrid_junction)), "hybrid_seq"] = 1

        return df

    @property
    def hybrid_seqs_supported(self) -> Set[str]:
        if self.hybrid_seq is None:
            return None
        tmp = self.product_ion_seqs_with_matching_peaks
        hybrid_seqs_supported = set(
            tmp[(tmp["hybrid_seq"] == 1) & (tmp["num_matching_peaks"] > 0)]["seq"]
        )
        return hybrid_seqs_supported

    @property
    def num_hybrid_seqs(self):
        if self.hybrid_seq is None:
            return None
        return sum(self.product_ion_seqs_with_matching_peaks["hybrid_seq"] == 1)

    @property
    def num_hybrid_seqs_supported(self):
        return len(self.hybrid_seqs_supported)

    @property
    def prop_hybrid_seqs_supported(self):
        return self.num_hybrid_seqs_supported / self.num_hybrid_seqs

    @property
    def prop_of_intensity_supported_by_seqs(self):
        return self._get_proportion_of_intensity_supported_by_seqs(ion_type="b+y")

    @property
    def prop_of_intensity_supported_by_b_seqs(self):
        return self._get_proportion_of_intensity_supported_by_seqs(ion_type="b")

    @property
    def prop_of_intensity_supported_by_y_seqs(self):
        return self._get_proportion_of_intensity_supported_by_seqs(ion_type="y")

    @property
    def prop_of_intensity_supported_by_hybrid_seqs(self):
        return self._get_proportion_of_intensity_supported_by_seqs(ion_type="hybrid")

    @property
    def num_b_ion_seqs(self):
        return len(self.peptide.seq)

    @property
    def num_y_ion_seqs(self):
        return len(self.peptide.seq)

    @property
    def b_ion_seqs_supported(self):
        return self.get_ion_seqs_supported(ion_type="b")

    @property
    def num_b_ion_seqs_supported(self):
        return len(self.get_ion_seqs_supported(ion_type="b"))

    @property
    def y_ion_seqs_supported(self):
        return self.get_ion_seqs_supported(ion_type="y")

    @property
    def num_y_ion_seqs_supported(self):
        return len(self.get_ion_seqs_supported(ion_type="y"))

    @property
    def num_ion_seqs_supported(self):
        return self.num_b_ion_seqs_supported + self.num_y_ion_seqs_supported

    @property
    def proportion_of_seqs_supported(self):
        return self.num_ion_seqs_supported / (self.num_b_ion_seqs + self.num_y_ion_seqs)

    @property
    def num_peaks_with_match(self):
        return len(
            set(
                [
                    peak.id
                    for peak in flatten_list_of_lists(
                        self.product_ions_with_matching_peaks["matching_peaks"]
                    )
                ]
            )
        )

    @property
    def ppm_diff(self):
        peptide_mz = compute_peptide_precursor_mz(
            aa_seq=self.peptide.seq, charge=self.spectrum.precursor_charge
        )
        return mass_difference_in_ppm(
            mass1=peptide_mz, mass2=self.spectrum.precursor_mz
        )

    def plot_ions(
        self,
        ax: Axes,
        intensity: Optional[float] = None,
        # log_intensity: bool = False,
        # txt_size: float = 5,
    ):
        ion_types = ["b", "y"]
        # Get peaks matching product ions
        peak_product_ion_matches = flatten_list_of_lists(
            [
                list(
                    filter(
                        lambda x: x.ion_type == ion_type, self.peak_product_ion_matches
                    )
                )
                for ion_type in ion_types
            ]
        )

        colors = sns.color_palette("husl", len(peak_product_ion_matches))
        idx = 0
        for match in sort_peak_product_ions(peak_product_ion_matches):
            if intensity is not None:
                match.peak.intensity = intensity
            label = f"{match.ion_type} (z={match.charge}): {match.seq}"
            plot_peaks(
                ax=ax,
                peaks=[match.peak],
                label=label,
                color=colors[idx],
                # log_intensity=log_intensity,
            )
            idx += 1

        #     # Label them
        #     max_intensity = max([peak.intensity for peak in self.spectrum.peaks])
        #     for peak_product_ion_match in peak_product_ion_matches:
        #         peak = peak_product_ion_match.peak
        #         _ = ax.text(
        #             x=peak.mz,
        #             # peak.intensity + 0.02 * max_intensity,  # y value
        #             y=max_intensity + 0.01 * max_intensity,
        #             s=peak_product_ion_match.seq,  # label
        #             rotation=90,
        #             ha="center",
        #             va="bottom",
        #             fontsize=txt_size,
        #             color="black",
        #         )

    def plot(
        self,
        ax: Optional[Axes] = None,
        intensity: Optional[float] = None,
        title: Optional[str] = None,
        # title_pad: int = 20,
        log_intensity: bool = False,
        # ion_txt_size: float = 5,
        # title_txt_size: float = 12,
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        self.spectrum.plot_spectrum(ax=ax, log_intensity=log_intensity)
        self.plot_ions(
            ax=ax,
            intensity=intensity,
            # txt_size=ion_txt_size,
        )
        if title is None:
            title = self.peptide.seq
        else:
            title += f"\n({self.peptide.seq})"
        _ = ax.set_title(title)
        # pad=title_pad, loc="left")
        # _ = ax.text(
        #     -0.4,
        #     1.2,
        #     title,  # adjust -0.1 to move farther left
        #     transform=ax.transAxes,
        #     ha="center",
        #     va="bottom",
        #     fontsize=title_txt_size,
        # )
        finalize(ax)
        _ = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

        return ax


@dataclass
class PsmSpectrumComparison:
    spectrum: Spectrum
    psm: CometPSM
    comp: PeptideSpectrumComparison = field(init=False)
    peak_ppm_tol: float = DEFAULT_PPM_TOLERANCE

    def __post_init__(self):
        self.comp = PeptideSpectrumComparison.from_str(
            spectrum=self.spectrum,
            peptide=self.psm.seq,
        )

    @property
    def comet_ions_matched(self):
        return self.psm.ions_matched

    @property
    def prop_intensity_supported(self):
        return self.comp.proportion_of_intensity_supported_by_peptide

    @property
    def prop_seqs_supported(self):
        return self.comp.proportion_of_seqs_supported

    @property
    def seq(self):
        return self.psm.seq

    @property
    def scan(self):
        assert (
            self.spectrum.scan == self.psm.scan
        ), f"{self.spectrum.scan} != {self.psm.scan}"
        return self.spectrum.scan

    @property
    def sample(self):
        assert (
            self.spectrum.sample == self.psm.sample
        ), f"{self.spectrum.sample} != {self.psm.sample}"
        return self.spectrum.sample

    @property
    def xcorr(self):
        return self.psm.xcorr

    @property
    def eval(self):
        return self.psm.eval

    @property
    def proteins(self):
        return self.psm.proteins


@dataclass
class HSDecision:
    decision: int
    reason: str


@dataclass
class ScanPsms:
    native: List[CometPSM]
    hybrid: List[CometPSM]

    @property
    def scan(self):
        scan = set(psm.scan for psm in self.native).union(
            psm.scan for psm in self.hybrid
        )
        assert len(scan) == 1
        return list(scan)[0]

    @property
    def sample(self):
        sample = set(psm.sample for psm in self.native).union(
            psm.sample for psm in self.hybrid
        )
        assert len(sample) == 1
        return list(sample)[0]


@dataclass
class CompareNativePSMsToHybridPSMs:
    native_psms: List[CometPSM]
    hybrid_psms: List[CometPSM]
    # true_hybrid: bool = False
    # top_ranking_native_psms: List[CometPSM] = field(init=False)
    # top_ranking_hybrid_psms: List[CometPSM] = field(init=False)
    # xcorr_rel_change: float = field(init=False)
    # ions_matched_rel_change: float = field(init=False)
    # ions_matched_change: float = field(init=False)
    # decision: int = field(init=False)
    # top_native_seqs: List[str] = field(init=False)
    # top_hybrid_seqs: List[str] = field(init=False)

    def set_top_psms(self):
        top_nums = [1]
        self.top_ranking_native_psms = list(
            filter(lambda psm: psm.num in top_nums, self.native_psms)
        )
        self.top_ranking_hybrid_psms = list(
            filter(lambda psm: psm.num in top_nums, self.hybrid_psms)
        )
        self.top_native_seqs = [
            psm.seq_with_hyphen for psm in self.top_ranking_native_psms
        ]
        self.top_hybrid_seqs = flatten_list_of_lists(
            [psm.seq_with_hyphen for psm in self.top_ranking_hybrid_psms]
        )

    def set_xcorr_fields(self):
        self.top_native_xcorr = max([getattr(psm, XCORR) for psm in self.native_psms])
        self.top_hybrid_xcorr = max([getattr(psm, XCORR) for psm in self.hybrid_psms])
        self.xcorr_change = self.top_hybrid_xcorr - self.top_native_xcorr
        self.xcorr_rel_change = self.xcorr_change / self.top_native_xcorr

    def set_eval_fields(self):
        self.top_native_eval = max([getattr(psm, "eval") for psm in self.native_psms])

    def set_ions_fields(self):
        self.top_native_ions_matched = max(
            [getattr(psm, IONS_MATCHED) for psm in self.native_psms]
        )
        self.top_hybrid_ions_matched = max(
            [getattr(psm, IONS_MATCHED) for psm in self.hybrid_psms]
        )
        self.ions_matched_change = (
            self.top_hybrid_ions_matched - self.top_native_ions_matched
        )

    def __post_init__(self):
        self.set_top_psms()
        self.set_xcorr_fields()
        self.set_eval_fields()
        self.set_ions_fields()

        # # Hypedsearch decision
        # if self.xcorr_rel_change >= 0:
        #     # the run 1 xcorr is greater than run 2 -> not hybrid
        #     self.decision = HSDecision(decision=NOT_HYBRID, reason="xcorr")

        # elif self.ions_matched_rel_change >= 0:
        #     self.decision = HSDecision(decision=NOT_HYBRID, reason="ions_matched")

        # else:
        #     self.decision = MAYBE_HYBRID

    def get_rel_change(self, column: str):
        tol = 1e-5
        # Relative change in xcorr from run 1 to run 2
        run_1_max = max([getattr(psm, column) for psm in self.top_ranking_native_psms])
        run_2_max = max([getattr(psm, column) for psm in self.top_ranking_hybrid_psms])
        if math.isclose(run_1_max, 0, rel_tol=tol):
            run_1_max = run_1_max + tol
        rel_change = (run_1_max - run_2_max) / run_1_max
        return rel_change

    def get_change(self, column: str):
        run_1_max = max([getattr(psm, column) for psm in self.top_ranking_native_psms])
        run_2_max = max([getattr(psm, column) for psm in self.top_ranking_hybrid_psms])
        change = run_1_max - run_2_max
        return change

    @cached_property
    def scan(self):
        scan = set(psm.scan for psm in self.native_psms).union(
            psm.scan for psm in self.hybrid_psms
        )
        assert len(scan) == 1
        return list(scan)[0]

    @cached_property
    def sample(self):
        sample = set(psm.sample for psm in self.native_psms).union(
            psm.sample for psm in self.hybrid_psms
        )
        assert len(sample) == 1
        return list(sample)[0]


@dataclass
class HybridSpectrumComparison:
    num_hybrid_supporting_ions: int
    num_hybrid_supporting_ions_with_peak_match: int
    df: pd.DataFrame


def get_peaks_near_mz(
    query_mz: float, peaks: List[Peak], ppm_tolerance: float
) -> List[Peak]:
    """
    Given a list of mass spectrum peaks and a query mass-to-charge ratio (m/z),
    find the peaks that are within the given PPM tolerance of the query m/z.
    """
    matching_peaks = []
    for peak in peaks:
        if mass_difference_in_ppm(mass1=peak.mz, mass2=query_mz) <= ppm_tolerance:
            matching_peaks.append(peak)
    return matching_peaks


# def sort_peak_product_ions(
#     ions: List[PeakProductIonMatch],
# ) -> List[PeakProductIonMatch]:
#     # Define desired ion_type order
#     ion_type_order = {"b": 0, "y": 1}  # Add more types if needed

#     return sorted(
#         ions,
#         key=lambda ion: (
#             ion_type_order.get(ion.ion_type, 99),  # Default to 99 if not in dict
#             len(ion.seq),
#         ),
#     )


@dataclass
class PeakIonMatch:
    ion_mz: float
    ion_charge: int
    ion_seq: str
    ion_type: str
    peak_mz: float
    peak_intensity: float
    sample: Optional[str]
    scan: Optional[int]

    def mz_diff(self, type: Literal["rel", "rel_ppm"] = "rel_ppm"):
        """
        Mass (more precisely, m/z) difference between the theoretical ion and the peak.
        Let x_i = theoretical ion mass, x_p = peak mass,
        then returns
            - (x_i - x_t) / x_i when type='rel'
            - ((x_i - x_t) / x_i) * (10**6) when type='rel_ppm'
        """
        if type == "rel":
            return (self.ion_mz - self.peak_mz) / self.ion_mz
        elif type == "rel_ppm":
            return ((self.ion_mz - self.peak_mz) / self.ion_mz) * (10**6)


def get_peaks_that_match_peptide_product_ions(
    spectrum: Spectrum,
    peptide: Union[Peptide, str],
    ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
    peak_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
    return_dataclasses: bool = False,
) -> Union[List[PeakIonMatch], pd.DataFrame]:
    """
    Compare the given spectrum to the given peptide. This method helps evaluate
    how strong the evidence is for a peptide-spectrum match (PSM).
    """
    if isinstance(peptide, str):
        peptide = Peptide(seq=peptide)
    # Product ions will have charge <= precursor's charge
    charges = list(range(1, spectrum.precursor_charge + 1))

    # Get product ions of the proposed peptide
    product_ions = peptide.product_ions(ion_types=ion_types, charges=charges)
    assert 2 * len(peptide.seq) * len(charges) == len(product_ions)

    # Get peaks that match a product ion
    product_ions_with_matching_peaks = []
    peak_ion_matches = []
    for ion in product_ions:
        matching_peaks = get_peaks_near_mz(
            query_mz=ion.mz,
            peaks=spectrum.peaks,
            ppm_tolerance=peak_ppm_tolerance,
        )

        if return_dataclasses:
            peak_ion_matches.extend(
                [
                    PeakIonMatch(
                        ion_mz=ion.mz,
                        ion_charge=ion.charge,
                        ion_type=ion.ion_type_as_str,
                        ion_seq=ion.seq,
                        peak_mz=peak.mz,
                        peak_intensity=peak.intensity,
                        sample=spectrum.sample,
                        scan=spectrum.scan,
                    )
                    for peak in matching_peaks
                ]
            )
        else:
            product_ions_with_matching_peaks.append(
                [
                    ion.charge,
                    ion.ion_type_as_str,
                    ion.mz,
                    ion.seq,
                    matching_peaks,
                ]
            )
    if return_dataclasses:
        return peak_ion_matches
    else:
        product_ions_with_matching_peaks = pd.DataFrame(
            product_ions_with_matching_peaks,
            columns=["charge", "ion_type", "m/z", "seq", "matching_peaks"],
        )
        return product_ions_with_matching_peaks


def group_product_ions_and_matching_peaks_by_charge_and_ion_type(
    product_ions_with_matching_peaks: pd.DataFrame,
):
    product_ion_seqs_with_matching_peaks = {"b": [], "y": []}
    for name, group in product_ions_with_matching_peaks.groupby(by=["seq", "ion_type"]):
        matching_peaks = flatten_list_of_lists(
            [matching_peaks for matching_peaks in group["matching_peaks"]]
        )
        num_matching_peaks = len(matching_peaks)
        product_ion_seqs_with_matching_peaks[name[1]].append(
            [
                name[0],
                name[1],
                num_matching_peaks,
                matching_peaks,
            ]
        )
    for key, tmp_data in product_ion_seqs_with_matching_peaks.items():
        df = pd.DataFrame(
            tmp_data,
            columns=["seq", "ion_type", "num_matching_peaks", "matching_peaks"],
        )
        df.sort_values(
            by=["seq"], key=lambda x: x.str.len(), inplace=True, ignore_index=True
        )
        product_ion_seqs_with_matching_peaks[key] = df

    product_ion_seqs_with_matching_peaks = pd.concat(
        [
            product_ion_seqs_with_matching_peaks["b"],
            product_ion_seqs_with_matching_peaks["y"],
        ],
        ignore_index=True,
    )
    return product_ion_seqs_with_matching_peaks


def ions_as_df(ions: List[DbKmer]):
    data = [
        [
            ion.protein_id,
            ion.inclusive_start,
            ion.exclusive_end,
            ion.charge,
            ion.neutral_mass,
            ION_INT_TO_TYPE[ion.ion_type],
            ion.aa_seq,
        ]
        for ion in ions
    ]
    df = pd.DataFrame(
        data, columns=["p_id", "start", "end", "charge", "m/z", "type", "seq"]
    )
    return df


def create_df_from_ions(ions: List[DbKmer]) -> pd.DataFrame:
    rows = []
    for ion in ions:
        rows.append(
            [
                ion.charge,
                ION_INT_TO_TYPE[ion.ion_type],
                ion.neutral_mass,
                ion.aa_seq,
                ion.protein_id,
                ion.inclusive_start,
                ion.exclusive_end,
            ]
        )
    df = pd.DataFrame(
        rows,
        columns=["charge", "ion_type", "m/z", "seq", "p_id", "i_start", "e_end"],
    )
    return df


# def compare_hybrid_to_spectrum(
#     hybrid: HybridPeptide, spectrum: Spectrum, peak_mz_ppm_tol: int
# ):
#     product_ions_with_matching_peaks = get_peaks_that_match_peptide_product_ions(
#         spectrum=spectrum,
#         peptide=Peptide(seq=hybrid.seq),
#         peak_ppm_tolerance=peak_mz_ppm_tol,
#     )
#     df = group_product_ions_and_matching_peaks_by_charge_and_ion_type(
#         product_ions_with_matching_peaks=product_ions_with_matching_peaks
#     )
#     df["hybrid"] = 1

#     df.loc[
#         (df["ion_type"] == "b") & (df["seq"].apply(lambda seq: seq in hybrid.b_seq)),
#         "hybrid",
#     ] = 0
#     df.loc[
#         (df["ion_type"] == "y") & (df["seq"].apply(lambda seq: seq in hybrid.y_seq)),
#         "hybrid",
#     ] = 0

#     # num_non_hybrid_b_ions = df[(df["hybrid"] == 0) & (df["ion_type"] == "b")].shape[0]
#     # num_non_hybrid_y_ions = df[(df["hybrid"] == 0) & (df["ion_type"] == "y")].shape[0]
#     num_hybrid_supporting_ions = sum(df["hybrid"])
#     num_hybrid_supporting_ions_with_peak_match = df[
#         (df["hybrid"] == 1) & (df["num_matching_peaks"] != 0)
#     ].shape[0]

#     return HybridSpectrumComparison(
#         num_hybrid_supporting_ions=num_hybrid_supporting_ions,
#         num_hybrid_supporting_ions_with_peak_match=num_hybrid_supporting_ions_with_peak_match,
#         df=df,
#     )


def get_scan_results(
    native_txt: Path, hybrid_txt: Path, sample: str, scan: Optional[int] = None
):
    native_run_psms = CometPSM.from_txt(txt=native_txt, sample=sample)
    hybrid_run_psms = CometPSM.from_txt(txt=hybrid_txt, sample=sample)

    if scan is not None:
        native_run_psms = list(filter(lambda psm: psm.scan == scan, native_run_psms))
        hybrid_run_psms = list(filter(lambda psm: psm.scan == scan, hybrid_run_psms))
    return native_run_psms, hybrid_run_psms


# def get_hybrid_peptide_from_comet_psm(psm: CometPSM):
#     if len(psm.proteins) > 1:
#         print("There is more than one protein for this hybrid!")
#         protein = psm.proteins[0]

#         # assert (
#         #     0 == 1
#         # ), f"This PSM doesn't seem to correspond to a hybrid. psm.proteins = {psm.proteins}"

#     protein = psm.proteins[0]
#     assert (
#         protein[:7] == "hybrid_"
#     ), f"This PSM doesn't seem to correspond to a hybrid. psm.proteins = {psm.proteins}"

#     protein = protein[7:]
#     b_seq, y_seq = protein.split("-")
#     return HybridPeptide(b_seq=b_seq, y_seq=y_seq)


@log_time(level=logging.DEBUG)
def compare_mzml_spectra_to_comet_psms(
    psms: List[CometPSM],
    mzml: Path,
    peak_ppm_tolerance: float,
    verbose: bool = False,
) -> List[PsmSpectrumComparison]:
    """
    Compare the spectra in the mzML file to the PSMs in the CometPSM list.
    """
    # Read in the mzML file
    spectra = Spectrum.parse_ms2_from_mzml(mzml)

    # Compare each spectrum to the PSMs
    results = []
    for idx, spectrum in enumerate(spectra):
        if verbose:
            if (idx + 1) % 500 == 0:
                logger.info(f"Comparing spectrum {idx + 1} of {len(spectra)}")
        # Get the PSMs that correspond to this spectrum
        scan_psms = list(
            filter(
                lambda psm: (psm.scan == spectrum.scan)
                and (psm.sample == spectrum.sample),
                psms,
            )
        )
        if len(scan_psms) == 0:
            continue
        psm = scan_psms[0]
        comp = PsmSpectrumComparison(
            spectrum=spectrum, psm=psm, peak_ppm_tol=peak_ppm_tolerance
        )
        results.append(comp)
    return results


# def process_one_mzml(
#     mzml_path, psms, peak_ppm_tolerance=10, verbose=False
# ) -> List[PsmSpectrumComparison]:
#     return compare_mzml_spectra_to_comet_psms(
#         psms=psms,
#         mzml=mzml_path,
#         peak_ppm_tolerance=peak_ppm_tolerance,
#         verbose=verbose,
#     )
