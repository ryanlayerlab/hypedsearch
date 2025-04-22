import logging
import math
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import List, Literal, Optional, Union

import pandas as pd
import seaborn as sns
from matplotlib.pyplot import Axes

from src.comet_utils import CometPSM
from src.constants import (
    DEFAULT_PPM_TOLERANCE,
    ION_INT_TO_TYPE,
    IONS_MATCHED,
    MAYBE_HYBRID,
    NOT_HYBRID,
    XCORR,
    IonTypes,
)
from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import Peak, Spectrum, plot_peaks
from src.peptides_and_ions import Peptide, ProductIon, compute_peptide_mz
from src.plot_utils import fig_setup, finalize
from src.protein_product_ion_database import (
    DbKmer,
    PeakWithMatchingProductIons,
    ProteinProductIonDb,
)
from src.utils import Position, flatten_list_of_lists, log_time, mass_difference_in_ppm

logger = logging.getLogger(__name__)


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: ProductIon
    peaks: List[Peak]


@dataclass
class PeakProductIonMatch:
    seq: str
    peak: Peak
    ion_type: str
    charge: int


@dataclass
class PeptideSpectrumComparison:
    spectrum: Spectrum
    peptide: Union[Peptide, str]
    peak_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE
    ion_types: List[IonTypes] = field(
        default_factory=lambda: [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]
    )

    def __post_init__(self):
        if isinstance(self.peptide, str):
            self.peptide = Peptide(seq=self.peptide)

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
            peak_ppm_tolerance=self.peak_ppm_tolerance,
        )

    @cached_property
    def product_ion_seqs_with_matching_peaks(self):
        return group_product_ions_and_matching_peaks_by_charge_and_ion_type(
            product_ions_with_matching_peaks=self.product_ions_with_matching_peaks
        )

    @property
    def proportion_of_intensity_supported_by_peptide(self):
        peaks = flatten_list_of_lists(
            self.product_ion_seqs_with_matching_peaks["matching_peaks"].to_list()
        )
        uniq_peak_ids = set(peak.id for peak in peaks)
        uniq_matching_peaks = list(filter(lambda peak: peak.id in uniq_peak_ids, peaks))

        prop_intensity_supported_by_peptide = (
            sum(peak.intensity for peak in uniq_matching_peaks)
            / self.spectrum.total_intensity
        )
        return prop_intensity_supported_by_peptide

    def compute_num_ion_seqs_matched(self, ion_type: Literal["b", "y"]):
        num_ion_seqs_matched = sum(
            self.product_ion_seqs_with_matching_peaks[
                self.product_ion_seqs_with_matching_peaks["ion_type"] == ion_type
            ]["num_matching_peaks"]
            > 0
        )
        return num_ion_seqs_matched

    @property
    def num_b_ion_seqs_matched(self):
        return self.compute_num_ion_seqs_matched(ion_type="b")

    @property
    def num_y_ion_seqs_matched(self):
        return self.compute_num_ion_seqs_matched(ion_type="y")

    @property
    def num_ion_seqs_matched(self):
        return self.num_b_ion_seqs_matched + self.num_y_ion_seqs_matched

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
        self.comp = PeptideSpectrumComparison(
            spectrum=self.spectrum,
            peptide=self.psm.seq,
            peak_ppm_tolerance=self.peak_ppm_tol,
        )

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


def get_peaks_near_mz(
    query_mz: float, peaks: List[Peak], ppm_tolerance: float
) -> List[Peak]:
    """
    Given a list of mass spectrum peaks and a query mass-to-charge ratio (m/z),
    find the peaks that are within the given PPM tolerance of the query m/z.
    """
    matching_peaks = []
    for peak in peaks:
        if (
            mass_difference_in_ppm(ref_mass=peak.mz, query_mass=query_mz)
            <= ppm_tolerance
        ):
            matching_peaks.append(peak)
    return matching_peaks


def sort_peak_product_ions(
    ions: List[PeakProductIonMatch],
) -> List[PeakProductIonMatch]:
    # Define desired ion_type order
    ion_type_order = {"b": 0, "y": 1}  # Add more types if needed

    return sorted(
        ions,
        key=lambda ion: (
            ion_type_order.get(ion.ion_type, 99),  # Default to 99 if not in dict
            len(ion.seq),
        ),
    )


def get_peaks_that_match_peptide_product_ions(
    spectrum: Spectrum,
    peptide: Peptide,
    ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
    peak_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
) -> pd.DataFrame:
    """
    Compare the given spectrum to the given peptide. This method helps evaluate
    how strong the evidence is for a peptide-spectrum match (PSM).
    """
    # Product ions will have charge <= precursor's charge
    charges = list(range(1, spectrum.precursor_charge + 1))

    # Get product ions of the proposed peptide
    product_ions = peptide.product_ions(ion_types=ion_types, charges=charges)
    assert 2 * len(peptide.seq) * len(charges) == len(product_ions)

    # Get peaks that match a product ion
    product_ions_with_matching_peaks = []
    for ion in product_ions:
        matching_peaks = get_peaks_near_mz(
            query_mz=ion.neutral_mass,
            peaks=spectrum.peaks,
            ppm_tolerance=peak_ppm_tolerance,
        )
        product_ions_with_matching_peaks.append(
            [ion.charge, ion.ion_type_as_str, ion.neutral_mass, ion.seq, matching_peaks]
        )
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


def get_ions_matching_peak(
    peak: Peak,
    precursor_charge: int,
    precursor_mz: float,
    ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
    db_path: Optional[str] = None,
    db: Optional[ProteinProductIonDb] = None,
) -> PeakWithMatchingProductIons:
    # Load database if it's not provided
    if db is None:
        db = ProteinProductIonDb(db_path=db_path)

    # Get product ions from database that are within the given PPM of the peak
    matching_ions = db.get_ions_within_mass_tolerance(
        query_mass=peak.mz, ppm_tolerance=ppm_tolerance
    )

    # Filter out ions with (1) charge > precursor charge or
    # (2) the charge=precursor charge m/z of the ion's AA seq is > precursor m/z
    charge_filtered_ions = []
    for ion in matching_ions:
        if ion.charge <= precursor_charge:
            aa_seq = ion.set_aa_seq(db=db)
            mz = compute_peptide_mz(aa_seq=aa_seq, charge=precursor_charge)
            if mz <= precursor_mz:
                charge_filtered_ions.append(ion)

    return PeakWithMatchingProductIons(peak=peak, ions=charge_filtered_ions)


def peak_to_product_ion_mapping(
    spectrum: Spectrum,
    db_path: str,
    num_cpus: Optional[int] = None,
    ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
) -> List[PeakWithMatchingProductIons]:

    process_peak_fcn = lambda peak: get_ions_matching_peak(
        db_path=db_path,
        peak=peak,
        precursor_charge=spectrum.precursor_charge,
        ppm_tolerance=ppm_tolerance,
        precursor_mz=spectrum.precursor_mz,
    )

    if num_cpus is None:
        num_cpus = multiprocessing.cpu_count()

    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        peaks_with_matches = list(executor.map(process_peak_fcn, spectrum.peaks))

    return peaks_with_matches


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


def get_start_and_end_positions_from_ions(ions: List[DbKmer]) -> Position:
    start = min([ion.inclusive_start for ion in ions])
    end = max([ion.exclusive_end for ion in ions])
    return Position(inclusive_start=start, exclusive_end=end)


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


@dataclass
class HybridSpectrumComparison:
    num_hybrid_supporting_ions: int
    num_hybrid_supporting_ions_with_peak_match: int
    df: pd.DataFrame


def compare_hybrid_to_spectrum(
    hybrid: HybridPeptide, spectrum: Spectrum, peak_mz_ppm_tol: int
):

    product_ions_with_matching_peaks = get_peaks_that_match_peptide_product_ions(
        spectrum=spectrum,
        peptide=Peptide(seq=hybrid.seq),
        peak_ppm_tolerance=peak_mz_ppm_tol,
    )
    df = group_product_ions_and_matching_peaks_by_charge_and_ion_type(
        product_ions_with_matching_peaks=product_ions_with_matching_peaks
    )
    df["hybrid"] = 1

    df.loc[
        (df["ion_type"] == "b") & (df["seq"].apply(lambda seq: seq in hybrid.b_seq)),
        "hybrid",
    ] = 0
    df.loc[
        (df["ion_type"] == "y") & (df["seq"].apply(lambda seq: seq in hybrid.y_seq)),
        "hybrid",
    ] = 0

    # num_non_hybrid_b_ions = df[(df["hybrid"] == 0) & (df["ion_type"] == "b")].shape[0]
    # num_non_hybrid_y_ions = df[(df["hybrid"] == 0) & (df["ion_type"] == "y")].shape[0]
    num_hybrid_supporting_ions = sum(df["hybrid"])
    num_hybrid_supporting_ions_with_peak_match = df[
        (df["hybrid"] == 1) & (df["num_matching_peaks"] != 0)
    ].shape[0]

    return HybridSpectrumComparison(
        num_hybrid_supporting_ions=num_hybrid_supporting_ions,
        num_hybrid_supporting_ions_with_peak_match=num_hybrid_supporting_ions_with_peak_match,
        df=df,
    )


def get_scan_results(
    native_txt: Path, hybrid_txt: Path, sample: str, scan: Optional[int] = None
):
    native_run_psms = CometPSM.from_txt(file_path=native_txt, sample=sample)
    hybrid_run_psms = CometPSM.from_txt(file_path=hybrid_txt, sample=sample)

    if scan is not None:
        native_run_psms = list(filter(lambda psm: psm.scan == scan, native_run_psms))
        hybrid_run_psms = list(filter(lambda psm: psm.scan == scan, hybrid_run_psms))
    return native_run_psms, hybrid_run_psms


@dataclass
class HSDecision:
    decision: int
    reason: str


@dataclass
class CompareNativeRunToHybridRun:
    native_psms: List[CometPSM]
    hybrid_psms: List[CometPSM]
    true_hybrid: bool = False
    top_ranking_native_psms: List[CometPSM] = field(init=False)
    top_ranking_hybrid_psms: List[CometPSM] = field(init=False)
    xcorr_rel_change: float = field(init=False)
    ions_matched_rel_change: float = field(init=False)
    ions_matched_change: float = field(init=False)
    decision: int = field(init=False)

    def __post_init__(self):
        self.set_top_psms()
        self.xcorr_rel_change = self.get_rel_change(column=XCORR)
        self.ions_matched_rel_change = self.get_rel_change(column=IONS_MATCHED)
        self.ions_matched_change = self.get_change(column=IONS_MATCHED)

        # Hypedsearch decision
        if self.xcorr_rel_change >= 0:
            # the run 1 xcorr is greater than run 2 -> not hybrid
            self.decision = HSDecision(decision=NOT_HYBRID, reason="xcorr")

        elif self.ions_matched_rel_change >= 0:
            self.decision = HSDecision(decision=NOT_HYBRID, reason="ions_matched")

        else:
            self.decision = MAYBE_HYBRID

    def set_top_psms(self):
        top_nums = [1]
        self.top_ranking_native_psms = list(
            filter(lambda psm: psm.num in top_nums, self.native_psms)
        )
        self.top_ranking_hybrid_psms = list(
            filter(lambda psm: psm.num in top_nums, self.hybrid_psms)
        )

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

    @property
    def scan(self):
        scan = set(psm.scan for psm in self.native_psms).union(
            psm.scan for psm in self.hybrid_psms
        )
        assert len(scan) == 1
        return list(scan)[0]

    @property
    def sample(self):
        sample = set(psm.sample for psm in self.native_psms).union(
            psm.sample for psm in self.hybrid_psms
        )
        assert len(sample) == 1
        return list(sample)[0]


def get_hybrid_peptide_from_comet_psm(psm: CometPSM):
    if len(psm.proteins) > 1:
        print(f"There is more than one protein for this hybrid!")
        protein = psm.proteins[0]

        # assert (
        #     0 == 1
        # ), f"This PSM doesn't seem to correspond to a hybrid. psm.proteins = {psm.proteins}"

    protein = psm.proteins[0]
    assert (
        protein[:7] == "hybrid_"
    ), f"This PSM doesn't seem to correspond to a hybrid. psm.proteins = {psm.proteins}"

    protein = protein[7:]
    b_seq, y_seq = protein.split("-")
    return HybridPeptide(b_seq=b_seq, y_seq=y_seq)


@log_time
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
    spectra = Spectrum.from_mzml(mzml)

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


def process_one_mzml(mzml_path, psms, peak_ppm_tolerance=10, verbose=False):
    return compare_mzml_spectra_to_comet_psms(
        psms=psms,
        mzml=mzml_path,
        peak_ppm_tolerance=peak_ppm_tolerance,
        verbose=verbose,
    )
