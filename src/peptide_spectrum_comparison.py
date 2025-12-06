import logging
import math
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import List, Literal, Optional, Set, Union

import pandas as pd
import seaborn as sns
from matplotlib.pyplot import Axes
from pydantic import BaseModel, model_validator

from src.comet_utils import CometPSM
from src.constants import (
    B_ION_TYPE,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PPM_TOLERANCE,
    ION_INT_TO_TYPE,
    IONS_MATCHED,
    XCORR,
    Y_ION_TYPE,
    IonTypes,
)

# from src.hypedsearch_utils import HybridPeptide
from src.kmer_database import DbKmer
from src.mass_spectra import (
    Peak,
    Spectrum,
    create_sample_scan_to_spectrum_map,
    plot_peaks,
)
from src.peptides_and_ions import (
    Peptide,
    UnpositionedProductIon,
    compute_peptide_precursor_mz,
)
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
from src.utils import (
    flatten_list_of_lists,
    get_b_ion_prefixes,
    get_y_ion_suffixes,
    list_to_df,
    log_time,
    mass_difference_in_ppm,
)

logger = logging.getLogger(__name__)


class PeakIonMatch(BaseModel):
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

    @property
    def ion_id(self):
        return f"{self.ion_type}-{self.ion_seq}-z{self.ion_charge}"


class PSM(BaseModel):
    peptide: str
    peak_ion_matches: List[PeakIonMatch]
    peak_to_ion_ppm_tolerance: float
    prop_intensity_supported: float
    prop_ions_matched: float
    spectrum: Optional[Spectrum] = None
    spectrum_id: Optional[str] = None
    xcorr: Optional[float] = None
    q_value: Optional[float] = None

    @classmethod
    def from_spectrum_and_seq(
        cls,
        spectrum: Spectrum,
        seq: str,
        peak_to_ion_ppm_tolerance: float,
    ):
        peak_to_ion_matches = get_peak_product_ion_matches(
            spectrum=spectrum,
            peptide=seq,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
        )
        intensity_supported = sum(
            [peak_ion_match.peak_intensity for peak_ion_match in peak_to_ion_matches]
        )
        # Compute the proportion of ions matched
        num_product_ions = len(
            Peptide(seq=seq).product_ions(
                charges=list(range(1, spectrum.precursor_charge + 1))
            )
        )
        num_found_ions = len(set(x.ion_id for x in peak_to_ion_matches))
        return cls(
            peptide=seq,
            peak_ion_matches=peak_to_ion_matches,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
            spectrum_id=spectrum.uid,
            prop_intensity_supported=intensity_supported / spectrum.total_intensity,
            prop_ions_matched=num_found_ions / num_product_ions,
        )

    @classmethod
    def from_spectrum_and_comet_psm(
        cls,
        comet_psm: Optional[CometPSM],
        spectrum: Spectrum,
        peak_to_ion_ppm_tolerance: float,
    ) -> Optional["PSM"]:
        if comet_psm is None:
            return None
        peptide = comet_psm.seq
        peak_to_ion_matches = get_peak_product_ion_matches(
            spectrum=spectrum,
            peptide=peptide,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
        )
        intensity_supported = sum(
            [peak_ion_match.peak_intensity for peak_ion_match in peak_to_ion_matches]
        )
        return cls(
            peptide=peptide,
            peak_ion_matches=peak_to_ion_matches,
            peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance,
            spectrum_id=spectrum.uid,
            xcorr=comet_psm.xcorr,
            q_value=comet_psm.q_value,
            prop_intensity_supported=intensity_supported / spectrum.total_intensity,
            prop_ions_matched=comet_psm.prop_ions_matched,
        )

    @classmethod
    def from_comet_psms(
        cls,
        comet_psms: List[CometPSM],
        spectra_dir: Union[str, Path],
        peak_to_ion_ppm_threshold: float,
    ) -> List["PSM"]:
        psms = []
        mzml_names = set(psm.sample for psm in comet_psms)
        mzmls = [Path(spectra_dir) / f"{name}.mzML" for name in mzml_names]
        sample_scan_to_spectrum_map = create_sample_scan_to_spectrum_map(mzmls=mzmls)
        for comet_psm in comet_psms:
            psms.append(
                cls(
                    spectrum=sample_scan_to_spectrum_map[
                        (comet_psm.sample, comet_psm.scan)
                    ],
                    peptide=comet_psm.seq,
                    peak_to_ion_ppm_tolerance=peak_to_ion_ppm_threshold,
                    comet_psm=comet_psm,
                )
            )
        return psms

    @property
    def seq(self):
        return self.peptide

    @property
    def df(self):
        return list_to_df(pydantic_list=self.peak_ion_matches)

    @property
    def num_b_ions_supported(self):
        return len(
            self.df[self.df.ion_type == B_ION_TYPE].groupby(
                by=["ion_charge", "ion_seq"]
            )
        )

    @property
    def num_ions_supported(self):
        return self.num_b_ions_supported + self.num_y_ions_supported

    @property
    def num_y_ions_supported(self):
        return len(
            self.df[self.df.ion_type == Y_ION_TYPE].groupby(
                by=["ion_charge", "ion_seq"]
            )
        )

    def sequences_supported(
        self,
        ion_type: Literal[B_ION_TYPE, Y_ION_TYPE],
    ):
        if len(self.peak_ion_matches) == 0:
            return set()
        return set(self.df.loc[self.df.ion_type == ion_type, "ion_seq"].unique())

    @property
    def prefixes_supported(self):
        return self.sequences_supported(ion_type=B_ION_TYPE)

    @property
    def suffixes_supported(self):
        return self.sequences_supported(ion_type=Y_ION_TYPE)

    def left_seq_support(self, left_seq: str):
        assert (
            self.peptide[: len(left_seq)] == left_seq
        ), f"{self.peptide} does not start with {left_seq}"
        prefixes = get_b_ion_prefixes(seq=left_seq)
        return len(
            self.df[
                (self.df.ion_type == B_ION_TYPE) & (self.df.ion_seq.isin(prefixes))
            ].groupby(by=["ion_seq"])
        )

    def right_seq_support(self, right_seq: str):
        assert (
            self.peptide[-len(right_seq) :] == right_seq
        ), f"{self.peptide} does not end with {right_seq}"
        suffixes = get_y_ion_suffixes(seq=right_seq)
        return len(
            self.df[
                (self.df.ion_type == Y_ION_TYPE) & (self.df.ion_seq.isin(suffixes))
            ].groupby(by=["ion_seq"])
        )

    @property
    def intensity_supported(self):
        return sum(
            [peak_ion_match.peak_intensity for peak_ion_match in self.peak_ion_matches]
        )

    @property
    def prop_prefixes_supported(self):
        return len(self.prefixes_supported) / len(self.peptide)

    @property
    def prop_suffixes_supported(self):
        return len(self.suffixes_supported) / len(self.peptide)


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: UnpositionedProductIon
    peaks: List[Peak]


def spectrum_peptide_plot(
    spectrum: Spectrum,
    seq: str,
    peak_to_ion_ppm_tolerance: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
):
    ion_intensity = max(peak.intensity for peak in spectrum.peaks) / 2
    psm = PSM.from_spectrum_and_seq(
        spectrum=spectrum, seq=seq, peak_to_ion_ppm_tolerance=peak_to_ion_ppm_tolerance
    )
    peptide = Peptide(seq=seq)
    product_ions = peptide.product_ions(
        charges=list(range(1, spectrum.precursor_charge + 1)),
    )
    _, axs = fig_setup()
    ax = axs[0]

    # Plot spectrum
    spectrum.plot_spectrum(
        ax=ax,
    )
    # Plot all product ions below spectrum
    plot_peaks(
        ax=ax,
        peaks=[Peak(mz=ion.mz, intensity=-1 * ion_intensity) for ion in product_ions],
    )
    # Plot ions that match peaks in a different color
    plot_peaks(
        ax=ax,
        peaks=[
            Peak(mz=match.ion_mz, intensity=-1 * ion_intensity)
            for match in psm.peak_ion_matches
        ],
        color="red",
    )
    # Plot peaks that match ions in a different color
    plot_peaks(
        ax=ax,
        peaks=[
            Peak(mz=match.peak_mz, intensity=match.peak_intensity)
            for match in psm.peak_ion_matches
        ],
        color="red",
        label="peak-to-ion matches",
    )
    # Plot y=0 line
    ax.axhline(0, color="black", linestyle="-", linewidth=0.5)

    # Finishing touches
    ax.set_ylim(bottom=-ion_intensity * 1.2)
    set_title_axes_labels(
        ax=ax,
        xlabel="m/z",
        ylabel="Intensity",
        title=f"peptide: {seq}\nspectrum: {spectrum.uid}\nPPM tol: {peak_to_ion_ppm_tolerance}",
    )
    finalize(ax)


def get_peak_product_ion_matches(
    spectrum: Spectrum,
    peptide: Union[Peptide, str],
    ion_types: Set[Literal[B_ION_TYPE, Y_ION_TYPE]] = {B_ION_TYPE, Y_ION_TYPE},
    peak_to_ion_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
) -> List[PeakIonMatch]:
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
    assert 2 * (len(peptide.seq) - 1) * len(charges) == len(product_ions)

    # Get peaks that match a product ion
    peak_ion_matches = []
    for ion in product_ions:
        matching_peaks = get_peaks_near_mz(
            query_mz=ion.mz,
            peaks=spectrum.peaks,
            ppm_tolerance=peak_to_ion_ppm_tolerance,
        )

        peak_ion_matches.extend(
            [
                PeakIonMatch(
                    ion_mz=ion.mz,
                    ion_charge=ion.charge,
                    ion_type=ion.ion_type,
                    ion_seq=ion.seq,
                    peak_mz=peak.mz,
                    peak_intensity=peak.intensity,
                    sample=spectrum.sample,
                    scan=spectrum.scan,
                )
                for peak in matching_peaks
            ]
        )

    return peak_ion_matches


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
