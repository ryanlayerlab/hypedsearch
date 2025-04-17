import logging
import multiprocessing
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from copy import deepcopy
from dataclasses import dataclass, field
from itertools import groupby, product
from pathlib import Path
from time import time
from typing import DefaultDict, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib.pyplot import Axes

from src.constants import (
    B_ION_AS_INT,
    DEFAULT_PPM_TOLERANCE,
    ION_INT_TO_TYPE,
    MOUSE_PROTEOME,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_AS_INT,
    IonTypes,
)
from src.mass_spectra import (
    Peak,
    Spectrum,
    get_specific_spectrum_by_sample_and_scan_num,
    plot_peaks,
)
from src.peptides_and_ions import (
    Peptide,
    ProductIon,
    compute_peptide_mz,
    get_product_ion_creator,
    write_fasta,
)
from src.plot_utils import fig_setup, finalize, set_title_axes_labels
from src.protein_product_ion_database import (
    DbKmer,
    IonWithSeq,
    PeakWithMatchingProductIons,
    PositionedIon,
    ProteinProductIonDb,
    get_aa_seq_from_db,
    get_positions_in_proteins_of_peak_matching_ions,
    get_product_ions_matching_spectrum,
)
from src.sql_database import Sqlite3Database, SqlTableRow
from src.utils import (
    Position,
    flatten_list_of_lists,
    get_time_in_diff_units,
    mass_difference_in_ppm,
    relative_ppm_tolerance_in_daltons,
    run_in_parallel,
)

logger = logging.getLogger(__name__)


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: ProductIon
    peaks: List[Peak]


class PeptideSpectrumComparison:
    def __init__(self, spectrum: Spectrum, peptide: Union[Peptide, str]):
        self.spectrum = spectrum
        if isinstance(peptide, str):
            self.peptide = Peptide(seq=peptide)
        elif isinstance(peptide, Peptide):
            self.peptide = peptide
        else:
            raise RuntimeError(f"Peptide isn't the correct type!")
        self.product_ions_with_matching_peaks = None
        self.product_ion_seqs_with_matching_peaks = None
        self.total_intensity = sum([peak.intensity for peak in spectrum.peaks])
        self.prop_intensity_supported_by_peptide = None
        self.uniq_matching_peaks = None
        self.compare()

    def compare(
        self,
        ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
        peak_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
    ):
        product_ions_with_matching_peaks = get_peaks_that_match_peptide_product_ions(
            spectrum=self.spectrum,
            peptide=self.peptide,
            ion_types=ion_types,
            peak_ppm_tolerance=peak_ppm_tolerance,
        )
        product_ion_seqs_with_matching_peaks = (
            group_product_ions_and_matching_peaks_by_charge_and_ion_type(
                product_ions_with_matching_peaks=product_ions_with_matching_peaks
            )
        )
        self.product_ions_with_matching_peaks = product_ions_with_matching_peaks
        self.product_ion_seqs_with_matching_peaks = product_ion_seqs_with_matching_peaks

        # Get the unique peaks that match a product ion
        peaks = flatten_list_of_lists(
            self.product_ion_seqs_with_matching_peaks["matching_peaks"].to_list()
        )
        uniq_peak_ids = set(peak.id for peak in peaks)
        self.uniq_matching_peaks = list(
            filter(lambda peak: peak.id in uniq_peak_ids, peaks)
        )

        self.prop_intensity_supported_by_peptide = (
            sum(peak.intensity for peak in self.uniq_matching_peaks)
            / self.total_intensity
        )

        # Convert DF to list of dataframes:
        peak_product_ion_matches = []
        for _, row in self.product_ions_with_matching_peaks.iterrows():
            for peak in row.matching_peaks:

                peak_product_ion_matches.append(
                    PeakProductIonMatch(
                        seq=row.seq, peak=peak, ion_type=row.ion_type, charge=row.charge
                    )
                )
        self.peak_product_ion_matches = peak_product_ion_matches

    def plot_ions(
        self,
        ax: Axes,
        ion_types=["b", "y"],
        intensity: Optional[float] = None,
        log_intensity: bool = False,
    ):
        colors = {"b": "blue", "y": "red"}
        for ion_type in ion_types:
            # Get peaks of given ion type
            peak_product_ion_matches = list(
                filter(lambda x: x.ion_type == ion_type, self.peak_product_ion_matches)
            )
            peaks = [
                peak_product_ion_match.peak
                for peak_product_ion_match in peak_product_ion_matches
            ]
            if intensity is not None:
                peaks = [Peak(mz=peak.mz, intensity=intensity) for peak in peaks]

            # Plot them
            label = f"{ion_type} (n={len(peaks)})"
            plot_peaks(
                ax=ax,
                peaks=peaks,
                label=label,
                color=colors[ion_type],
                log_intensity=log_intensity,
            )

            # Label them
            max_intensity = max([peak.intensity for peak in self.spectrum.peaks])
            for peak_product_ion_match in peak_product_ion_matches:
                peak = peak_product_ion_match.peak
                _ = ax.text(
                    x=peak.mz,
                    # peak.intensity + 0.02 * max_intensity,  # y value
                    y=max_intensity + 0.01 * max_intensity,
                    s=peak_product_ion_match.seq,  # label
                    rotation=90,
                    ha="center",
                    va="bottom",
                    fontsize=5,
                    color="black",
                )

    def plot(
        self,
        ax: Optional[Axes] = None,
        intensity: Optional[float] = None,
        title: Optional[str] = None,
        title_pad: int = 20,
        log_intensity: bool = False,
    ):
        if ax is None:
            _, axs = fig_setup()
            ax = axs[0]
        self.spectrum.plot_spectrum(ax=ax, log_intensity=log_intensity)
        self.plot_ions(ax=ax, intensity=intensity, log_intensity=log_intensity)
        if title is None:
            _ = ax.set_title(self.peptide.seq, pad=title_pad)
        else:
            _ = ax.set_title(title, pad=title_pad)

        return ax


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


@dataclass
class PeakProductIonMatch:
    seq: str
    peak: Peak
    ion_type: str
    charge: int


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
