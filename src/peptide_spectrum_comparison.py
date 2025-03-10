from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from src.constants import IonTypes
from src.mass_spectra import Peak, Spectrum
from src.peptides_and_ions import Peptide, ProductIon, get_product_ion_creator
from src.protein_product_ion_database import (
    DbProductIon,
    ProteinProductIonDb,
    get_aa_seq_from_db,
)
from src.utils import Position, mass_difference_in_ppm


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: ProductIon
    peaks: List[Peak]


@dataclass
class PeptideSpectrumComparison:
    num_peaks: int
    num_peaks_with_a_product_ion_match: int
    num_product_ions: int
    num_product_ions_with_match: int
    num_peaks_matching_product_ion: Dict[Tuple[str, str], int]


def get_peaks_near_mz(
    query_mz: float, peaks: List[Peak], ppm_tolerance: float
) -> List[Peak]:
    matching_peaks = []
    for peak in peaks:
        if (
            mass_difference_in_ppm(ref_mass=peak.mz, query_mass=query_mz)
            <= ppm_tolerance
        ):
            matching_peaks.append(peak)
    return matching_peaks


def compare_peptide_to_spectrum(
    peptide: Peptide,
    spectrum: Spectrum,
    ppm_tolerance: int,
    ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
):
    """ """
    # The product ions will have charge <= precursor's charge
    charges_to_consider = list(range(1, spectrum.precursor_charge + 1))

    # For each theoretical product ion, get the peaks that match product ion
    product_ions = peptide.product_ions(
        ion_types=ion_types, charges=charges_to_consider
    )
    product_ion_peak_matches = []
    for _, product_ion in enumerate(product_ions):
        product_ion_peak_matches.append(
            ProductIonWithMatchingPeaks(
                product_ion=product_ion,
                peaks=get_peaks_near_mz(
                    query_mz=product_ion.neutral_mass,
                    peaks=spectrum.peaks,
                    ppm_tolerance=ppm_tolerance,
                ),
            )
        )

    # Get the peaks that have a product ion match
    num_peaks = len(spectrum.peaks)
    peaks_with_match = set()
    for product_ion_with_matching_peaks in product_ion_peak_matches:
        matching_peaks = product_ion_with_matching_peaks.peaks
        if len(matching_peaks) > 0:
            peaks_with_match.update([peak.id for peak in matching_peaks])
    num_peaks_with_match = len(peaks_with_match)

    # Group theoretical product ions (and their matching peaks) by the product ion's AA seq
    # ion type: (seq, ion_type)
    product_ions_by_seq_ion_type = defaultdict(list)
    for peaks_matching_product_ion in product_ion_peak_matches:
        key = (
            peaks_matching_product_ion.product_ion.seq,
            peaks_matching_product_ion.product_ion.ion_type_as_str,
        )  # Tuple of attributes to group by
        product_ions_by_seq_ion_type[key].append(peaks_matching_product_ion)

    num_peaks_matching_product_ion = {}
    for key, product_ion_matches in product_ions_by_seq_ion_type.items():
        num_peaks_matching_product_ion[key] = sum(
            [len(x.peaks) for x in product_ion_matches]
        )
    num_product_ions_with_match = sum(
        [num_peaks > 0 for num_peaks in num_peaks_matching_product_ion.values()]
    )

    # Verify that the number of (seq, ion_type) groups is the same as the number of
    # product ion sequences (ignoring charge)
    # TODO: THIS VALIDATION SHOULDN'T BE HERE; IT SHOULD BE IN A TEST
    num_product_ions = 0
    for ion_type in ion_types:
        seq_generator = get_product_ion_creator(ion_type=ion_type)
        num_product_ions += len(
            seq_generator.generate_product_ion_seqs(seq=peptide.seq)
        )
    assert (
        len(product_ions_by_seq_ion_type) == num_product_ions
    ), "Something weird is going on with the number of product ions"

    return PeptideSpectrumComparison(
        num_peaks=num_peaks,
        num_peaks_with_a_product_ion_match=num_peaks_with_match,
        num_product_ions=num_product_ions,
        num_peaks_matching_product_ion=num_peaks_matching_product_ion,
        num_product_ions_with_match=num_product_ions_with_match,
    )


@dataclass
class Cluster:
    # protein_id: int
    ions: List[DbProductIon]
    ion_type: IonTypes
    # region: Position = field(init=False)

    def __post_init__(self):
        if self.ion_type == IonTypes.B_ION_TYPE:
            # Start positions should be the same
            end_pos = self.ions[0].inclusive_start
            assert all([ion.inclusive_start == end_pos for ion in self.ions])
        else:
            # End positions should be the same
            end_pos = self.ions[0].exclusive_end
            assert all([ion.exclusive_end == end_pos for ion in self.ions])

    @property
    def protein_id(self):
        uniq_protein_ids = np.unique([ion.protein_id for ion in self.ions])
        assert (
            len(uniq_protein_ids) == 1
        ), f"There should only be one protein ID in a cluster. Found {len(uniq_protein_ids)}"
        return uniq_protein_ids[0]

    @property
    def region(self) -> Position:
        min_start = min([ion.inclusive_start for ion in self.ions])
        max_end = max([ion.exclusive_end for ion in self.ions])
        return Position(inclusive_start=min_start, exclusive_end=max_end)

    @property
    def num_supporting_ions(self):
        return len(self.ions)

    def get_aa_seq(self, db: ProteinProductIonDb):
        region = self.region
        return get_aa_seq_from_db(
            protein_id=self.protein_id,
            inclusive_start=region.inclusive_start,
            exclusive_end=region.exclusive_end,
            db=db,
        )
