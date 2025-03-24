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
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from src.constants import (
    B_ION_AS_INT,
    DEFAULT_PPM_TOLERANCE,
    ION_INT_TO_TYPE,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_AS_INT,
    IonTypes,
)
from src.mass_spectra import Peak, Spectrum
from src.peptides_and_ions import (
    Peptide,
    ProductIon,
    compute_peptide_mz,
    get_product_ion_creator,
    write_fasta,
)
from src.protein_product_ion_database import (
    DbProductIon,
    IonWithSeq,
    PeakWithMatchingProductIons,
    PositionedIon,
    ProteinProductIonDb,
    get_aa_seq_from_db,
)
from src.utils import (
    Position,
    flatten_list_of_lists,
    mass_difference_in_ppm,
    relative_ppm_tolerance_in_daltons,
)

logger = logging.getLogger(__name__)


@dataclass
class ProductIonWithMatchingPeaks:
    product_ion: ProductIon
    peaks: List[Peak]


class PeptideSpectrumComparison:
    def __init__(self, spectrum: Spectrum, peptide: Peptide):
        self.spectrum = spectrum
        self.peptide = peptide
        self.product_ions_with_matching_peaks = None
        self.product_ion_seqs_with_matching_peaks = None
        self.total_intensity = sum([peak.intensity for peak in spectrum.peaks])
        self.prop_intensity_supported_by_peptide = None
        self.uniq_matching_peaks = None

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


# def compare_peptide_to_spectrum(
#     spectrum: Spectrum,
#     peptide: Peptide,
#     ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
#     peak_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
# ) -> List[pd.DataFrame]:
#     """
#     Compare the given spectrum to the given peptide. This method helps evaluate
#     how strong the evidence is for a peptide-spectrum match (PSM).
#     """
#     # Product ions will have charge <= precursor's charge
#     charges = list(range(1, spectrum.precursor_charge + 1))

#     # Get product ions of the proposed peptide
#     product_ions = peptide.product_ions(ion_types=ion_types, charges=charges)
#     assert 2 * len(peptide.seq) * len(charges) == len(product_ions)

#     # Get peaks that match a product ion
#     product_ions_with_matching_peaks = []
#     for ion in product_ions:
#         matching_peaks = get_peaks_near_mz(
#             query_mz=ion.neutral_mass,
#             peaks=spectrum.peaks,
#             ppm_tolerance=peak_ppm_tolerance,
#         )
#         product_ions_with_matching_peaks.append(
#             [ion.charge, ion.ion_type_as_str, ion.neutral_mass, ion.seq, matching_peaks]
#         )
#     product_ions_with_matching_peaks = pd.DataFrame(
#         product_ions_with_matching_peaks,
#         columns=["charge", "ion_type", "m/z", "seq", "matching_peaks"],
#     )

#     # Add metadata to dataframe
#     product_ions_with_matching_peaks.attrs["peak_ppm_tolerance"] = peak_ppm_tolerance
#     product_ions_with_matching_peaks.attrs["ion_types"] = [
#         ion_type.value for ion_type in ion_types
#     ]
#     # data.head(3)

#     # Group by (seq, ion_type) ignoring charge
#     product_ion_seqs_with_matching_peaks = {"b": [], "y": []}
#     for name, group in product_ions_with_matching_peaks.groupby(by=["seq", "ion_type"]):
#         matching_peaks = flatten_list_of_lists(
#             [matching_peaks for matching_peaks in group["matching_peaks"]]
#         )
#         num_matching_peaks = len(matching_peaks)
#         product_ion_seqs_with_matching_peaks[name[1]].append(
#             [
#                 name[0],
#                 name[1],
#                 num_matching_peaks,
#                 matching_peaks,
#             ]
#         )
#     for key, tmp_data in product_ion_seqs_with_matching_peaks.items():
#         df = pd.DataFrame(
#             tmp_data,
#             columns=["seq", "ion_type", "num_matching_peaks", "matching_peaks"],
#         )
#         df.sort_values(
#             by=["seq"], key=lambda x: x.str.len(), inplace=True, ignore_index=True
#         )
#         product_ion_seqs_with_matching_peaks[key] = df

#     product_ion_seqs_with_matching_peaks = pd.concat(
#         [
#             product_ion_seqs_with_matching_peaks["b"],
#             product_ion_seqs_with_matching_peaks["y"],
#         ],
#         ignore_index=True,
#     )

#     return product_ions_with_matching_peaks, product_ion_seqs_with_matching_peaks


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


# def compare_peptide_to_spectrum(
#     peptide: Peptide,
#     spectrum: Spectrum,
#     ppm_tolerance: int,
#     ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
# ):
#     """ """
#     # The product ions will have charge <= precursor's charge
#     charges_to_consider = list(range(1, spectrum.precursor_charge + 1))

#     # For each theoretical product ion, get the peaks that match product ion
#     product_ions = peptide.product_ions(
#         ion_types=ion_types, charges=charges_to_consider
#     )
#     product_ion_peak_matches = []
#     for _, product_ion in enumerate(product_ions):
#         product_ion_peak_matches.append(
#             ProductIonWithMatchingPeaks(
#                 product_ion=product_ion,
#                 peaks=get_peaks_near_mz(
#                     query_mz=product_ion.neutral_mass,
#                     peaks=spectrum.peaks,
#                     ppm_tolerance=ppm_tolerance,
#                 ),
#             )
#         )

#     # Get the peaks that have a product ion match
#     num_peaks = len(spectrum.peaks)
#     peaks_with_match = set()
#     for product_ion_with_matching_peaks in product_ion_peak_matches:
#         matching_peaks = product_ion_with_matching_peaks.peaks
#         if len(matching_peaks) > 0:
#             peaks_with_match.update([peak.id for peak in matching_peaks])
#     num_peaks_with_match = len(peaks_with_match)

#     # Group theoretical product ions (and their matching peaks) by the product ion's AA seq
#     # ion type: (seq, ion_type)
#     product_ions_by_seq_ion_type = defaultdict(list)
#     for peaks_matching_product_ion in product_ion_peak_matches:
#         key = (
#             peaks_matching_product_ion.product_ion.seq,
#             peaks_matching_product_ion.product_ion.ion_type_as_str,
#         )  # Tuple of attributes to group by
#         product_ions_by_seq_ion_type[key].append(peaks_matching_product_ion)

#     num_peaks_matching_product_ion = {}
#     for key, product_ion_matches in product_ions_by_seq_ion_type.items():
#         num_peaks_matching_product_ion[key] = sum(
#             [len(x.peaks) for x in product_ion_matches]
#         )
#     num_product_ions_with_match = sum(
#         [num_peaks > 0 for num_peaks in num_peaks_matching_product_ion.values()]
#     )

#     # Verify that the number of (seq, ion_type) groups is the same as the number of
#     # product ion sequences (ignoring charge)
#     # TODO: THIS VALIDATION SHOULDN'T BE HERE; IT SHOULD BE IN A TEST
#     num_product_ions = 0
#     for ion_type in ion_types:
#         seq_generator = get_product_ion_creator(ion_type=ion_type)
#         num_product_ions += len(
#             seq_generator.generate_product_ion_seqs(seq=peptide.seq)
#         )
#     assert (
#         len(product_ions_by_seq_ion_type) == num_product_ions
#     ), "Something weird is going on with the number of product ions"

#     return PeptideSpectrumComparison(
#         num_peaks=num_peaks,
#         num_peaks_with_a_product_ion_match=num_peaks_with_match,
#         num_product_ions=num_product_ions,
#         num_peaks_matching_product_ion=num_peaks_matching_product_ion,
#         num_product_ions_with_match=num_product_ions_with_match,
#     )


def ions_as_df(ions: List[DbProductIon]):
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


@dataclass
class Cluster:
    # ions: List[DbProductIon]
    ions: List[PositionedIon]
    ion_type: str = field(init=False)
    protein_id: int = field(init=False)
    inclusive_start: str = field(init=False)
    exclusive_end: int = field(init=False)

    def __post_init__(self):
        # Get ion type
        ion_type = self.ions[0].ion_type
        assert all(
            [ion.ion_type == ion_type for ion in self.ions]
        ), "All ions should have the same ion type!"
        self.ion_type = ion_type

        # Get protein ID
        p_id = self.ions[0].protein_id
        assert all(
            [ion.protein_id == p_id for ion in self.ions]
        ), "All ions should have same protein ID"
        self.protein_id = p_id

        # Make sure that
        # - for b-ions, start positions should all be the same
        # - for y-ions, end positions should all be the same
        if self.ion_type == "b":
            start = self.ions[0].inclusive_start
            assert all(
                [ion.inclusive_start == start for ion in self.ions]
            ), "All b-ions should have same start position"
        else:
            end = self.ions[0].exclusive_end
            assert all(
                [ion.exclusive_end == end for ion in self.ions]
            ), "All y-ions should have same end position"

        # Set start and end position
        position = get_start_and_end_positions_from_ions(ions=self.ions)
        self.inclusive_start = position.inclusive_start
        self.exclusive_end = position.exclusive_end

    def get_aa_seq(self, db: ProteinProductIonDb):
        return get_aa_seq_from_db(
            protein_id=self.protein_id,
            inclusive_start=self.inclusive_start,
            exclusive_end=self.exclusive_end,
            db=db,
        )

    @property
    def ions_as_df(self):
        return ions_as_df(ions=self.ions)

    def extend_cluster_to_precursor_mass(
        self, precursor_mz: float, precursor_charge: float
    ):
        pass

    @property
    def support(self):
        return len(self.ions)

    @property
    def length(self):
        return self.exclusive_end - self.inclusive_start


def get_b_clusters(
    ions: List[PositionedIon],
    # List[DbProductIon]
) -> List[Cluster]:
    """
    b-clusters are all those ions with the same protein ID AND start position
    """
    b_ions = list(filter(lambda ion: ion.ion_type == "b", ions))
    b_ions.sort(key=lambda ion: (ion.protein_id, ion.inclusive_start))

    b_clusters = []
    for (p_id, loc), group in groupby(
        b_ions, key=lambda ion: (ion.protein_id, ion.inclusive_start)
    ):
        b_clusters.append(Cluster(ions=list(group)))
    return b_clusters


def get_y_clusters(ions: List[PositionedIon]) -> List[Cluster]:
    """
    y-clusters are all those ions with the same protein ID AND end position
    """
    y_ions = list(filter(lambda ion: ion.ion_type == "y", ions))
    y_ions.sort(key=lambda ion: (ion.protein_id, ion.exclusive_end))

    y_clusters = []
    for (p_id, loc), group in groupby(
        y_ions, key=lambda ion: (ion.protein_id, ion.exclusive_end)
    ):
        y_clusters.append(Cluster(ions=list(group)))
    return y_clusters


def get_start_and_end_positions_from_ions(ions: List[DbProductIon]) -> Position:
    start = min([ion.inclusive_start for ion in ions])
    end = max([ion.exclusive_end for ion in ions])
    return Position(inclusive_start=start, exclusive_end=end)


def create_df_from_ions(ions: List[DbProductIon]) -> pd.DataFrame:
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
class SeqWithMz:
    seq: str
    mz: float


@dataclass
class ExtendedCluster:
    cluster: Cluster
    num_extended_aa: int
    seqs: List[SeqWithMz]


def extend_cluster(
    cluster: Cluster,
    precursor_charge: int,
    precursor_mz: float,
    db: Optional[ProteinProductIonDb] = None,
    db_path: Optional[str] = None,
) -> ExtendedCluster:

    if db is None:
        db = ProteinProductIonDb(db_path=db_path)

    aa_seq = cluster.get_aa_seq(db=db)
    cluster_mz = compute_peptide_mz(aa_seq=aa_seq, charge=precursor_charge)
    protein = db.get_protein_by_id(protein_id=cluster.protein_id).seq

    extended_cluster_mz = cluster_mz
    num_extended_aa = 0
    seqs = []
    while extended_cluster_mz < precursor_mz:
        seqs.append(
            SeqWithMz(
                seq=aa_seq,
                mz=extended_cluster_mz,
            )
        )
        num_extended_aa += 1
        if cluster.ion_type == "b":
            start = cluster.inclusive_start
            end = cluster.exclusive_end + num_extended_aa
        else:
            start = cluster.inclusive_start - num_extended_aa
            end = cluster.exclusive_end

        # Exit loop if we've reached the end of the protein
        if (start < 0) or (end > len(protein)):
            break

        aa_seq = protein[start:end]
        extended_cluster_mz = compute_peptide_mz(
            aa_seq=aa_seq,
            charge=precursor_charge,
        )

    extended_cluster = ExtendedCluster(
        cluster=cluster, num_extended_aa=num_extended_aa - 1, seqs=seqs
    )
    return extended_cluster


def filter_clusters(
    clusters: List[Cluster], min_len: int = 3, min_support: int = 2
) -> List[Cluster]:
    filtered_clusters = list(
        filter(
            lambda cluster: (cluster.support >= min_support)
            and (cluster.length >= min_len),
            clusters,
        )
    )
    return filtered_clusters


@dataclass
class SpectrumClusters:
    b: List[ExtendedCluster]
    y: List[ExtendedCluster]


@dataclass
class SpectrumExtendedClusters:
    b: List[ExtendedCluster]
    y: List[ExtendedCluster]


@dataclass
class SpectrumClusters:
    b: List[Cluster]
    y: List[Cluster]


def get_clusters_from_ions(ions: List[PositionedIon]) -> SpectrumClusters:
    # Group product ions to form clusters
    clusters = {}
    clusters["b"] = get_b_clusters(ions=ions)
    clusters["y"] = get_y_clusters(ions=ions)
    msg = f"Before filtering, number b-clusters = {len(clusters['b'])}; number y-clusters = {len(clusters['y'])}"
    logger.info(msg)

    # Filter the clusters by some criteria in filter_clusters
    for cluster_type in ["b", "y"]:
        clusters[cluster_type] = filter_clusters(clusters=clusters[cluster_type])
    msg = f"After filtering, number b-clusters = {len(clusters['b'])}; number y-clusters = {len(clusters['y'])}"
    logger.info(msg)

    return SpectrumClusters(b=clusters["b"], y=clusters["y"])


def get_extended_clusters(
    spectrum_clusters: SpectrumClusters, spectrum: Spectrum, db: ProteinProductIonDb
) -> SpectrumExtendedClusters:
    # Extend the clusters
    t0 = time()
    # b-clusters
    b_extended_clusters = [
        extend_cluster(
            cluster=cluster,
            db=db,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz=spectrum.precursor_mz,
        )
        for cluster in spectrum_clusters.b
    ]

    # y-clusters
    y_extended_clusters = [
        extend_cluster(
            cluster=cluster,
            db=db,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz=spectrum.precursor_mz,
        )
        for cluster in spectrum_clusters.y
    ]
    logger.info(f"Getting extended clusters took {round(time()-t0, 2)} seconds")

    return SpectrumExtendedClusters(b=b_extended_clusters, y=y_extended_clusters)


@dataclass
class HybridPeptide:
    b_seq: SeqWithMz
    b_protein_id: int
    y_seq: SeqWithMz
    y_protein_id: int


def get_hybrids_from_b_and_y_cluster(
    b_cluster: ExtendedCluster, y_cluster: ExtendedCluster, min_mz: float, max_mz: float
) -> List[HybridPeptide]:
    possible_hybrids = []
    b_prot_id = b_cluster.cluster.protein_id
    y_prot_id = y_cluster.cluster.protein_id
    # Case 1: b-cluster and y-cluster are from DIFFERENT proteins
    # if b_prot_id == y_prot_id:

    # if b_prot_id != y_prot_id:
    for b_seq, y_seq in product(b_cluster.seqs, y_cluster.seqs):
        if min_mz <= (b_seq.mz + y_seq.mz) <= max_mz:
            hybrid = HybridPeptide(
                b_seq=b_seq,
                y_seq=y_seq,
                b_protein_id=b_prot_id,
                y_protein_id=y_prot_id,
            )
            possible_hybrids.append(hybrid)

    # Case 2: b-cluster and y-cluster are from the same protein
    # else:

    return possible_hybrids


def get_possible_hybrids(
    extended_clusters: SpectrumExtendedClusters,
    spectrum: Spectrum,
    precursor_mz_ppm_tolerance: int = DEFAULT_PPM_TOLERANCE,
) -> List[HybridPeptide]:
    logger.info("Creating possible hybrids...")
    start_time = time()
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=precursor_mz_ppm_tolerance, ref_mass=spectrum.precursor_mz
    )
    adjusted_precursor_mz = (
        spectrum.precursor_mz + (WATER_MASS / spectrum.precursor_charge) + PROTON_MASS
    )
    min_mz = adjusted_precursor_mz - mz_tolerance
    max_mz = adjusted_precursor_mz + mz_tolerance

    possible_hybrids = []
    for b_cluster, y_cluster in product(extended_clusters.b, extended_clusters.y):
        possible_hybrids.extend(
            get_hybrids_from_b_and_y_cluster(
                b_cluster=b_cluster,
                y_cluster=y_cluster,
                min_mz=min_mz,
                max_mz=max_mz,
            )
        )

    # def process_pair(pair):
    #     b_cluster, y_cluster = pair
    #     return get_hybrids_from_b_and_y_cluster(
    #         b_cluster=b_cluster,
    #         y_cluster=y_cluster,
    #         min_mz=min_mz,
    #         max_mz=max_mz,
    #     )

    # with ThreadPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
    #     possible_hybrids = executor.map(
    #         process_pair, product(extended_clusters.b, extended_clusters.y)
    #     )
    logger.info(
        f"Getting possible hybrids took {round(time() - start_time, 2)} seconds"
    )
    return possible_hybrids


def create_new_fasta_including_hybrids(
    db_proteins: List[Peptide],
    hybrids: List[HybridPeptide],
    fasta_path: Path,
    protein_id_to_name_map: Dict[int, str],
):
    """ """
    new_proteins = deepcopy(db_proteins)
    for idx, hybrid in enumerate(hybrids):
        b_prot = protein_id_to_name_map[hybrid.b_protein_id]
        y_prot = protein_id_to_name_map[hybrid.y_protein_id]
        b_seq, y_seq = hybrid.b_seq.seq, hybrid.y_seq.seq
        name = f"hybrid-{idx}-{b_prot}_{b_seq}-{y_prot}_{y_seq}"
        seq = b_seq + y_seq
        new_peptide = Peptide(seq=seq, name=name)
        new_proteins.append(new_peptide)

    write_fasta(peptides=new_proteins, output_path=fasta_path)
