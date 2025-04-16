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

    @property
    def support(self):
        return len(self.ions)

    @property
    def length(self):
        return self.exclusive_end - self.inclusive_start

    @property
    def smallest_ion(self):
        smallest_ion = min(self.ions, key=lambda ion: ion.exclusive_end)
        return smallest_ion.seq


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
class SeqWithMz:
    seq: str
    mz: float


@dataclass
class ExtendedCluster:
    cluster: Cluster
    # num_extended_aa: int
    # seqs: List[SeqWithMz]
    extended_seq: Optional[str]

    @property
    def protein_id(self):
        return self.cluster.protein_id

    @property
    def start(self):
        return self.cluster.inclusive_start

    @property
    def end(self):
        return self.cluster.exclusive_end

    @property
    def smallest_ion(self):
        return self.cluster.smallest_ion


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

    if len(seqs) == 0:
        extended_seq = None
    else:
        extended_seq = seqs[-1].seq
    extended_cluster = ExtendedCluster(cluster=cluster, extended_seq=extended_seq)
    # extended_cluster = ExtendedCluster(
    #     cluster=cluster, num_extended_aa=num_extended_aa - 1, seqs=seqs
    # )
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


def extend_clusters(
    spectrum_clusters: SpectrumClusters,
    db: ProteinProductIonDb,
    precursor_mz: float,
    precursor_charge: int,
) -> SpectrumExtendedClusters:
    # Extend the clusters
    t0 = time()
    # b-clusters
    b_extended_clusters = [
        extend_cluster(
            cluster=cluster,
            db=db,
            precursor_charge=precursor_charge,
            precursor_mz=precursor_mz,
        )
        for cluster in spectrum_clusters.b
    ]

    # y-clusters
    y_extended_clusters = [
        extend_cluster(
            cluster=cluster,
            db=db,
            precursor_charge=precursor_charge,
            precursor_mz=precursor_mz,
        )
        for cluster in spectrum_clusters.y
    ]
    logger.info(f"Getting extended clusters took {round(time()-t0, 2)} seconds")

    return SpectrumExtendedClusters(b=b_extended_clusters, y=y_extended_clusters)


@dataclass
class HybridPeptide:
    b_seq: str
    y_seq: str
    b_prot_ids: Set[int] = field(default_factory=set)
    b_prot_names: Set[str] = field(default_factory=set)
    y_prot_ids: Set[int] = field(default_factory=set)
    y_prot_names: Set[str] = field(default_factory=set)

    @property
    def seq(self):
        return self.b_seq + self.y_seq

    @property
    def fasta_name(self):
        return f"hybrid_{self.b_seq}-{self.y_seq}"

    def set_protein_names(self, prot_id_to_name_map: Dict[int, str]):
        self.b_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.b_prot_ids
        )
        self.y_prot_names = set(
            prot_id_to_name_map[prot_id] for prot_id in self.y_prot_ids
        )

    def set_fasta_description(
        self,
        prot_id_to_name_map: Dict[int, str],
    ):
        self.set_protein_names(prot_id_to_name_map=prot_id_to_name_map)
        self.fasta_description = f"b-prots:{','.join(self.b_prot_names)} y-prots:{','.join(self.y_prot_names)}"

    def mz(self, charge: int):
        return compute_peptide_mz(aa_seq=self.b_seq + self.y_seq, charge=charge)


def create_hybrids_fasta(
    hybrids: List[HybridPeptide],
    fasta_path: Path,
    prot_id_to_name_map: Dict[int, str],
    other_prots: Optional[Union[List[Peptide], Path]] = None,
) -> List[Peptide]:
    prots = []
    if other_prots is not None:
        if isinstance(other_prots, list):
            prots = other_prots
        elif isinstance(other_prots, Path):
            prots = Peptide.from_fasta(fasta_path=other_prots)
    for idx, hybrid in enumerate(hybrids):
        hybrid.set_fasta_description(prot_id_to_name_map=prot_id_to_name_map)
        new_peptide = Peptide(
            seq=hybrid.seq, name=hybrid.fasta_name, desc=hybrid.fasta_description
        )

        prots.append(new_peptide)

    write_fasta(peptides=prots, output_path=fasta_path)
    return prots


def get_sequences_from_extended_cluster(cluster: ExtendedCluster):
    return [
        cluster.extended_seq[:i]
        for i in range(len(cluster.smallest_ion), len(cluster.extended_seq) + 1)
    ]


def get_sequences_to_form_hybrids_out_of(
    b_ext_clusters: List[ExtendedCluster],
    y_ext_clusters: List[ExtendedCluster],
) -> List[HybridPeptide]:
    t0 = time()
    logger.info("Making hybrids from extended clusters...")

    b_seqs = defaultdict(set)
    for cluster in b_ext_clusters:
        prot_id = cluster.cluster.protein_id
        for seq in get_sequences_from_extended_cluster(cluster=cluster):
            b_seqs[seq].add(prot_id)

    y_seqs = defaultdict(set)
    for cluster in y_ext_clusters:
        prot_id = cluster.cluster.protein_id
        y_seq = cluster.extended_seq
        y_kmers = [
            kmer.seq for kmer in Peptide(seq=y_seq).kmers(min_k=1, max_k=len(y_seq))
        ]
        for kmer in y_kmers:
            y_seqs[kmer].add(prot_id)

    return b_seqs, y_seqs


@dataclass
class SeqWithMass(SqlTableRow):
    seq: str
    mz: float

    @classmethod
    def from_seq(cls, seq: str, charge: int):
        mz = compute_peptide_mz(aa_seq=seq, charge=charge)
        return cls(seq=seq, mz=mz)


def get_hybrids_from_b_and_y_seqs(
    b_seqs: DefaultDict[str, Set[int]],
    y_seqs: DefaultDict[str, Set[int]],
    precursor_charge: int,
    precursor_mz: float,
    precursor_mz_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
):
    t0 = time()

    # Create database of y-sequences
    y_rows = [SeqWithMass.from_seq(seq=seq, charge=precursor_charge) for seq in y_seqs]
    db = Sqlite3Database()
    table_name = "ys"
    db.create_table_from_dataclass(table_name=table_name, obj=SeqWithMass)
    db.insert_dataclasses(table_name=table_name, data_classes=y_rows)
    db.add_index(table_name=table_name, index_name="mass", colms_to_index=["mz"])

    # For each b-sequence search the y-sequence database to find the hybrids
    # that would produce a peptide within the given PPM tolerance of the precursor m/z
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
    )
    adjusted_precursor_mz = precursor_mz + (WATER_MASS / precursor_charge) + PROTON_MASS
    potench_hybrids = []
    for b_seq in b_seqs:
        b_seq_mz = compute_peptide_mz(aa_seq=b_seq, charge=precursor_charge)
        lower_bdd = adjusted_precursor_mz - mz_tolerance - b_seq_mz
        upper_bdd = adjusted_precursor_mz + mz_tolerance - b_seq_mz
        query = f"""
            SELECT
                *
            FROM {table_name} as ion
            WHERE ion.mz BETWEEN {lower_bdd} AND {upper_bdd}
        """
        matches = [match["seq"] for match in db.read_query(query=query)]
        for y_seq in matches:
            potench_hybrids.append(
                HybridPeptide(
                    b_seq=b_seq,
                    y_seq=y_seq,
                    b_prot_ids=b_seqs[b_seq],
                    y_prot_ids=y_seqs[y_seq],
                )
            )
    t1 = time()
    logger.info(f"Creating hybrids took {get_time_in_diff_units(t1 - t0)}")
    return potench_hybrids


@dataclass
class HSResults:
    spectrum: Spectrum
    clusters: SpectrumExtendedClusters
    b_seqs: DefaultDict[str, Set[int]]
    y_seqs: DefaultDict[str, Set[int]]
    hybrids: List[HybridPeptide]


def get_extended_clusters(
    spectrum: Spectrum,
    db_path: Path,
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
) -> SpectrumExtendedClusters:
    # Connect to the database
    db = ProteinProductIonDb(db_path=db_path, overwrite=False)

    # Find ions that match spectrum peaks
    peaks_with_matches = get_product_ions_matching_spectrum(
        spectrum=spectrum,
        db=db,
        peak_product_ion_ppm_tolerance=peak_to_ion_ppm_tol,
    )
    positioned_ions = get_positions_in_proteins_of_peak_matching_ions(
        peaks_with_matches=peaks_with_matches,
        db=db,
    )

    # Get clusters
    logger.info("Getting clusters...")
    t0 = time()
    clusters = get_clusters_from_ions(ions=positioned_ions)
    logger.info(f"Getting clusters took {get_time_in_diff_units(time() - t0)}")

    # Extend clusters
    logger.info("Extending clusters...")
    t0 = time()
    extended_clusters = extend_clusters(
        spectrum_clusters=clusters,
        db=db,
        precursor_charge=spectrum.precursor_charge,
        precursor_mz=spectrum.precursor_mz,
    )
    logger.info(f"Extending clusters took {get_time_in_diff_units(time() - t0)}")

    return extended_clusters


def get_hybrids(
    db_path: Path,
    # sample: str,
    # scan_num: int,
    spectrum: Spectrum,
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    precursor_mz_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
):
    # # Get spectrum
    # spectrum = get_specific_spectrum_by_sample_and_scan_num(
    #     sample=sample, scan_num=scan_num
    # )
    # spectrum.filter_to_top_n_peaks(n=100)
    clusters = get_extended_clusters(
        spectrum=spectrum, db_path=db_path, peak_to_ion_ppm_tol=peak_to_ion_ppm_tol
    )
    if (len(clusters.b) == 0) or (len(clusters.y) == 0):
        return HSResults(
            spectrum=spectrum,
            clusters=clusters,
            b_seqs=[],
            y_seqs=[],
            hybrids=[],
        )

    b_seqs, y_seqs = get_sequences_to_form_hybrids_out_of(
        b_ext_clusters=clusters.b, y_ext_clusters=clusters.y
    )

    hybrids = get_hybrids_from_b_and_y_seqs(
        b_seqs=b_seqs,
        y_seqs=y_seqs,
        precursor_charge=spectrum.precursor_charge,
        precursor_mz=spectrum.precursor_mz,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
    )
    return HSResults(
        spectrum=spectrum,
        clusters=clusters,
        b_seqs=b_seqs,
        y_seqs=y_seqs,
        hybrids=hybrids,
    )
