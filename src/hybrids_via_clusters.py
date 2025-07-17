import logging
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path
from time import time
from typing import DefaultDict, Dict, List, Literal, Optional, Set, Union
from venv import logger

from src.constants import (
    B_ION_TYPE,
    DEFAULT_PPM_TOLERANCE,
    MIN_CLUSTER_LENGTH,
    MIN_CLUSTER_SUPPORT,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
)
from src.create_db import KmerDatabase
from src.hypedsearch_utils import HybridPeptide, SeqWithMass
from src.mass_spectra import Spectrum

# from src.peptide_spectrum_comparison import (
#     get_start_and_end_positions_from_ions,
#     ions_as_df,
# )
from src.peptides_and_ions import Peptide, ProductIon, compute_peptide_mz
from src.sql_database import Sqlite3Database
from src.utils import (
    Position,
    get_positions_of_subseq_in_seq,
    get_time_in_diff_units,
    log_time,
    relative_ppm_tolerance_in_daltons,
)

logger = logging.getLogger(__name__)


@dataclass
class PositionedIon:
    seq: str
    charge: int
    ion_type: str
    protein: Union[int, str]
    inclusive_start: int
    exclusive_end: int


@dataclass
class Cluster:
    """Object to represent a b- or y-cluster of product ions."""

    ions: List[PositionedIon]
    ion_type: str = field(init=False)
    protein: int = field(init=False)
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
        protein = self.ions[0].protein
        assert all(
            [ion.protein == protein for ion in self.ions]
        ), "All ions should have same protein ID"
        self.protein = protein

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
        position = Position(
            inclusive_start=min([ion.inclusive_start for ion in self.ions]),
            exclusive_end=max([ion.exclusive_end for ion in self.ions]),
        )
        self.inclusive_start = position.inclusive_start
        self.exclusive_end = position.exclusive_end

    @staticmethod
    def get_cluster_start_and_end_position_from_constituent_ions(
        ions: List[PositionedIon],
    ) -> Position:
        start = min([ion.inclusive_start for ion in ions])
        end = max([ion.exclusive_end for ion in ions])
        return Position(inclusive_start=start, exclusive_end=end)

    # @property
    # def ions_as_df(self):
    #     pass
    #     # return ions_as_df(ions=self.ions)

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


@dataclass
class ExtendedCluster:
    cluster: Cluster
    extended_seq: Optional[str]

    @property
    def protein_id(self):
        return self.cluster.protein

    @property
    def start(self):
        return self.cluster.inclusive_start

    @property
    def end(self):
        return self.cluster.exclusive_end

    @property
    def smallest_ion(self):
        return self.cluster.smallest_ion

    @classmethod
    def from_cluster(
        cls,
        cluster: Cluster,
        protein_seq: str,
        precursor_mz_ppm_tol: float,
        precursor_charge: int,
        precursor_mz: float,
    ):
        """
        Extend a cluster until the precursor m/z is reached.
        """
        # A cluster has a start and end position and corresponds to some peptide.
        # Get the amino acid sequence of the peptide the cluster corresponds to, get the peptide's
        # m/z
        aa_seq = protein_seq[cluster.inclusive_start : cluster.exclusive_end]
        cluster_mz = compute_peptide_mz(aa_seq=aa_seq, charge=precursor_charge)

        extended_cluster_mz = cluster_mz
        num_extended_aa = 0
        seqs = []
        da_tol = relative_ppm_tolerance_in_daltons(
            ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
        )
        while extended_cluster_mz < precursor_mz + da_tol:
            # Extend the cluster by one amino acid
            seqs.append(
                SeqWithMass(
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
            if (start < 0) or (end > len(protein_seq)):
                break

            aa_seq = protein_seq[start:end]
            extended_cluster_mz = compute_peptide_mz(
                aa_seq=aa_seq,
                charge=precursor_charge,
            )

        # Grab the last sequence of the extended cluster
        if len(seqs) == 0:
            extended_seq = None
        else:
            extended_seq = seqs[-1].seq

        return cls(cluster=cluster, extended_seq=extended_seq)


@dataclass
class SpectrumClusters:
    b_clusters: List[Cluster]
    y_clusters: List[Cluster]

    @staticmethod
    def _filter_clusters(
        clusters: List[Cluster], min_len: int, min_support: int
    ) -> List[Cluster]:
        """
        Filter clusters by minimum length and support.
        """
        filtered_clusters = list(
            filter(
                lambda cluster: (cluster.support >= min_support)
                and (cluster.length >= min_len),
                clusters,
            )
        )
        return filtered_clusters

    def filter_clusters(
        self,
        min_len: int = MIN_CLUSTER_LENGTH,
        min_support: int = MIN_CLUSTER_SUPPORT,
    ):
        """
        Filter b- and y-clusters.
        """
        self.b_clusters = self._filter_clusters(
            clusters=self.b_clusters,
            min_len=min_len,
            min_support=min_support,
        )
        self.y_clusters = self._filter_clusters(
            clusters=self.y_clusters,
            min_len=min_len,
            min_support=min_support,
        )

    @staticmethod
    def get_b_clusters(
        positioned_ions: List[PositionedIon],
    ) -> List[Cluster]:
        """
        b-clusters are all those ions with the same protein ID AND start position
        """
        b_ions = list(filter(lambda ion: ion.ion_type == "b", positioned_ions))
        b_ions.sort(key=lambda ion: (ion.protein, ion.inclusive_start))
        b_clusters = []
        for _, group in groupby(
            b_ions, key=lambda ion: (ion.protein, ion.inclusive_start)
        ):
            b_clusters.append(Cluster(ions=list(group)))

        return b_clusters

    @staticmethod
    def get_clusters(
        positioned_ions: List[PositionedIon], ion_type: Literal["b", "y"]
    ) -> List[Cluster]:
        """
        b-clusters are all those ions with the same protein ID AND start position
        y-clusters are all those ions with the same protein ID AND end position
        """
        if ion_type == B_ION_TYPE:
            grouping_attr = "inclusive_start"
        elif ion_type == Y_ION_TYPE:
            grouping_attr = "exclusive_end"
        else:
            raise ValueError(f"Unknown ion type: {ion_type}")

        ions = list(filter(lambda ion: ion.ion_type == ion_type, positioned_ions))
        ions.sort(key=lambda ion: (ion.protein, getattr(ion, grouping_attr)))
        clusters = [
            Cluster(ions=list(group))
            for _, group in groupby(
                ions, key=lambda ion: (ion.protein, getattr(ion, grouping_attr))
            )
        ]
        return clusters

    @classmethod
    def from_positioned_ions(
        cls, positioned_ions: List[PositionedIon]
    ) -> "SpectrumClusters":
        """
        This functions groups the product-ions into b- and y-clusters by (protein ID, start position)
        and (protein ID, end position), respectively.
        """
        return cls(
            b_clusters=cls.get_clusters(
                positioned_ions=positioned_ions, ion_type=B_ION_TYPE
            ),
            y_clusters=cls.get_clusters(
                positioned_ions=positioned_ions, ion_type=Y_ION_TYPE
            ),
        )


@dataclass
class SpectrumExtendedClusters:
    """
    Object to represent the extended b- and y-clusters. Extended clusters are clusters
    that have been extended until the precursor m/z is reached.
    """

    b_ext_clusters: List[ExtendedCluster]
    y_ext_clusters: List[ExtendedCluster]

    @staticmethod
    def get_seqs_from_clusters(
        ext_clusters: List[ExtendedCluster], ion_type: Literal["b", "y"]
    ):
        """
        Given a an extended cluster corresponding to sequence ABCDEFG (where ABCDEFG is a peptide
        within X PPM of the precursor m/z), the subsequences that a hybrid could form out of.
        If ABC is a b-extended cluster, then the subsequences are A, AB, ABC.
        If ABC is a y-extended cluster, then the subsequences are C, BC, ABC
        """
        if ion_type == B_ION_TYPE:
            ion_seq_creator = ProductIon.get_b_ion_seqs
        elif ion_type == Y_ION_TYPE:
            ion_seq_creator = ProductIon.get_y_ion_seqs
        else:
            raise ValueError(f"Unknown ion type: {ion_type}")
        seq_to_prot_id_map = defaultdict(set)
        for cluster in ext_clusters:
            prot_id = cluster.cluster.protein
            seq = cluster.extended_seq
            for ion_seq in ion_seq_creator(seq=seq):
                seq_to_prot_id_map[ion_seq].add(prot_id)

        return seq_to_prot_id_map

    def form_hybrids(
        self,
        precursor_charge: int,
        precursor_mz: float,
        precursor_mz_ppm_tol: float,
    ) -> List[HybridPeptide]:
        t0 = time()
        # Get the sequences to form hybrids out of from the extended clusters
        seq_to_prot_id_map = {}
        seq_to_prot_id_map[B_ION_TYPE] = self.get_seqs_from_clusters(
            ext_clusters=self.b_ext_clusters, ion_type=B_ION_TYPE
        )
        seq_to_prot_id_map[Y_ION_TYPE] = self.get_seqs_from_clusters(
            ext_clusters=self.y_ext_clusters, ion_type=Y_ION_TYPE
        )

        # Create database of y-sequences
        y_rows = [
            SeqWithMass.from_seq(seq=y_seq, charge=precursor_charge)
            for y_seq in seq_to_prot_id_map["y"]
        ]
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
        adjusted_precursor_mz = (
            precursor_mz + (WATER_MASS / precursor_charge) + PROTON_MASS
        )
        potential_hybrids = []
        for b_seq in seq_to_prot_id_map[B_ION_TYPE]:
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
                potential_hybrids.append(
                    HybridPeptide(
                        b_seq=b_seq,
                        y_seq=y_seq,
                        b_prot_ids=seq_to_prot_id_map["b"][b_seq],
                        y_prot_ids=seq_to_prot_id_map["y"][y_seq],
                    )
                )
        t1 = time()
        logger.debug(f"Creating hybrids took {get_time_in_diff_units(t1 - t0)}")
        return potential_hybrids

    @classmethod
    def from_spectrum(
        cls,
        spectrum: Spectrum,
        db_path: Path,
        proteins: Union[str, Path],
        peak_to_ion_ppm_tol: float,
        precursor_mz_ppm_tol: float,
        min_cluster_len: int = MIN_CLUSTER_LENGTH,
        min_cluster_support: int = MIN_CLUSTER_SUPPORT,
    ) -> "SpectrumExtendedClusters":
        """
        Workhorse function to get the extended b- and y-clusters corresponding to the given
        spectrum and database of product-ions
        -
        """
        # Load the database and a protein name-to-sequence map
        db = KmerDatabase(db_path=db_path)
        protein_name_to_seq_map = {
            pep.name: pep.seq for pep in Peptide.from_fasta(proteins)
        }
        # Get peak-ion matches for the spectrum
        peak_ion_matches = db.get_spectrum_peak_ion_matches(
            spectrum=spectrum,
            ppm_tolerance=peak_to_ion_ppm_tol,
        )
        # For each peak-ion match, find all its potential locations in the proteins
        positioned_ions = []
        for peak_ion_match in peak_ion_matches:
            for protein in peak_ion_match.ion.proteins:
                locations_in_protein = get_positions_of_subseq_in_seq(
                    subseq=peak_ion_match.seq,
                    seq=protein_name_to_seq_map[protein],
                )
                positioned_ions.extend(
                    [
                        PositionedIon(
                            seq=peak_ion_match.seq,
                            charge=peak_ion_match.ion.charge,
                            ion_type=peak_ion_match.ion.ion_type,
                            protein=protein,
                            inclusive_start=loc.inclusive_start,
                            exclusive_end=loc.exclusive_end,
                        )
                        for loc in locations_in_protein
                    ]
                )

        # Get clusters
        logger.debug("Getting clusters...")
        t0 = time()
        clusters = SpectrumClusters.from_positioned_ions(
            positioned_ions=positioned_ions
        )
        logger.debug(f"Getting clusters took {get_time_in_diff_units(time() - t0)}")
        msg = f"Before filtering, number b-clusters = {len(clusters.b_clusters)}; number y-clusters = {len(clusters.y_clusters)}"
        logger.debug(msg)

        # Filter the clusters
        clusters.filter_clusters(
            min_len=min_cluster_len, min_support=min_cluster_support
        )
        msg = f"After filtering, number b-clusters = {len(clusters.b_clusters)}; number y-clusters = {len(clusters.y_clusters)}"
        logger.debug(msg)

        # Extend clusters
        b_ext_clusters, y_ext_clusters = [], []
        for cluster in clusters.b_clusters + clusters.y_clusters:
            extended_cluster = ExtendedCluster.from_cluster(
                cluster=cluster,
                protein_seq=protein_name_to_seq_map[cluster.protein],
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
                precursor_charge=spectrum.precursor_charge,
                precursor_mz=spectrum.precursor_mz,
            )
            if cluster.ion_type == "b":
                b_ext_clusters.append(extended_cluster)
            elif cluster.ion_type == "y":
                y_ext_clusters.append(extended_cluster)

        return cls(
            b_ext_clusters=b_ext_clusters,
            y_ext_clusters=y_ext_clusters,
        )


def get_hybrids_via_clusters(
    db_path: Path,
    spectrum: Spectrum,
    peak_to_ion_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    precursor_mz_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
) -> List[HybridPeptide]:
    ext_clusters = SpectrumExtendedClusters.from_spectrum(
        spectrum=spectrum,
        db_path=db_path,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
    )

    if (len(ext_clusters.b_ext_clusters) == 0) or (
        len(ext_clusters.y_ext_clusters) == 0
    ):
        return []
    potential_hybrids = ext_clusters.form_hybrids(
        precursor_charge=spectrum.precursor_charge,
        precursor_mz=spectrum.precursor_mz,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
    )
    return potential_hybrids
