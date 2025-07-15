import logging
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path
from time import time
from typing import DefaultDict, Dict, List, Literal, Optional, Set
from venv import logger

from src.constants import DEFAULT_PPM_TOLERANCE, PROTON_MASS, WATER_MASS
from src.hypedsearch_utils import HybridPeptide, SeqWithMass
from src.mass_spectra import Spectrum
from src.peptide_spectrum_comparison import (
    get_start_and_end_positions_from_ions,
    ions_as_df,
)
from src.peptides_and_ions import BIonCreator, YIonCreator, compute_peptide_mz
from src.protein_product_ion_database import (
    PositionedIon,
    ProteinProductIonDb,
    get_aa_seq_from_db,
    get_positions_in_proteins_of_peak_matching_ions,
    get_product_ions_matching_spectrum,
)
from src.sql_database import Sqlite3Database
from src.utils import get_time_in_diff_units, relative_ppm_tolerance_in_daltons

logger = logging.getLogger(__name__)


@dataclass
class Cluster:
    """Object to represent a b- or y-cluster of product ions."""

    ions: List[PositionedIon]
    ion_type: str = field(init=False)
    protein_id: int = field(init=False)
    inclusive_start: str = field(init=False)
    exclusive_end: int = field(init=False)

    def __post_init__(self):
        # Get ion type
        ion_type = self.ions[0].ion_type
        assert all([ion.ion_type == ion_type for ion in self.ions]), (
            "All ions should have the same ion type!"
        )
        self.ion_type = ion_type

        # Get protein ID
        p_id = self.ions[0].protein_id
        assert all([ion.protein_id == p_id for ion in self.ions]), (
            "All ions should have same protein ID"
        )
        self.protein_id = p_id

        # Make sure that
        # - for b-ions, start positions should all be the same
        # - for y-ions, end positions should all be the same
        if self.ion_type == "b":
            start = self.ions[0].inclusive_start
            assert all([ion.inclusive_start == start for ion in self.ions]), (
                "All b-ions should have same start position"
            )
        else:
            end = self.ions[0].exclusive_end
            assert all([ion.exclusive_end == end for ion in self.ions]), (
                "All y-ions should have same end position"
            )

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


@dataclass
class ExtendedCluster:
    cluster: Cluster
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

    @classmethod
    def from_cluster(
        cls,
        cluster: Cluster,
        precursor_mz_ppm_tol: float,
        precursor_charge: int,
        precursor_mz: float,
        db: Optional[ProteinProductIonDb] = None,
        db_path: Optional[str] = None,
    ):
        """
        Extend a cluster until the precursor m/z is reached.
        """
        # Connect to the database if not already connected
        if db is None:
            db = ProteinProductIonDb(db_path=db_path)

        # A cluster has a start and end position and corresponds to some peptide.
        # Get the amino acid sequence of the peptide the cluster corresponds to, get the peptide's
        # m/z
        aa_seq = cluster.get_aa_seq(db=db)
        cluster_mz = compute_peptide_mz(aa_seq=aa_seq, charge=precursor_charge)
        protein_seq = db.get_protein_by_id(protein_id=cluster.protein_id).seq

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
    def filter_clusters(
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

    def filter_b_and_y_clusters(
        self,
        min_len: int = 3,
        min_support: int = 2,
    ):
        """
        Filter b- and y-clusters.
        """
        self.b_clusters = self.filter_clusters(
            clusters=self.b_clusters,
            min_len=min_len,
            min_support=min_support,
        )
        self.y_clusters = self.filter_clusters(
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
        b_ions.sort(key=lambda ion: (ion.protein_id, ion.inclusive_start))
        b_clusters = [
            Cluster(ions=list(group))
            for (p_id, loc), group in groupby(
                b_ions, key=lambda ion: (ion.protein_id, ion.inclusive_start)
            )
        ]
        return b_clusters

    @staticmethod
    def get_y_clusters(positioned_ions: List[PositionedIon]) -> List[Cluster]:
        """
        y-clusters are all those ions with the same protein ID AND end position
        """
        y_ions = list(filter(lambda ion: ion.ion_type == "y", positioned_ions))
        y_ions.sort(key=lambda ion: (ion.protein_id, ion.exclusive_end))
        y_clusters = [
            Cluster(ions=list(group))
            for (p_id, loc), group in groupby(
                y_ions, key=lambda ion: (ion.protein_id, ion.exclusive_end)
            )
        ]
        return y_clusters

    @classmethod
    def from_positioned_ions(
        cls, positioned_ions: List[PositionedIon]
    ) -> "SpectrumClusters":
        """
        This functions groups the product-ions into b- and y-clusters by (protein ID, start position)
        and (protein ID, end position), respectively.
        """
        return cls(
            b_clusters=cls.get_b_clusters(positioned_ions=positioned_ions),
            y_clusters=cls.get_y_clusters(positioned_ions=positioned_ions),
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
        if ion_type == "b":
            ion_seq_creator = BIonCreator().generate_product_ion_seqs
        elif ion_type == "y":
            ion_seq_creator = YIonCreator().generate_product_ion_seqs
        else:
            raise ValueError(f"Invalid ion type: {ion_type}. Must be 'b' or 'y'.")

        seq_to_prot_id_map = defaultdict(set)
        for cluster in ext_clusters:
            prot_id = cluster.cluster.protein_id
            seq = cluster.extended_seq
            for ion_seq in ion_seq_creator(seq=seq):
                seq_to_prot_id_map[ion_seq].add(prot_id)

        return seq_to_prot_id_map

    def get_b_and_y_seqs_to_form_hybrids_out_of(
        self,
    ) -> Dict[str, DefaultDict[str, Set[int]]]:
        """ """
        seq_to_prot_id_map = {
            "b": defaultdict(set),
            "y": defaultdict(set),
        }

    @classmethod
    def from_spectrum_clusters(
        cls,
        clusters: SpectrumClusters,
        db: ProteinProductIonDb,
        precursor_mz: float,
        precursor_charge: int,
        precursor_mz_ppm_tol: float,
    ) -> "SpectrumExtendedClusters":
        # Extend the clusters
        t0 = time()
        # b-clusters
        b_extended_clusters = [
            ExtendedCluster.from_cluster(
                cluster=cluster,
                precursor_charge=precursor_charge,
                precursor_mz=precursor_mz,
                db=db,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            )
            for cluster in clusters.b_clusters
        ]

        # y-clusters
        y_extended_clusters = [
            ExtendedCluster.from_cluster(
                cluster=cluster,
                precursor_charge=precursor_charge,
                precursor_mz=precursor_mz,
                db=db,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            )
            for cluster in clusters.y_clusters
        ]
        logger.debug(f"Getting extended clusters took {round(time() - t0, 2)} seconds")

        return cls(
            b_ext_clusters=b_extended_clusters, y_ext_clusters=y_extended_clusters
        )

    def form_hybrids(
        self,
        precursor_charge: int,
        precursor_mz: float,
        precursor_mz_ppm_tol: float,
    ) -> List[HybridPeptide]:
        t0 = time()
        # Get the sequences to form hybrids out of from the extended clusters
        seq_to_prot_id_map = {}
        seq_to_prot_id_map["b"] = self.get_seqs_from_clusters(
            ext_clusters=self.b_ext_clusters, ion_type="b"
        )
        seq_to_prot_id_map["y"] = self.get_seqs_from_clusters(
            ext_clusters=self.y_ext_clusters, ion_type="y"
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
        for b_seq in seq_to_prot_id_map["b"]:
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
        peak_to_ion_ppm_tol: float,
        precursor_mz_ppm_tol: float,
    ):
        """
        Workhorse function to get the extended b- and y-clusters corresponding to the given
        spectrum and database of product-ions
        -
        """
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
        logger.debug("Getting clusters...")
        t0 = time()
        clusters = SpectrumClusters.from_positioned_ions(
            positioned_ions=positioned_ions
        )
        logger.debug(f"Getting clusters took {get_time_in_diff_units(time() - t0)}")
        msg = f"Before filtering, number b-clusters = {len(clusters.b_clusters)}; number y-clusters = {len(clusters.y_clusters)}"
        logger.debug(msg)

        # Filter the clusters
        clusters.filter_b_and_y_clusters()
        msg = f"After filtering, number b-clusters = {len(clusters.b_clusters)}; number y-clusters = {len(clusters.y_clusters)}"
        logger.debug(msg)

        # Extend clusters
        logger.debug("Extending clusters...")
        t0 = time()
        extended_clusters = cls.from_spectrum_clusters(
            clusters=clusters,
            db=db,
            precursor_mz=spectrum.precursor_mz,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        )
        logger.debug(f"Extending clusters took {get_time_in_diff_units(time() - t0)}")

        return extended_clusters


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
