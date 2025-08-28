import logging
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union
from venv import logger

import click
import numpy as np

from src.constants import (
    B_ION_TYPE,
    DEFAULT_MIN_CLUSTER_LENGTH,
    DEFAULT_MIN_CLUSTER_SUPPORT,
    DEFAULT_PEAK_TO_ION_PPM_TOL,
    DEFAULT_PRECURSOR_MZ_PPM_TOL,
    HS_PREFIX,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Spectrum
from src.peptides_and_ions import (
    Fasta,
    Peptide,
    UnpositionedProductIon,
    compute_peptide_precursor_mz,
)
from src.sql_database import Sqlite3Database, SqlTableRow
from src.utils import (
    Position,
    get_positions_of_subseq_in_seq,
    load_json,
    log_time,
    relative_ppm_tolerance_in_daltons,
    setup_logger,
    to_json,
)

logger = logging.getLogger(__name__)


@dataclass
class SeqWithMass(SqlTableRow):
    """
    Dataclass to represent a sequence treated as a product ion for when we search via an
    SQL database (which is why this dataclass is a subclass of SqlTableRow)
    for hybrids, <b-seq>-<y-seq>, that are within X PPM of the precursor m/z.
    """

    seq: str
    mz: float

    @classmethod
    def from_seq_and_charge(cls, seq: str, charge: int):
        mz = compute_peptide_precursor_mz(seq=seq, charge=charge)
        return cls(seq=seq, mz=mz)


@dataclass
class PositionedProductIon:
    seq: str
    charge: int
    ion_type: str
    protein: Union[int, str]
    inclusive_start: int
    exclusive_end: int

    @classmethod
    def from_unpositioned_product_ion(
        cls,
        unpositioned_product_ion: UnpositionedProductIon,
        protein_name_to_seq_map: Dict[str, str],
    ) -> List["PositionedProductIon"]:
        """
        Create a PositionedProductIon from an UnpositionedProductIon.
        This function will find all the positions of the product ion's sequence in the protein
        sequences and create a PositionedProductIon for each position.
        """
        positioned_ions = []
        for protein in unpositioned_product_ion.proteins:
            locations_in_protein = get_positions_of_subseq_in_seq(
                subseq=unpositioned_product_ion.seq,
                seq=protein_name_to_seq_map[protein],
            )
            positioned_ions.extend(
                [
                    PositionedProductIon(
                        seq=unpositioned_product_ion.seq,
                        charge=unpositioned_product_ion.charge,
                        ion_type=unpositioned_product_ion.ion_type,
                        protein=protein,
                        inclusive_start=loc.inclusive_start,
                        exclusive_end=loc.exclusive_end,
                    )
                    for loc in locations_in_protein
                ]
            )
        return positioned_ions


@dataclass
class Cluster:
    """Object to represent a b- or y-cluster of product ions."""

    ions: List[PositionedProductIon]
    ion_type: Literal[B_ION_TYPE, Y_ION_TYPE] = field(init=False)
    protein: Union[str, int] = field(init=False)
    inclusive_start: int = field(init=False)
    exclusive_end: int = field(init=False)
    extended_seq: str = field(default=None, init=False, repr=False)

    def __post_init__(self):
        # Get ion type
        ion_type = self.ions[0].ion_type
        assert all(
            [ion.ion_type == ion_type for ion in self.ions]
        ), "All ions should have the same ion type!"
        self.ion_type = ion_type

        # Get protein identifier
        protein = self.ions[0].protein
        assert all(
            [ion.protein == protein for ion in self.ions]
        ), "All ions should have same protein"
        self.protein = protein

        # Make sure that
        # - for b-ions, start positions should all be the same
        # - for y-ions, end positions should all be the same
        if self.ion_type == B_ION_TYPE:
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

    @property
    def support(self):
        return len(self.ions)

    @property
    def length(self):
        return self.exclusive_end - self.inclusive_start

    def set_extended_seq(
        self, protein_seq, precursor_mz_ppm_tol, precursor_charge, precursor_mz
    ):
        """
        Extend a cluster until the precursor m/z is reached.
        Say ABCD has m/z 99 and ABCDE has m/z 105 and the precursor's m/z + the PPM tolerance in
        daltons is 100, then we want to return ABCD--the sequence that's closest but less than
        the precursor m/z + PPM tolerance in daltons.
        """
        # A cluster has a start and end position and corresponds to some peptide.
        # Get the amino acid sequence of the peptide the cluster corresponds to, get the peptide's
        # m/z
        aa_seq = protein_seq[self.inclusive_start : self.exclusive_end]
        cluster_mz = compute_peptide_precursor_mz(seq=aa_seq, charge=precursor_charge)

        extended_cluster_mz = cluster_mz
        num_extended_aa = 0
        seqs = []
        da_tol = relative_ppm_tolerance_in_daltons(
            ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
        )
        while extended_cluster_mz < precursor_mz + da_tol:
            # Extend the cluster by one amino acid
            seqs.append(aa_seq)
            num_extended_aa += 1
            if self.ion_type == B_ION_TYPE:
                start = self.inclusive_start
                end = self.exclusive_end + num_extended_aa
            elif self.ion_type == Y_ION_TYPE:
                start = self.inclusive_start - num_extended_aa
                end = self.exclusive_end
            else:
                raise ValueError(f"Unknown ion type: {self.ion_type}")

            # Exit loop if we've reached the end of the protein
            if (start < 0) or (end > len(protein_seq)):
                break

            aa_seq = protein_seq[start:end]
            extended_cluster_mz = compute_peptide_precursor_mz(
                seq=aa_seq,
                charge=precursor_charge,
            )

        # Grab the last sequence of the extended cluster
        if len(seqs) == 0:
            extended_seq = None
        else:
            extended_seq = seqs[-1]

        self.extended_seq = extended_seq

    def get_seqs_for_hybrids(self) -> List[str]:
        """
        These are the sequences that we'll use for forming hybrids.
        For b-ions,
        """
        shortest_seq = self.ions[np.argmin([len(ion.seq) for ion in self.ions])].seq
        if self.ion_type == B_ION_TYPE:
            # For b-ions, we take all prefixes of the extended sequence starting from the
            # sequence of the shortest supported ion
            return [
                self.extended_seq[:idx]
                for idx in range(len(shortest_seq), len(self.extended_seq) + 1)
            ]
        elif self.ion_type == Y_ION_TYPE:
            # For y-ions, we take all suffixes of the extended sequence starting from the
            # sequence of the shortest supported ion
            return [
                self.extended_seq[idx:]
                for idx in range(0, len(self.extended_seq) - len(shortest_seq) + 1)
            ]
        else:
            raise ValueError(f"Unknown ion type: {self.ion_type}")


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

    @log_time(level=logging.INFO)
    def filter_clusters(
        self,
        min_len: int = DEFAULT_MIN_CLUSTER_LENGTH,
        min_support: int = DEFAULT_MIN_CLUSTER_SUPPORT,
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
    def _get_clusters(
        positioned_ions: List[PositionedProductIon],
        ion_type: Literal[B_ION_TYPE, Y_ION_TYPE],
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
    @log_time(level=logging.INFO)
    def from_positioned_ions(
        cls, positioned_ions: List[PositionedProductIon]
    ) -> "SpectrumClusters":
        """
        This functions groups the product-ions into b- and y-clusters by (protein ID, start position)
        and (protein ID, end position), respectively.
        """
        return cls(
            b_clusters=cls._get_clusters(
                positioned_ions=positioned_ions, ion_type=B_ION_TYPE
            ),
            y_clusters=cls._get_clusters(
                positioned_ions=positioned_ions, ion_type=Y_ION_TYPE
            ),
        )


def form_extended_clusters_for_spectrum(
    kmer_db: KmerDatabase,
    spectrum: Spectrum,
    protein_name_to_seq_map: Dict[str, str],
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH,
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
) -> SpectrumClusters:
    # Get peak-ion matches for the spectrum
    peak_ion_matches = kmer_db.get_peak_ion_matches_for_spectrum(
        spectrum=spectrum,
        ppm_tolerance=peak_to_ion_ppm_tol,
    )
    # For each peak-ion match, find all its potential locations in the proteins
    positioned_ions = []
    for peak_ion_mach in peak_ion_matches:
        positioned_ions.extend(
            PositionedProductIon.from_unpositioned_product_ion(
                unpositioned_product_ion=peak_ion_mach.ion,
                protein_name_to_seq_map=protein_name_to_seq_map,
            )
        )

    # Get clusters for spectrum
    clusters = SpectrumClusters.from_positioned_ions(positioned_ions=positioned_ions)
    logger.info(
        f"Before filtering, there are {len(clusters.b_clusters)} b-clusters and {len(clusters.y_clusters)} y-clusters."
    )
    clusters.filter_clusters(min_len=min_cluster_len, min_support=min_cluster_support)
    logger.info(
        f"After filtering, there are {len(clusters.b_clusters)} b-clusters and {len(clusters.y_clusters)} y-clusters."
    )

    # Extend the sequence
    for cluster in clusters.b_clusters + clusters.y_clusters:
        cluster.set_extended_seq(
            protein_seq=protein_name_to_seq_map[cluster.protein],
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz=spectrum.precursor_mz,
        )

    return clusters


@dataclass
class HybridPeptide:
    left_seq: str
    right_seq: str
    left_proteins: Set[str]
    right_proteins: Set[str]
    fasta_description: str = field(init=False)
    scan: Optional[int] = None
    sample: Optional[str] = None

    def __post_init__(self):
        self.fasta_description = f"left-prots:{','.join(self.left_proteins)};right-prots:{','.join(self.right_proteins)}"

    def serialize(self):
        return {
            "left_seq": self.left_seq,
            "right_seq": self.right_seq,
            "left_proteins": list(self.left_proteins),
            "right_proteins": list(self.right_proteins),
            "scan": self.scan,
            "sample": self.sample,
        }

    @property
    def seq(self):
        return self.left_seq + self.right_seq

    @property
    def seq_with_hyphen(self):
        return f"{self.left_seq}-{self.right_seq}"

    @property
    def fasta_name(self):
        return f"{HS_PREFIX}{self.left_seq}-{self.right_seq}"

    @classmethod
    def from_name(cls, name: str):
        b_seq, y_seq = cls.parse_name_to_b_and_y_seqs(name=name)
        return cls(
            b_seq=b_seq,
            y_seq=y_seq,
        )

    @staticmethod
    def parse_name_to_b_and_y_seqs(name: str):
        """
        Parse the name of a hybrid peptide to extract the b and y sequences.
        Depends on how .fasta_name works
        """
        if name.startswith(HS_PREFIX):
            name = name[len(HS_PREFIX) :]
        # For back compatibility
        elif name.startswith("hybrid_"):
            name = name[len("hybrid_") :]
        else:
            raise ValueError(
                f"Invalid hybrid peptide name: {name}. Expected prefix '{HS_PREFIX}' or 'hybrid_'."
            )

        b_seq, y_seq = name.split("-")
        return b_seq, y_seq

    def set_protein_names(self, prot_id_to_name_map: Dict[int, str]):
        self.left_proteins = set(
            prot_id_to_name_map[prot_id] for prot_id in self.b_prot_ids
        )
        self.right_proteins = set(
            prot_id_to_name_map[prot_id] for prot_id in self.y_prot_ids
        )

    def set_fasta_info(
        self,
        prot_id_to_name_map: Optional[Dict[int, str]],
    ):
        self.set_protein_names(prot_id_to_name_map=prot_id_to_name_map)
        self.fasta_description = f"b-prots:{','.join(self.left_proteins)} y-prots:{','.join(self.right_proteins)}"

    def mz(self, charge: int):
        return compute_peptide_precursor_mz(
            seq=self.left_seq + self.right_seq, charge=charge
        )


@log_time(level=logging.INFO)
def form_hybrids_from_clusters(
    b_clusters: List[Cluster],
    y_clusters: List[Cluster],
    precursor_charge: int,
    precursor_mz: float,
    precursor_mz_ppm_tol: float,
    scan: int,
    sample: str,
) -> List[HybridPeptide]:
    """
    Given b- and y-clusters, form hybrids.
    """
    seq_protein_map = {
        B_ION_TYPE: defaultdict(set),
        Y_ION_TYPE: defaultdict(set),
    }
    for cluster in b_clusters:
        for seq in cluster.get_seqs_for_hybrids():
            seq_protein_map[B_ION_TYPE][seq].add(cluster.protein)
    for cluster in y_clusters:
        for seq in cluster.get_seqs_for_hybrids():
            seq_protein_map[Y_ION_TYPE][seq].add(cluster.protein)
    left_seqs = list(seq_protein_map[B_ION_TYPE].keys())
    right_seqs = list(seq_protein_map[Y_ION_TYPE].keys())
    hybrid_seqs = form_hybrids_from_left_and_right_seqs(
        left_seqs=left_seqs,
        right_seqs=right_seqs,
        precursor_charge=precursor_charge,
        precursor_mz=precursor_mz,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
    )
    hybrids = []
    for hybrid_seq in hybrid_seqs:
        left_seq, right_seq = hybrid_seq.split("-")
        hybrids.append(
            HybridPeptide(
                left_seq=left_seq,
                right_seq=right_seq,
                left_proteins=seq_protein_map[B_ION_TYPE].get(left_seq),
                right_proteins=seq_protein_map[Y_ION_TYPE].get(right_seq),
                scan=scan,
                sample=sample,
            )
        )
    return hybrids


@log_time(level=logging.INFO)
def form_hybrids_from_left_and_right_seqs(
    left_seqs: List[str],
    right_seqs: List[str],
    precursor_charge: int,
    precursor_mz: float,
    precursor_mz_ppm_tol: float,
) -> List[str]:
    """
    Consider all possible hybrids of form <left seq>-<right seq>.
    To do so, we
    - create a database of right sequences
    - for each left sequence, search the right sequence database to find the hybrids
      that would produce a peptide within the given PPM tolerance of the precursor m/z
    """
    # Constants
    table_name = "right_seqs"
    index_name = "precursor_mz"

    # Create database of y-sequences
    right_seq_rows = [
        SeqWithMass.from_seq_and_charge(seq=right_seq, charge=precursor_charge)
        for right_seq in right_seqs
    ]
    db = Sqlite3Database()
    db.create_table_from_dataclass(table_name=table_name, obj=SeqWithMass)
    db.insert_dataclasses(table_name=table_name, data_classes=right_seq_rows)
    db.add_index(table_name=table_name, index_name=index_name, colms_to_index=["mz"])

    # For each left-sequence search the right-sequence database to find the hybrids
    # that would produce a peptide within the given PPM tolerance of the precursor m/z
    da_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
    )
    adjusted_precursor_mz = precursor_mz + (WATER_MASS / precursor_charge) + PROTON_MASS
    hybrids = []
    for left_seq in left_seqs:
        left_seq_mz = compute_peptide_precursor_mz(
            seq=left_seq, charge=precursor_charge
        )
        lower_bdd = adjusted_precursor_mz - da_tolerance - left_seq_mz
        upper_bdd = adjusted_precursor_mz + da_tolerance - left_seq_mz
        query = f"""
            SELECT
                *
            FROM {table_name} as ion
            WHERE ion.mz BETWEEN {lower_bdd} AND {upper_bdd}
        """
        hybrids.extend(
            [f"{left_seq}-{match['seq']}" for match in db.read_query(query=query)]
        )
    return hybrids


@log_time(level=logging.INFO)
def form_spectrum_hybrids_via_clustering(
    spectrum: Spectrum,
    kmer_db: KmerDatabase,
    protein_name_to_seq_map: Dict[str, str],
    kmer_to_proteins_map: Dict[str, List[str]],
    precursor_mz_ppm_tol: float = DEFAULT_PRECURSOR_MZ_PPM_TOL,
    peak_to_ion_ppm_tol: float = DEFAULT_PEAK_TO_ION_PPM_TOL,
    min_cluster_len: int = DEFAULT_MIN_CLUSTER_LENGTH,
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT,
) -> Dict[str, List[HybridPeptide]]:
    """
    This function will
    1. form hybrids for the given spectrum via clustering
    2. remove native hybrids
    3. return a dictionary mapping hybrid sequences to lists of HybridPeptide objects
    """
    # Form clusters
    clusters = form_extended_clusters_for_spectrum(
        kmer_db=kmer_db,
        spectrum=spectrum,
        protein_name_to_seq_map=protein_name_to_seq_map,
        peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        min_cluster_len=min_cluster_len,
        min_cluster_support=min_cluster_support,
    )
    hybrids = form_hybrids_from_clusters(
        b_clusters=clusters.b_clusters,
        y_clusters=clusters.y_clusters,
        precursor_charge=spectrum.precursor_charge,
        precursor_mz=spectrum.precursor_mz,
        precursor_mz_ppm_tol=precursor_mz_ppm_tol,
        scan=spectrum.scan,
        sample=spectrum.sample,
    )
    # Remove hybrids that are native sequences and
    # group hybrids by sequence (e.g., group A-BC with AB-C)
    logger.info(
        "Removing hybrids that correspond to native sequences. Then grouping the non-native "
        "hybrid peptides by sequence"
    )
    seq_to_hybrids = defaultdict(list)
    for hybrid in hybrids:
        if hybrid.seq in kmer_to_proteins_map:
            continue
        seq_to_hybrids[hybrid.seq].append(hybrid)

    return seq_to_hybrids


def serialize_hybrids(seq_to_hybrids: Dict[str, List[HybridPeptide]]) -> Dict:
    return {
        seq: [hybrid.serialize() for hybrid in hybrids]
        for seq, hybrids in seq_to_hybrids.items()
    }


@click.command(
    name="form-hybrids",
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
    help=(
        "Form hybrids for the given spectra. If no `--scan` (`-s`) is given, then hybrids "
        "will be formed for all scans in the MZML file. The hybrids will be saved as "
        "'<scan>.json' in the specified output directory."
    ),
)
@click.option(
    "--mzml",
    "-m",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the MZML file.",
)
@click.option(
    "--scan",
    "-s",
    type=int,
    required=False,
    help="Scan number of the spectrum to form hybrids for. If not provided, all scans will be processed.",
)
@click.option(
    "--database",
    "-d",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the kmer database",
)
@click.option(
    "--precursor_mz_ppm_tol",
    "-pmpt",
    type=float,
    default=DEFAULT_PRECURSOR_MZ_PPM_TOL,
    show_default=True,
    help="Precursor m/z PPM tolerance. Hybrids will be within this PPM of the precursor m/z.",
)
@click.option(
    "--peak_to_ion_ppm_tol",
    "-pipt",
    type=float,
    default=DEFAULT_PEAK_TO_ION_PPM_TOL,
    show_default=True,
    help="The PPM tolerance within which spectra peaks and matched fragment ions must be",
)
@click.option(
    "--fasta",
    "-f",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to the FASTA file",
)
@click.option(
    "--out_dir",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Where hybrid .json files will be saved as '<scan>.json'",
)
@log_time(level=logging.INFO)
def cli_form_hybrids(
    mzml: Path,
    database: Path,
    fasta: Path,
    out_dir: Path,
    precursor_mz_ppm_tol: float,
    peak_to_ion_ppm_tol: float,
    scan: Optional[int],
):
    if scan is None:
        logger.info(f"Processing all scans in {mzml}")
        for spectrum in Spectrum.parse_ms2_from_mzml(mzml):
            logger.info(f"Processing scan {spectrum.scan} in {mzml}")
            seq_to_hybrids = form_spectrum_hybrids_via_clustering(
                kmer_db=database,
                fasta=fasta,
                precursor_mz_ppm_tol=precursor_mz_ppm_tol,
                peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
                spectrum=spectrum,
            )
            # Save the hybrids to a JSON file
            to_json(
                data=serialize_hybrids(seq_to_hybrids=seq_to_hybrids),
                path=out_dir / f"{spectrum.scan}.json",
            )
    else:
        logger.info(f"Processing scan {scan} in {mzml}")
        spectrum = Spectrum.get_spectrum(scan=scan, mzml=mzml)
        seq_to_hybrids = form_spectrum_hybrids_via_clustering(
            kmer_db=database,
            fasta=fasta,
            precursor_mz_ppm_tol=precursor_mz_ppm_tol,
            peak_to_ion_ppm_tol=peak_to_ion_ppm_tol,
            spectrum=spectrum,
        )
        # Save the hybrids to a JSON file
        to_json(
            data=serialize_hybrids(seq_to_hybrids=seq_to_hybrids),
            path=out_dir / f"{spectrum.scan}.json",
        )


if __name__ == "__main__":
    logger = setup_logger()
    cli_form_hybrids()
