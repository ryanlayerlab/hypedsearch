import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union

import click

from src.constants import (
    AMINO_ACID_MASSES,
    B_ION_TYPE,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    MEMORY,
    PRODUCT_ION_TABLE,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
    IonTypes,
)
from src.mass_spectra import Peak, Spectrum
from src.peptides_and_ions import (
    Peptide,
    ProductIon,
    compute_peptide_mz,
    get_proteins_by_name,
)
from src.sql_database import Sqlite3Database, SqlTableRow
from src.utils import (
    PathType,
    decompress_and_depickle,
    load_json,
    log_time,
    pickle_and_compress,
    relative_ppm_tolerance_in_daltons,
    setup_logger,
    to_json,
)

logger = logging.getLogger(__name__)


@dataclass
class KmerToProteinMap:
    kmer_to_protein_map: Dict[str, List[Union[int, str]]]

    @classmethod
    def create(
        cls,
        fasta: Optional[Path] = None,
        proteins: Optional[List[Peptide]] = None,
        min_k: int = DEFAULT_MIN_KMER_LEN,
        max_k: int = DEFAULT_MAX_KMER_LEN,
        protein_attr: Literal["id", "name"] = "id",
        protein_names: Optional[Union[List[str], Path]] = None,
    ) -> "KmerToProteinMap":
        """
        Create a kmer to protein map from a FASTA file or a list of Peptides
        """
        # Load proteins
        logger.info("Loading proteins from FASTA file...")
        if proteins is None:
            proteins = Peptide.from_fasta(fasta_path=fasta)
        if protein_names is not None:
            proteins = get_proteins_by_name(
                proteins=proteins, protein_names=protein_names
            )
        # Create kmer-to-protein-id map
        logger.info("Creating kmer-to-protein map...")
        return cls(
            kmer_to_protein_map=cls.get_uniq_kmer_to_protein_map(
                proteins=proteins, min_k=min_k, max_k=max_k, protein_attr=protein_attr
            )
        )

    # Don't use @log_params because 'proteins' can be a list of many, many peptides
    @staticmethod
    @log_time(level=logging.INFO)
    def get_uniq_kmer_to_protein_map(
        proteins: List[Peptide],
        min_k: int = DEFAULT_MIN_KMER_LEN,
        max_k: int = DEFAULT_MAX_KMER_LEN,
        protein_attr: str = "id",
        verbose: bool = True,
    ) -> Dict[str, List[Union[int, str]]]:
        """ """
        uniq_kmer_to_protein_map = defaultdict(list)
        num_proteins = len(proteins)
        for p_idx, protein in enumerate(proteins):
            if verbose and (p_idx % 10 == 0):
                logger.info(f"Processing protein {p_idx + 1} of {num_proteins}")
            uniq_kmers = set(
                kmer.seq for kmer in protein.kmers(min_k=min_k, max_k=max_k)
            )
            for kmer in uniq_kmers:
                uniq_kmer_to_protein_map[kmer].append(getattr(protein, protein_attr))
        logger.info(f"Number of unique kmers {len(uniq_kmer_to_protein_map)}")
        return dict(uniq_kmer_to_protein_map)

    def save(self, out_path: Union[str, Path]) -> None:
        out_path = Path(out_path)
        if out_path.suffix == ".pklz":
            pickle_and_compress(obj=self.kmer_to_protein_map, file_path=out_path)
        elif out_path.suffix == ".json":
            to_json(data=self.kmer_to_protein_map, out_path=out_path)
        else:
            raise ValueError("Output path must be a .pklz or .json file.")

    @staticmethod
    def load(path: Union[str, Path]) -> "KmerToProteinMap":
        path = Path(path)
        if path.suffix == ".pklz":
            kmer_to_protein_map = decompress_and_depickle(path)
        elif path.suffix == ".json":
            kmer_to_protein_map = load_json(path)
        else:
            raise ValueError("Input path must be a .pklz or .json file.")
        return KmerToProteinMap(kmer_to_protein_map=kmer_to_protein_map)


@dataclass
class DbKmer(SqlTableRow):
    seq: str
    aa_mass: float
    proteins: str

    @classmethod
    def from_seq_and_proteins(
        cls,
        seq: str,
        proteins: List[Union[int, str]],
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ) -> "DbKmer":
        aa_mass = sum([amino_acid_mass_lookup[aa] for aa in seq])
        return cls(seq=seq, aa_mass=aa_mass, proteins=",".join(map(str, proteins)))

    def set_charge(self, charge: int):
        self._charge = charge

    @property
    def charge(self):
        return self._charge

    @property
    def z(self):
        return self._charge

    @property
    def ion_type(self):
        return self._ion_type

    @property
    def ion(self):
        return self._ion_type

    @property
    def protein_ids_as_ints(self):
        return [int(p_id) for p_id in self.proteins.split(",")]

    def ionize(self, charge: int, ion_type: Literal["b", "y"]):
        self._ion_type = ion_type
        self._charge = charge
        return self

    def __str__(self):
        return f"{self.ion}-{self.z}-{self.seq}"


@dataclass
class PeakIonMatch:
    ion: ProductIon
    peak: Peak

    @property
    def seq(self):
        return self.ion.seq

    @property
    def charge(self):
        return self.ion.charge


@dataclass
class KmerDatabase:
    db_path: Path
    db: Sqlite3Database = field(init=False)
    table_name = "kmers"
    index_name = "mass"
    row_object = DbKmer
    index_colm = "aa_mass"

    def __post_init__(self):
        self.db_path = Path(self.db_path)
        if self.db_path.exists():
            self.db = Sqlite3Database(db_path=self.db_path)
        else:
            raise RuntimeError(
                "Database path does not exist. Please create the database first."
            )

    @classmethod
    def create_db(
        cls,
        db_path: Path,
        kmer_to_protein_map: Union[Path, KmerToProteinMap, str],
        overwrite: bool = True,
    ) -> "KmerDatabase":
        if isinstance(kmer_to_protein_map, (Path, str)):
            kmer_to_protein_map = KmerToProteinMap.load(kmer_to_protein_map)

        db = Sqlite3Database(db_path=db_path, overwrite=overwrite)
        db.create_table_from_dataclass(table_name=cls.table_name, obj=cls.row_object)
        for seq, proteins in kmer_to_protein_map.kmer_to_protein_map.items():
            kmer = DbKmer.from_seq_and_proteins(seq=seq, proteins=proteins)
            db.insert_dataclasses(table_name=cls.table_name, data_classes=[kmer])
        db.add_index(
            table_name=cls.table_name,
            index_name=cls.index_name,
            colms_to_index=[cls.index_colm],
        )
        return cls(
            db_path=db_path,
        )

    def get_matching_product_ions(
        self,
        query_mz: float,
        charge: int,
        ppm_tolerance: float,
        ion_type: Literal["b", "y"],
    ) -> List[ProductIon]:
        # Constant--only used in this function
        bounds_fcns = {
            "b": b_ion_bounds_on_amino_acid_mass,
            "y": y_ion_bounds_on_amino_acid_mass,
        }

        # Get the product ions that match the given m/z
        bounds_fcn = bounds_fcns[ion_type]
        lower_bound, upper_bound = bounds_fcn(
            ppm_tolerance=ppm_tolerance,
            ref_mass=query_mz,
            charge=charge,
        )
        query = f"""
        SELECT
            *
        FROM {self.table_name} as kmer
        WHERE kmer.aa_mass BETWEEN {lower_bound} AND {upper_bound}
        """
        product_ions = []
        for kmer in self.db.read_query(query=query):
            kmer = DbKmer(**kmer)
            proteins = [int(p) if p.isdigit() else p for p in kmer.proteins.split(",")]
            product_ions.append(
                ProductIon(
                    seq=kmer.seq,
                    charge=charge,
                    ion_type=ion_type,
                    proteins=proteins,
                )
            )
        return product_ions

    def get_matching_product_ions_to_spectrum_peak(
        self,
        peak: Peak,
        precursor_charge: int,
        precursor_mz: float,
        ppm_tolerance: float,
        ion_types: Set[Literal["b", "y"]] = {"b", "y"},
    ) -> List[PeakIonMatch]:
        peak_ion_matches = []
        for charge in range(1, precursor_charge + 1):
            for ion_type in ion_types:

                charge_ion_matches = self.get_matching_product_ions(
                    query_mz=peak.mz,
                    charge=charge,
                    ppm_tolerance=ppm_tolerance,
                    ion_type=ion_type,
                )
                # Remove the ions whose m/z at the given precursor charge would be greater
                # than the precursor's m/z
                tmp = list(
                    filter(
                        lambda ion: compute_peptide_mz(
                            aa_seq=ion.seq, charge=precursor_charge
                        )
                        <= precursor_mz,
                        charge_ion_matches,
                    ),
                )

                peak_ion_matches.extend(
                    [PeakIonMatch(ion=ion_match, peak=peak) for ion_match in tmp]
                )
        return peak_ion_matches

    def get_spectrum_peak_ion_matches(
        self,
        spectrum: Spectrum,
        ppm_tolerance: float,
        ion_types: Set[Literal["b", "y"]] = {"b", "y"},
    ) -> List[PeakIonMatch]:
        peak_ion_matches = []
        for peak in spectrum.peaks:
            peak_ion_matches.extend(
                self.get_matching_product_ions_to_spectrum_peak(
                    peak=peak,
                    precursor_charge=spectrum.precursor_charge,
                    precursor_mz=spectrum.precursor_mz,
                    ppm_tolerance=ppm_tolerance,
                    ion_types=ion_types,
                )
            )
        return peak_ion_matches


def b_ion_bounds_on_amino_acid_mass(ppm_tolerance, ref_mass, charge):
    """
    Let A = the sum of the masses of the amino acids in the peptide, Bp = the mass of a proton,
    m_{zb} = the mass of a b-ion with charge z is (A/z) + Bp, and
    Da(M, X) = X PPM of M in daltons.
    Then you can show that if you want all b-ions for which m_{zb} is within within X ppm of a
    given mass M (called the "reference mass"), then that means that you want
    (M - Da(M, X)) <= m_{zb} = (A/z) + Bp <= (M + Da(M, X)) which means
    z*((M-Da(M, X)) - Bp) <= A <= z*((M + Da(M, X)) - Bp)

    If you want the b-ion m/z, m_{zb}, to be within X ppm of a given mass M (called the "reference mass")
    e1 <= m_{zb} <= e2, then you want want
    """
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=ppm_tolerance, ref_mass=ref_mass
    )

    lower_mz_bound = ref_mass - mz_tolerance
    upper_mz_bound = ref_mass + mz_tolerance

    lower_bound = charge * (lower_mz_bound - PROTON_MASS)
    upper_bound = charge * (upper_mz_bound - PROTON_MASS)

    return lower_bound, upper_bound


def y_ion_bounds_on_amino_acid_mass(ppm_tolerance, ref_mass, charge):
    """
    Let A = the sum of the masses of the amino acids in the peptide, Bp = the mass of a proton,
    Bw = the mass of water,
    m_{zy} = the mass of a y-ion with charge z is ((A+Bw)/z) + Bp, and
    Da(M, X) = X PPM of M in daltons.
    Then you can show that if you want all y-ions for which m_{zy} is within within X ppm of a
    given mass M (called the "reference mass"), then that means that you want
    (M - Da(M, X)) <= m_{zy} <= (M + Da(M, X)) which means
    z*((M-Da(M, X)) - Bp) - Bw <= A <= z*((M + Da(M, X)) - Bp) + Bw
    """
    b_lower, b_upper = b_ion_bounds_on_amino_acid_mass(
        ppm_tolerance=ppm_tolerance, ref_mass=ref_mass, charge=charge
    )
    y_lower, y_upper = b_lower - WATER_MASS, b_upper - WATER_MASS
    return y_lower, y_upper


# @log_time(level=logging.INFO)
# def create_kmer_database(
#     kmer_to_protein_map: Union[Dict, Path],
#     min_k: int = DEFAULT_MIN_KMER_LEN,
#     max_k: int = DEFAULT_MAX_KMER_LEN,
# ):
#     logger.info("Getting proteins whose kmers will be stored in the database...")
#     db_proteins = get_proteins_for_db(proteins=proteins, fasta_path=fasta_path)

#     db = Sqlite3Database(db_path=db_path)
#     create_protein_table(db=db, db_proteins=db_proteins)

#     # Create product/fragment ion table
#     create_product_ion_table(
#         db=db,
#         db_proteins=db_proteins,
#         min_k=min_k,
#         max_k=max_k,
#     )


# def create_product_ion_table(
#     db: ProteinProductIonDb,
#     db_proteins: List[Peptide],
#     min_k: int = DEFAULT_MIN_KMER_LEN,
#     max_k: int = DEFAULT_MAX_KMER_LEN,
# ) -> None:
#     uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(
#         proteins=db_proteins, min_k=min_k, max_k=max_k
#     )
#     logger.info("Creating 'DbKmer' objects for the database")
#     ions = [
#         DbKmer.seq_to_ion(seq=seq, protein_ids=protein_ids)
#         for seq, protein_ids in uniq_kmer_to_protein_map.items()
#     ]
#     logger.info(
#         "Creating the table for the 'DbKmer' object and inserting the object instances"
#     )
#     db.create_table_from_dataclass(table_name=db.product_ion_table_name, obj=DbKmer)
#     db.insert_dataclasses(table_name=db.product_ion_table_name, data_classes=ions)
#     logger.info("Creating 'mass' index for the table")
#     db.add_index(
#         table_name=db.product_ion_table_name,
#         index_name="mass",
#         colms_to_index=["aa_mass"],
#     )


def create_db(
    kmer_to_protein_path: Path,
    fasta: Path = None,
    # proteins: Optional[Union[List[str], Path]] = None,
    db_path: Union[Path, str] = MEMORY,
    min_k: int = DEFAULT_MIN_KMER_LEN,
    max_k: int = DEFAULT_MAX_KMER_LEN,
) -> KmerDatabase:
    if not kmer_to_protein_path.exists():
        logger.info(
            f"kmer-to-protein map doesn't exist at {kmer_to_protein_path}. "
            "Creating it now..."
        )
        kmer_to_protein_map = KmerToProteinMap.create(
            fasta=fasta,
            min_k=min_k,
            max_k=max_k,
            protein_attr="name",
        )
        kmer_to_protein_path.parent.mkdir(parents=True, exist_ok=True)
        kmer_to_protein_map.save(kmer_to_protein_path)
    if not db_path.exists():
        logger.info(f"Database does not exist at {db_path}. Creating it now...")
        KmerDatabase.create_db(
            db_path=db_path,
            kmer_to_protein_map=kmer_to_protein_map,
        )


@click.command(
    name="create-db",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="Create a database of product ions",
)
@click.option(
    "--db_path",
    "-d",
    type=PathType(),
    required=True,
    help="Path to the database.",
)
@click.option(
    "--fasta",
    "-f",
    type=PathType(),
    required=True,
    help="Path to the FASTA file.",
)
@click.option(
    "--kmer_to_protein",
    "-ktp",
    type=PathType(),
    required=True,
    help=("Path to the kmer-to-protein map."),
)
@click.option(
    "--min_k",
    "-mk",
    type=int,
    default=DEFAULT_MIN_KMER_LEN,
    show_default=True,
    help="Minimum kmer length to consider.",
)
@click.option(
    "--max_k",
    "-Mk",
    type=int,
    default=DEFAULT_MAX_KMER_LEN,
    show_default=True,
    help="Maximum kmer length to consider.",
)
def cli_create_db(
    kmer_to_protein: Path,
    fasta: Path,
    db_path: Path,
    min_k: int,
    max_k: int,
):
    create_db(
        fasta=fasta,
        kmer_to_protein_path=kmer_to_protein,
        db_path=db_path,
        min_k=min_k,
        max_k=max_k,
    )


if __name__ == "__main__":
    setup_logger()
    cli_create_db()
