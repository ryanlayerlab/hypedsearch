import logging
from dataclasses import dataclass, fields
from typing import Any, List, Literal, Optional

from src.constants import (
    B_ION_TYPE,
    DEFAULT_MAX_K,
    DEFAULT_MIN_K,
    EXCLUSIVE_END,
    INCLUSIVE_START,
    ION_TYPE_TO_INT,
    MEMORY,
    NEUTRAL_MASS,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    Y_ION_TYPE,
    IonTypes,
)
from src.peptides_and_ions import Peptide, ProductIon
from src.sql_database import (
    PrimaryKey,
    SqlColumn,
    Sqlite3Database,
    SqlTableRow,
    add_data_classes_to_table,
)
from src.utils import Kmer, Position, log_params, relative_ppm_tolerance_in_daltons

logger = logging.getLogger(__name__)


@dataclass
class DbProtein(SqlTableRow):
    id: PrimaryKey[int]
    seq: str

    @classmethod
    def from_peptide(cls, peptide: Peptide) -> "DbProtein":
        assert peptide.id is not None
        return cls(id=peptide.id, seq=peptide.seq)


@dataclass
class DbProductIon(SqlTableRow):
    protein_id: int
    inclusive_start: int
    exclusive_end: int
    charge: int
    neutral_mass: float
    ion_type: int

    @classmethod
    def kmer_to_ion(
        cls, kmer: Kmer, charge: int, ion_type: IonTypes, protein_id: int
    ) -> "DbProductIon":
        product_ion = ProductIon(seq=kmer.seq, charge=charge, ion_type=ion_type)
        return cls(
            protein_id=protein_id,
            # position=kmer.position,
            inclusive_start=kmer.position.inclusive_start,
            exclusive_end=kmer.position.exclusive_end,
            charge=product_ion.charge,
            neutral_mass=product_ion.neutral_mass,
            ion_type=ION_TYPE_TO_INT[product_ion.ion_type.value],
        )

    @classmethod
    def from_peptide(
        cls,
        peptide: Peptide,
        charges: List[int],
        ion_types: List[IonTypes],
        max_k: int = DEFAULT_MAX_K,
        min_k: int = DEFAULT_MIN_K,
    ) -> List["DbProductIon"]:
        kmers = peptide.kmers(min_k=min_k, max_k=max_k)
        ions = []
        for kmer in kmers:
            for charge in charges:
                for ion_type in ion_types:
                    ions.append(
                        cls.kmer_to_ion(
                            kmer=kmer,
                            charge=charge,
                            ion_type=ion_type,
                            protein_id=peptide.id,
                        )
                    )
        return ions


class ProteinProductIonDb(Sqlite3Database):
    def __init__(
        self,
        charges: List[int],
        ion_types: List[IonTypes],
        db_path: Optional[str] = MEMORY,
        protein_table_name: str = PROTEIN_TABLE,
        product_ion_table_name: str = PRODUCT_ION_TABLE,
        min_k: int = DEFAULT_MIN_K,
        max_k: int = DEFAULT_MAX_K,
        protein_obj: type[DbProtein] = DbProtein,
        product_ion_obj: type[DbProductIon] = DbProductIon,
    ):
        super().__init__(db_path=db_path)
        self.protein_table_name = protein_table_name
        self.product_ion_table_name = product_ion_table_name
        self.min_k = min_k
        self.max_k = max_k
        self.charges = charges
        self.ion_types = ion_types
        self.protein_obj = protein_obj
        self.product_ion_obj = product_ion_obj

    def create_protein_table(self):
        self.create_table_from_dataclass(
            table_name=self.protein_table_name, obj=self.protein_obj
        )

    def create_product_ion_table(self):
        self.create_table_from_dataclass(
            table_name=self.product_ion_table_name, obj=self.product_ion_obj
        )

    def add_peptide_and_product_ions(
        self,
        peptide: Peptide,
    ) -> None:
        # Add peptide to DB
        db_protein = self.protein_obj.from_peptide(peptide=peptide)
        self.insert_dataclasses(
            table_name=self.protein_table_name, data_classes=[db_protein]
        )

        # Add product ions to DB
        db_ions = self.product_ion_obj.from_peptide(
            peptide=peptide,
            charges=self.charges,
            ion_types=self.ion_types,
            min_k=self.min_k,
            max_k=self.max_k,
        )
        self.insert_dataclasses(
            table_name=self.product_ion_table_name, data_classes=db_ions
        )

    def get_ions_within_mass_tolerance(
        self,
        query_mass: float,
        mz_tolerance: Optional[float] = None,
        ppm_tolerance: Optional[float] = None,
    ) -> DbProductIon:

        if ppm_tolerance is not None:
            mz_tolerance = relative_ppm_tolerance_in_daltons(
                ppm=ppm_tolerance, ref_mass=query_mass
            )

        upper_bound = query_mass + mz_tolerance
        lower_bound = query_mass - mz_tolerance

        query = f"""
            SELECT
                *
            FROM {self.product_ion_table_name} as product_ion
            WHERE product_ion.{NEUTRAL_MASS} BETWEEN {lower_bound} AND {upper_bound}
            ORDER BY {PROTEIN_ID}, {INCLUSIVE_START}, {EXCLUSIVE_END};
        """
        matching_ions = self.read_query(query=query)
        return [DbProductIon(**dict(ion)) for ion in matching_ions]

    def get_proteins(self):
        rows = self.all_table_rows(table_name=self.protein_table_name)
        return [self.protein_obj(**row) for row in rows]

    def get_product_ions(self):
        rows = self.all_table_rows(table_name=self.product_ion_table_name)
        return [self.product_ion_obj(**row) for row in rows]


@log_params
def create_and_populate_protein_and_product_ion_database(
    charges: List[int],
    ion_types: List[IonTypes],
    protein_obj: DbProtein = DbProtein,
    product_ion_obj: DbProductIon = DbProductIon,
    protein_table_name: str = PROTEIN_TABLE,
    product_ion_table_name: str = PRODUCT_ION_TABLE,
    db_path: str = MEMORY,
    min_k: int = DEFAULT_MIN_K,
    max_k: int = DEFAULT_MAX_K,
    peptides: Optional[List[Peptide]] = None,
    fasta_path: Optional[str] = None,
) -> ProteinProductIonDb:
    logger.info("Creating protein-product ion database...")
    if fasta_path is not None:
        peptides = Peptide.from_fasta(fasta_path=fasta_path)

    db = ProteinProductIonDb(
        charges=charges,
        ion_types=ion_types,
        db_path=db_path,
        protein_table_name=protein_table_name,
        product_ion_table_name=product_ion_table_name,
        min_k=min_k,
        max_k=max_k,
        protein_obj=protein_obj,
        product_ion_obj=product_ion_obj,
    )

    db.create_protein_table()
    db.create_product_ion_table()
    num_peptides = len(peptides)
    for p_idx, peptide in enumerate(peptides):
        logger.info(f"Adding protein {p_idx + 1} / {num_peptides}...")
        db.add_peptide_and_product_ions(peptide=peptide)

    return db


def protein_product_ion_database_file_name(
    charges: List[int],
    ion_types: List[IonTypes],
    file_name_prefix: str,
    min_k: int,
    max_k: int,
):
    ion_types = ",".join([ion_type.value for ion_type in ion_types])
    charges = ",".join([str(charge) for charge in charges])
    return f"{file_name_prefix}_ionTypes={ion_types}_charges={charges}_minK={min_k}_maxK={max_k}.db"
