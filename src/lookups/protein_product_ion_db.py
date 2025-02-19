import logging
import os
import shutil
import sqlite3
import sys
import time
from copy import deepcopy
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

from src.erik_constants import (
    B_ION_AS_INT,
    CHARGE,
    EXCLUSIVE_END,
    INCLUSIVE_START,
    ION,
    ION_CHARGES_TO_CONSIDER,
    MASS,
    MAX_KMER_LEN,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SEQ,
    SUBSEQ,
    Y_ION_AS_INT,
)
from src.erik_utils import get_time_in_diff_units, relative_ppm_tolerance_in_daltons
from src.fasta_utils import get_proteins_from_fasta
from src.lookups.constants import (
    AMINO_ACID_MASSES,
    DOUBLY_CHARGED_B_BASE,
    DOUBLY_CHARGED_Y_BASE,
    PROTON_MASS,
    SINGLY_CHARGED_B_BASE,
    SINGLY_CHARGED_Y_BASE,
    WATER_MASS,
)
from src.lookups.data_classes import (
    IonWithProteinInfo,
    Kmer,
    KmerIons,
    ProductIonTableRow,
    Protein,
)

logger = logging.getLogger(__name__)


def compute_b_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + (charge * PROTON_MASS)) / charge
    return neutral_mass


def compute_y_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + WATER_MASS + (charge * PROTON_MASS)) / charge
    return neutral_mass


def load_existing_protein_product_ion_db(db_path: str):
    db = ProteinProductIonDb(
        db_path=db_path,
        reset=False,
    )
    return db


class ProteinProductIonDb:
    def __init__(
        self,
        db_path: str,
        reset: bool = True,
        max_kmer_len: Optional[int] = None,
    ):
        self.connection = sqlite3.connect(db_path)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        self.max_kmer_len = max_kmer_len

        if reset:
            # Create table to store proteins
            self.cursor.execute(f"DROP TABLE IF EXISTS {PROTEIN_TABLE}")
            self.cursor.execute(
                f"""
                CREATE TABLE {PROTEIN_TABLE} (
                    {PROTEIN_ID} INTEGER PRIMARY KEY,
                    {SEQ} TEXT
                )
                """
            )

            self.cursor.execute(f"DROP TABLE IF EXISTS {PRODUCT_ION_TABLE}")
            self.cursor.execute(
                f"""
                CREATE TABLE {PRODUCT_ION_TABLE} (
                    {MASS} REAL,
                    {INCLUSIVE_START} INTEGER,
                    {EXCLUSIVE_END} INTEGER,
                    {ION} INTEGER,
                    {CHARGE} INTEGER,
                    {PROTEIN_ID} INTEGER
                )
                """
            )

    def insert_product_ions(self, product_ions: List[ProductIonTableRow]):
        product_ions = [asdict(kmer) for kmer in product_ions]
        self.cursor.executemany(
            f"INSERT INTO {PRODUCT_ION_TABLE} VALUES(:{MASS}, :{INCLUSIVE_START}, :{EXCLUSIVE_END}, :{ION}, :{CHARGE}, :{PROTEIN_ID})",
            product_ions,
        )
        self.connection.commit()

    def insert_proteins(self, proteins: List[Protein]):
        protein_data = []
        for protein in proteins:
            protein_data.append({PROTEIN_ID: protein.protein_id, SEQ: protein.sequence})
        self.cursor.executemany(
            f"INSERT INTO {PROTEIN_TABLE} VALUES(:{PROTEIN_ID}, :{SEQ})", protein_data
        )
        self.connection.commit()

    def get_ions_within_mass_tolerance(
        self,
        query_mass: float,
        mz_tolerance: Optional[float] = None,
        ppm_tolerance: Optional[float] = None,
        # precursor_charge: Optional[int] = 1000,
    ) -> List[IonWithProteinInfo]:
        # Constants
        table_name = "filtered_ions"
        if ppm_tolerance is not None:
            mz_tolerance = relative_ppm_tolerance_in_daltons(
                ppm=ppm_tolerance, ref_mass=query_mass
            )

        upper_bound = query_mass + mz_tolerance
        lower_bound = query_mass - mz_tolerance

        query = f"""
            WITH {table_name} AS (
                SELECT
                    k.*,
                    SUBSTR(p.{SEQ}, k.{INCLUSIVE_START} + 1, k.{EXCLUSIVE_END} - k.{INCLUSIVE_START}) AS {SUBSEQ}
                FROM {PRODUCT_ION_TABLE} AS k
                INNER JOIN {PROTEIN_TABLE} AS p ON k.{PROTEIN_ID} = p.{PROTEIN_ID}
                WHERE k.{MASS} BETWEEN {lower_bound} AND {upper_bound}
            )
            SELECT
                {PROTEIN_ID},
                {MASS},
                {INCLUSIVE_START},
                {EXCLUSIVE_END},
                {ION},
                {CHARGE},
                {SUBSEQ}
            FROM filtered_ions
            ORDER BY {PROTEIN_ID}, {INCLUSIVE_START}, {EXCLUSIVE_END};
        """
        matching_ions = self.run_query(query=query)
        return [IonWithProteinInfo(**dict(ion)) for ion in matching_ions]

    def create_index_on_product_ion_mass(self):
        start_time = time.time()
        self.cursor.execute(f"CREATE INDEX {MASS} ON {PRODUCT_ION_TABLE}({MASS})")
        logger.info(
            f"Took {get_time_in_diff_units(time_sec=(time.time() - start_time))}"
        )

    def index_ion_mass_b_kmers(self):
        ctime = time.time()
        self.cursor.execute(
            "CREATE INDEX b_mass_ion_idx ON kmers(protein, location_start)"
        )
        print("time to index by protein", time.time() - ctime)

    # Helper methods for querying and investigating the database
    def run_query(self, query):
        results = self.cursor.execute(query).fetchall()
        return results

    def get_tables_in_db(self):
        name = "name"
        query = f"SELECT {name} FROM sqlite_master WHERE type = 'table';"
        rows = self.run_query(query=query)
        return [dict(row)[name] for row in rows]

    def get_indices_in_db(self):
        name = "name"
        query = f"SELECT {name} FROM sqlite_master WHERE type='index';"
        rows = self.run_query(query=query)
        return [dict(row)[name] for row in rows]

    def read_all_table_rows(self, table_name: str):
        query = f"SELECT * FROM {table_name}"
        return self.run_query(query=query)

    def get_number_rows(self, table_name: str):
        return len(self.read_all_table_rows(table_name=table_name))

    def get_protein_by_id(self, protein_id):
        query = f"SELECT * FROM {PROTEIN_TABLE} WHERE {PROTEIN_ID} = {protein_id}"
        return self.run_query(query=query)

    def get_random_sample_of_rows(self, table_name: str, sample_size: int = 100):
        query = f"SELECT * FROM {table_name} ORDER BY RANDOM() LIMIT {sample_size};"
        rows = self.run_query(query=query)
        return rows

    def get_unique_values_in_column(self, table_name: str, colm_name: str):
        query = f"SELECT DISTINCT {colm_name} FROM {table_name}"
        rows = self.run_query(query=query)
        rows = [dict(row) for row in rows]
        return rows

    def get_charges_in_db(self):
        charges = self.get_unique_values_in_column(
            table_name=PRODUCT_ION_TABLE, colm_name=CHARGE
        )
        return [x[CHARGE] for x in charges]

    def count_ion_mass_kmers(self, mass, tol, ion):
        upper = mass + tol
        lower = mass - tol
        qtime = time.time()
        rows = self.cursor.execute(
            "SELECT count(*) FROM kmers where ion = ? and mass between ? and ?",
            (ion, lower, upper),
        ).fetchall()
        print("Query time:", time.time() - qtime)
        print(mass, ion, rows)
        return rows

    def print_db_info(self):
        print(
            f"Number product ions = {self.get_number_rows(table_name=PRODUCT_ION_TABLE)}"
        )
        print(f"Tables in DB: {self.get_tables_in_db()}")
        print(f"Charges in DB: {self.get_charges_in_db()}")


def filter_ions_by_charge(
    objects_with_charge: List[Any], max_charge: int, min_charge: int = 1
) -> List[Any]:
    return list(
        filter(
            lambda x: filter(
                lambda x: min_charge <= x.charge <= max_charge, objects_with_charge
            )
        )
    )


def format_rows_of_product_ion_table(rows):
    return [ProductIonTableRow(**dict(row)) for row in rows]


def create_protein_product_ion_db(
    db_path: str,
    fasta_path: Optional[str] = None,
    protein_seqs: Optional[List[str]] = None,
    max_kmer_len: int = MAX_KMER_LEN,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> ProteinProductIonDb:
    """ """
    fcn_start_time = time.time()
    # Need to either provide a FASTA from which to
    if fasta_path is not None:
        logger.info(f"Adding proteins from fasta {fasta_path}")
        proteins = list(get_proteins_from_fasta(fasta_path=fasta_path))
    else:
        assert protein_seqs is not None
        proteins = [
            Protein(sequence=seq, protein_id=protein_idx)
            for protein_idx, seq in enumerate(protein_seqs)
        ]
    num_proteins = len(proteins)
    db = ProteinProductIonDb(db_path=db_path, max_kmer_len=max_kmer_len)
    for protein_idx, protein in enumerate(proteins):
        iteration_start_time = time.time()
        logger.info(
            f"Adding product ions for protein {protein_idx+1} of {num_proteins}"
        )
        add_protein_and_its_product_ions_to_db(
            protein=protein,
            db=db,
            charges_to_consider=charges_to_consider,
            amino_acid_mass_lookup=amino_acid_mass_lookup,
        )
        logger.info(
            f"Adding protein's product ions took {round(time.time()-iteration_start_time, 2)} seconds"
        )
        logger.info(
            f"Total run-time thus far is {round(time.time()-fcn_start_time, 2)} seconds"
        )

    return db


def get_b_ion_and_y_ion_corresponding_to_kmer(
    kmer: Kmer,
    protein_id: int,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> KmerIons:

    b_mass = compute_b_ion_neutral_mass(
        aa_seq=kmer.seq,
        charge=charge,
        amino_acid_mass_lookup=amino_acid_mass_lookup,
    )
    y_mass = compute_y_ion_neutral_mass(
        aa_seq=kmer.seq,
        charge=charge,
        amino_acid_mass_lookup=amino_acid_mass_lookup,
    )

    return KmerIons(
        b_ion=ProductIonTableRow(
            mass=b_mass,
            inclusive_start=kmer.inclusive_start,
            exclusive_end=kmer.exclusive_end,
            ion=B_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        ),
        y_ion=ProductIonTableRow(
            mass=y_mass,
            inclusive_start=kmer.inclusive_start,
            exclusive_end=kmer.exclusive_end,
            ion=Y_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        ),
    )


def add_protein_and_its_product_ions_to_db(
    protein: Protein,
    db: ProteinProductIonDb,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> None:
    from src.erik import generate_kmers

    # Add protein to protein table
    db.insert_proteins(proteins=[protein])

    # Add protein's kmers to kmer table
    kmers = generate_kmers(peptide=protein.sequence, max_k=db.max_kmer_len)
    ions_for_table = []
    for kmer in kmers:
        for charge in charges_to_consider:
            kmer_ions = get_b_ion_and_y_ion_corresponding_to_kmer(
                kmer=kmer,
                protein_id=protein.protein_id,
                charge=charge,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
            )
            ions_for_table.extend([kmer_ions.b_ion, kmer_ions.y_ion])

    db.insert_product_ions(product_ions=ions_for_table)


def get_average_mass_search_time(db: ProductIonTableRow, sample_size: int = 100):
    # query = f"SELECT * FROM {PRODUCT_ION_TABLE} ORDER BY RANDOM() LIMIT {sample_size};"
    # rows = db.run_query(query=query)
    rows = db.get_random_sample_of_rows(
        table_name=PRODUCT_ION_TABLE, sample_size=sample_size
    )
    rows = format_rows_of_product_ion_table(rows=rows)

    # for row in rows:

    return rows
