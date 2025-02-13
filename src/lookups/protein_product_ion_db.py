import logging
import os
import shutil
import sqlite3
import sys
import time
from copy import deepcopy
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

from src.erik_constants import (
    B_ION_AS_INT,
    CHARGE,
    END,
    ION,
    ION_CHARGES_TO_CONSIDER,
    MASS,
    MAX_KMER_LEN,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SEQ,
    START,
    SUBSEQ,
    Y_ION_AS_INT,
)
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


class ProteinProductIonDb:
    def __init__(
        self, db_path: str, max_kmer_len: int = MAX_KMER_LEN, reset: bool = True
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
                    {START} INTEGER,
                    {END} INTEGER,
                    {ION} INTEGER,
                    {CHARGE} INTEGER,
                    {PROTEIN_ID} INTEGER
                )
                """
            )

    def insert_product_ions(self, kmers: List[ProductIonTableRow]):
        kmers = [asdict(kmer) for kmer in kmers]
        self.cursor.executemany(
            f"INSERT INTO {PRODUCT_ION_TABLE} VALUES(:{MASS}, :{START}, :{END}, :{ION}, :{CHARGE}, :{PROTEIN_ID})",
            kmers,
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
        self, query_mass: float, tolerance: float
    ) -> List[IonWithProteinInfo]:
        # Constants
        table_name = "temp.mass"
        upper_bound = query_mass + tolerance
        lower_bound = query_mass - tolerance

        # This query:
        # 1. selects the product ions within the mass bounds
        # 2. gets the protein information and peptide within the protein that the product
        # ion corresponds to;
        # 3. combine the selected product ions and the corresponding protein info into a new table
        self.cursor.execute(
            f"""
            CREATE TABLE {table_name} AS
            SELECT
                k.*,
                SUBSTR(p.{SEQ}, k.{START}, k.{END} - k.{START} + 1) AS {SUBSEQ}
            FROM {PRODUCT_ION_TABLE} AS k
            INNER JOIN {PROTEIN_TABLE} AS p ON k.{PROTEIN_ID} = p.{PROTEIN_ID}
            WHERE k.{MASS} BETWEEN ? AND ?
            ORDER BY k.{PROTEIN_ID}, k.{START};
            """,
            (lower_bound, upper_bound),
        )

        # The command above creates a new table but doesn't return the rows in that table.
        # This query returns the rows created above ordered by (protein, start, end)
        matching_ions = self.cursor.execute(
            f"""
            SELECT
                {PROTEIN_ID},
                {MASS},
                {START},
                {END},
                {ION},
                {CHARGE},
                {SUBSEQ}
            FROM {table_name}
            ORDER BY {PROTEIN_ID}, {START}, {END};
            """
        ).fetchall()

        # Remove temporary table
        self.cursor.execute(f"DROP TABLE {table_name}")
        return [IonWithProteinInfo(**dict(ion)) for ion in matching_ions]

    def database_info(self):
        query = "SELECT name FROM sqlite_master WHERE type = 'table';"
        return self.run_query(query=query)

    def read_all_product_ions(self):
        query = f"SELECT * FROM {PRODUCT_ION_TABLE}"
        return self.run_query(query=query)

    def read_all_proteins(self):
        query = f"SELECT * FROM {PROTEIN_TABLE}"
        return self.run_query(query=query)

    def query_sequence_kmers(self, pid, start, end):
        rows = self.cursor.execute(
            "SELECT * FROM kmers where protein = ? and location_start = ? and location_end = ?",
            (pid, start, end),
        ).fetchall()
        return rows

    def query_extensions_to_length_kmers(self, target_mass, pid, start, end, dist):
        query_start = time.time()
        rows_cursor = self.cursor.execute(
            "SELECT * FROM kmers where protein = ? and location_start = ? and location_end = ? and mass <= ? and ion = 0 and charge = 2",
            (pid, start, end + dist - 1, target_mass),
        )  # and ion = 0 and charge = 1
        query_time = time.time() - query_start
        fetchall_start = time.time()
        rows = rows_cursor.fetchall()
        fetchall_time = time.time() - fetchall_start
        self.query_protein_average = (
            self.query_protein_average * self.protein_count + query_time
        ) / (self.protein_count + 1)
        self.fetchall_protein_average = (
            self.fetchall_protein_average * self.protein_count + fetchall_time
        ) / (self.protein_count + 1)
        self.protein_count += 1
        return rows

    def query_extensions_b_kmers(self, target_mass, pid, start, end, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute(
            "SELECT * FROM kmers where protein = ? and location_start = ? and location_end >= ? and mass < ? and ion = ? and charge = 2 order by location_end",
            (pid, start, end, target_mass, ion),
        )
        query_time = time.time() - query_start
        self.query_protein_average = (
            self.query_protein_average * self.protein_count + query_time
        ) / (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor

    def query_extensions_y_kmers(self, target_mass, pid, end, start, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute(
            "SELECT * FROM kmers where protein = ? and location_end = ? and location_start <= ? and mass < ? and ion = ? and charge = 2 order by location_start desc",
            (pid, end, start, target_mass, ion),
        )
        query_time = time.time() - query_start
        self.query_protein_average = (
            self.query_protein_average * self.protein_count + query_time
        ) / (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor

    def query_fetchall(self):
        return self.cursor.fetchall()

    def index_ion_mass_b_kmers(self):
        ctime = time.time()
        self.cursor.execute(
            "CREATE INDEX b_mass_ion_idx ON kmers(protein, location_start)"
        )
        print("time to index by protein", time.time() - ctime)

    def index_ion_mass_y_kmers(self):
        ctime = time.time()
        self.cursor.execute(
            "CREATE INDEX y_mass_ion_idx ON kmers(protein, location_end)"
        )
        print("time to index by protein", time.time() - ctime)

    def index_ion_mass_kmers(self):
        ctime = time.time()
        self.cursor.execute(
            "CREATE INDEX mass_ion_idx ON kmers(mass, protein, location_start)"
        )
        print("time to index", time.time() - ctime)

    def check_sizes_kmers(self):
        size = self.cursor.execute("SELECT count(*) FROM kmers").fetchall()
        print(size)

    def run_query(self, query):
        results = self.cursor.execute(query).fetchall()
        return results

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

    def get_all_rows_from_table(self, table_name: str):
        rows = self.cursor.execute(f"SELECT * FROM {table_name}").fetchall()
        return rows

    def get_protein_by_id(self, protein_id):
        row = self.cursor.execute(
            f"SELECT * FROM {PROTEIN_TABLE} WHERE {PROTEIN_ID} = ?", (protein_id,)
        ).fetchone()
        return row

    # def get_max_mass(self, sequence, ion, charge):
    #     if ion == "y":
    #         total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
    #         total += sum([AMINO_ACIDS[aa] for aa in sequence])
    #         mz = total / charge
    #         return mz
    #     else:
    #         total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
    #         total += sum([AMINO_ACIDS[aa] for aa in sequence])
    #     mz = total / charge
    #     return mz

    # def get_kmers_for_protein(self, kmer, start, end, protein_id, ion):
    #     data_list = []
    #     for charge in [1, 2]:
    #         mass = self.get_max_mass(kmer, ion=ion, charge=charge)
    #         ion_int = 0 if ion == "b" else 1
    #         input_tuple = (mass, start, end, ion_int, charge, protein_id)
    #         data_list.append(input_tuple)
    #     return data_list

    def check_for_enough_disk_space(self, protein_id, plen, last_percent):
        percent = int((protein_id + 1) * 100 / plen)
        print(f"\rOn protein {protein_id + 1}/{plen} [{percent}%]", end="")
        if percent != last_percent:
            last_percent = percent
            free = shutil.disk_usage("/")[2]
            free = free / (1024**3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB")
                return False, last_percent
            else:
                return True, last_percent
        return True, last_percent

    def populate_database(
        self,
        kv_proteins,
        max_peptide_length,
        digest_left,
        digest_right,
        number_decimal_places,
    ):
        self.insert_prepped_proteins(kv_proteins)
        self.insert_prepped_kmers(
            kv_proteins,
            max_peptide_length,
            digest_left,
            digest_right,
            number_decimal_places,
        )
        self.index_ion_mass_kmers()
        self.index_ion_mass_b_kmers()
        self.index_ion_mass_y_kmers()


def create_protein_product_ion_db(
    db_path: str,
    fasta_path: str,
    max_kmer_len: int = MAX_KMER_LEN,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    verbose: Optional[bool] = False,
) -> ProteinProductIonDb:
    """ """
    fcn_start_time = time.time()
    logger.info(f"Adding proteins from fasta {fasta_path}")
    proteins = list(get_proteins_from_fasta(fasta_path=fasta_path))
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

    template_kmer_for_table = {
        START: kmer.inclusive_start,
        END: kmer.exclusive_end,
        CHARGE: charge,
        PROTEIN_ID: protein_id,
    }
    kmer_as_b_ion = deepcopy(template_kmer_for_table)
    kmer_as_b_ion[MASS] = b_mass
    kmer_as_y_ion = deepcopy(template_kmer_for_table)
    kmer_as_y_ion[MASS] = y_mass

    return KmerIons(
        b_ion=ProductIonTableRow(
            mass=b_mass,
            start=kmer.inclusive_start,
            end=kmer.exclusive_end,
            ion=B_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        ),
        y_ion=ProductIonTableRow(
            mass=y_mass,
            start=kmer.inclusive_start,
            end=kmer.exclusive_end,
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

    db.insert_product_ions(kmers=ions_for_table)
