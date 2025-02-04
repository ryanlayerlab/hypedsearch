import os
import shutil
import sqlite3
import sys
import time
from copy import deepcopy
from dataclasses import asdict, dataclass
from typing import Dict, List, Tuple

from src.erik_constants import (
    B_ION_AS_INT,
    CHARGE,
    END,
    ID,
    ION,
    ION_CHARGES_TO_CONSIDER,
    KMER_TABLE,
    MASS,
    MAX_KMER_LEN,
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
from src.lookups.data_classes import Kmer, KmerIons, KmerTableRow, Protein


def b_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + (charge * PROTON_MASS)) / charge
    return neutral_mass


def y_ion_neutral_mass(
    aa_seq: str,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> float:
    aa_mass_sum = sum([amino_acid_mass_lookup[aa] for aa in aa_seq])
    neutral_mass = (aa_mass_sum + WATER_MASS + (charge * PROTON_MASS)) / charge
    return neutral_mass


class ProteinKmerDb:
    def __init__(
        self, db_path: str, max_kmer_len: int = MAX_KMER_LEN, reset: bool = True
    ):
        self.connection = sqlite3.connect(db_path)
        self.cursor = self.connection.cursor()
        self.max_kmer_len = max_kmer_len
        # self.query_protein_average = 0
        # self.query_mass_average = 0
        # self.protein_count = 0
        # self.mass_count = 0
        # self.fetchall_protein_average = 0
        # self.fetchall_mass_average = 0
        if reset:
            # Create table to store proteins
            self.cursor.execute(f"DROP TABLE IF EXISTS {PROTEIN_TABLE}")
            self.cursor.execute(
                f"""
                CREATE TABLE {PROTEIN_TABLE} (
                    {ID} INTEGER PRIMARY KEY,
                    {SEQ} TEXT
                )
                """
            )

            self.cursor.execute(f"DROP TABLE IF EXISTS {KMER_TABLE}")
            self.cursor.execute(
                f"""
                CREATE TABLE {KMER_TABLE} (
                    {MASS} REAL,
                    {START} INTEGER,
                    {END} INTEGER,
                    {ION} INTEGER,
                    {CHARGE} INTEGER,
                    {PROTEIN_ID} INTEGER
                )
                """
            )

    def insert_kmers(self, kmers: List[KmerTableRow]):
        kmers = [asdict(kmer) for kmer in kmers]
        self.cursor.executemany(
            f"INSERT INTO {KMER_TABLE} VALUES(:{MASS}, :{START}, :{END}, :{ION}, :{CHARGE}, :{PROTEIN_ID})",
            kmers,
        )
        self.connection.commit()

    def insert_proteins(self, proteins: List[Protein]):
        protein_data = []
        for protein in proteins:
            protein_data.append({ID: protein.id, SEQ: protein.seq})
        self.cursor.executemany(
            f"INSERT INTO {PROTEIN_TABLE} VALUES(:{ID}, :{SEQ})", protein_data
        )
        self.connection.commit()

    def read_kmers(self):
        rows = self.cursor.execute(f"SELECT * FROM {KMER_TABLE}").fetchall()
        print(rows)

    def query_mass_kmers(
        self,
        fragment_id,
        precursor_mass,
        precursor_charge,
        mass,
        tolerance,
        number_decimal_places,
    ):
        """
        Gets all the kmers within the tolerance of a given mass.
        """
        upper_bound = mass + tolerance
        lower_bound = mass - tolerance
        # ERIK: the kmer table doesn't actual have the kmer amino acid sequences so this
        # query gets each kmer's amino acid sequence from the database
        # ERIK: this won't be very fast
        self.cursor.execute(
            f"""
            CREATE TABLE temp.mass AS
            SELECT
                ? AS fragment_id,
                ? AS precursor_mass,
                ? AS precursor_charge,
                k.*,
                SUBSTR(p.{SEQ}, k.{START}, k.{END} - k.{START} + 1) AS subsequence
            FROM {KMER_TABLE} AS k
            INNER JOIN {PROTEIN_TABLE} AS p ON k.{PROTEIN_ID} = p.{ID}
            WHERE k.{MASS} BETWEEN ? AND ?
            ORDER BY k.{PROTEIN_ID}, k.{START};
            """,
            (fragment_id, precursor_mass, precursor_charge, lower_bound, upper_bound),
        )
        b_rows = self.cursor.execute(
            f"""
            SELECT
                fragment_id,
                precursor_mass,
                precursor_charge,
                {PROTEIN_ID},
                ROUND(mass, ?),
                {START},
                {END},
                {ION},
                {CHARGE},
                subsequence,
                'N'
            FROM temp.mass
            WHERE ion = {B_ION_AS_INT}
            ORDER BY {PROTEIN_ID}, {START}, {END};
            """,
            (number_decimal_places,),
        ).fetchall()
        y_rows = self.cursor.execute(
            f"""
            SELECT
                fragment_id,
                precursor_mass,
                precursor_charge,
                {PROTEIN_ID},
                ROUND(mass, ?),
                {START},
                {END},
                {ION},
                {CHARGE},
                subsequence,
                'N'
            FROM temp.mass
            WHERE ion = {Y_ION_AS_INT}
            ORDER BY {PROTEIN_ID}, {START}, {END};
            """,
            (number_decimal_places,),
        ).fetchall()
        self.cursor.execute("DROP TABLE temp.mass")
        return b_rows, y_rows

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
        print(results)

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

    def get_protein(self, pid):
        row = self.cursor.execute(
            "SELECT * FROM proteins WHERE id = ?", (pid,)
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

    def db_make_set_for_protein_digest(
        self,
        protein_id,
        protein,
        max_peptide_length,
        digest_left,
        digest_right,
        number_decimal_places,
    ):
        data = []
        seq_len = len(protein)
        count_max = 1000000
        for size in range(2, max_peptide_length + 1):  # k-mer size
            for start in range(0, seq_len - size + 1):
                end = start + size
                kmer = protein[start:end]  # "ABCD"
                if (
                    kmer[0] in digest_left  # [""]
                    or digest_left == ["-"]
                    or (start > 0 and protein[start - 1] in digest_right)
                ):
                    bad_chars = ["B", "X", "U", "Z", "O", "J"]
                    if not any(x in bad_chars for x in kmer):
                        data_list = self.get_kmers_for_protein(
                            kmer, start, end, protein_id, "b"
                        )
                        data.extend(data_list)
                        if len(data) > count_max:
                            rounded_data = round(data, number_decimal_places)
                            self.insert(rounded_data)
                            data.clear()
                if (
                    kmer[-1] in digest_right  # ["K", "R"]
                    or digest_right == ["-"]
                    or (end < seq_len and protein[end] in digest_left)
                ):
                    bad_chars = ["B", "X", "U", "Z", "O", "J"]
                    if not any(x in bad_chars for x in kmer):
                        data_list = self.get_kmers_for_protein(
                            kmer, start, end, protein_id, "y"
                        )
                        data.extend(data_list)
                        if len(data) > count_max:
                            rounded_data = round(data, number_decimal_places)
                            self.insert(rounded_data)
                            data.clear()
        return data

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

    def insert_prepped_kmers(
        self,
        kv_proteins,
        max_peptide_length,
        digest_left,
        digest_right,
        number_decimal_places,
    ):
        plen = len(kv_proteins)
        last_percent = 0
        all_data = []

        for protein_id, (_, protein) in enumerate(kv_proteins):
            enough_space, last_percent = self.check_for_enough_disk_space(
                protein_id, plen, last_percent
            )
            if enough_space:
                data = self.db_make_set_for_protein_digest(
                    protein_id,
                    protein,
                    max_peptide_length,
                    digest_left,
                    digest_right,
                    number_decimal_places,
                )
                all_data.extend(data)
            else:
                sys.exit(1)

        if len(all_data) != 0:
            self.insert_kmers(all_data)

    def insert_prepped_proteins(self, kv_proteins):
        all_data = []
        protein_id = 0
        for kv_protein in kv_proteins:
            (description, amino_acids) = kv_protein
            input_tuple = (protein_id, description, amino_acids)
            all_data.append(input_tuple)
            protein_id = protein_id + 1
        self.insert_proteins(all_data)

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


def prepare_protein_kmer_database(
    db_path: str,
    fasta_path: str,
    max_kmer_len: int = MAX_KMER_LEN,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> ProteinKmerDb:
    """ """
    proteins = get_proteins_from_fasta(fasta_path=fasta_path)
    protein_kmer_db = ProteinKmerDb(db_path=db_path, max_kmer_len=max_kmer_len)
    for protein in proteins:
        add_protein_and_its_kmers_to_db(
            protein=protein,
            db=protein_kmer_db,
            charges_to_consider=charges_to_consider,
            amino_acid_mass_lookup=amino_acid_mass_lookup,
        )
    return protein_kmer_db


def kmer_rows_to_structured_object(kmer_rows: List) -> KmerTableRow:
    return [KmerTableRow(*kmer) for kmer in kmer_rows]


def add_protein_and_its_kmers_to_db(
    protein: Protein,
    db: ProteinKmerDb,
    charges_to_consider: List[int] = ION_CHARGES_TO_CONSIDER,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> None:
    from src.erik import generate_kmers

    # Add protein to protein table
    db.insert_proteins(proteins=[protein])

    # Add protein's kmers to kmer table
    kmers = generate_kmers(peptide=protein.seq, max_k=db.max_kmer_len)
    kmers_for_table = []
    for kmer in kmers:
        for charge in charges_to_consider:
            kmer_ions = get_b_ion_and_y_ion_corresponding_to_kmer(
                kmer=kmer,
                protein_id=protein.id,
                charge=charge,
                amino_acid_mass_lookup=amino_acid_mass_lookup,
            )
            kmers_for_table.extend([kmer_ions.b_ion, kmer_ions.y_ion])

    db.insert_kmers(kmers=kmers_for_table)


def get_b_ion_and_y_ion_corresponding_to_kmer(
    kmer: Kmer,
    protein_id: int,
    charge: int,
    amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
) -> KmerIons:

    b_mass = b_ion_neutral_mass(
        aa_seq=kmer.seq,
        charge=charge,
        amino_acid_mass_lookup=amino_acid_mass_lookup,
    )
    y_mass = y_ion_neutral_mass(
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
        b_ion=KmerTableRow(
            mass=b_mass,
            start=kmer.inclusive_start,
            end=kmer.exclusive_end,
            ion=B_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        ),
        y_ion=KmerTableRow(
            mass=y_mass,
            start=kmer.inclusive_start,
            end=kmer.exclusive_end,
            ion=Y_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        ),
    )


@dataclass
class MatchedRow:
    protein_id: int
    mass: float
    start: int
    stop: int
    ion: int
    charge: int
    seq: str


def create_table_of_ions_within_ppm_tolerance(
    db: ProteinKmerDb, mass: float, table_name: str, ppm_tol: int
):
    from src.erik import relative_ppm_tolerance_in_daltons

    # frag_id, p_mass, p_charge = "fragment_id", "precursor_mass", "precursor_charge"
    da_tol = relative_ppm_tolerance_in_daltons(ppm=ppm_tol, ref_mass=mass)
    lower_bound, upper_bound = mass - da_tol, mass + da_tol
    db.cursor.execute(
        f"""
        CREATE TABLE {table_name} AS
        SELECT
            k.*,
            SUBSTR(p.{SEQ}, k.{START}, k.{END} - k.{START} + 1) AS {SUBSEQ}
        FROM {KMER_TABLE} AS k
        INNER JOIN {PROTEIN_TABLE} AS p ON k.{PROTEIN_ID} = p.{ID}
        WHERE k.mass BETWEEN ? AND ?
        ORDER BY k.{PROTEIN_ID}, k.{START}; 
        """,
        (lower_bound, upper_bound),
    )


def query_mass_kmers(
    db: ProteinKmerDb,
    # fragment_id: int,
    # precursor_mass: float,
    # precursor_charge: int,
    mass: float,
    ppm_tol: int,
    number_decimal_places: int,
):
    from src.erik import relative_ppm_tolerance_in_daltons

    # frag_id, p_mass, p_charge = "fragment_id", "precursor_mass", "precursor_charge"
    table_name, subseq = "temp_mass", "subsequence"
    da_tol = relative_ppm_tolerance_in_daltons(ppm=ppm_tol, ref_mass=mass)
    lower_bound, upper_bound = mass - da_tol, mass + da_tol
    # ERIK: the kmer table doesn't actual have the kmer amino acid sequences so this
    # query gets each kmer's amino acid sequence from the database
    # ERIK: this won't be very fast
    db.cursor.execute(
        f"""
        CREATE TABLE {table_name} AS
        SELECT
            k.*,
            SUBSTR(p.{SEQ}, k.{START}, k.{END} - k.{START} + 1) AS {subseq}
        FROM {KMER_TABLE} AS k
        INNER JOIN {PROTEIN_TABLE} AS p ON k.{PROTEIN_ID} = p.{ID}
        WHERE k.mass BETWEEN ? AND ?
        ORDER BY k.{PROTEIN_ID}, k.{START}; 
        """,
        (lower_bound, upper_bound),
    )

    def ion_rows_query(ion: str):
        return f"""
            SELECT
                {PROTEIN_ID},
                ROUND({MASS}, ?),
                {START},
                {END},
                {ION},
                {CHARGE},
                {subseq}
            FROM {table_name}
            WHERE ion = {ion}
            ORDER BY {PROTEIN_ID}, {START}, {END};
            """

    b_rows = db.cursor.execute(
        ion_rows_query(ion=B_ION_AS_INT),
        (number_decimal_places,),
    ).fetchall()
    y_rows = db.cursor.execute(
        ion_rows_query(ion=Y_ION_AS_INT),
        (number_decimal_places,),
    ).fetchall()
    db.cursor.execute(f"DROP TABLE {table_name}")
    return [MatchedRow(*ion) for ion in b_rows] + [MatchedRow(*ion) for ion in y_rows]
