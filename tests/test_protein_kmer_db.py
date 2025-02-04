import random
import time
from collections import Counter
from pathlib import Path
from sqlite3 import IntegrityError
from typing import List
from unittest.mock import patch

import numpy as np
import pytest
from pyteomics import mzml

from src.erik import (  # Peptide,
    Spectrum,
    generate_kmers,
    get_b_ion_sequences,
    get_b_y_ion_sequences,
    get_data_for_spectrum,
    get_indices_of_largest_elements,
    get_y_ion_sequences,
    parse_spectrum,
    query_database,
)
from src.erik_constants import (
    B_ION_AS_INT,
    CHARGE,
    END,
    ID,
    ION,
    KMER_TABLE,
    MASS,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SEQ,
    SPECTRA_DIR,
    START,
    TEST_DIR,
    Y_ION_AS_INT,
)
from src.fasta_utils import get_proteins_from_fasta
from src.lookups.constants import (
    DOUBLY_CHARGED_B_BASE,
    DOUBLY_CHARGED_Y_BASE,
    SINGLY_CHARGED_B_BASE,
    SINGLY_CHARGED_Y_BASE,
)
from src.lookups.data_classes import Kmer, Protein
from src.lookups.protein_kmer_db import (
    KmerIons,
    KmerTableRow,
    ProteinKmerDb,
    add_protein_and_its_kmers_to_db,
    get_b_ion_and_y_ion_corresponding_to_kmer,
    kmer_rows_to_structured_object,
    prepare_protein_kmer_database,
)
from tests.fixtures_and_helpers import create_fasta


def create_db_with_given_kmer_masses(
    b_ion_masses: List[float],
    y_ion_masses: List[float],
    db_path: str = ":memory:",
    max_kmer_len: int = 10,
    protein_seqs: List[str] = ["ABC"],
):
    # Constants
    charge = 1
    start = 0

    # Create kmers for table
    kmers = []
    for _, mass in enumerate(b_ion_masses):
        # Randomly select one of the proteins
        protein_idx = random.randint(0, len(protein_seqs) - 1)
        end = random.randint(1, len(protein_seqs[protein_idx]) - 1)
        kmers.append(
            KmerTableRow(
                mass=mass,
                start=start,
                end=end,
                ion=B_ION_AS_INT,
                charge=charge,
                protein_id=protein_idx,
            )
        )

    for _, mass in enumerate(y_ion_masses):
        # Randomly select one of the proteins
        protein_idx = random.randint(0, len(protein_seqs) - 1)
        end = random.randint(1, len(protein_seqs[protein_idx]) - 1)
        kmers.append(
            KmerTableRow(
                mass=mass,
                start=start,
                end=end,
                ion=Y_ION_AS_INT,
                charge=charge,
                protein_id=protein_idx,
            )
        )

    # Create database
    db = ProteinKmerDb(db_path=db_path, max_kmer_len=max_kmer_len)
    db.insert_proteins(
        proteins=[
            Protein(id=p_id, seq=p_seq) for p_id, p_seq in enumerate(protein_seqs)
        ]
    )
    db.insert_kmers(kmers=kmers)

    return db


class TestCreateDbWithGivenKmerMasses:
    @staticmethod
    def test_smoke():
        b_ion_masses = [1, 2]
        y_ion_masses = [3, 4]
        proteins = ["ABC", "DEF"]
        db = create_db_with_given_kmer_masses(
            b_ion_masses=b_ion_masses, y_ion_masses=y_ion_masses, protein_seqs=proteins
        )
        kmer_rows = db.get_all_rows_from_table(table_name=KMER_TABLE)
        protein_rows = db.get_all_rows_from_table(table_name=PROTEIN_TABLE)
        assert len(kmer_rows) == len(b_ion_masses) + len(y_ion_masses)
        assert len(proteins) == len(protein_rows)


class TestProteinKmerDbInit:
    @staticmethod
    def test_initialization(tmp_path):
        # Arrange
        db_path = tmp_path / "test.db"
        # Act
        db = ProteinKmerDb(db_path=db_path, max_kmer_len=10)
        # Assert

        # Check that database has correct tables
        result = query_database(
            query="SELECT name FROM sqlite_master WHERE type = 'table'", db_path=db_path
        )
        table_names = [table[0] for table in result]

        assert len(table_names) == 2
        assert PROTEIN_TABLE in table_names
        assert KMER_TABLE in table_names

        # Check that tables have correct columns
        # kmer table
        expected_colms = set([MASS, START, END, ION, CHARGE, PROTEIN_ID])
        result = query_database(
            query=f"PRAGMA table_info({KMER_TABLE})", db_path=db_path
        )
        column_names = set(colm[1] for colm in result)
        assert column_names == expected_colms

        # protein table
        expected_colms = set([ID, SEQ])
        result = query_database(
            query=f"PRAGMA table_info({PROTEIN_TABLE})", db_path=db_path
        )
        column_names = set(colm[1] for colm in result)
        assert column_names == expected_colms


class TestInsertProteins:
    @staticmethod
    def test_multiple_proteins():
        db = ProteinKmerDb(db_path=":memory:", max_kmer_len=3)
        expected_out = [(0, "ATG"), (1, "CGT")]
        proteins = [Protein(seq="ATG", id=0), Protein(seq="CGT", id=1)]
        db.insert_proteins(proteins=proteins)
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert actual == expected_out

    @staticmethod
    def test_one_protein():
        """Test adding one protein at a time"""
        db = ProteinKmerDb(db_path=":memory:", max_kmer_len=3)
        db.insert_proteins(proteins=[Protein(seq="ATG", id=1)])
        expected_out = [(1, "ATG")]
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert actual == expected_out

        db.insert_proteins(proteins=[Protein(seq="BC", id=3)])
        expected_out = [(1, "ATG"), (3, "BC")]
        actual = db.cursor.execute("SELECT * FROM proteins").fetchall()
        assert actual == expected_out

    @staticmethod
    def test_non_unique_ids():
        """Test to make sure protein IDs added to database are unique"""
        db = ProteinKmerDb(db_path=":memory:", max_kmer_len=3)
        with pytest.raises(IntegrityError):
            db.insert_proteins(
                proteins=[Protein(seq="A", id=1), Protein(seq="B", id=1)]
            )


class TestInsertKmers:
    @staticmethod
    def test_basic():
        # Arrange
        db = ProteinKmerDb(db_path=":memory:", max_kmer_len=3)
        kmer1, kmer2 = (1, 0, 3, 4, 5, 6), (7, 8, 9, 1, 10, 11)
        kmers = [
            KmerTableRow(
                mass=kmer1[0],
                start=kmer1[1],
                end=kmer1[2],
                ion=kmer1[3],
                charge=kmer1[4],
                protein_id=kmer1[5],
            ),
            KmerTableRow(
                mass=kmer2[0],
                start=kmer2[1],
                end=kmer2[2],
                ion=kmer2[3],
                charge=kmer2[4],
                protein_id=kmer2[5],
            ),
        ]
        # Act
        db.insert_kmers(kmers=kmers)

        # Assert
        actual = db.cursor.execute(f"SELECT * FROM {KMER_TABLE}").fetchall()
        assert actual == [kmer1, kmer2]


class TestGetBIonAndYIonCorrespondingToKmer:
    @staticmethod
    def test_basic():
        # Arrange
        seq, start, end = "ABC", 0, 3
        kmer = Kmer(seq=seq, inclusive_start=start, exclusive_end=end)
        b_mass, y_mass = 10, 20
        charge, protein_id = 3, 4
        expected_b_ion = KmerTableRow(
            mass=b_mass,
            start=start,
            end=end,
            ion=B_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        )
        expected_y_ion = KmerTableRow(
            mass=y_mass,
            start=start,
            end=end,
            ion=Y_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        )

        # Act
        with patch(
            "src.lookups.protein_kmer_db.b_ion_neutral_mass", return_value=b_mass
        ), patch("src.lookups.protein_kmer_db.y_ion_neutral_mass", return_value=y_mass):
            actual = get_b_ion_and_y_ion_corresponding_to_kmer(
                kmer=kmer, protein_id=protein_id, charge=charge
            )

        # Assert
        assert actual == KmerIons(b_ion=expected_b_ion, y_ion=expected_y_ion)


class TestAddProteinAndItsKmersToDb:
    @staticmethod
    def test_1():
        # Arrange
        max_k, charges = 3, [1]
        db = ProteinKmerDb(db_path=":memory:", max_kmer_len=max_k)
        seq, protein_id = "ACDE", 3
        kmers = ["A", "C", "D", "E", "AC", "CD", "DE", "ACD", "CDE"]
        protein = Protein(seq=seq, id=protein_id)

        # Act
        add_protein_and_its_kmers_to_db(
            protein=protein, db=db, charges_to_consider=charges
        )

        # Assert
        # Check that proteins are added
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        expected = [(protein_id, seq)]
        assert actual == expected
        actual = db.cursor.execute(f"SELECT COUNT(*) FROM {KMER_TABLE}").fetchall()
        assert (
            actual[0][0] == len(kmers) * 2
        )  # multiply by 2 because we're adding b- and y-ions


class TestPrepareProteinKmerDatabase:
    @staticmethod
    def test_basic(tmp_path):
        # Arrange
        db_path = ":memory:"
        fasta_name = "test.fasta"
        max_k = 2
        seqs = [
            "ACD",
            "EFG",
        ]  # avoid using disallowed amino acid characters in generate_kmers
        create_fasta(folder=tmp_path, file_name=fasta_name, seqs=seqs)
        expected_proteins = [(seq_num, seq) for seq_num, seq in enumerate(seqs)]
        b_mass, y_mass = 1, 2
        charges = [1, 2]

        # Expected kmers
        expected_kmers = [
            # A
            KmerTableRow(
                mass=1, start=0, end=1, ion=B_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=0, end=1, ion=Y_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=1, start=0, end=1, ion=B_ION_AS_INT, charge=2, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=0, end=1, ion=Y_ION_AS_INT, charge=2, protein_id=0
            ),
            # C
            KmerTableRow(
                mass=1, start=1, end=2, ion=B_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=1, end=2, ion=Y_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=1, start=1, end=2, ion=B_ION_AS_INT, charge=2, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=1, end=2, ion=Y_ION_AS_INT, charge=2, protein_id=0
            ),
            # D
            KmerTableRow(
                mass=1, start=2, end=3, ion=B_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=2, end=3, ion=Y_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=1, start=2, end=3, ion=B_ION_AS_INT, charge=2, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=2, end=3, ion=Y_ION_AS_INT, charge=2, protein_id=0
            ),
            # AC
            KmerTableRow(
                mass=1, start=0, end=2, ion=B_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=0, end=2, ion=Y_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=1, start=0, end=2, ion=B_ION_AS_INT, charge=2, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=0, end=2, ion=Y_ION_AS_INT, charge=2, protein_id=0
            ),
            # CD
            KmerTableRow(
                mass=1, start=1, end=3, ion=B_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=1, end=3, ion=Y_ION_AS_INT, charge=1, protein_id=0
            ),
            KmerTableRow(
                mass=1, start=1, end=3, ion=B_ION_AS_INT, charge=2, protein_id=0
            ),
            KmerTableRow(
                mass=2, start=1, end=3, ion=Y_ION_AS_INT, charge=2, protein_id=0
            ),
            # E
            KmerTableRow(
                mass=1, start=0, end=1, ion=B_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=0, end=1, ion=Y_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=1, start=0, end=1, ion=B_ION_AS_INT, charge=2, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=0, end=1, ion=Y_ION_AS_INT, charge=2, protein_id=1
            ),
            # F
            KmerTableRow(
                mass=1, start=1, end=2, ion=B_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=1, end=2, ion=Y_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=1, start=1, end=2, ion=B_ION_AS_INT, charge=2, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=1, end=2, ion=Y_ION_AS_INT, charge=2, protein_id=1
            ),
            # G
            KmerTableRow(
                mass=1, start=2, end=3, ion=B_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=2, end=3, ion=Y_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=1, start=2, end=3, ion=B_ION_AS_INT, charge=2, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=2, end=3, ion=Y_ION_AS_INT, charge=2, protein_id=1
            ),
            # EF
            KmerTableRow(
                mass=1, start=0, end=2, ion=B_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=0, end=2, ion=Y_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=1, start=0, end=2, ion=B_ION_AS_INT, charge=2, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=0, end=2, ion=Y_ION_AS_INT, charge=2, protein_id=1
            ),
            # FG
            KmerTableRow(
                mass=1, start=1, end=3, ion=B_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=1, end=3, ion=Y_ION_AS_INT, charge=1, protein_id=1
            ),
            KmerTableRow(
                mass=1, start=1, end=3, ion=B_ION_AS_INT, charge=2, protein_id=1
            ),
            KmerTableRow(
                mass=2, start=1, end=3, ion=Y_ION_AS_INT, charge=2, protein_id=1
            ),
        ]

        # Act
        with patch(
            "src.lookups.protein_kmer_db.b_ion_neutral_mass", return_value=b_mass
        ), patch("src.lookups.protein_kmer_db.y_ion_neutral_mass", return_value=y_mass):
            db = prepare_protein_kmer_database(
                db_path=db_path,
                fasta_path=tmp_path / fasta_name,
                max_kmer_len=max_k,
                charges_to_consider=charges,
            )

        # Assert
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert actual == expected_proteins

        kmer_rows = db.get_all_rows_from_table(table_name=KMER_TABLE)
        kmer_rows = kmer_rows_to_structured_object(kmer_rows)
        assert Counter(expected_kmers) == Counter(kmer_rows)


class TestQueryMassKmers:
    @staticmethod
    def test_smoke():
        # Arrange
        b_ion_masses = [1, 1.8, 1.96]
        y_ion_masses = [2.02, 2.04, 3]
        query_mass, tol = 2, 0.05
        fragment_id = 1
        precursor_mass, precursor_charge = 10, 11
        proteins = ["ABC", "DEF"]
        db = create_db_with_given_kmer_masses(
            b_ion_masses=b_ion_masses, y_ion_masses=y_ion_masses, protein_seqs=proteins
        )

        # Act
        b_rows, y_rows = db.query_mass_kmers(
            fragment_id=fragment_id,
            precursor_mass=precursor_mass,
            precursor_charge=precursor_charge,
            mass=query_mass,
            tolerance=tol,
            number_decimal_places=10,
        )

        # Assert
        assert len(b_rows) + len(y_rows) == 3
