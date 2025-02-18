import random
import sqlite3
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from sqlite3 import IntegrityError
from typing import List, Optional
from unittest.mock import patch

import numpy as np
import pandas as pd
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
    EXCLUSIVE_END,
    INCLUSIVE_START,
    ION,
    MASS,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    SEQ,
    SPECTRA_DIR,
    SUBSEQ,
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
from src.lookups.data_classes import IonWithProteinInfo, Kmer, Protein
from src.lookups.protein_product_ion_db import (
    KmerIons,
    ProductIonTableRow,
    ProteinProductIonDb,
    add_protein_and_its_product_ions_to_db,
    create_protein_product_ion_db,
    get_b_ion_and_y_ion_corresponding_to_kmer,
)
from tests.fixtures_and_helpers import create_fasta


def create_db_with_given_kmer_masses(
    b_ion_masses: List[float],
    y_ion_masses: List[float],
    protein_indices: Optional[List[int]] = None,
    db_path: str = ":memory:",
    max_kmer_len: int = 10,
    protein_seqs: List[str] = ["ABC"],
):
    if protein_indices is not None:
        assert len(protein_indices) == len(b_ion_masses) + len(y_ion_masses)
    else:
        protein_indices = [
            random.randint(0, len(protein_seqs) - 1)
            for _ in range(len(b_ion_masses) + len(y_ion_masses))
        ]
    # Constants
    charge = 1
    start = 0

    # Create kmers for table
    kmers = []
    idx = 0
    for _, mass in enumerate(b_ion_masses):
        # Randomly select one of the proteins
        protein_idx = protein_indices[idx]
        idx += 1
        end = random.randint(1, len(protein_seqs[protein_idx]) - 1)
        kmers.append(
            ProductIonTableRow(
                mass=mass,
                inclusive_start=start,
                exclusive_end=end,
                ion=B_ION_AS_INT,
                charge=charge,
                protein_id=protein_idx,
            )
        )
        if protein_indices is not None:
            protein_idx += 1

    for _, mass in enumerate(y_ion_masses):
        # Randomly select one of the proteins
        protein_idx = protein_indices[idx]
        idx += 1
        end = random.randint(1, len(protein_seqs[protein_idx]) - 1)
        kmers.append(
            ProductIonTableRow(
                mass=mass,
                inclusive_start=start,
                exclusive_end=end,
                ion=Y_ION_AS_INT,
                charge=charge,
                protein_id=protein_idx,
            )
        )

    # Create database
    db = ProteinProductIonDb(db_path=db_path, max_kmer_len=max_kmer_len)
    db.insert_proteins(
        proteins=[
            Protein(protein_id=p_id, sequence=p_seq)
            for p_id, p_seq in enumerate(protein_seqs)
        ]
    )
    db.insert_product_ions(product_ions=kmers)

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
        product_ion_rows = db.run_query(query=f"SELECT * FROM {PRODUCT_ION_TABLE}")
        protein_rows = db.run_query(query=f"SELECT * FROM {PROTEIN_TABLE}")
        assert len(product_ion_rows) == len(b_ion_masses) + len(y_ion_masses)
        assert len(proteins) == len(protein_rows)


class TestProteinKmerDbInit:
    @staticmethod
    def test_initialization(tmp_path):
        # Arrange
        db_path = tmp_path / "test.db"
        # Act
        db = ProteinProductIonDb(db_path=db_path, max_kmer_len=10)
        # Assert

        # Check that database has correct tables
        result = query_database(
            query="SELECT name FROM sqlite_master WHERE type = 'table'", db_path=db_path
        )
        table_names = [table[0] for table in result]

        assert len(table_names) == 2
        assert PROTEIN_TABLE in table_names
        assert PRODUCT_ION_TABLE in table_names

        # Check that tables have correct columns
        # kmer table
        expected_colms = set(
            [MASS, INCLUSIVE_START, EXCLUSIVE_END, ION, CHARGE, PROTEIN_ID]
        )
        result = query_database(
            query=f"PRAGMA table_info({PRODUCT_ION_TABLE})", db_path=db_path
        )
        column_names = set(colm[1] for colm in result)
        assert column_names == expected_colms

        # protein table
        expected_colms = set([PROTEIN_ID, SEQ])
        result = query_database(
            query=f"PRAGMA table_info({PROTEIN_TABLE})", db_path=db_path
        )
        column_names = set(colm[1] for colm in result)
        assert column_names == expected_colms


class TestInsertProteins:
    @staticmethod
    def test_multiple_proteins():
        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=3)
        proteins = [
            Protein(sequence="ATG", protein_id=0),
            Protein(sequence="CGT", protein_id=1),
        ]
        db.insert_proteins(proteins=proteins)
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert [Protein(**dict(protein)) for protein in actual] == proteins

    @staticmethod
    def test_one_protein():
        """Test adding one protein at a time"""
        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=3)
        protein_1 = Protein(sequence="ATG", protein_id=1)
        protein_2 = Protein(sequence="BC", protein_id=3)
        db.insert_proteins(proteins=[protein_1])
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert [Protein(**dict(protein)) for protein in actual] == [protein_1]

        db.insert_proteins(proteins=[protein_2])
        actual = db.cursor.execute("SELECT * FROM proteins").fetchall()
        assert [Protein(**dict(protein)) for protein in actual] == [
            protein_1,
            protein_2,
        ]

    @staticmethod
    def test_non_unique_ids():
        """Test to make sure protein IDs added to database are unique"""
        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=3)
        with pytest.raises(IntegrityError):
            db.insert_proteins(
                proteins=[
                    Protein(sequence="A", protein_id=1),
                    Protein(sequence="B", protein_id=1),
                ]
            )


class TestInsertKmers:
    @staticmethod
    def test_basic():
        # Arrange
        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=3)
        kmers = [
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=3,
                ion=4,
                charge=5,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=7,
                inclusive_start=8,
                exclusive_end=9,
                ion=10,
                charge=11,
                protein_id=1,
            ),
        ]
        # Act
        db.insert_product_ions(product_ions=kmers)

        # Assert
        actual = db.cursor.execute(f"SELECT * FROM {PRODUCT_ION_TABLE}").fetchall()
        assert [ProductIonTableRow(**dict(kmer)) for kmer in actual] == kmers


class TestGetBIonAndYIonCorrespondingToKmer:
    @staticmethod
    def test_basic():
        # Arrange
        seq, start, end = "ABC", 0, 3
        kmer = Kmer(seq=seq, inclusive_start=start, exclusive_end=end)
        b_mass, y_mass = 10, 20
        charge, protein_id = 3, 4
        expected_b_ion = ProductIonTableRow(
            mass=b_mass,
            inclusive_start=start,
            exclusive_end=end,
            ion=B_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        )
        expected_y_ion = ProductIonTableRow(
            mass=y_mass,
            inclusive_start=start,
            exclusive_end=end,
            ion=Y_ION_AS_INT,
            charge=charge,
            protein_id=protein_id,
        )

        # Act
        with patch(
            "src.lookups.protein_product_ion_db.compute_b_ion_neutral_mass",
            return_value=b_mass,
        ), patch(
            "src.lookups.protein_product_ion_db.compute_y_ion_neutral_mass",
            return_value=y_mass,
        ):
            actual = get_b_ion_and_y_ion_corresponding_to_kmer(
                kmer=kmer, protein_id=protein_id, charge=charge
            )

        # Assert
        assert actual == KmerIons(b_ion=expected_b_ion, y_ion=expected_y_ion)


class TestAddProteinAndItsProductIonsToDb:
    @staticmethod
    def test_add_one_protein():
        # Arrange
        max_k, charges = 3, [1]
        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=max_k)
        seq, protein_id = "ACD", 3
        kmers = [
            "A",
            "C",
            "D",
            "AC",
            "CD",
            "ACD",
        ]

        # Act
        add_protein_and_its_product_ions_to_db(
            protein=Protein(sequence=seq, protein_id=protein_id),
            db=db,
            charges_to_consider=charges,
        )

        # Assert
        # Check that proteins are added
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert len(actual) == 1
        assert dict(actual[0]) == {PROTEIN_ID: protein_id, SEQ: seq}

        # Check ions are added
        actual = pd.read_sql(f"SELECT * FROM {PRODUCT_ION_TABLE}", db.connection)
        # actual = db.cursor.execute(f"SELECT COUNT(*) FROM {KMER_TABLE}").fetchall()
        assert (
            actual.shape[0] == len(kmers) * 2
        )  # multiply by 2 because we're adding b- and y-ions


class TestPrepareProteinKmerDatabase:
    @staticmethod
    def test_smoke(tmp_path):
        # Arrange
        db_path = ":memory:"
        fasta_name = "test.fasta"
        max_k = 2
        seqs = [
            "ACD",
            "EFG",
        ]  # avoid using disallowed amino acid characters in generate_kmers
        create_fasta(folder=tmp_path, file_name=fasta_name, seqs=seqs)
        expected_proteins = [
            {PROTEIN_ID: seq_num, SEQ: seq} for seq_num, seq in enumerate(seqs)
        ]
        b_mass, y_mass = 1, 2
        charges = [1, 2]

        # Expected kmers
        expected_kmers = [
            # A
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=1,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=1,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=1,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=1,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            # C
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            # D
            ProductIonTableRow(
                mass=1,
                inclusive_start=2,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=2,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=2,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=2,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            # AC
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            # CD
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            # E
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=1,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=1,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=1,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=1,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            # F
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            # G
            ProductIonTableRow(
                mass=1,
                inclusive_start=2,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=2,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=2,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=2,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            # EF
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=0,
                exclusive_end=2,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            # FG
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=1,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=1,
                inclusive_start=1,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=2,
                inclusive_start=1,
                exclusive_end=3,
                ion=Y_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
        ]

        # Act
        with patch(
            "src.lookups.protein_product_ion_db.compute_b_ion_neutral_mass",
            return_value=b_mass,
        ), patch(
            "src.lookups.protein_product_ion_db.compute_y_ion_neutral_mass",
            return_value=y_mass,
        ):
            db = create_protein_product_ion_db(
                db_path=db_path,
                fasta_path=tmp_path / fasta_name,
                max_kmer_len=max_k,
                charges_to_consider=charges,
            )

        # Assert
        actual = db.cursor.execute(f"SELECT * FROM {PROTEIN_TABLE}").fetchall()
        assert [dict(protein) for protein in actual] == expected_proteins

        actual = db.cursor.execute(f"SELECT * FROM {PRODUCT_ION_TABLE}").fetchall()
        kmer_rows = [ProductIonTableRow(**dict(kmer)) for kmer in actual]
        assert Counter(expected_kmers) == Counter(kmer_rows)


class TestGetIonsWithinMassTolerance:
    @staticmethod
    def test_smoke():
        # Arrange
        b_ion_masses = [1, 1.8, 1.96]
        y_ion_masses = [2.02, 3]
        query_mass, tol = 2, 0.05
        proteins = ["ABC", "DEF"]
        db = create_db_with_given_kmer_masses(
            b_ion_masses=b_ion_masses, y_ion_masses=y_ion_masses, protein_seqs=proteins
        )

        # Act
        matching_ions = db.get_ions_within_mass_tolerance(
            query_mass=query_mass, mz_tolerance=tol
        )

        # Assert
        assert len(matching_ions) == 2
        assert (matching_ions[0].mass == 2.02) or (matching_ions[1].mass == 2.02)
        assert (matching_ions[0].mass == 1.96) or (matching_ions[1].mass == 1.96)

    @staticmethod
    def test_blah():
        # Define proteins and product ions to put in database
        proteins = [
            Protein(protein_id=0, sequence="ACDEFG"),
            Protein(protein_id=1, sequence="HIKLMN"),
        ]
        product_ions = [
            ProductIonTableRow(
                mass=1,
                inclusive_start=0,
                exclusive_end=2,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=1.96,
                inclusive_start=0,
                exclusive_end=1,
                ion=B_ION_AS_INT,
                charge=1,
                protein_id=0,
            ),
            ProductIonTableRow(
                mass=2.02,
                inclusive_start=1,
                exclusive_end=3,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
            ProductIonTableRow(
                mass=3,
                inclusive_start=2,
                exclusive_end=5,
                ion=B_ION_AS_INT,
                charge=2,
                protein_id=1,
            ),
        ]

        expected = [
            IonWithProteinInfo(**asdict(product_ions[1]) | {SUBSEQ: "A"}),
            IonWithProteinInfo(**asdict(product_ions[2]) | {SUBSEQ: "IK"}),
        ]

        db = ProteinProductIonDb(db_path=":memory:", max_kmer_len=10)
        db.insert_proteins(proteins=proteins)
        db.insert_product_ions(product_ions=product_ions)

        # Act
        query_mass, tol = 2, 0.05
        result = db.get_ions_within_mass_tolerance(
            query_mass=query_mass, mz_tolerance=tol
        )

        # Assert
        assert result == expected


class TestSQLite3:
    @staticmethod
    def test_one_based_and_inclusive_of_endpoints():
        """
        Test to verify that the SUBSTR command is 1-based
        """
        # Arrange
        seq = "ABCD"
        start, length = 2, 3
        expected = "BCD"

        # Act
        conn = sqlite3.connect(":memory:")
        cursor = conn.cursor()
        result = cursor.execute(
            f"SELECT SUBSTR('{seq}', {start}, {length});"
        ).fetchall()[0][0]
        assert result == expected

    @staticmethod
    def test_convert_python_zero_based_indexing_to_sql_one_based():
        # Arrange
        seq = "ABCDEF"
        inclusive_start, exclusive_end = 1, 5
        length = exclusive_end - inclusive_start
        expected = "BCDE"

        conn = sqlite3.connect(":memory:")
        cursor = conn.cursor()
        result = cursor.execute(
            f"SELECT SUBSTR('{seq}', {inclusive_start + 1}, {length});"
        ).fetchall()[0][0]
        assert result == expected
