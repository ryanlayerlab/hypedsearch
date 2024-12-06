import sqlite3

import pytest
from click.testing import CliRunner

from src.create_kmer_db import create_kmer_mass_db, main
from src.erik_constants import MASS, SEQ, TEST_DIR
from src.erik_utils import Protein


class TestCreateKmerMassDb:
    @staticmethod
    def test_basic_test(tmp_path):
        # Arrange
        db_file = tmp_path / "kmers.db"
        table_name = "test"
        # not using "B" because it's not an amino acid symbol
        aa_masses = {
            "A": 1,
            "H": 2,
            "C": 3,
            "D": 4,
            "E": 5,
            "F": 6,
            "G": 7,
        }
        proteins = [
            Protein(
                seq="AHC",
                desc="AHC",
            ),
            Protein(seq="DEFG", desc="DEFG"),
            Protein(seq="AHCDEFG", desc="AHCDEFG"),
        ]
        expected = [
            # 1-mers
            ("A", 1),
            ("H", 2),
            ("C", 3),
            ("D", 4),
            ("E", 5),
            ("F", 6),
            ("G", 7),
            # 2-mers
            ("AH", 1 + 2),
            ("HC", 2 + 3),
            ("DE", 4 + 5),
            ("EF", 5 + 6),
            ("FG", 6 + 7),
            ("CD", 3 + 4),
            # 3-merrs
            ("AHC", 1 + 2 + 3),
            ("DEF", 4 + 5 + 6),
            ("EFG", 5 + 6 + 7),
            ("HCD", 2 + 3 + 4),
            ("CDE", 3 + 4 + 5),  # 3-mers
        ]

        # Act
        create_kmer_mass_db(
            proteins=proteins,
            db_file=db_file,
            table_name=table_name,
            amino_acid_mass_lookup=aa_masses,
            max_peptide_len=3,
        )

        # Assert
        connection = sqlite3.connect(db_file, timeout=10)
        cursor = connection.cursor()
        cursor.execute(f"SELECT {SEQ}, {MASS} FROM {table_name}")
        actual = cursor.fetchall()
        assert len(actual) == len(expected)
        for kmer in expected:
            assert kmer in actual


class TestMain:
    @staticmethod
    def test_main(tmp_path, caplog):
        fasta_path = TEST_DIR / "test.fasta"
        db_file = tmp_path / "kmers.db"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta_path", fasta_path, "--db_file", db_file])
        assert db_file.exists()
        assert result.exit_code == 0
