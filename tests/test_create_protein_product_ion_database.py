from pathlib import Path

import pytest
from click.testing import CliRunner

from create_protein_product_ion_database import (
    create_protein_and_product_ion_database,
    database_protein_product_ion_db_filename,
)
from src.lookups.protein_product_ion_db import ProteinProductIonDb


class TestDatabaseName:
    @staticmethod
    def test_smoke(tmp_path):
        # Arrange
        fasta_path = tmp_path / "test.fasta"
        charges = [1, 2]
        max_k = 10
        expected_name = f"{fasta_path.name}_max_k={max_k}_charges={charges}.db"

        # Act
        db_name = database_protein_product_ion_db_filename(
            fasta_path=fasta_path, max_k=max_k, charges_to_consider=charges
        )

        # Assert
        assert db_name == expected_name


class TestCreateProteinAndProductIonDatabase:
    @staticmethod
    def test_smoke(tmp_path):
        # Arrange
        fasta_path = Path("tests/test.fasta")
        charges = [1, 2]
        max_k = 5
        runner = CliRunner()
        result = runner.invoke(
            create_protein_and_product_ion_database,
            [
                "-f",
                f"{fasta_path}",
                "-d",
                f"{tmp_path}",
                "-k",
                f"{max_k}",
                "-c",
                "1",
                "-c",
                "2",
            ],
        )

        # Assert
        # Load database
        db_name = database_protein_product_ion_db_filename(
            fasta_path=fasta_path, max_k=max_k, charges_to_consider=charges
        )
        db_path = tmp_path / db_name
        db = ProteinProductIonDb(db_path=db_path, max_kmer_len=max_k, reset=False)
        db.database_info
        pass
