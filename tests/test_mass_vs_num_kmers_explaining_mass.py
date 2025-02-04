import sqlite3
from pathlib import Path

import pytest

from scripts.mass_vs_num_kmers_explaining_mass import (
    get_num_of_explanatory_kmers,
    plot_explanatory_kmers,
    sample_from_kmer_db,
    sample_in_parallel,
    sample_masses_and_get_num_explanatory_kmers,
)
from src.erik_constants import GIT_REPO_DIR, KMER_TABLE, PLOTS_DIR

TEST_DB_PATH = GIT_REPO_DIR / "./dbs/test.fasta.db"
assert TEST_DB_PATH.exists()


class TestSampleFromKmerDb:
    @staticmethod
    def test_basic():
        db_path = TEST_DB_PATH
        num_rows = 3
        rows = sample_from_kmer_db(num_rows=num_rows, db_path=db_path)
        assert len(rows) == num_rows


class TestGetSumOfExplanatoryKmers:
    @staticmethod
    def test_basic():
        mass = 321.1324707
        db_path = TEST_DB_PATH
        ppm_tol = 10
        actual = get_num_of_explanatory_kmers(
            mass=mass, ppm_tol=ppm_tol, db_path=db_path
        )
        assert actual == 3


class TestSampleMassesAndGetNumExplanatoryKmers:
    @staticmethod
    def test_basic():
        num_samples = 10
        ppm_tol = 10
        explanatory_kmers = sample_masses_and_get_num_explanatory_kmers(
            num_samples=num_samples, kmer_db_path=TEST_DB_PATH, ppm_tol=ppm_tol
        )


class TestSampleInParallel:
    @staticmethod
    def test_basic():
        num_samples = 10
        num_folds = 5
        ppm_tol = 10
        actual = sample_in_parallel(
            num_samples=num_samples,
            num_folds=num_folds,
            ppm_tol=ppm_tol,
            kmer_db_path=TEST_DB_PATH,
        )

        assert len(actual) == (num_samples * num_folds)

        file_name = TEST_DB_PATH.name
        save_path = PLOTS_DIR / f"{file_name}.ppm_tol={ppm_tol}.png"
        plot_explanatory_kmers(
            explanatory_kmers=actual, save_path=save_path, ppm_tol=ppm_tol
        )
