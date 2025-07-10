import os
import sys
from pathlib import Path

import pytest

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir / "src"))

from src.constants import GIT_REPO_DIR


@pytest.fixture
def test_data_dir():
    return GIT_REPO_DIR / "tests/data"


@pytest.fixture
def comet_params(test_data_dir):
    return test_data_dir / "comet.params"


@pytest.fixture
def mouse_mzml(test_data_dir):
    return test_data_dir / "spectra/10_mouse_spectra.mzML"


@pytest.fixture
def mouse_fasta(test_data_dir):
    return test_data_dir / "mouse_proteome.fasta"


@pytest.fixture
def mouse_db(test_data_dir):
    return test_data_dir / "mouse_top_10_proteins.db"


@pytest.fixture
def unit_tests_dir():
    return GIT_REPO_DIR / "tests/unit"


@pytest.fixture
def integration_tests_dir():
    return GIT_REPO_DIR / "tests/integration"


@pytest.fixture
def snapshot_dir():
    return GIT_REPO_DIR / "tests/snapshots"


@pytest.fixture
def comet_txt(test_data_dir):
    return test_data_dir / "test_comet_result.txt"


@pytest.fixture
def crux_txt(test_data_dir):
    return test_data_dir / "crux.comet.1-10.txt"
