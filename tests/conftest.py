import sys
from pathlib import Path

import pytest

from src.mass_spectra import Spectrum
from src.psm import CometPSM

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir / "src"))


from src.constants import GIT_REPO_DIR


@pytest.fixture
def test_data_dir():
    return Path("tests/data")


@pytest.fixture
def comet_params(test_data_dir):
    return test_data_dir / "comet.params"


@pytest.fixture
def crux_comet_params(test_data_dir):
    return test_data_dir / "crux.comet.params"


@pytest.fixture
def test_hs_config_path(test_data_dir):
    return test_data_dir / "test.hs.config.json"


@pytest.fixture
def mouse_mzml_path(test_data_dir):
    return test_data_dir / "BMEM_AspN_Fxn4_scans1-20.mzML"


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


@pytest.fixture
def assign_confidence_txt(test_data_dir):
    return test_data_dir / "BMEM_AspN_Fxn4.assign-confidence.txt"


@pytest.fixture
def comet_psm(test_data_dir):
    path = test_data_dir / "example_comet_psm.json"
    psms = CometPSM.from_txt(
        txt="tests/data/BMEM_AspN_Fxn4/assign-confidence.target.txt"
    )
    psms = [psm for psm in psms if "sp|P99027|RLA2_MOUSE" in psm.proteins]
    psms[0].save_to_json(path=path)
    return CometPSM.load(path=path)


@pytest.fixture
def mouse_spectrum(mouse_mzml_path) -> Spectrum:
    return Spectrum.get_spectrum(scan=7, mzml=mouse_mzml_path)
