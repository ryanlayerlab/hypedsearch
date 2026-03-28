import sys
from pathlib import Path
from typing import Optional

import pytest

from src.hypedsearch import HypedsearchRunConfig
from src.mass_spectra import Mzml, Spectrum
from src.psm import CometPSM
from src.utils import load_json, to_json

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
def hs_config_path(test_data_dir):
    return test_data_dir / "test.hs.config.json"


@pytest.fixture
def hs_config(hs_config_path):
    return HypedsearchRunConfig.from_json(path=hs_config_path)


@pytest.fixture
def mouse_mzml_path(test_data_dir):
    return test_data_dir / "BMEM_AspN_Fxn4_scans1-20.mzML"


@pytest.fixture
def mouse_spectrum(mouse_mzml_path) -> Spectrum:
    return Spectrum.get_spectrum(scan=7, mzml=mouse_mzml_path)


@pytest.fixture
def mouse_mzml(mouse_mzml_path):
    return Mzml(path=mouse_mzml_path)


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


def update_hs_config_out_dir_and_save_json(
    old_config_path: str | Path,
    out_dir: str | Path,
    new_config_path: Optional[str | Path] = None,
) -> Path:
    data = load_json(path=old_config_path)
    data["parent_output_dir"] = str(out_dir)
    if new_config_path is None:
        new_config_path = Path(out_dir) / "hs.config.json"
    to_json(data=data, path=new_config_path)
    return new_config_path
