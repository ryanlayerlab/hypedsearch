import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from src.psm import CometPSM

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir / "src"))

from click.testing import CliRunner

from src.constants import GIT_REPO_DIR, MOUSE_PROTEOME
from src.hypedsearch import HybridFormer, HybridPSMScorer, HypedsearchRunConfig


@pytest.fixture
def test_data_dir():
    return GIT_REPO_DIR / "tests/data"


@pytest.fixture
def comet_params(test_data_dir):
    return test_data_dir / "comet.params"


@pytest.fixture
def crux_comet_params(test_data_dir):
    return test_data_dir / "crux.comet.params"


@pytest.fixture
def mouse_mzml(test_data_dir):
    return test_data_dir / "spectra/10_mouse_spectra.mzML"


@pytest.fixture
def mouse_fasta(test_data_dir):
    return test_data_dir / "mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta"


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
def comet_psm(test_data_dir):
    path = test_data_dir / "example_comet_psm.json"
    psms = CometPSM.from_txt(
        txt="tests/data/BMEM_AspN_Fxn4/assign-confidence.target.txt"
    )
    psms = [psm for psm in psms if "sp|P99027|RLA2_MOUSE" in psm.proteins]
    psms[0].save(path=path)
    return CometPSM.load(path=path)


def default_test_config(test_data_dir: Path, out_dir: Path) -> HypedsearchRunConfig:
    db_path = (
        test_data_dir / "sp-P99027-RLA2_MOUSE_mzml=BMEM_AspN_Fxn4;scan=7_kmer_db.db"
    )
    return HypedsearchRunConfig(
        mzml_to_scans={test_data_dir / "spectra/BMEM_AspN_Fxn4_scans1-20.mzML": "all"},
        parent_out_dir=out_dir,
        hybrid_former=HybridFormer(
            kmer_db=db_path,
            fasta=MOUSE_PROTEOME,
        ),
        psm_scorer=HybridPSMScorer(
            comet_params=test_data_dir / "comet.params",
        ),
    )


def default_hs_run(test_data_dir: Path, out_dir: Path):
    # Arrange
    hs_config = default_test_config(
        test_data_dir=test_data_dir,
        out_dir=out_dir,
    )
    config_path = out_dir / "config.json"
    hs_config.to_json(path=config_path)
    cmd_parts = ["./src/run_hypedsearch.sh", f"--config {config_path}", f"--cores 8"]
    cmd = " ".join(cmd_parts)
    _ = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        shell=True,
    )
    return hs_config


def create_comet_psm(test_data_dir):
    psms = CometPSM.from_txt(
        txt="results/hs_mouse_samples/native_run/assign-confidence.txt"
    )
    best_psm = max(psms, key=lambda obj: obj.xcorr)
    best_psm.save()
