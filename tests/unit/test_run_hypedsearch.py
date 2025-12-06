"""
Tests for run_hypedsearch.sh
"""

from src.utils import CmdLineRunner
from tests.conftest import default_config


@staticmethod
def test_smoke(tmp_path, test_data_dir):
    # Arrange
    config_path = tmp_path / "hs_config.json"
    hs_config = default_config(test_data_dir=test_data_dir, out_dir=tmp_path)
    expected_native_txt = (
        hs_config.native_run_dir / "BMEM_AspN_Fxn4_scans1-20.comet.target.txt"
    )
    expected_hybrid_txt = (
        hs_config.hybrid_run_dir / "BMEM_AspN_Fxn4_scans1-20.comet.target.txt"
    )
    # # Act
    hs_config.to_json(path=config_path)
    cmd = [
        "./src/run_hypedsearch.sh",
        f"--config {config_path}",
        f"--cores 1",
    ]
    result = CmdLineRunner.run_cmd(cmd=cmd)
    # Assert
    assert result.returncode == 0
    assert expected_native_txt.exists()
    assert expected_hybrid_txt.exists()
