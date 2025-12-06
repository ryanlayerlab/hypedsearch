import subprocess
import sys

from src.constants import RUN_HYPEDSEARCH_SMK
from tests.unit.test_hypedsearch import default_config


class Test_run_via_snakemake:
    @staticmethod
    def test_smoke(tmp_path, test_data_dir):
        hs_config = default_config(
            test_data_dir=test_data_dir,
            out_dir=tmp_path,
        )
        config_path = tmp_path / "hs_config.json"
        hs_config.to_json(path=config_path)
        cmd = f"snakemake -s {RUN_HYPEDSEARCH_SMK} --config hs_config={config_path} --cores 8 --scheduler greedy"
        result = subprocess.run(
            cmd,
            stdout=sys.stdout,
            stderr=sys.stderr,
            text=True,
            shell=True,
        )
        assert result.returncode == 0
        mzml = list(hs_config._mzml_to_scans.keys())[0]
        for scan in hs_config._mzml_to_scans[mzml]:
            comet_txt = (
                hs_config.hybrid_run_scan_results_dir
                / f"{mzml.stem}.comet.{scan}-{scan}.target.txt"
            )
            assert comet_txt.exists()
