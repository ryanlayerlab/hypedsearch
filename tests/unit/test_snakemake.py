import subprocess
import sys

from src.constants import RUN_HYPEDSEARCH_SMK
from src.hypedsearch import HypedsearchRunConfig
from src.utils import load_json, to_json


class Test_run_via_snakemake:
    @staticmethod
    def test_smoke(tmp_path, test_hs_config_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config_path = tmp_path / "hs.config.json"
        config = HypedsearchRunConfig(**data)
        expected_outputs = list(config.missing_hybrid_run_scan_target_txts)
        to_json(data=data, path=config_path)
        cmd = f"snakemake -s {RUN_HYPEDSEARCH_SMK} --config hs_config={config_path} --cores 8 --scheduler greedy"
        result = subprocess.run(
            cmd,
            stdout=sys.stdout,
            stderr=sys.stderr,
            text=True,
            shell=True,
        )
        assert result.returncode == 0
        for output in expected_outputs:
            assert output.exists()
