import os
import shutil
import sys
from pathlib import Path

from slurm.fiji_utils import SbatchConfig, create_sbatch_script_to_run_hypedsearch
from src.constants import GIT_REPO_DIR, MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.crux import CometRun
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
    run_in_parallel,
)
from src.mass_spectra import Spectrum
from src.utils import from_pickle, load_json, setup_logger

logger = setup_logger()

n_cores = 80
test_data_dir = Path("tests/data")
mouse_mzml_path = test_data_dir / "BMEM_AspN_Fxn4_scans1-20.mzML"
test_hs_config_path = test_data_dir / "test.hs.config.json"
data = load_json(path=test_hs_config_path)
out_dir = Path(data["parent_output_dir"])
if out_dir.exists():
    logger.info(f"Deleting parent_output_dir: {out_dir}")
    shutil.rmtree(out_dir)
hs_config = HypedsearchRunConfig(**data)

sbatch_config = SbatchConfig(
    name=hs_config.name,
)
create_sbatch_script_to_run_hypedsearch(
    sbatch_config=sbatch_config,
    config=test_hs_config_path,
)
