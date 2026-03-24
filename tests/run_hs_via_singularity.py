import os
import shutil
import sys
from pathlib import Path

from src.constants import GIT_REPO_DIR, MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.crux import CometRun
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
    run_hypedsearch,
)
from src.mass_spectra import Spectrum
from src.utils import from_pickle, load_json, setup_logger

logger = setup_logger()
if sys.platform == "darwin":
    logger.info("Running on Mac...")
    crux_path = MAC_CRUX_EXECUTABLE
    on_singularity = False
elif sys.platform == "linux":
    logger.info("Running on Linux...")
    crux_path = None
    on_singularity = True
else:
    raise ValueError(f"Unsupported platform: {sys.platform}")

test_data_dir = Path("tests/data")
mouse_mzml_path = test_data_dir / "BMEM_AspN_Fxn4_scans1-20.mzML"
test_hs_config_path = test_data_dir / "test.hs.config.json"
scan = 7
data = load_json(path=test_hs_config_path)

out_dir = Path(data["parent_output_dir"])
if out_dir.exists():
    logger.info(f"Deleting parent_output_dir: {out_dir}")
    shutil.rmtree(out_dir)
config = HypedsearchRunConfig(**data)

logger.info("Running CometRun.run_comet_and_keep_only_results")
run = CometRun(
    fasta=config.fasta,
    mzml=mouse_mzml_path,
    crux_comet_params=config.crux_comet_params,
    out_dir=config.native_run_dir,
    decoy_search=2,
    scan_min=scan,
    scan_max=scan,
)
process = run.run_comet_and_keep_only_results(
    on_singularity=on_singularity, crux_path=crux_path
)
assert process.returncode == 0
assert run.nonstandardized_comet_outputs.target.exists()
assert run.nonstandardized_comet_outputs.decoy.exists()

logger.info("Running hybrid_run_on_spectrum")
params = HybridRunParams(
    kmer_db=config.kmer_db_path,
    fasta=config.fasta,
    fasta_fm_index=from_pickle(path=config.fasta_fm_index),
    crux_comet_params=config.crux_comet_params,
    out_dir=config.hybrid_run_dir,
)
cmd_result, run = hybrid_run_on_spectrum(
    spectrum=Spectrum.get_spectrum(scan=scan, mzml=mouse_mzml_path),
    params=params,
    on_singularity=on_singularity,
    crux_path=crux_path,
)

logger.info("Running run_in_parallel...")
run_hypedsearch(
    config=test_hs_config_path,
    n_cores=4,
    on_singularity=on_singularity,
    crux_path=crux_path,
)
