import sys
from pathlib import Path

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from src.comet_utils import CometExe, run_comet_and_save_params_file
from src.constants import (
    COMET_EXECUTABLE,
    COMET_PARAMS,
    COMET_RUN_1_DIR,
    SPECTRA_DIR,
    THOMAS_SAMPLES,
)
from src.utils import setup_logger

setup_logger()

comet = CometExe(exe=COMET_EXECUTABLE, params=COMET_PARAMS)
parent_output_dir = COMET_RUN_1_DIR


mzml_paths = [SPECTRA_DIR / f"{sample}.mzML" for sample in THOMAS_SAMPLES]

for mzml_path in mzml_paths:
    run_comet_and_save_params_file(
        mzml_path=mzml_path, comet=comet, parent_output_dir=parent_output_dir
    )
