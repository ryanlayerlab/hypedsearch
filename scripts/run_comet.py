import sys
from pathlib import Path

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from src.comet_utils import run_comet_on_one_mzml
from src.constants import (
    COMET_EXECUTABLE,
    COMET_PARAMS,
    COMET_RUN_1_DIR,
    SPECTRA_DIR,
    THOMAS_SAMPLES,
)
from src.utils import setup_logger

logger = setup_logger()


mzmls = list(SPECTRA_DIR.glob("*.mzML"))

logger.info("About to run Comet on ")
for mzml_path in mzml_paths:
    run_comet_and_save_params_file(
        mzml_path=mzml_path, comet=comet, parent_output_dir=parent_output_dir
    )
