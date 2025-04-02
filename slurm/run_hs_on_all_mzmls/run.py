import sys
from pathlib import Path

repo_dir = Path(__file__).parents[2]
sys.path.append(str(repo_dir))

import subprocess

from src.constants import SPECTRA_DIR

mzml_paths = list(SPECTRA_DIR.glob("*.mzML"))
mzml_paths = mzml_paths[:2]
for mzml in mzml_paths:
    cmd = [
        "sbatch",
        "slurm/run_hs_on_all_mzmls/run_on_one_mzml.sbatch",
        mzml,
    ]
    subprocess.Popen(cmd)  # Runs the command without waiting for it to finish
