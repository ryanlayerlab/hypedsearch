import sys
from pathlib import Path

repo_dir = Path(__file__).parents[2]
assert repo_dir.name == "hypedsearch"
sys.path.append(str(repo_dir))

import subprocess

from src.comet_utils import CometPSM
from src.constants import SPECTRA_DIR

# mzml_paths = list(SPECTRA_DIR.glob("*.mzML"))

# mzml_paths = [SPECTRA_DIR / "BMEM_AspN_Fxn5.mzML"]  # for testing

df = CometPSM.from_txt(
    txt_path=Path("results/new_fasta/BMEM_AspN_Fxn5.txt"), as_df=True
)
scans = list(df["scan"].unique())
print(f"Number of total scans = {len(scans)}")

mzml = (SPECTRA_DIR / "BMEM_AspN_Fxn5.mzML").absolute()
for scan in scans[0:100]:
    # for scan in scans[1000:2000]:
    # for scan in scans[2000:3000]:
    # for scan in scans[3000:]:
    cmd = [
        "sbatch",
        "slurm/run_hs_on_all_mzmls/run_on_one_mzml.sbatch",
        mzml,
        scan,
    ]
    subprocess.Popen(cmd)  # Runs the command without waiting for it to finish
