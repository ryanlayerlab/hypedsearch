import shlex
import subprocess

from src.constants import GIT_REPO_DIR

mzml = (GIT_REPO_DIR / "data/spectra/BMEM_AspN_Fxn5.mzML").absolute()
protein_names = (GIT_REPO_DIR / "tests/unit/data/top_10_proteins.txt").absolute()
fasta_path = (GIT_REPO_DIR / "fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta").absolute()
output_dir = (GIT_REPO_DIR / "/tmp/test").absolute()

scans = 9

cmd = f"sbatch slurm/form_all_hybrids.sbatch --array=9 {fasta_path} {protein_names} {mzml} {output_dir}"
run_args = shlex.split(cmd)

result = subprocess.run(run_args, capture_output=True, text=True)
print(f"stdout = {result.stdout}")
print(f"stderr = {result.stderr}")
print(f"Return code = {result.returncode}")
