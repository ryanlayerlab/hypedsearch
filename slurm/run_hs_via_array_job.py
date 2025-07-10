import os
import shlex
import subprocess
import sys
from pathlib import Path

repo_dir = Path(__file__).parents[1]
# print(repo_dir)
sys.path.append(str(repo_dir))

from src.comet_utils import CometPSM
from src.constants import GIT_REPO_DIR, RESULTS_DIR
from src.utils import make_directory

ARRAY_BATCH_SIZE = 800

# Parse command line input
assert (
    len(sys.argv) == 2
), f"len(sys.argv) should be 2. Just passing the sample number. It's {len(sys.argv)}"
sample_num = int(sys.argv[1])

# Set up paths
mzml = (GIT_REPO_DIR / f"data/spectra/BMEM_AspN_Fxn{sample_num}.mzML").absolute()
protein_names = (GIT_REPO_DIR / "tests/unit/data/top_10_proteins.txt").absolute()
fasta_path = (GIT_REPO_DIR / "fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta").absolute()
output_dir = (GIT_REPO_DIR / "results/new_fasta/all_hybrids").absolute()
mzmls_comet_run_txt = (
    GIT_REPO_DIR / f"results/new_fasta/BMEM_AspN_Fxn{sample_num}.txt"
).absolute()

make_directory(output_dir)

for path in [mzml, protein_names, fasta_path, output_dir, mzmls_comet_run_txt]:
    assert (
        path.exists()
    ), f"Expecting the following path to exist but it doesn't:\n{path}"

# Get the scans to run
mzml_comet_run = CometPSM.from_txt(
    txt_path=mzmls_comet_run_txt, sample=mzml.stem, as_df=True
)
scans = list(mzml_comet_run["scan"].unique())
scans_formatted_as_strings = ",".join([str(scan) for scan in scans])
# print(scans)
# print(scans_formatted_as_strings)

# Run the sbatch script
cmd = f"sbatch --array={scans_formatted_as_strings}%{ARRAY_BATCH_SIZE} slurm/form_all_hybrids.sbatch {fasta_path} {protein_names} {mzml} {output_dir}"
run_args = shlex.split(cmd)
print(f"Running cmd:\n{cmd}")
result = subprocess.run(run_args, capture_output=True, text=True)
print(f"\tstdout = {result.stdout}")
print(f"\tstderr = {result.stderr}")
print(f"\treturn code = {result.returncode}")
