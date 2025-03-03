import os
import subprocess

import numpy as np

# Print current working directory
# cwd = os.getcwd()
# print(f"Current working directory: {cwd}")

# for num_peptides in [10, 20]:
for num_peptides in [2000, 4000, 8000, 16000]:
    # for num_peptides in np.arange(, 21001, 500):
    cmd = ["sbatch", "slurm/create_protein_product_ion_dbs.sbatch", str(num_peptides)]
    subprocess.Popen(cmd)  # Runs the command without waiting
