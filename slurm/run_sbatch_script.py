import os
import subprocess

import numpy as np

# Print current working directory
# cwd = os.getcwd()
# print(f"Current working directory: {cwd}")

# Run a command-line program (example: `ls` on Linux/macOS, `dir` on Windows)
command = "ls"  # Change this based on your OS and needs

import os
import subprocess

# Print current working directory
cwd = os.getcwd()
print(f"Current working directory: {cwd}")

# Run a command-line program (example: `ls` on Linux/macOS, `dir` on Windows)
command = "ls"  # Change this based on your OS and needs

for num_peptides in [10, 20]:
    # for num_peptides in [10, 50, 100, 250, 500, 1000, 1500, 2000, 4000, 8000, 16000]:
    # for num_peptides in np.arange(, 21001, 500):
    cmd = ["sbatch", "slurm/argument_accepting.sbatch", str(num_peptides)]
    subprocess.Popen(cmd)  # Runs the command without waiting
