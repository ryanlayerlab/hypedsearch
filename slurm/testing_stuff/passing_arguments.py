import os
import subprocess

# Print current working directory
cwd = os.getcwd()
print(f"Current working directory: {cwd}")

# Run a command-line program (example: `ls` on Linux/macOS, `dir` on Windows)
command = "ls"  # Change this based on your OS and needs

for thing in [1, "apple", 1.567]:
    cmd = f"sbatch slurm/argument_accepting.sbatch {thing}"
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Print command output
    print(process.stdout)
