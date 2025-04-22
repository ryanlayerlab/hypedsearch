#!/bin/bash
#SBATCH --partition=highmem
#SBATCH --node=1
#SBATCH --ntasks=192
#SBATCH --mem=1536000mb
#SBATCH --time=23:00:00
#SBATCH --job-name=scan_batch

python results/new_fasta/slurm/run_hs_on_scans.py