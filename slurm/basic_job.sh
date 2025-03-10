#!/bin/bash    
#SBATCH -p long              # Partition or queue. In this case, short!
#SBATCH --job-name=get-unique-peptides    # Job name
#SBATCH --mail-type=NONE               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOU@YOUREMAIL.COM
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=50gb                   # Memory limit
#SBATCH --time=23:59:59               # Time limit hrs:min:sec
# Standard output and error log
# %j inserts job ID, %x inserts job name
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.out   # Standard output and error log
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.err   # %j inserts job number


python scripts/get_uniq_peptides.py