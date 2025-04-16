#!/bin/bash    
#SBATCH -p short              # Partition or queue. In this case, short!
#SBATCH --job-name=create-kmer-db    # Job name
#SBATCH --mail-type=NONE               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOU@YOUREMAIL.COM
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=50gb                   # Memory limit
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/slurm_%j.out   # Standard output and error log
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/slurm_%j.err   # %j inserts job number
 
# Setup
pwd
# source ./constants.sh
# cd $REPO_DIR
# K=50
# FASTA="Uniprot_mouse.fasta"
# # Testing
# # K=50
# # FASTA="test.fasta"

# python -m src.create_kmer_db -f "fastas/${FASTA}" -k $K -c 1 -d "dbs/kmer_db_fasta=${FASTA}_max_k=${K}_charges=[1].db"