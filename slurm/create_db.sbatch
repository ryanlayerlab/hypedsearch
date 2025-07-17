#!/bin/bash    
#SBATCH --partition=short
#SBATCH --job-name=create-db
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=75gb               
#SBATCH --time=23:00:00
#SBATCH --output=logs/create_db/%j.out
#SBATCH --error=logs/create_db/%j.err           

python -m src.create_db -d tmp/create_db_cli/test.fasta.db -f tests/data/test.fasta -ktp tmp/create_db_cli/test.fasta.json