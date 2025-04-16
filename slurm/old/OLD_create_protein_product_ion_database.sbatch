#!/bin/bash    
#SBATCH -p long              # Partition or queue. In this case, short!
#SBATCH --job-name=create-protein-product-ion-db    # Job name
#SBATCH --mail-type=NONE               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOU@YOUREMAIL.COM
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=50gb                   # Memory limit
#SBATCH --time=23:59:59               # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/slurm_%j.out   # Standard output and error log
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/slurm_%j.err   # %j inserts job number
 
# TESTING
# python create_protein_product_ion_database.py -f tests/test.fasta -d tmp -k 30 -c 1 -c 2 -c 3

# FOR REALS
python create_protein_product_ion_database.py -f fastas/Uniprot_mouse.fasta -d dbs -k 30 -c 1 -c 2 -c 3