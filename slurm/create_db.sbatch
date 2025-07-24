#!/bin/bash    
#SBATCH --partition=short
#SBATCH --job-name=create-db
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=75gb               
#SBATCH --time=23:00:00
#SBATCH --output=logs/create_db/%j.out
#SBATCH --error=logs/create_db/%j.err           


# For testing
# python -m src.create_db \
#     -f tests/data/test.fasta \
#     -ktp tmp/create_db_cli/test.fasta.json \
#     -ktp tmp/create_db_cli/test.fasta.pklz \
#     -d tmp/create_db_cli/test.fasta.db 

# For actuals
# python -m src.create_db \
#     -f fastas/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fasta \
#     -ktp fastas/db/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fasta.pklz \
#     -d fastas/db/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fastas.db 


python -m src.create_db \
    -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
    -ktp fastas/db/SwissProt.TAW_mouse_w_NOD_IAPP.fasta.pklz \
    -d fastas/db/SwissProt.TAW_mouse_w_NOD_IAPP.fasta.db 
