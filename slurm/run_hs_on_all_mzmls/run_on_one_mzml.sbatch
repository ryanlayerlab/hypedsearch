#!/bin/bash    

#SBATCH -p short              # Partition or queue. In this case, short!
#SBATCH --job-name=run-hs   # Job name
#SBATCH --mail-type=NONE               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOU@YOUREMAIL.COM
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=5gb                   # Memory limit
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
# Standard output and error log
# %j inserts job ID, %x inserts job name
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.out 
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.err 

# python hypedsearch.py \
# --mzml_dir data/spectra \
# --mzml_path $1 \
# --output_dir results/april-1-run \
# --db_path results/april-1-run/test.db \
# --top_n_proteins 50 \
# --num_peaks 100 \
# --comet_exe_path comet/comet.linux.exe \
# --comet_params_path comet/comet.params

python -m src.run_hypedsearch \
    -d results/new_fasta/top_10_prots.db \
    -o results/new_fasta/5 \
    -m $1 \
    -s $2 \
    -n 0 -p 10 -P 20 -ow True -kf True 