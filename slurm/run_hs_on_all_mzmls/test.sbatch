#!/bin/bash    

#SBATCH -p long              # Partition or queue. In this case, short!
#SBATCH --job-name=run-hs   # Job name
#SBATCH --mail-type=NONE               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOU@YOUREMAIL.COM
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=50gb                   # Memory limit
#SBATCH --time=23:00:00               # Time limit hrs:min:sec
# Standard output and error log
# %j inserts job ID, %x inserts job name
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.out   # Standard output and error log
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.err   # %j inserts job number

python hypedsearch.py \
--mzml_dir data/spectra \
--mzml_path data/spectra/BMEM_AspN_Fxn4.mzML \
--scan_num 7 \
--output_dir results/april-1-run \
--db_path test.db \
--top_n_proteins 50 \
--num_peaks 100 \
--comet_exe_path comet/comet.linux.exe \
--comet_params_path comet/comet.params