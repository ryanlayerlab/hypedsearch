#!/bin/bash
#SBATCH --job-name=run-hs
#SBATCH --output=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.out
#SBATCH --error=/scratch/Shares/layer/hypedsearch/hypedsearch/slurm/logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --mail-type=NONE
#SBATCH --mail-user=YOU@YOUREMAIL.COM

python 