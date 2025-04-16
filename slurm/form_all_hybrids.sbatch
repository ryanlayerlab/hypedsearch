#!/bin/bash    
#SBATCH --partition=short
#SBATCH --job-name=form-all-hybrids   
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=5gb               
#SBATCH --time=00:05:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null            

source slurm/slurm_fcns.sh
out_dir="slurm/logs/$SLURM_JOB_NAME"
mkdir -p $out_dir

redirect_output $out_dir

if is_slurm_array_job; then
    python -m src.form_all_hybrids \
        -pn results/new_fasta/top_10_proteins.txt \
        -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
        -m data/spectra/BMEM_AspN_Fxn5.mzML \
        -pt 20 \
        -s $SLURM_ARRAY_TASK_ID \
        -o tmp
else
    python -m src.form_all_hybrids \
        -pn results/new_fasta/top_10_proteins.txt \
        -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
        -m data/spectra/BMEM_AspN_Fxn5.mzML \
        -pt 20 \
        -s 9 \
        -o tmp
fi

remove_empty_file $ERRFILE
