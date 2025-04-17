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

# Parse command-line arguments
FASTA=$1
PROTEINS=$2
MZML=$3
HS_OUTPUTDIR=$4
echo -e "fasta=${FASTA}\nproteins=${PROTEINS}\nmzml=${MZML}\nout dir=${HS_OUTPUTDIR}"

# Redirect Slurm output to a specific folder
out_dir="slurm/logs/$SLURM_JOB_NAME"
mkdir -p $out_dir
redirect_output $out_dir
echo "Redirected SLURM output to $out_dir"

# Copy needed files to local machine
TMP_DIR="/tmp"
FASTA_TMP="$TMP_DIR/$(basename "$FASTA")"
PROTEINS_TMP="$TMP_DIR/$(basename "$PROTEINS")"
MZML_TMP="$TMP_DIR/$(basename "$MZML")"

cp "$FASTA" "$FASTA_TMP"
cp "$PROTEINS" "$PROTEINS_TMP"
cp "$MZML" "$MZML_TMP"
echo "Done copying files to local machine!"

# Run python script
if is_slurm_array_job; then
    print_job_array_stuff
    python -m src.hybrids_via_all_hybrids \
        --protein_names $PROTEINS_TMP \
        --fasta_path $FASTA_TMP \
        --mzml $MZML_TMP \
        --precursor_mz_ppm_tol 20 \
        --scan $SLURM_ARRAY_TASK_ID \
        --output_dir $HS_OUTPUTDIR
else
    echo "Expected an array job! Nothing happened"
fi

remove_empty_file $ERRFILE
