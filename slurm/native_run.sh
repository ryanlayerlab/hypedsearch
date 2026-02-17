#!/bin/bash

# Default values
MEM="200GB"
PARTITION="highmem"
CORES=192

# Parse flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --name) NAME="$2"; shift ;;
        --config) CONFIG="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --partition) PARTITION="$2"; shift ;;
        --cores) CORES="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$NAME" || -z "$CONFIG" ]]; then
    echo "Error: You must provide both --name and --config"
    exit 1
fi

# mkdir -p logs/$NAME

sbatch <<EOT
#!/bin/bash

#SBATCH --job-name="$NAME"
#SBATCH --mem=$MEM
#SBATCH --ntasks=$CORES
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/native_run/$NAME.out
#SBATCH --error=logs/native_run/$NAME.err           

snakemake -s snakefiles/run_comet.smk \
    --configfile $CONFIG \
    --cores $CORES \
    --use-singularity \
    --keep-going \
    --rerun-incomplete \
    --nolock \
    --retries 5 || true

EOT