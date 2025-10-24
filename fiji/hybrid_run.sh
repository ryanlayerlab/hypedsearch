#!/bin/bash

# Default values
MEM="500GB"
CORES=192
PARTITION="sandbox"
FIJI_NODE="fijinode-62"
DATA_DIR="/localscratch"

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
#SBATCH --time=24:00:00
#SBATCH --nodelist=$FIJI_NODE
#SBATCH --output=logs/hybrid_run/$NAME.out
#SBATCH --error=logs/hybrid_run/$NAME.err           

# Move files to a directory that's on the node itself for faster I/O
python -m fiji.fiji_utils prep-hybrid-run -c $CONFIG -d $DATA_DIR

# Hybrid run
snakemake -s snakefiles/run_hypedsearch.smk \
    --configfile /localscratch/$NAME/hybrid_run/hybrid_run_config.json \
    --cores 80 \
    --scheduler greedy \
    --use-singularity \
    --singularity-args "-B $DATA_DIR:$DATA_DIR" \
    --keep-going \
    --rerun-incomplete \
    --nolock \
    --retries 5 || true

EOT