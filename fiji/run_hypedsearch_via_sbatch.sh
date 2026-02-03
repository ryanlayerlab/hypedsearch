#!/bin/bash


# Default values
MEM="500GB"
SLURM_CORES=180
SMK_CORES=80
DATA_DIR="/localscratch"
LOG_DIR="logs/hypedsearch"

show_help() {
  echo "Usage: $0 [--options]"
  echo
  echo "Options:"
  echo "  --config    Path to Hypedsearch JSON config (required)"
  echo "  --mem    Memory allocation for SLURM job (default: $MEM)"
#   echo "  --data   Path to the directory on the node where files will be temporarily copied for I/O optimization (default: $DATA_DIR)"
  echo "  --part   Partition to use (required; recommend 'highmem')"
  echo "  --slurm-cores   Number of cores to ask SLURM for (default $SLURM_CORES)"
  echo "  --smk-cores   Number of cores to provide to Snakemake (default $SMK_CORES)"
  echo "  -h, --help    Show this help message"
}

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --config) HS_CONFIG="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --part) PARTITION="$2"; shift ;;
        --slurm-cores) SLURM_CORES="$2"; shift ;;
        --smk-cores) SMK_CORES="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        # --node) FIJI_NODE="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 0 ;;
    esac
    shift
done

NAME=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['name'])" "$HS_CONFIG")

# Validation
if [ -z "$HS_CONFIG" ]; then
    echo "Error: You must provide both a Hypedsearch JSON config via '--config'"
    exit 0
fi
if [ -z "$PARTITION" ]; then
    echo "Error: You must provide a partition for the SLURM job via '--part'. Recommend 'highmem'"
    exit 0
fi

# Parse the name from Hypedsearch config

cat <<EOF
Running Hypedsearch via SLURM with:
- NAME: $NAME
- HS_CONFIG: $HS_CONFIG
- PARTITION: $PARTITION
- SLURM_CORES: $SLURM_CORES
- SMK_CORES: $SMK_CORES
- MEM: $MEM
- LOG_FILES: $LOG_DIR/$NAME.out and $LOG_DIR/$NAME.err
EOF


sbatch <<EOT
#!/bin/bash

#SBATCH --job-name="$NAME"
#SBATCH --mem=$MEM
#SBATCH --ntasks=$SLURM_CORES
#SBATCH --partition="$PARTITION"
#SBATCH --nodes=1
#SBATCH --time=24:00:00
# #SBATCH --nodelist=$FIJI_NODE
#SBATCH --output=$LOG_DIR/$NAME.out
#SBATCH --error=$LOG_DIR/$NAME.err 
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=erjo3868@colorado.edu

./src/run_hypedsearch.sh \
    --config $HS_CONFIG \
    --cores $SMK_CORES \
    --data $DATA_DIR \
    --singularity

EOT