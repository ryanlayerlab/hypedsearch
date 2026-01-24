#!/bin/bash


# Default values
MEM="500GB"
CORES=180
# FIJI_NODE="fijinode-68"
DATA_DIR="/localscratch"
LOG_DIR="logs/hypedsearch"

show_help() {
  echo "Usage: $0 [--options]"
  echo
  echo "Options:"
#   echo "  --name   Number of cores to use (required)"
  echo "  --config    Path to Hypedsearch JSON config (required)"
  echo "  --mem    Memory allocation for SLURM job (default: $MEM)"
#   echo "  --data   Path to the directory on the node where files will be temporarily copied for I/O optimization (default: $DATA_DIR)"
  echo "  --part   Partition to use (required; recommend 'highmem')"
  echo "  --cores   Number of cores to use (default $CORES)"
  echo "  -h, --help    Show this help message"
}

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        # --name) NAME="$2"; shift ;;
        --config) HS_CONFIG="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --part) PARTITION="$2"; shift ;;
        --cores) CORES="$2"; shift ;;
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
# if [ -z "$NAME" ]; then
#     echo "Error: You must provide a name for the SLURM job name via '--name'"
#     exit 0
# fi
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
- CORES: $CORES
- MEM: $MEM
- LOG_FILES: $LOG_DIR/$NAME.out and $LOG_DIR/$NAME.err
EOF


sbatch <<EOT
#!/bin/bash

#SBATCH --job-name="$NAME"
#SBATCH --mem=$MEM
#SBATCH --ntasks=$CORES
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
    --cores 80 \
    --data $DATA_DIR \
    --singularity

EOT