#!/bin/bash


# Default values
MEM="500GB"
CORES=180
LOG_DIR="logs/run_commands_via_slurm"
TIME="24:00:00"
PARTITION="highmem"

show_help() {
  echo "Usage: $0 [--options]"
  echo
  echo "Options:"
  echo "  --cmd    Command to run via SLURM (required)"
  echo "  --name    Name (required)"
  echo "  --time    Time to ask SLURM for (default: $TIME)"
  echo "  --mem    Memory allocation for SLURM job (default: $MEM)"
  echo "  --part   Partition to use (required; recommend 'highmem'; default: $PARTITION)"
  echo "  --cores   Number of cores to ask SLURM for (default: $CORES)"
  echo "  -h, --help    Show this help message"
}

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cmd) CMD="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --name) NAME="$2"; shift ;;
        --time) TIME="$2"; shift ;;
        --part) PARTITION="$2"; shift ;;
        --cores) CORES="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown parameter: $1"; exit 0 ;;
    esac
    shift
done

cat <<EOF
Running Hypedsearch via SLURM with:
- NAME: $NAME
- PARTITION: $PARTITION
- CORES: $CORES
- TIME: $TIME
- MEM: $MEM
- LOG_FILES: $LOG_DIR/$NAME.out and $LOG_DIR/$NAME.err
- CMD: $CMD
EOF


sbatch <<EOT
#!/bin/bash

#SBATCH --job-name="$NAME"
#SBATCH --mem=$MEM
#SBATCH --ntasks=$CORES
#SBATCH --partition="$PARTITION"
#SBATCH --nodes=1
#SBATCH --time=$TIME
#SBATCH --output=$LOG_DIR/$NAME.out
#SBATCH --error=$LOG_DIR/$NAME.err 
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=erjo3868@colorado.edu

$CMD

EOT