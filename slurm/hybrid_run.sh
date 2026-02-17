#!/bin/bash

# Default values
RUN_TYPE="hybrid"
MEM="500GB"
CORES=192
PARTITION="highmem" 
# FIJI_NODE="fijinode-68"
DATA_DIR="/localscratch"

# Parse flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --name) NAME="$2"; shift ;;
        --config) HS_CONFIG="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --partition) PARTITION="$2"; shift ;;
        --cores) CORES="$2"; shift ;;
        --run) RUN_TYPE="$2"; shift ;;
        # --node) FIJI_NODE="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$NAME" || -z "$HS_CONFIG" ]]; then
    echo "Error: You must provide both a name (via '--name') and a Hypedsearch JSON config (via '--config')"
    return
fi

sbatch <<EOT
#!/bin/bash

#SBATCH --job-name="$NAME"
#SBATCH --mem=$MEM
#SBATCH --ntasks=$CORES
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --time=24:00:00
# #SBATCH --nodelist=$FIJI_NODE
#SBATCH --output=logs/$RUN_TYPE/$NAME.out
#SBATCH --error=logs/$RUN_TYPE/$NAME.err 

if [[ "$RUN_TYPE" == "native" ]]; then
    # Create native run snakemake config
    NATIVE_CONFIG="tmp/${NAME}.native_run.smk.json"
    python -m src.hypedsearch create-native-run-snakemake-config \
        --config $HS_CONFIG --out_path=$NATIVE_CONFIG

    # Native run via snakemake
    snakemake -s snakefiles/run_comet.smk \
        --configfile $NATIVE_CONFIG \
        --cores 80 \
        --scheduler greedy \
        --use-singularity \
        --keep-going \
        --rerun-incomplete \
        --nolock \
        --retries 5 || true

elif [[ "$RUN_TYPE" == "hybrid" ]]; then
    # Move files to a directory that's on the node itself for faster I/O
    python -m fiji.fiji_utils prep-files-on-fiji \
        --config $HS_CONFIG --data_dir $DATA_DIR

    # Hybrid run
    snakemake -s snakefiles/run_hypedsearch.smk \
        --config hs_config=$DATA_DIR/$NAME/hs.fiji.config.json \
        --cores 80 \
        --scheduler greedy \
        --use-singularity \
        --singularity-args "-B $DATA_DIR:$DATA_DIR" \
        --keep-going \
        --rerun-incomplete \
        --nolock \
        --retries 5 || true
else:
    echo "Unrecognized run type: $RUN_TYPE. Expected 'native' or 'hybrid'."
    return
fi

EOT