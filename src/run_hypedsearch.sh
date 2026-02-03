#!/bin/bash

# Default values
USE_SINGULARITY=false

show_help() {
  echo "Usage: $0 [--options]"
  echo
  echo "Options:"
  echo "  --config    Path to Hypedsearch JSON config. Required"
  echo "  --data    Path to data directory. Defaults to a temporary directory created by 'mktemp -d'."
  echo "  --cores    Number of cores to give to snakemake. Required"
  echo "  --singularity    Flag that sets whether or not to use sinuglarity. Default: $USE_SINGULARITY"
  echo "  -h, --help    Show this help message"
}

# Argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            HS_CONFIG="$2"
            shift 2
            ;;
        --data)
            DATA_DIR="$2"
            shift 2
            ;;
        --singularity)
            USE_SINGULARITY=true
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1"
            show_help
            exit 0
            ;;
    esac
done

# Validation
if [ -z "$HS_CONFIG"  ]; then
    echo "Error: You must provide both a Hypedsearch JSON config via '--config'"
    exit 0
fi
if [ -z "$CORES"  ]; then
    echo "Error: You must provide the number of cores to use via '--cores'"
    exit 0
fi
# Parse the name from Hypedsearch config
NAME=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['name'])" "$HS_CONFIG")

if [ -z "$DATA_DIR" ]; then
    DATA_DIR="$(mktemp -d)"
    echo "No data directory provided. Using temporary directory: $DATA_DIR"
fi

cat <<EOF
Running Hypedsearch with:
- NAME: $NAME
- HS_CONFIG: $HS_CONFIG
- DATA_DIR: $DATA_DIR
- CORES: $CORES
- USE_SINGULARITY: $USE_SINGULARITY
EOF


# Move files to a directory that's on the node itself for faster I/O
if $USE_SINGULARITY; then
    python -m fiji.fiji_utils prep-files-on-fiji \
        --config $HS_CONFIG \
        --data_dir $DATA_DIR
fi

# Run Hypedsearch
if $USE_SINGULARITY; then
    snakemake -s snakefiles/run_hypedsearch.smk \
        --config hs_config=$DATA_DIR/$NAME/hs.fiji.config.json \
        --cores $CORES \
        --scheduler greedy \
        --use-singularity \
        --singularity-args "-B $DATA_DIR:$DATA_DIR" \
        --keep-going \
        --rerun-incomplete \
        --nolock \
        --retries 5 || true
else
    snakemake -s snakefiles/run_hypedsearch.smk \
        --config hs_config=$HS_CONFIG \
        --cores $CORES \
        --scheduler greedy \
        --keep-going \
        --rerun-incomplete \
        --nolock \
        --retries 5 || true
fi
if $USE_SINGULARITY; then
    # Move scan results back to persistent storage
    echo "Moving scan results back to persistent storage"
    python -m fiji.fiji_utils move-scan-results \
        --config $HS_CONFIG \
        --data_dir $DATA_DIR

    # Clean up data directory on the node
    echo "Removing data directory $DATA_DIR/$NAME"
    rm -rf $DATA_DIR/$NAME
fi

echo "Finished running Hypedsearch!"