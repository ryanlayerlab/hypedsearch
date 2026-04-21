#!/bin/bash

CONFIG="tests/data/ex1.hypedsearch.config.json"
CORES=50
# Remove parent output directoyr that HS creates
rm -rf tests/outputs

# Run HS
echo "Native run!" 
python cli.py native-run -c $CONFIG -os
echo "Finished native run!"

echo "Hybrid run!"
python cli.py run-hypedsearch -c $CONFIG -n $CORES -os -p
echo "Finished hybrid run!"
