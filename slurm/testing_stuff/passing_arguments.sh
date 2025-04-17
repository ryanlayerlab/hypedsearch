#!/bin/bash    

for item in apple banana cherry; do
    sbatch slurm/argument_accepting.sbatch "$item"
done