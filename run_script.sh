#!/bin/bash -l

# modified for maintenance window
# #SBATCH --time=34:00:00

# original:
#SBATCH --time=96:00:00

#SBATCH --ntasks=1
#SBATCH --mem=50g
#SBATCH --tmp=1g


conda activate snakemake


set -euo pipefail

mkdir -p logs

export PATH=$PWD/bin:$PATH



snakemake \
    --use-conda \
    -k \
    --configfile config/config.yaml \
    --jobs 100 \
    --profile profile \
    --latency-wait 60 \
    --restart-times 3 \
    "$@"


#    --local-cores 100 \
#    --cores 100 \
#    --max-inventory-time 10000 \
#    --resources load=1000 \
#    --scheduler greedy \
#    --latency-wait 60 \
#    --restart-times 1 \
#    --rerun-incomplete \
