#!/bin/bash -l

# modified for maintenance window
#SBATCH --time=24:00:00

# original:
# #SBATCH --time=96:00:00

#SBATCH --ntasks=1
#SBATCH --mem=50g
#SBATCH --tmp=1g

source /projects/standard/hsiehph/shared/bin/initialize_conda.sh
conda activate fastcn2


set -euo pipefail

mkdir -p logs

export PATH=$PWD/bin:$PATH



snakemake \
    --use-conda \
    --configfile config/config.yaml \
    --local-cores 20 \
    --cores 20 \
    --max-inventory-time 10000 \
    --resources load=1000 \
    --scheduler greedy \
    --latency-wait 60 \
    --restart-times 3 \
    --jobs 100 \
    --profile profile \
    "$@"


#    --local-cores 100 \
#    --cores 100 \
#    --max-inventory-time 10000 \
#    --resources load=1000 \
#    --scheduler greedy \
#    --latency-wait 60 \
#    --restart-times 1 \
#    --rerun-incomplete \
