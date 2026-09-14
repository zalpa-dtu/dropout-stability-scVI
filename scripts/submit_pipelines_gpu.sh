#!/bin/bash

set -e

# Usage:
# bash scripts/submit_pipelines_gpu.sh DATAPATH OUTPUTPATH SEED_LIST UNIQUE_GROUPS RESULTS_DIR

DATAPATH="$1"
OUTPUTPATH="$2"
SEED_LIST="$3"
UNIQUE_GROUPS="$4"
RESULTS_DIR="$5"

# Convert the comma-separated list into an array.
IFS=',' read -r -a SEEDS <<< "${SEED_LIST}"

mkdir -p "${OUTPUTPATH}"

PROJECT_DIR="${PROJECT_DIR:-/dtu-compute/digitstem/zuzanna/zuzanna_dropout_project}"
PYTHON="${PROJECT_DIR}/.conda-env/bin/python"

for seed in "${SEEDS[@]}"
do
  bsub -q gpuv100 \
       -o ${OUTPUTPATH}/Output_DropoutMan_${seed}_%J.out \
       -e ${OUTPUTPATH}/Error_DropoutMan_${seed}_%J.err \
       -J DropoutMan_${seed} \
       -n 4 \
       -R "span[hosts=1]" \
       -R "rusage[mem=8GB]" \
       -gpu "num=1:mode=exclusive_process" \
       -W 21:00 \
       "cd \"${PROJECT_DIR}\" && \
        \"${PYTHON}\" run_exps.py \
        \"path=${DATAPATH}\" \
        seed=${seed} \
        unique_groups=${UNIQUE_GROUPS} \
        \"hydra.run.dir=${RESULTS_DIR%/}/seed_${seed}_job_\${LSB_JOBID}\""
done