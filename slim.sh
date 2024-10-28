#!/bin/bash
#SBATCH --job-name=gaia_temporal_testing
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --array=1-10
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/processing_logs/slim_sh_logs/%x-%A_%a.log

set -e

CWD=$(pwd)

module load Bioinformatics
module load slim

echo "Starting simulation REP=${SLURM_ARRAY_TASK_ID} on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"

slim -d "PWD='$CWD'" -d "REP=${SLURM_ARRAY_TASK_ID}" run-pareto.slim > "${CWD}/slim_sh_logs/simulation_rep_${SLURM_ARRAY_TASK_ID}.log" 2>&1

echo "Completed simulation REP=${SLURM_ARRAY_TASK_ID} at $(date)"
