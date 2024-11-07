#!/bin/bash
#SBATCH --job-name=simplify_trees
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=24G
#SBATCH --time=3:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/processing_logs/simplify_trees_sh_logs/%x-%A_%a.log

set -e

# Load required modules
module load python/3.11

# Activate virtual environment
source ~/venv/tskit_env/bin/activate

# Set working directory
CWD=$(pwd)

# Create log directory if it doesn't exist
mkdir -p .logs/processing_logs/simplify_trees_sh_logs

echo "Starting tree simplification on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Number of tree files (CPUs): $SLURM_CPUS_PER_TASK"

# Run the Python script
python simplify_trees_cont_unif.py \
    --trees-dir ./trees \
    --logs-dir ./logs/SLiM_logs \
    --num-processes $SLURM_CPUS_PER_TASK

echo "Completed tree simplification at $(date)"
