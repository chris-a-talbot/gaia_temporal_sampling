#!/bin/bash
#SBATCH --job-name=sample_loc_matrix
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/processing_logs/sample_loc_mat_logs/%x-%A_%a.log

set -e

# Load required modules
module load python/3.11

# Activate virtual environment
source ~/venv/tskit_env/bin/activate

# Set working directory
CWD=$(pwd)

# Create log directory if it doesn't exist
mkdir -p ./logs/processing_logs/sample_loc_mat_logs

# Create output directory for sample locations if it doesn't exist
mkdir -p ./sample_locations

echo "Starting sample location matrix generation on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"

# Install required Python packages if not already installed
pip install --quiet tskit pandas tqdm

# Run the Python script
python sample_loc_mat_from_trees.py

echo "Completed sample location matrix generation at $(date)"
