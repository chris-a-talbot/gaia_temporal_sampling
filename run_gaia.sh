#!/bin/bash
#SBATCH --job-name=run_gaia
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --time=3:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/processing_logs/run_gaia_sh_logs/%x-%A_%a.log

set -e

# Load required modules
module load R/4.4.0
module load gcc/13.2.0

# Set working directory
CWD=$(pwd)

# Create log directories if they don't exist
mkdir -p ./logs/processing_logs/run_gaia_sh_logs

echo "Starting tree processing on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "R user library path: $R_LIBS_USER"

# Now run the main R script
Rscript run_gaia.R "$CWD"

echo "Completed tree processing at $(date)"
