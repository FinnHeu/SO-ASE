#!/bin/bash
#SBATCH --job-name=fesom_postproc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --account=ab0995
#SBATCH --output=logs/postprocessing/cav_icb/postproc_%a.out
#SBATCH --error=logs/postprocessing/cav_icb/postproc_%a.err

# ============================================================================
#                    FESOM2 Postprocessing Job Script
# ============================================================================
# Array job where each task processes one year.
# The array task ID IS the year to process.
#
# Usage:
#   sbatch --array=1850-1900 job_fesom_postprocessing.sh
#   sbatch --array=1850-1900%20 job_fesom_postprocessing.sh  # max 20 concurrent
#
# ============================================================================

# ============================================================================
#                         CONFIGURATION
# ============================================================================

# Default year (used if not running as array job)
YEAR_START=1850

# Conda environment path
CONDA_ENV="/home/a/a270186/.conda/envs/so_ase/"

# Python script path
PYTHON_SCRIPT="/home/a/a270186/python_modules/SO-ASE/workflows_fesom/workflows_output/workflows_postprocessing/python/fesom_postprocessing_FLUX-OFF.py"

# Log directory
LOG_DIR="./logs/postprocessing/FLUX-OFF/"

# ============================================================================
#                         SETUP
# ============================================================================

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Get year directly from array task ID (use --array=1850-1900 etc.)
YEAR=${SLURM_ARRAY_TASK_ID:-$YEAR_START}

echo "============================================================================"
echo "FESOM2 Postprocessing Job"
echo "============================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing Year: ${YEAR}"
echo "Node: $(hostname)"
echo "Start Time: $(date)"
echo "============================================================================"

# ============================================================================
#                         RUN POSTPROCESSING
# ============================================================================

${CONDA_ENV}/bin/python ${PYTHON_SCRIPT} ${YEAR}

EXIT_CODE=$?

echo "============================================================================"
echo "End Time: $(date)"
echo "Exit Code: ${EXIT_CODE}"
echo "============================================================================"

exit ${EXIT_CODE}
