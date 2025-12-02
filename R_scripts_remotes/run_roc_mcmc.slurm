#!/bin/bash
#SBATCH --job-name=ROC_MCMC
#SBATCH --output=roc_mcmc_%j.out
#SBATCH --error=roc_mcmc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#
# ==============================================================================
# SLURM Job Template for Pure R ROC Model MCMC
#
# Modify the parameters below according to your analysis needs.
# The C-accelerated version typically needs:
#   - 1-2 hours for 10,000 samples with ~20,000 genes
#   - 8-16 GB RAM
#   - Single CPU is sufficient with C acceleration
#
# Usage:
#   1. Copy this file and edit parameters below
#   2. Submit with: sbatch run_roc_mcmc.slurm
# ==============================================================================

# Load required modules (adjust for your cluster)
module load R/4.3.1  # or your R version

# ==============================================================================
# USER CONFIGURATION - MODIFY THESE
# ==============================================================================

# Base directory (where the R scripts are)
SCRIPT_DIR="/path/to/Codon_bias_analysis/R_scripts_remotes"

# Input files
FASTA_FILE="/path/to/sequences.fa"
PHI_FILE="/path/to/observed_expression_multitissue.csv"  # Optional
DM_FILE="/path/to/initial_dM.csv"                        # Optional

# Output directory
OUTPUT_DIR="/path/to/output/run_$(date +%Y%m%d_%H%M%S)"

# MCMC parameters
N_SAMPLES=10000      # Number of posterior samples (after thinning)
THINNING=10          # Save every Nth iteration
WITH_PHI=true        # Include expression in likelihood (true/false)
FIX_DM=true          # Fix mutation parameters (true/false)

# ==============================================================================
# END OF USER CONFIGURATION
# ==============================================================================

echo "============================================"
echo "  ROC Model MCMC - SLURM Job"
echo "============================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE} MB"
echo "Started: $(date)"
echo "============================================"

# Change to script directory
cd "${SCRIPT_DIR}"

# Ensure C library is compiled (only on first run)
if [ ! -f "src/roc_likelihood.so" ]; then
    echo "C library not found, compiling..."
    bash setup_cluster.sh
fi

# Build command
CMD="Rscript ROC_model_pure_R.R \
    -i ${FASTA_FILE} \
    -o ${OUTPUT_DIR} \
    -s ${N_SAMPLES} \
    -t ${THINNING}"

# Add optional arguments
if [ -n "${PHI_FILE}" ] && [ -f "${PHI_FILE}" ]; then
    CMD="${CMD} --phi ${PHI_FILE}"
    
    if [ "${WITH_PHI}" = true ]; then
        CMD="${CMD} --with_phi"
    fi
fi

if [ -n "${DM_FILE}" ] && [ -f "${DM_FILE}" ]; then
    CMD="${CMD} --dM ${DM_FILE}"
    
    if [ "${FIX_DM}" = true ]; then
        CMD="${CMD} --fix_dM"
    fi
fi

echo ""
echo "Running command:"
echo "${CMD}"
echo ""

# Run the MCMC
eval ${CMD}

EXIT_CODE=$?

echo ""
echo "============================================"
echo "Job completed at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================"

# Copy SLURM logs to output directory
if [ -d "${OUTPUT_DIR}" ]; then
    cp "roc_mcmc_${SLURM_JOB_ID}.out" "${OUTPUT_DIR}/" 2>/dev/null || true
    cp "roc_mcmc_${SLURM_JOB_ID}.err" "${OUTPUT_DIR}/" 2>/dev/null || true
fi

exit ${EXIT_CODE}
