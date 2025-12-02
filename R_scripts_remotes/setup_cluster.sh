#!/bin/bash
#
# setup_cluster.sh - Set up the Pure R ROC model on an HPC cluster
#
# This script:
#   1. Compiles the C acceleration library
#   2. Tests that the library loads correctly
#   3. Runs a quick benchmark
#
# Usage:
#   cd /path/to/Codon_bias_analysis/R_scripts_remotes
#   bash setup_cluster.sh
#
# SLURM example:
#   sbatch -N 1 -c 4 --mem=8G -t 00:30:00 setup_cluster.sh
#

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/src"

echo "============================================"
echo "  Pure R ROC Model - Cluster Setup"
echo "============================================"
echo ""
echo "Script directory: ${SCRIPT_DIR}"
echo "Source directory: ${SRC_DIR}"
echo ""

# Check for R
if ! command -v R &> /dev/null; then
    echo "ERROR: R not found. Please load the R module first."
    echo "  Example: module load R/4.3.1"
    exit 1
fi

echo "R version: $(R --version | head -1)"
echo ""

# Create src directory if needed
mkdir -p "${SRC_DIR}"

# Check if C source exists
if [ ! -f "${SRC_DIR}/roc_likelihood.c" ]; then
    echo "ERROR: C source file not found at ${SRC_DIR}/roc_likelihood.c"
    exit 1
fi

# Compile the C code
echo "Compiling C acceleration library..."
cd "${SRC_DIR}"

# Clean old builds
rm -f roc_likelihood.o roc_likelihood.so roc_likelihood.dll

# Compile
R CMD SHLIB roc_likelihood.c -o roc_likelihood.so

if [ -f "roc_likelihood.so" ]; then
    echo "✓ Compilation successful!"
    echo "  Library: ${SRC_DIR}/roc_likelihood.so"
    ls -lh roc_likelihood.so
else
    echo "✗ Compilation failed!"
    exit 1
fi

echo ""
echo "Testing library loading..."

cd "${SCRIPT_DIR}"

# Test that the library loads
Rscript -e "
tryCatch({
  dyn.load('${SRC_DIR}/roc_likelihood.so')
  cat('✓ Library loads successfully!\n')
  
  # Check that functions are registered
  if (is.loaded('C_calc_log_lik_all_genes')) {
    cat('✓ C_calc_log_lik_all_genes is available\n')
  }
  if (is.loaded('C_batch_update_phi')) {
    cat('✓ C_batch_update_phi is available\n')
  }
  
  cat('\nC acceleration is ready for use!\n')
}, error = function(e) {
  cat('✗ Failed to load library:', e\$message, '\n')
  quit(status = 1)
})
"

echo ""
echo "============================================"
echo "  Quick Benchmark"
echo "============================================"

# Run a quick benchmark
if [ -f "${SCRIPT_DIR}/benchmark_roc.R" ]; then
    Rscript "${SCRIPT_DIR}/benchmark_roc.R"
else
    echo "Benchmark script not found, skipping..."
fi

echo ""
echo "============================================"
echo "  Setup Complete!"
echo "============================================"
echo ""
echo "You can now run the ROC MCMC with:"
echo ""
echo "  Rscript ROC_model_pure_R.R \\"
echo "    -i /path/to/sequences.fa \\"
echo "    --phi /path/to/expression.csv \\"
echo "    --with_phi \\"
echo "    --dM /path/to/initial_dM.csv \\"
echo "    --fix_dM \\"
echo "    -s 10000 -t 10 \\"
echo "    -o /path/to/output"
echo ""
echo "For SLURM jobs, use the provided SLURM template."
echo "============================================"
