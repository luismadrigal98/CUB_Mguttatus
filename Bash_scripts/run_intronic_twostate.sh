#!/bin/bash
# Usage: sbatch Bash_scripts/run_intronic_twostate.sbatch <vcf.gz> <gff3> <out.csv>
# Example: sbatch Bash_scripts/run_intronic_twostate.sbatch data/all.vcf.gz data/annotations.gff3 results/Two_allele_pi.csv

#SBATCH --job-name=intronic_pi
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --output=run_intronic_twostate.%j.out
#SBATCH --error=run_intronic_twostate.%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

VCF_GZ=${1:-}
GFF3=${2:-}
OUT=${3:-}

if [[ -z "$VCF_GZ" || -z "$GFF3" || -z "$OUT" ]]; then
  echo "Usage: $0 <vcf.gz> <gff3> <out.csv>"
  exit 1
fi

# Determine workers from SLURM, fallback to SBATCH default above
CPUS=${SLURM_CPUS_PER_TASK:-8}
WORKERS=${CPUS}

# Tunable defaults (can be overridden in environment before sbatch)
BUFFER_MB_DEFAULT=64
BATCH_SIZE_DEFAULT=20000
TRIM_BP_DEFAULT=30
MIN_WIDTH_DEFAULT=86

BUFFER_MB=${BUFFER_MB:-$BUFFER_MB_DEFAULT}
BATCH_SIZE=${BATCH_SIZE:-$BATCH_SIZE_DEFAULT}
TRIM_BP=${TRIM_BP:-$TRIM_BP_DEFAULT}
MIN_WIDTH=${MIN_WIDTH:-$MIN_WIDTH_DEFAULT}

echo "Job: ${SLURM_JOB_ID:-local}"
echo "VCF: $VCF_GZ"
echo "GFF: $GFF3"
echo "OUT: $OUT"
echo "CPUS: $CPUS, workers: $WORKERS, buffer_mb: $BUFFER_MB, batch_size: $BATCH_SIZE"

# Load modules if your cluster requires it (uncomment & adjust)

if command -v module >/dev/null 2>&1; then
  module load conda
fi

if command -v conda >/dev/null 2>&1; then
  set +u  # Temporarily allow unbound variables for conda setup
  eval "$(conda shell.bash hook)"
  if conda env list | awk '{print $1}' | grep -qx 'PyR'; then
    conda activate PyR
  else
    echo "Warning: conda environment 'PyR' not found; using current Python." >&2
  fi
  set -u  # Re-enable unbound variable check
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=${WORKDIR:-$(cd "$SCRIPT_DIR/.." && pwd)}
cd "$REPO_ROOT"

# Stream the VCF into the Python script. Auto-detect gzip-compressed vs plain text input.
if file -b --mime-type "$VCF_GZ" | grep -q '^application/gzip$'; then
  STREAM_CMD=(zcat "$VCF_GZ")
else
  STREAM_CMD=(cat "$VCF_GZ")
fi

"${STREAM_CMD[@]}" | python3 miscellanea_code/calculate_intronic_twostatepi.py /dev/stdin "$GFF3" "$OUT" \
  --buffer-mb "$BUFFER_MB" \
  --workers "$WORKERS" \
  --batch-size "$BATCH_SIZE" \
  --trim-bp "$TRIM_BP" \
  --min-width "$MIN_WIDTH"

echo "Finished: $(date)"

exit 0