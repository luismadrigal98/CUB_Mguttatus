#!/bin/bash
# Usage: sbatch Bash_scripts/run_intronic_twostate.sh <all_sites.vcf[.gz]> <gff3> <out.csv> [sfs_out.csv]
#
# Cluster invocation (from the repo root):
#   sbatch Bash_scripts/run_intronic_twostate.sh \
#       /home/l338m483/scratch/VCF_lines/Lines187.to.IM767v2.1.vcf.gz \
#       data/Mguttatusvar_IM767_887_v2.1.gene.gff3 \
#       data/intron_2_allele.csv \
#       data/intron_sfs_by_gene.csv
#
# The VCF MUST contain invariant sites. Lines187.to.IM767v2.1_ready*.vcf do not
# (vcf_cleaner.py removed 277 M invariant lines); the May 2026 run used one of
# them and produced C frequencies of ~0.32 over ~19% of intronic positions. The
# Python script now exits with status 2, writing nothing, if it is given one.
#
# Outputs (both read by main_CUB.R, Section 12):
#   <out.csv>      per-gene two-allele statistics (C and G systems)
#   [sfs_out.csv]  per-gene SFS of segregating intronic sites
#
# Expected log lines on success: invariant share ~0.95, intronic C ~0.15, G ~0.17.

#SBATCH --job-name=intronic_pi
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --output=run_intronic_twostate.%j.out
#SBATCH --error=run_intronic_twostate.%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

VCF=${1:-}
GFF3=${2:-}
OUT=${3:-}
SFS_OUT=${4:-}

if [[ -z "$VCF" || -z "$GFF3" || -z "$OUT" ]]; then
  echo "Usage: $0 <all_sites.vcf[.gz]> <gff3> <out.csv> [sfs_out.csv]"
  exit 1
fi
if [[ "$(basename "$VCF")" == *_ready* ]]; then
  echo "ERROR: $VCF is a variants-only VCF (invariant sites removed). Use the all-sites VCF." >&2
  exit 1
fi
if [[ -z "$SFS_OUT" ]]; then
  SFS_OUT="$(dirname "$OUT")/intron_sfs_by_gene.csv"
fi

CPUS=${SLURM_CPUS_PER_TASK:-8}
WORKERS=${CPUS}

# Tunable defaults (override in the environment before sbatch)
BUFFER_MB=${BUFFER_MB:-64}
BATCH_SIZE=${BATCH_SIZE:-20000}
TRIM_BP=${TRIM_BP:-30}
MIN_WIDTH=${MIN_WIDTH:-86}
MIN_CALLED=${MIN_CALLED:-10}

echo "Job: ${SLURM_JOB_ID:-local}"
echo "VCF: $VCF"
echo "GFF: $GFF3"
echo "OUT: $OUT"
echo "SFS_OUT: $SFS_OUT"
echo "CPUS: $CPUS, workers: $WORKERS, buffer_mb: $BUFFER_MB, batch_size: $BATCH_SIZE"
echo "trim_bp: $TRIM_BP, min_width: $MIN_WIDTH, min_called: $MIN_CALLED"

if command -v module >/dev/null 2>&1; then
  module load conda
fi

if command -v conda >/dev/null 2>&1; then
  set +u  # conda's activation scripts reference unset variables
  eval "$(conda shell.bash hook)"
  if conda env list | awk '{print $1}' | grep -qx 'PyR'; then
    conda activate PyR
  else
    echo "Warning: conda environment 'PyR' not found; using current Python." >&2
  fi
  set -u
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEFAULT_REPO_ROOT=/home/l338m483/scratch/CUB/CUB_Mguttatus
if [[ -n "${WORKDIR:-}" ]]; then
  REPO_ROOT=$WORKDIR
elif [[ -f "$SCRIPT_DIR/../miscellanea_code/calculate_intronic_twostatepi.py" ]]; then
  REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
elif [[ -f "$DEFAULT_REPO_ROOT/miscellanea_code/calculate_intronic_twostatepi.py" ]]; then
  REPO_ROOT=$DEFAULT_REPO_ROOT
else
  echo "Could not locate miscellanea_code/calculate_intronic_twostatepi.py." >&2
  echo "Set WORKDIR to your repo root or update DEFAULT_REPO_ROOT in the script." >&2
  exit 1
fi

PYTHON_SCRIPT="$REPO_ROOT/miscellanea_code/calculate_intronic_twostatepi.py"
cd "$REPO_ROOT"

# Decide by extension: `file --mime-type` reports gzip differently across versions.
if [[ "$VCF" == *.gz ]]; then
  STREAM_CMD=(zcat "$VCF")
else
  STREAM_CMD=(cat "$VCF")
fi

"${STREAM_CMD[@]}" | python3 "$PYTHON_SCRIPT" /dev/stdin "$GFF3" "$OUT" \
  --sfs-out "$SFS_OUT" \
  --buffer-mb "$BUFFER_MB" \
  --workers "$WORKERS" \
  --batch-size "$BATCH_SIZE" \
  --trim-bp "$TRIM_BP" \
  --min-width "$MIN_WIDTH" \
  --min-called "$MIN_CALLED"

echo "Finished: $(date)"

exit 0
