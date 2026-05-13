#!/bin/bash
# Usage: sbatch Bash_scripts/run_sfs_introns.sh <vcf.gz> <gff3> <out_dir>
# Example: sbatch Bash_scripts/run_sfs_introns.sh data/all.vcf.gz data/annotations.gff3 data/
#
# Produces two files in <out_dir>:
#   sfs_introns_G.csv   (n, k, count) for G-system intronic sites
#   sfs_introns_C.csv   (n, k, count) for C-system intronic sites
#
# Both files feed main.R's load_and_estimate_neutral_params() to recover the
# 4N-scaled Beta shape parameters (alpha, beta) that calibrate the Wright
# two-allele model in Branch A.

#SBATCH --job-name=sfs_introns
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --output=run_sfs_introns.%j.out
#SBATCH --error=run_sfs_introns.%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

VCF_GZ=${1:-}
GFF3=${2:-}
OUT_DIR=${3:-}

if [[ -z "$VCF_GZ" || -z "$GFF3" || -z "$OUT_DIR" ]]; then
  echo "Usage: $0 <vcf.gz> <gff3> <out_dir>"
  exit 1
fi

CPUS=${SLURM_CPUS_PER_TASK:-8}
WORKERS=${CPUS}

BATCH_SIZE_DEFAULT=20000
BATCH_SIZE=${BATCH_SIZE:-$BATCH_SIZE_DEFAULT}

echo "Job: ${SLURM_JOB_ID:-local}"
echo "VCF: $VCF_GZ"
echo "GFF: $GFF3"
echo "OUT_DIR: $OUT_DIR"
echo "CPUS: $CPUS, workers: $WORKERS, batch_size: $BATCH_SIZE"

if command -v module >/dev/null 2>&1; then
  module load conda
fi

if command -v conda >/dev/null 2>&1; then
  set +u
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
elif [[ -f "$SCRIPT_DIR/../miscellanea_code/filter_vcf_for_introns.py" ]]; then
  REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
elif [[ -f "$DEFAULT_REPO_ROOT/miscellanea_code/filter_vcf_for_introns.py" ]]; then
  REPO_ROOT=$DEFAULT_REPO_ROOT
else
  echo "Could not locate miscellanea_code/filter_vcf_for_introns.py." >&2
  echo "Set WORKDIR to your repo root or update DEFAULT_REPO_ROOT in the script." >&2
  exit 1
fi

PYTHON_SCRIPT="$REPO_ROOT/miscellanea_code/filter_vcf_for_introns.py"

# Resolve VCF and GFF to absolute paths before changing directory.
VCF_ABS=$(readlink -f "$VCF_GZ")
GFF_ABS=$(readlink -f "$GFF3")
OUT_DIR_ABS=$(readlink -f "$OUT_DIR")
mkdir -p "$OUT_DIR_ABS"

# The Python script writes sfs_introns_{G,C}.csv to the current working
# directory; run it from OUT_DIR so outputs land in the requested location.
cd "$OUT_DIR_ABS"

# Stream the VCF into the Python script. Auto-detect gzip-compressed vs plain text.
if file -b --mime-type "$VCF_ABS" | grep -q '^application/gzip$'; then
  STREAM_CMD=(zcat "$VCF_ABS")
else
  STREAM_CMD=(cat "$VCF_ABS")
fi

"${STREAM_CMD[@]}" | python3 "$PYTHON_SCRIPT" \
  --stream \
  --gff "$GFF_ABS" \
  --workers "$WORKERS" \
  --batch_size "$BATCH_SIZE"

echo "Wrote: $OUT_DIR_ABS/sfs_introns_G.csv"
echo "Wrote: $OUT_DIR_ABS/sfs_introns_C.csv"
echo "Finished: $(date)"

exit 0
