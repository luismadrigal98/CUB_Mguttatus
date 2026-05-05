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

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /home/l338m483/scratch/CUB/CUB_Mguttatus

# Stream-compressed VCF into the Python script. Use /dev/stdin as VCF path.
zcat "$VCF_GZ" | python3 miscellanea_code/calculate_intronic_twostatepi.py /dev/stdin "$GFF3" "$OUT" \
  --buffer-mb "$BUFFER_MB" \
  --workers "$WORKERS" \
  --batch-size "$BATCH_SIZE" \
  --trim-bp "$TRIM_BP" \
  --min-width "$MIN_WIDTH"

echo "Finished: $(date)"

exit 0