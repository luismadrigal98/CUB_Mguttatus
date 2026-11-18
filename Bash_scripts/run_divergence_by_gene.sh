#!/bin/bash
# Usage: sbatch Bash_scripts/run_divergence_by_gene.sh [aln_dir] [out_dir]
#
# Per-gene divergence from IM767 to JUNG1 (outgroup) and KONZA (check) at 4-fold
# and trimmed intronic sites, from the AnchorWave proAli alignments of 2026-11-04.
# Writes into data/divergence/ (default):
#   fourfold_sites.tsv.gz, intron_divergence_by_gene.tsv, zerofold_divergence_by_gene.tsv,
#   intron_outgroup_{JUNG1,KONZA}.fa.gz (outgroup base per aligned intronic site),
#   intron_aligned_{JUNG1,KONZA}.bed, summary_by_chrom.tsv
# The base tracks are the --outgroup inputs of the next step (fixed differences):
#   OUTGROUPS="JUNG1=data/divergence/intron_outgroup_JUNG1.fa.gz KONZA=data/divergence/intron_outgroup_KONZA.fa.gz" \
#     sbatch Bash_scripts/run_intronic_twostate.sh <all-sites vcf.gz> <gff3> data/intron_2_allele.csv data/intron_sfs_by_gene.csv
#
# ~10 s per chromosome; memory is dominated by the genome (~0.5 GB).

#SBATCH --job-name=div_by_gene
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=run_divergence_by_gene.%j.out
#SBATCH --error=run_divergence_by_gene.%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

ALN_DIR=${1:-/home/l338m483/scratch/767_SF_KONZA_ANCHORWAVE/anchorwave_proali_out/chrom}
OUT_DIR=${2:-data/divergence}
TAXA=${TAXA:-JUNG1,KONZA}

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
REPO_ROOT=${WORKDIR:-$(cd "$SCRIPT_DIR/.." && pwd)}
[[ -f "$REPO_ROOT/miscellanea_code/divergence_by_gene.py" ]] || REPO_ROOT=/home/l338m483/scratch/CUB/CUB_Mguttatus
cd "$REPO_ROOT"

echo "Job: ${SLURM_JOB_ID:-local}  ALN_DIR: $ALN_DIR  OUT_DIR: $OUT_DIR  TAXA: $TAXA"
python3 -c "import numpy" || { echo "numpy is required" >&2; exit 1; }

python3 miscellanea_code/divergence_by_gene.py \
  --genome data/Mguttatusvar_IM767_887_v2.0.fa \
  --gff data/Mguttatusvar_IM767_887_v2.1.gene.gff3 \
  --aln-dir "$ALN_DIR" --taxa "$TAXA" --chroms all --out-dir "$OUT_DIR"

echo "Finished: $(date)"
