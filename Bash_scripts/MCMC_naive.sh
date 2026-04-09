#!/bin/bash
#SBATCH --job-name=MCMC_naive
#SBATCH --output=MCMC_naive_output
#SBATCH --partition=kelly,kucg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2g
#SBATCH --time=20-00:00:00
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /home/l338m483/scratch/CUB/CUB_Mguttatus

Rscript R_scripts_remotes/AnaCoDa_pipeline.R \
  -i ./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnlyClean.fa \
  -o ./results/MCMC_results/results_naive_2 \
  -s 10000 \
  --est_csp \
  --est_phi \
  --est_hyp \
  -n 10 \
  -d 4000 \
  -a 25 \
  --max_num_runs 6

echo "Job finished on $(date)"
