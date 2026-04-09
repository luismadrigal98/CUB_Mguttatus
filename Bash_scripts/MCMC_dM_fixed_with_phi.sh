#!/bin/bash
#SBATCH --job-name=MCMC_dM_fixed_with_phi
#SBATCH --output=MCMC_dM_fixed_with_phi_output
#SBATCH --partition=kelly,kucg,eeb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50g
#SBATCH --time=20-00:00:00
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /home/l338m483/scratch/CUB/CUB_Mguttatus

Rscript R_scripts_remotes/AnaCoDa_pipeline.R \
  -i ./data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa \
  -o ./results/MCMC_results/results_dM_fixed_with_phi_final \
  --phi ./data/compiled_expression_IM767.txt \
  -s 10000 \
  --est_csp \
  --est_phi \
  --est_hyp \
  -n 1 \
  -d 4000 \
  -a 25 \
  --max_num_runs 6 \
  --fix_dM \
  --dM ./data/Mguttatus_intron_derived_dM.csv \
  --sphi_init 1.4583 \
  --sepsilon_init '0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5'

Rscript R_scripts_remotes/AnaCoDa_pipeline.R \
  -i ./data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa \
  -o ./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final \
  --phi ./data/compiled_expression_IM767.txt \
  -s 10000 \
  --est_csp \
  --est_phi \
  --est_hyp \
  -n 1 \
  -d 4000 \
  -a 25 \
  --max_num_runs 6 \
  --fix_dM \
  --dM ./data/Mguttatus_intergenic_derived_dM.csv \
  --sphi_init 1.4583 \
  --sepsilon_init '0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5'

echo "Job finished on $(date)"
