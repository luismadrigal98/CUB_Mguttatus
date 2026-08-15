## ============================================================================
## supplementary_anacoda.R — ROC-SEMPPR / AnaCoDa supplementary analyses
##
## PROVENANCE
##   Extracted from `main.R` / `main_CUB.R` on 2026-08-14, when the preferred-
##   codon call was moved off AnaCoDa onto the expression-regime estimator in
##   `src/detect_preferred_codons.R`.  The pre-split pipeline is preserved at
##   `archive/main_full_analysis_2026-08-14.R`.
##
## WHY THIS IS A SEPARATE SCRIPT
##   AnaCoDa is archived upstream, so `main_CUB.R` must run without it — and
##   does: it no longer loads the package or any ROC-derived quantity.  What
##   remains here is the material that genuinely needs the MCMC posterior.
##
## WHAT IT IS FOR — AND WHAT IT IS NOT
##   KEEP: the correspondence between the AnaCoDa load/efficacy axes and the
##   polymorphism-based S_Wright.  S_ROC tracks S_Wright closely, and that
##   agreement between two very differently-derived quantities is worth showing
##   (it is the density figure in the paper).
##
##   DO NOT present the preferred-codon agreement (18/19 families) as
##   independent validation.  The selection coefficients used here come from
##   `results_dM_fixed_with_phi_final` — the run that was GIVEN the empirical
##   expression data.  Both that fit and the expression-regime estimator are
##   reading the same expression-codon association, so their agreement is a
##   reproducibility statement, not corroboration.
##
##   Independent validation of the preferred-codon call comes from tRNA gene
##   content instead — see `src/trna_supply_validation.R`, which uses no
##   expression data at all.
##
##   Known diagnostic worth reporting honestly: the expression-FREE run
##   (`results_dM_fixed`, --est_phi with no observed expression) fits the codon
##   data well, but its phi estimates ANTICORRELATE with empirical expression.
##   That is consistent with the two-regime structure documented in
##   src/detect_preferred_codons.R: across the drift-dominated bulk of genes,
##   GC3-ending codon usage DECREASES with expression, so a model inferring phi
##   from codon composition alone learns the bulk axis and assigns high phi to
##   genes that are in fact lowly expressed.  Feeding empirical expression in
##   (`_with_phi`) pins the axis the right way round, after which the model
##   recovers the C-ending preferred codons from the upper tail.  This is also
##   the most likely reason the trajectory figure looked like extrapolation to
##   Referee 1: the genes that identify the selection regime are a tiny
##   minority of the data.
##
## HOW TO RUN
##   Run `main_CUB.R` first.  It writes `results/cub_handoff_for_anacoda.rds`
##   containing msd_data, the SW_group classification, bin_sw, S_BARRIER and
##   the U/V calibration.  Then source this file.  Requires the AnaCoDa package
##   and the MCMC output tree under results/MCMC_results/.
##
## OUTPUTS
##   results/Drift_barrier_overview.pdf          three-panel A/B/C figure
##   results/Go_enrichment_load_ROC_eff.csv      GO for the top-50 load group
##   results/GO_dotplot_load_ROC_eff.pdf
##   results/Top_genes_strong_selection_load.csv
##   results/Wright_per_gene_ROC_eff_S_Wright.csv
##   results/preferred_codon_concordance_vs_ROC.csv
##   results/ROC_codon_trajectories.pdf, convergence + phi diagnostics
## ============================================================================

source("./src/set_environment.R")
required_libraries <- c('data.table', 'Biostrings', 'ggplot2', 'dplyr', 'tidyr',
                        'AnaCoDa', 'coda', 'patchwork', 'gprofiler2', 'scales')
set_environment(required_pckgs = required_libraries, personal_seed = 1998,
                parallel_backend = FALSE)

handoff <- readRDS("./results/cub_handoff_for_anacoda.rds")
msd_data              <- handoff$msd_data
plot_barrier          <- handoff$plot_barrier   # carries SW_group from main_CUB.R
bin_sw                <- handoff$bin_sw
S_BARRIER             <- handoff$S_BARRIER
U_bin_calib           <- handoff$U_bin_calib
V_bin_calib           <- handoff$V_bin_calib
custom_bag            <- handoff$custom_bag
integrated_data       <- handoff$integrated_data
codon_usage           <- handoff$codon_usage
exp_complete          <- handoff$exp_complete
genetic_code_dna_long <- handoff$genetic_code
preferred_detection   <- list(preferred = handoff$preferred)
barrier_colors        <- handoff$barrier_colors

## ****************************************************************************
## 1) Neutral mutation baseline (dM) and MCMC convergence ----
##    (was main.R Section 8 preamble and 8.1)
## ____________________________________________________________________________

## Estimate mutation rates ----

# Use the wrapper function to generate dM files from both introns and intergenic regions
# This replaces ~130 lines of duplicated code with a single function call

dM_results <- estimate_dM_from_neutral_regions(
 fasta_file = "./data/Mguttatusvar_IM767_887_v2.0.hardmasked.fa",
  ann_file = "./data/Mguttatusvar_IM767_887_v2.1.gene.gff3",
  output_dir = "./data",
  output_prefix = "Mguttatus",
  source = "both",  # Generate dM from BOTH introns and intergenic regions
  window_size = 100000,
  min_bp = 1000,
  max_N_freq = 0.25,
  organism = "Mimulus guttatus",
  return_intermediates = TRUE  # Keep intermediate data for further analysis if needed
)
## Results from AnaCoDa framework can be obtained by running:

# Rscript R_scripts_remotes/AnaCoDa_pipeline.R \
# -i ./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnlyClean.fa \
# -o ./MCMC_results/results_dM_fixed \
# -s 10000 \
# --est_csp \
# --est_phi \
# --est_hyp \
# -n 10 \
# -d 4000 \
# -a 25 \
# --max_num_runs 6 \
# --fix_dM \
# --dM ./data/Mguttatus_intron_derived_dM.csv
# 
# echo "Job finished on $(date)"

# 8.1) Retrieving AnaCoDa results to analyze congruence between runs ----

# 8.1.1) Naive model ----

# Setup paths for the 6 runs
run_dirs <- c(
  "./results/MCMC_results/results_naive/run_1",
  "./results/MCMC_results/results_naive/run_2",
  "./results/MCMC_results/results_naive/run_3",
  "./results/MCMC_results/results_naive/run_4",
  "./results/MCMC_results/results_naive/run_5",
  "./results/MCMC_results/results_naive/run_6"
)

Naive_conv <- GR_convergence(run_dirs)

# Convergence: FALSE

# 8.1.2) dM-fixed model ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed/run_1",
  "./results/MCMC_results/results_dM_fixed/run_2",
  "./results/MCMC_results/results_dM_fixed/run_3",
  "./results/MCMC_results/results_dM_fixed/run_4",
  "./results/MCMC_results/results_dM_fixed/run_5",
  "./results/MCMC_results/results_dM_fixed/run_6"
)

dM_fixed_conv <- GR_convergence(run_dirs, parameter = 'selection') # Mutation is fixed

# Convergence: TRUE

# 8.1.2.1) Checking the correlation between estimates of phi and the expression data ----

# From now on, we will work with chain 1, as an example

phi_hat_dM_fixed <- read.csv(file = "results/MCMC_results/results_dM_fixed/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed <- exp_complete |>
  left_join(phi_hat_dM_fixed, by = join_by("Gene_name" == "GeneID")) |>
  na.exclude()

cor.test(phi_dM_fixed$Mean.log10.Phi, phi_dM_fixed$Mean_Log10_Exp)

# We would expect a positive correlation. A negative rho suggest a model 
# misspecification

# Visualization

ggplot(data = phi_dM_fixed, aes(x = Mean.log10.Phi,
                                y = Mean_Log10_Exp)) +
  geom_point() +
  geom_smooth() +
  theme_custom() +
  xlab("Estimated phi (log10)") +
  ylab("Empirical Max Expresion (log10)")

ggsave("./results/phi_estimates_vs_expression_dM_fixed.pdf",
       width = 6, height = 5)

# There is no good correspondence with empirical data
# Next step is to pass expression data to the AnaCoDa

# 8.1.3) Preparing the expression data ----

# 1. Filter for complete cases (Intersection of expresion sources)
# We strictly remove genes with 0 counts in any tissue
multi_tissue_phi <- exp_complete |>
  dplyr::select(Gene_name, contains(c("IM62", "IM767"))) |>
  dplyr::rename(GeneID = Gene_name) # AnaCoDa expects "GeneID" as first col
  
multi_tissue_phi <- multi_tissue_phi |>
  dplyr::filter(rowSums(as.matrix(multi_tissue_phi[, -1])) > 0) |>
  dplyr::filter(GeneID %in% names(trans)) # Ensures correspondence with transcriptome file

# 2. Calculate sphi (Global Prior)
# We estimate the "True Phi" shape by taking the mean of the log-expressions
# This gives the model the "width" of the overall distribution.
log_means <- rowMeans(log(multi_tissue_phi[, -1] + 1))
sphi_init <- sd(log_means)

# 3. Calculate sepsilon (Noise per tissue)
# AnaCoDa needs a vector: c(noise_leaf, noise_bud)
# A good heuristic for initialization is 0.5.
# (The model will refine this during MCMC, but this puts it in the right ballpark)

num_tissues <- ncol(multi_tissue_phi) - 1
sepsilon_init <- rep(0.5, num_tissues)

sphi_str <- paste(round(sphi_init, 4), collapse = ",")
sepsilon_str <- paste(round(sepsilon_init, 4), collapse = ",")

message("\nUse these flags in your script:\n")
message("--sphi_init ", sphi_str, "\n")
message("--sepsilon_init ", shQuote(sepsilon_str), "\n")

# 4. Write empirical expression data
write.table(
  multi_tissue_phi, 
  file = "./data/observed_expression_multitissue.csv", 
  sep = ",", 
  row.names = FALSE, 
  quote = FALSE 
)

# Memory cleanup: phi estimation intermediates ---
# Keeping: phi_hat_dM_fixed, Naive_conv, dM_fixed_conv
rm(phi_dM_fixed, multi_tissue_phi,
   log_means, sphi_init, num_tissues, sepsilon_init,
   sphi_str, sepsilon_str)
gc()

# 8.1.3.1) dM-fixed-with_phi ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1",
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_2",
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_3",
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_4",
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_5",
  "./results/MCMC_results/results_dM_fixed_with_phi_final/run_6"
)

dM_fixed_with_phi_conv <- GR_convergence(run_dirs, 
                                         parameter = 'selection') # Mutation is fixed

# Checking the correlation between phi and empirical values

phi_hat_dM_fixed_with_phi <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed_with_phi <- exp_complete |>
  left_join(phi_hat_dM_fixed_with_phi, by = join_by("Gene_name" == "GeneID")) |>
  na.exclude()

cor.test(phi_dM_fixed_with_phi$Mean.log10.Phi, 
         phi_dM_fixed_with_phi$Mean_Log10_Exp)

# 8.1.3.2) Codon frequency trajectories across expression levels ----

# This section visualizes whether the ROC multinomial model:
#   P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
# correctly predicts how codon frequencies change with expression.

# 1. Prepare codon frequency data from the codon_usage already in environment
# codon_usage is a data.table with Gene_name column and codon count columns
codon_freq_long <- codon_usage |>
  as.data.frame() |>
  dplyr::rename(Gene = Gene_name) |>
  tidyr::pivot_longer(cols = -Gene, names_to = "Codon", values_to = "Count")

# Map codons to amino acids
codon_to_aa <- Biostrings::GENETIC_CODE
codon_to_aa_df <- data.frame(
  Codon = names(codon_to_aa),
  AA = as.character(codon_to_aa),
  stringsAsFactors = FALSE
)

codon_freq_long <- codon_freq_long |>
  dplyr::left_join(codon_to_aa_df, by = "Codon") |>
  dplyr::filter(AA != "*")  # Remove stop codons

# Calculate frequency within each gene's AA family
codon_freq_long <- codon_freq_long |>
  dplyr::group_by(Gene, AA) |>
  dplyr::mutate(
    AA_total = sum(Count, na.rm = TRUE),
    Observed_freq = ifelse(AA_total > 0, Count / AA_total, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(Observed_freq))

cat(sprintf("Codon frequencies: %d gene-codon observations\n", nrow(codon_freq_long)))

# 2. Prepare expression data from exp_complete already in environment
expr_data <- exp_complete |>
  dplyr::mutate(Exp_log10 = Max_Log10_Exp) |>
  dplyr::select(Gene_name, Exp_log10) |>
  dplyr::rename(Gene = Gene_name)

cat(sprintf("Expression data: %d genes\n", nrow(expr_data)))

# 3. Run the trajectory analysis using the convenience wrapper
trajectory_results <- run_trajectory_analysis(
  mutation_file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Mutation.csv",
  selection_file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv",
  codon_freq_df = codon_freq_long,
  expression_df = expr_data,
  output_file = "./results/ROC_codon_trajectories.pdf",
  n_bins = 10
)

# Memory cleanup: trajectory analysis intermediates ---
# Keeping: trajectory_results
rm(expr_data, codon_to_aa, codon_to_aa_df)
gc()

# 8.1.4) dM-fixed-intergenic ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed_intergenic/run_1",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_2",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_3",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_4",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_5",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_6"
)

dM_fixed_intergenic <- GR_convergence(run_dirs, 
                                       parameter = 'selection') # Mutation is fixed

# Convergence: FALSE

# 8.1.5) dM-fixed-with-phi-intergenic ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_1",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_2",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_3",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_4",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_5",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_6"
)

dM_fixed_with_phi_intergenic <- GR_convergence(run_dirs, 
                                       parameter = 'selection') # Mutation is fixed

# Convergence: TRUE

# Checking the correlation between phi and empirical values

phi_hat_dM_fixed_with_phi_intergenic <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed_with_phi_intergenic <- exp_complete |>
  left_join(phi_hat_dM_fixed_with_phi_intergenic, by = join_by("Gene_name" == "GeneID")) |>
  na.exclude()

cor.test(phi_dM_fixed_with_phi_intergenic$Mean.log10.Phi, 
         phi_dM_fixed_with_phi_intergenic$Mean_Log10_Exp)

# Memory cleanup: convergence diagnostics and phi comparisons ---
# Keeping: dM_fixed_with_phi_conv, dM_fixed_intergenic, dM_fixed_with_phi_intergenic,
#          phi_hat_dM_fixed_with_phi_intergenic, exp_complete
rm(phi_dM_fixed_with_phi, phi_dM_fixed_with_phi_intergenic, run_dirs)
gc()


## ****************************************************************************
## 2) Concordance with the expression-regime preferred-codon call ----
##    Reproducibility check, NOT independent validation — see header.
## ____________________________________________________________________________
# --- Supplementary: concordance with the archived AnaCoDa/ROC-SEMPPR fit -----
# ROC-SEMPPR is retained only as a cross-check (AnaCoDa is archived upstream
# and its trajectory figure drew a referee objection).  eta_data is still read
# because Section 8.3 uses it for the supplementary S_ROC load axis.

eta_data <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv")

roc_reference_codons <- eta_data |>
  dplyr::group_by(AA) |>
  dplyr::slice_min(Mean, n = 1, with_ties = FALSE) |>
  dplyr::pull(Codon)

preferred_codon_concordance <- compare_preferred_codon_sets(
  preferred_detection$preferred, roc_reference_codons, genetic_code_dna_long
)
data.table::fwrite(preferred_codon_concordance,
                   "./results/preferred_codon_concordance_vs_ROC.csv")

cat(sprintf(
  "[Preferred codons] Expression-regime vs ROC-SEMPPR: %d/%d families agree.\n",
  sum(preferred_codon_concordance$Agrees, na.rm = TRUE),
  sum(!is.na(preferred_codon_concordance$ROC_Codon))
))

# Regime figures — replace the ROC trajectory panel (former Figure 4).
# Bin edges are dense in the upper tail on purpose: a uniform decile grid hides
# the reversal inside its top bin.

## ****************************************************************************
## 3) Per-gene load (L_ROC) and efficacy (S_ROC / ROC_eff) ----
##    (was main.R 8.2 tail, 8.3.1, 8.3.2)
## ____________________________________________________________________________
genome <- initializeGenomeObject(file = 'data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa',
                                 match.expression.by.id = TRUE,
                                 observed.expression.file = 'data/compiled_expression_IM767.txt') # Warnings are expected if genes are missing from expression file

parameter_object <- loadParameterObject(file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/R_objects/parameter.Rda")

stopifnot(length(getNames(genome)) ==
          nrow(parameter_object$calculateSelectionCoefficients(1)))

# Visualizing cost per codons and confidence intervals

plot_data <- eta_data |>
  dplyr::mutate(
    # Check if 0 is inside the credible interval
    is_significant = (X2.5. > 0) | (X97.5. < 0),
    
    # Identify the Reference (Mean is exactly 0)
    is_reference = (Mean == 0),
    
    # Create a clean category factor for coloring
    Status = case_when(
      is_reference ~ "Reference (Fixed)",
      is_significant ~ "Significant Deviation",
      TRUE ~ "Not Significant"
    )
  )

p <- ggplot(plot_data, aes(x = Codon, y = Mean, color = Status)) +
  
  # A. The Reference Line (The Baseline)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", 
             linewidth = 0.5) +
  
  # B. The Estimates with Error Bars
  # geom_pointrange is perfect for Mean + Credible Intervals
  geom_pointrange(aes(ymin = X2.5., ymax = X97.5.), size = 0.3) +
  
  # C. Organization: Facet by Amino Acid
  # 'scales = "free_x"' ensures you only see relevant codons per AA panel
  facet_wrap(~AA, scales = "free", ncol = 6) +
  
  # D. Custom Colors to highlight the story
  scale_color_manual(values = c(
    "Significant Deviation" = "#E41A1C", # Red for strong signal
    "Not Significant" = "gray70",        # Faint gray for noise
    "Reference (Fixed)" = "black"        # Black anchor for the reference
  )) +
  
  # E. Aesthetics
  labs(
    y = "Relative Codon Costs, deta_eta (Mean ± 95% CI)",
    x = NULL # Codon labels are self-explanatory
  ) +
  theme_custom() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.background = element_rect(fill = "#f0f0f0"), # Light gray headers for AA
    strip.text = element_text(face = "bold")
  )

ggsave(file = "./results/Codon_Selection_Inefficiency_Estimates.pdf",
       plot = p, width = 12, 
       height = 10)

# Get selection coefficients which extracted as log(s)

selection_coeff <- getSelectionCoefficients(genome = genome, 
                                            parameter = parameter_object, 
                                            samples = 1000)

counts_df <- as.data.frame(codon_usage)
rownames(counts_df) <- counts_df$Gene_name
counts_df$Gene_name <- NULL

sel_mat <- as.matrix(selection_coeff)

common_genes <- intersect(rownames(counts_df), rownames(sel_mat))
common_codons <- intersect(colnames(counts_df), colnames(sel_mat))

counts_aligned <- as.matrix(counts_df[common_genes, common_codons])
sel_aligned <- sel_mat[common_genes, common_codons]

# Identify synonymous codons (AA families with >1 codon, i.e. excluding Met, Trp, STOP)
synonymous_aa <- names(which(table(genetic_code_dna_long) > 1))
synonymous_aa <- setdiff(synonymous_aa, c("Met", "Trp", "STOP"))
synonymous_codons <- names(genetic_code_dna_long)[genetic_code_dna_long %in% synonymous_aa]
synonymous_codons_aligned <- intersect(synonymous_codons, common_codons)

# n_synonymous_codons: per-gene count of synonymous codon sites
n_synonymous_codons <- rowSums(counts_aligned[, synonymous_codons_aligned], na.rm = TRUE)

aa_for_aligned <- genetic_code_dna_long[common_codons]

# Total selection intensity (phi-scaled; sel_aligned already incorporates φ)
total_selection_intensity <- rowSums(counts_aligned * abs(sel_aligned), na.rm = TRUE)

# L_ROC: per-gene mean |Δη × φ| over synonymous codons (translational load being paid).
# sel_aligned is phi-scaled with per-AA preferred codon = 0; |sel| = cost of each codon.
# High L_ROC = high-phi gene still paying selection cost (e.g. Rubisco/photosynthesis).
L_ROC <- total_selection_intensity / n_synonymous_codons

# eta_vec: unscaled AnaCoDa η posterior means per codon.
# Preferred codons η < 0; reference codon η = 0; disfavored η > 0.
eta_vec <- setNames(rep(0, length(common_codons)), common_codons)
m_eta   <- match(common_codons, eta_data$Codon)
eta_vec[!is.na(m_eta)] <- eta_data$Mean[m_eta[!is.na(m_eta)]]

is_4fold_codon <- sapply(names(genetic_code_dna_long), function(cdn) {
  prefix   <- substr(cdn, 1, 2)
  variants <- paste0(prefix, c("A", "T", "C", "G"))
  aa_set   <- unique(genetic_code_dna_long[variants])
  length(aa_set) == 1 && !("STOP" %in% aa_set)
})
fourfold_codons_syn <- intersect(names(is_4fold_codon)[is_4fold_codon], synonymous_codons_aligned)

syn_counts <- counts_aligned[, synonymous_codons_aligned]

n_4fold_syn_sites <- rowSums(syn_counts[, fourfold_codons_syn, drop = FALSE], na.rm = TRUE)

# ROC_eff: per-gene signed codon-usage efficacy at strictly 4-fold sites.
# Defined as −mean(η) weighted by 4-fold site counts. Preferred codons (η < 0)
# make this positive; well-optimized genes have ROC_eff > 0.
# Diagnostic shows Spearman rho ≈ +0.77 with S_Wright (vs −0.16 for old φ×Δη).
# φ is deliberately NOT used: AnaCoDa φ scale is incompatible with Wright 4Nes.
eta_4fold_vec <- eta_vec[fourfold_codons_syn]
ROC_eff <- ifelse(
  n_4fold_syn_sites > 0,
  -as.numeric(syn_counts[, fourfold_codons_syn, drop = FALSE] %*% eta_4fold_vec) /
    n_4fold_syn_sites,
  NA_real_
)
names(ROC_eff) <- common_genes
ROC_eff_4 <- ROC_eff   # backwards-compatible alias (same formula)
names(ROC_eff_4) <- common_genes

selection_metrics <- data.frame(
  Gene_name = common_genes,

  # Signed codon-usage efficacy at 4-fold sites: −mean(η_4fold).
  # Positive = gene uses preferred codons (well-optimized).
  ROC_eff = ROC_eff,

  # Backwards-compatible alias for ROC_eff (same formula).
  ROC_eff_4 = ROC_eff_4,

  # Phi-scaled translational load: mean |Δη × φ| per synonymous codon.
  L_ROC = L_ROC,

  n_codons = n_synonymous_codons,

  row.names = common_genes
)

selection_metrics <- selection_metrics |>
  left_join(phi_hat_dM_fixed_with_phi |> dplyr::select(GeneID, Mean.log10.Phi, MeanPhi),
            by = join_by(Gene_name == GeneID))

# Memory cleanup: AnaCoDa genome/parameter objects and selection matrices ---
rm(genome, parameter_object,
   counts_df, sel_mat, counts_aligned, sel_aligned,
   common_genes, common_codons, phi_hat_dM_fixed_with_phi, p,
   eta_vec, m_eta, eta_4fold_vec,
   is_4fold_codon, fourfold_codons_syn,
   syn_counts, n_4fold_syn_sites, total_selection_intensity)
gc()

# 8.3.1) Relationship between L_ROC and phi ----

final_analysis_data <- selection_metrics |>
  dplyr::filter(L_ROC > 0) |>
  dplyr::mutate(
    Intrinsic_Inefficiency = L_ROC / MeanPhi
  )

p_load <- ggplot(final_analysis_data, aes(x = Mean.log10.Phi,
                                          y = L_ROC)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       name = "Gene Count") +
  geom_smooth(method = "gam", color = "cyan", size = 1.2, se = TRUE) +
  labs(
    x = expression(bold(Log[10]("Expression" ~ (Phi)))),
    y = expression(bold("Translational Load" ~ (L[ROC])))
  ) +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/Load_vs_Expression_Plot.pdf", p_load, width = 6, height = 5)

p_optim <- ggplot(final_analysis_data, 
                  aes(x = Mean.log10.Phi, 
                      y = Intrinsic_Inefficiency)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Gene Count") +
  
  # Trend line
  geom_smooth(method = "gam", color = "green1", size = 1.2, se = TRUE) +
  
  # Log scale y-axis for Inefficiency to see the drop clearly
  scale_y_log10() +
  
  labs(x = expression(bold(Log[10]("Expression" ~ (Phi)))),
    y = expression(bold("Intrinsic Inefficiency" ~ (Delta~eta)))
  ) +
  theme_custom() +
  theme(legend.position = "right")

ggsave("./results/Intrinsic_Inefficiency_vs_Expression_Plot.pdf", 
       p_optim, width = 6, height = 5)

cor_load <- cor.test(final_analysis_data$Mean.log10.Phi,
                     final_analysis_data$L_ROC, method = "spearman",
                     exact = F)
cor_eff <- cor.test(final_analysis_data$Mean.log10.Phi, 
                    final_analysis_data$Intrinsic_Inefficiency, 
                    method = "spearman",
                    exact = F)
# Only in tail
selection_genes <- final_analysis_data |> dplyr::filter(Mean.log10.Phi >= 1)
cor_selection <- cor.test(selection_genes$Mean.log10.Phi, 
                          selection_genes$Intrinsic_Inefficiency, method = "spearman",
                          exact = F)

# Memory cleanup: section 8.3.1 plot objects ---
# Keeping: final_analysis_data, cor_load, cor_eff, cor_selection
rm(p_load, p_optim, selection_genes)
gc()

# 8.3.2) Analyzing the correlation between total selective pressure and CAI and CDC ----

n_pre_sel_join <- nrow(integrated_data)
integrated_data <- integrated_data |>
  dplyr::left_join(selection_metrics, by = "Gene_name") |>
  dplyr::filter(!is.na(L_ROC))
cat(sprintf("integrated_data: %d -> %d genes after left_join(selection_metrics) + filter(!is.na(L_ROC)) (dropped %d; these lack AnaCoDa estimates)\n",
            n_pre_sel_join, nrow(integrated_data), n_pre_sel_join - nrow(integrated_data)))
rm(n_pre_sel_join)

# Correlation between selection metrics and CUB metrics
cor_S_and_bias <- corrr::correlate(
  x = as.matrix(integrated_data[, c("L_ROC", "ROC_eff", "CAI", "CDC", "ENC")]),
  method = "spearman", use = "complete.obs"
)
cat("\n=== Spearman correlations among L_ROC / ROC_eff / CUB metrics ===\n")
print(cor_S_and_bias)

cat("\n--- L_ROC (load) vs CUB metrics ---\n")
print(cor.test(integrated_data$L_ROC, integrated_data$CAI, method = "spearman", exact = FALSE))
print(cor.test(integrated_data$L_ROC, integrated_data$CDC, method = "spearman", exact = FALSE))
print(cor.test(integrated_data$L_ROC[integrated_data$Max_Log10_Exp > 3.5],
               integrated_data$CDC[integrated_data$Max_Log10_Exp > 3.5],
               method = "spearman", exact = FALSE))

cat("\n--- ROC_eff (signed codon efficacy, 4-fold) vs CUB metrics ---\n")
print(cor.test(integrated_data$ROC_eff, integrated_data$CAI, method = "spearman", exact = FALSE))
print(cor.test(integrated_data$ROC_eff, integrated_data$CDC, method = "spearman", exact = FALSE))
print(cor.test(integrated_data$ROC_eff[integrated_data$Max_Log10_Exp > 3.5],
               integrated_data$CDC[integrated_data$Max_Log10_Exp > 3.5],
               method = "spearman", exact = FALSE))
## ============================================================================

## ****************************************************************************
## 4) Correspondence with S_Wright ----
##    This is the section the paper keeps: ROC-derived load and efficacy
##    against the polymorphism-derived drift-barrier classification.
## ____________________________________________________________________________

# Attach the ROC metrics just computed to the polymorphism-side objects handed
# over by main_CUB.R.  The SW_group classification is NOT recomputed here — it
# is carried across in `plot_barrier` so the figure splits genes exactly the
# way main_CUB.R did.
roc_metrics_join <- integrated_data |>
  dplyr::select(Gene_name, L_ROC, ROC_eff, ROC_eff_4)

msd_data <- msd_data |>
  dplyr::left_join(roc_metrics_join, by = "Gene_name")

plot_barrier <- plot_barrier |>
  dplyr::left_join(roc_metrics_join, by = "Gene_name") |>
  dplyr::filter(is.finite(L_ROC), is.finite(ROC_eff))

n_sel_barrier    <- sum(plot_barrier$SW_group == "Selection")
n_nearly_neutral <- sum(plot_barrier$SW_group == "Nearly neutral")
n_drift_barrier  <- sum(plot_barrier$SW_group == "Drift")

# Panel A is rebuilt here (identical to the one main_CUB.R saves on its own)
# so the three-panel figure can be assembled in this script.
p_sw_dist <- ggplot(plot_barrier, aes(x = S_Wright_signed, fill = SW_group)) +
  geom_histogram(binwidth = 0.05, boundary = 0, position = "stack") +
  geom_vline(xintercept = S_BARRIER,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = barrier_colors, name = NULL) +
  scale_y_continuous(trans  = "log1p",
                     breaks = c(0, 10, 100, 1000, 10000),
                     labels = scales::comma_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(
    x        = expression(S[Wright] ~ "(per-gene, non-negative)"),
    y        = "Gene count (log1p)",
    subtitle = sprintf(
      "Drift %d genes  |  Nearly neutral: %d genes  |  Selection: %d genes",
      n_drift_barrier, n_nearly_neutral, n_sel_barrier
    )
  ) +
  theme_custom() +
  theme(legend.position = "top")
bin_roc <- msd_data |>
  dplyr::filter(!is.na(ROC_eff), !is.na(Q_pref_base)) |>
  dplyr::arrange(ROC_eff) |>
  dplyr::mutate(ROC_eff_bin = ntile(ROC_eff, 30)) |>
  dplyr::group_by(ROC_eff_bin) |>
  dplyr::summarize(
    n_genes    = dplyr::n(),
    mean_ROC_eff = mean(ROC_eff),
    sites_total = sum(N_4fold_sites),
    Q_bin      = sum(N_preferred_base) / sum(N_4fold_sites),
    pi_bin     = sum(Pi_sum_4fold, na.rm = TRUE) /
                 sum(Sites_4fold,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    S_Wright_bin = vapply(Q_bin, function(q) {
      tryCatch(wright_invert_Q(q, U = U_bin_calib, V = V_bin_calib),
               error = function(e) NA_real_)
    }, numeric(1))
  )

bin_roc$pi_se <- sqrt(bin_roc$pi_bin * (1 - bin_roc$pi_bin / 2) /
                      pmax(bin_roc$sites_total, 1))
# thr_sel: L_ROC threshold = load of the 50th-highest L_ROC gene ----
# "Outlier load group" for GO / BGS isolation = top 50 by translational load.

thr_sel <- sort(integrated_data$L_ROC, decreasing = TRUE, na.last = NA)[50]
if (!is.finite(thr_sel)) stop("thr_sel could not be determined: fewer than 50 finite L_ROC values.")
attr(thr_sel, "criterion") <- "top50_L_ROC"
attr(thr_sel, "U_empirical") <- U_bin_calib
attr(thr_sel, "V_empirical") <- V_bin_calib

cat(sprintf(
  "\n>>> thr_sel = %.6f  (L_ROC of the 50th-highest gene; top-50 load group)\n",
  as.numeric(thr_sel)
))
cat(sprintf(
  "    %d / %d genes (%.1f%%) above thr_sel (= top 50 by L_ROC).\n",
  sum(integrated_data$L_ROC > thr_sel, na.rm = TRUE),
  nrow(integrated_data),
  100 * mean(integrated_data$L_ROC > thr_sel, na.rm = TRUE)
))
# Per-gene ROC_eff_4 vs S_Wright correlations (non-zero, well-covered genes) ----
per_gene_pool <- msd_data |>
  dplyr::filter(!is.na(S_Wright_signed), !is.na(ROC_eff_4),
                N_4fold_sites >= 50, ROC_eff_4 != 0)
cor_roc_eff_4_spearman <- cor(per_gene_pool$ROC_eff_4, per_gene_pool$S_Wright_signed,
                              method = "spearman")
cor_roc_eff_4_pearson  <- cor(per_gene_pool$ROC_eff_4, per_gene_pool$S_Wright_signed,
                              method = "pearson")
cat(sprintf(
  "[Validation] cor(ROC_eff_4, S_Wright_signed): Spearman = %+.3f, Pearson = %+.3f (n = %d)\n",
  cor_roc_eff_4_spearman, cor_roc_eff_4_pearson, nrow(per_gene_pool)
))
# Panel B: L_ROC density split by group (log scale)
p_lroc_split <- ggplot(plot_barrier, aes(x = L_ROC, fill = SW_group)) +
  geom_density(alpha = 0.55, color = NA) +
  scale_fill_manual(values = barrier_colors, name = NULL) +
  scale_x_log10() +
  coord_cartesian(xlim = c(quantile(plot_barrier$L_ROC[plot_barrier$L_ROC > 0], 0.01, na.rm = TRUE),
                           quantile(plot_barrier$L_ROC, 0.995, na.rm = TRUE))) +
  labs(
    x = expression(L[ROC] ~ "(translational load, log scale)"),
    y = "Density"
  ) +
  theme_custom() +
  theme(legend.position = "none")

# Panel C: ROC_eff density split by group
p_roc_eff_split <- ggplot(plot_barrier, aes(x = ROC_eff, fill = SW_group)) +
  geom_density(alpha = 0.55, color = NA) +
  scale_fill_manual(values = barrier_colors, name = NULL) +
  labs(
    x = expression(S[ROC] ~ "(signed AnaCoDa efficacy, 4-fold)"),
    y = "Density"
  ) +
  theme_custom() +
  theme(legend.position = "none")

p_barrier_overview <- (p_sw_dist / (p_lroc_split | p_roc_eff_split)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave("./results/Drift_barrier_overview.pdf",
       p_barrier_overview, width = 11, height = 9, device = cairo_pdf)

# Per-gene table pairing the ROC axes with S_Wright
write.csv(
  msd_data |> dplyr::select(Gene_name, ROC_eff, ROC_eff_4, L_ROC, Q_pref_base,
                            S_Wright_signed, S_Wright_raw, is_drift,
                            Mean_Log10_Exp, Max_Log10_Exp, N_4fold_sites),
  "./results/Wright_per_gene_ROC_eff_S_Wright.csv", row.names = FALSE
)

## ****************************************************************************
## 5) GO enrichment and top genes for the load-paying group ----
##    The S_Wright selection group is handled in main_CUB.R; this is its
##    AnaCoDa-defined counterpart, kept for the load-vs-efficacy contrast.
## ____________________________________________________________________________

# (a) Load-paying group: top 50 by L_ROC ---------------------------------
subset_load_paying <- integrated_data |>
  dplyr::filter(L_ROC > thr_sel) |>
  dplyr::pull(Gene_name)

GO_results_load <- gost(query = subset_load_paying,
                        organism = "gp__q7VP_EAck_dZk",
                        multi_query = FALSE, significant = TRUE,
                        correction_method = "fdr",
                        domain_scope = "custom", custom_bg = custom_bag,
                        user_threshold = 0.05)
write.csv(x = GO_results_load$result |> dplyr::select(-parents),
          file = "./results/Go_enrichment_load_ROC_eff.csv",
          quote = TRUE, row.names = FALSE)
cat(sprintf("[GO] Load-paying group (top 50 L_ROC; thr_sel = %.6f): n = %d genes\n",
            as.numeric(thr_sel), length(subset_load_paying)))

# (b) Top 50 genes by translational load ---------------------------------

detailed_annotation_full <- read.delim(
  "data/Mguttatusvar_IM767_887_v2.1.annotation_info.txt",
  header = TRUE, sep = "\t", comment.char = "", quote = "",
  fill = TRUE, na.strings = ""
) |>
  dplyr::select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline) |>
  dplyr::distinct()

top_L_ROC <- integrated_data |>
  dplyr::arrange(desc(L_ROC)) |>
  dplyr::select(Gene_name, L_ROC)
top_L_ROC <- top_L_ROC[1:50, ]
top_L_ROC <- top_L_ROC |>
  dplyr::left_join(detailed_annotation_full,
                   by = dplyr::join_by("Gene_name" == "locusName"))

write.csv(top_L_ROC,
          "./results/Top_genes_strong_selection_load.csv",
          quote = TRUE, row.names = FALSE)

# (c) Summary of the correspondence the paper reports --------------------

cat("\n=== ROC vs Wright correspondence (supplementary) ===\n")
cat(sprintf("  cor(ROC_eff_4, S_Wright_signed): Spearman = %+.3f, Pearson = %+.3f (n = %d)\n",
            cor_roc_eff_4_spearman, cor_roc_eff_4_pearson, nrow(per_gene_pool)))
cat(sprintf("  Load group (top 50 L_ROC) vs selection group (S_Wright >= %.2f): %d vs %d genes\n",
            S_BARRIER, length(subset_load_paying),
            sum(plot_barrier$SW_group == "Selection")))
cat(sprintf("  Overlap between the two groups: %d genes\n",
            length(intersect(subset_load_paying,
                             plot_barrier$Gene_name[plot_barrier$SW_group == "Selection"]))))

write.csv(
  data.frame(
    cor_ROC_eff_4_S_Wright_spearman = cor_roc_eff_4_spearman,
    cor_ROC_eff_4_S_Wright_pearson  = cor_roc_eff_4_pearson,
    n_per_gene_pool                 = nrow(per_gene_pool),
    thr_sel                         = as.numeric(thr_sel),
    n_load_group                    = length(subset_load_paying),
    n_selection_group               = sum(plot_barrier$SW_group == "Selection"),
    n_overlap                       = length(intersect(
      subset_load_paying,
      plot_barrier$Gene_name[plot_barrier$SW_group == "Selection"]))
  ),
  "./results/ROC_vs_Wright_correspondence.csv", row.names = FALSE
)
