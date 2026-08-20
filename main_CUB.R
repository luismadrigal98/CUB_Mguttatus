##' @title CUB in Mimulus guttatus — Paper-replication pipeline (CUB-only scope)
##'
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' @date   2026-08-14
##' ____________________________________________________________________________
##'
##' This script contains ONLY the code required to reproduce the figures,
##' tables, and quantitative claims in the codon-usage-bias manuscript
##' resubmitted to New Phytologist (successor to NPH-MS-2026-57041).
##'
##' ---------------------------------------------------------------------------
##' SCOPE CHANGE, 2026-08-14
##' ---------------------------------------------------------------------------
##' Built from `main.R` (preserved at `archive/main_full_analysis_2026-08-14.R`)
##' after the editor approved refocusing the paper on codon usage bias alone.
##' Two changes of substance:
##'
##'   1. WITHIN-GENE POSITIONAL ANALYSES REMOVED — π by distance from gene
##'      start, the first-300 bp decomposition, and the translational-ramp
##'      models moved verbatim to `paper2_linked_selection.R`.  Referee 1 did
##'      not dispute the pattern; the objection was that with no direct measure
##'      of recombination rate, background-selection amelioration cannot be
##'      separated from GC-biased gene conversion or from mutagenic
##'      recombination.  That needs polymorphism-vs-divergence, which needs the
##'      JUNG1 outgroup and the within-gene recombination map — both pending.
##'      The π-versus-expression / π-versus-selection result STAYS: it is the
##'      McVean & Charlesworth expectation and never used position.
##'
##'   2. PREFERRED CODONS NO LONGER COME FROM AnaCoDa — Section 8.2 now calls
##'      `detect_preferred_codons()` (src/detect_preferred_codons.R), which
##'      estimates the optimal codon from the expression-regime contrast in this
##'      study's own data.  AnaCoDa is archived upstream and its trajectory
##'      figure drew a referee objection.  ROC-SEMPPR is retained only as a
##'      supplementary cross-check (18/19 families agree) and for the S_ROC load
##'      axis.  Note S_Wright is unaffected: it is inverted from observed Q at
##'      4-fold sites and needs only the preferred base identity, so the drift
##'      barrier, the selection group and the GO enrichment carry over intact.
##'
##' Section headers below map onto the paper's Results subsections; each lists
##' the figures, tables, and cited values it produces.  Exploratory analyses,
##' alternative parameterizations, diagnostics and model-selection runs live in
##' `full_analysis.R` (verbatim copy of the historical pipeline).
##'
##' Pipeline order is preserved.  Inter-section dependencies (e.g.
##' Pi_mean_4fold pre-load in Section 5.5, preferred_codons from Section 8 used
##' by Sections 9, 11, 12, 14, 16) are unchanged.
##'
##' ____________________________________________________________________________

## ============================================================================
## SETUP — Working directory, libraries, helper functions, reference tables
##   Loads all source files in ./src and the comparative-plant codon table.
## ============================================================================

## *****************************************************************************
## 1) Set work directory ----
## _____________________________________________________________________________

setwd(".")

## *****************************************************************************
## 2) Load required libraries and set up environment ----
## _____________________________________________________________________________

# Source the set_environment function first
source("./src/set_environment.R")

required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel',
                        'doFuture', 'ggplot2', 'grid', 'gridExtra',
                        'ggseqlogo', 'FactoMineR',
                        'factoextra', 'dplyr', 'GenomicFeatures',
                        'ape', 'tidyr', 'caret', 'ggpointdensity',
                        'DescTools', 'mgcv', 'nnet', 'VGAM', 
                        'viridis', 'cubar', 'kohonen',
                        'rtracklayer', 'tidyverse',
                        'txdbmaker', 'Rsamtools', 'purrr',
                        'abind', 'scales', 'mclust', 'coda',
                        'admisc', 'corrr', 'patchwork', 'gprofiler2',
                        'ggnewscale', 'broom', 'reshape2',
                        'furrr', 'tidyr', 'gsl', 'rcompanion',
                        'FSA', 'matrixStats', 'ggpubr',
                        'boot', 'gratia', 'marginaleffects',
                        'corrr', 'nortest', 'patchwork',
                        'betareg')#, 'brms', 'cmdstanr')

set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = T)

# 1.1) Definition of globals ----
# Look-up table

genetic_code_dna_long <- c(
  "TTT"="Phe", "TTC"="Phe", "TTA"="Leu_2", "TTG"="Leu_2",
  "TCT"="Ser_4", "TCC"="Ser_4", "TCA"="Ser_4", "TCG"="Ser_4",
  "TAT"="Tyr", "TAC"="Tyr", "TAA"="STOP", "TAG"="STOP",
  "TGT"="Cys", "TGC"="Cys", "TGA"="STOP", "TGG"="Trp",
  "CTT"="Leu_4", "CTC"="Leu_4", "CTA"="Leu_4", "CTG"="Leu_4",
  "CCT"="Pro", "CCC"="Pro", "CCA"="Pro", "CCG"="Pro",
  "CAT"="His", "CAC"="His", "CAA"="Gln", "CAG"="Gln",
  "CGT"="Arg_4", "CGC"="Arg_4", "CGA"="Arg_4", "CGG"="Arg_4",
  "ATT"="Ile", "ATC"="Ile", "ATA"="Ile", "ATG"="Met",
  "ACT"="Thr", "ACC"="Thr", "ACA"="Thr", "ACG"="Thr",
  "AAT"="Asn", "AAC"="Asn", "AAA"="Lys", "AAG"="Lys",
  "AGT"="Ser_2", "AGC"="Ser_2", "AGA"="Arg_2", "AGG"="Arg_2",
  "GTT"="Val", "GTC"="Val", "GTA"="Val", "GTG"="Val",
  "GCT"="Ala", "GCC"="Ala", "GCA"="Ala", "GCG"="Ala",
  "GAT"="Asp", "GAC"="Asp", "GAA"="Glu", "GAG"="Glu",
  "GGT"="Gly", "GGC"="Gly", "GGA"="Gly", "GGG"="Gly"
)

# Define amino acid chemistry groups
aa_chemistry <- list(
  "Nonpolar_Aliphatic" = c("Ala", "Gly", "Ile", "Leu_2", "Leu_4", "Met", 
                           "Pro", "Val"),
  "Aromatic" = c("Phe", "Trp", "Tyr"),
  "Polar_Uncharged" = c("Asn", "Cys", "Gln", "Ser_2", "Ser_4", "Thr"),
  "Positively_Charged" = c("Arg_2", "Arg_4", "His", "Lys"),
  "Negatively_Charged" = c("Asp", "Glu")
)

aa_chemistry_df <- as.data.frame(stack(aa_chemistry))
colnames(aa_chemistry_df) <- c('AA', 'class')

# Preferred codons in three additional model plants

model_plants_PC <- read.table(file = "data/plant_preferred_codons.txt", 
                              header = T, sep = ',')
## ============================================================================
## DATA — Coding sequences (CDS) and multi-tissue gene expression
##   Builds `trans` (DNAStringSet of primary-transcript CDS), `codon_usage`,
##   and the multi-tissue expression matrix `exp_complete` (Max_Log10_Exp,
##   Exp_breadth).  All downstream sections depend on these objects.
## ============================================================================

## *****************************************************************************
# 3) Load the data ----
## _____________________________________________________________________________

# 3.1) Analysis from transcript ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnlyClean.fa", 
                                      format = 'fasta')

trans <- trans[check_canonical_start(trans)] |> check_cds()

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

# Loading expression data (multi-source) ----

exp_complete <- read.table(file = "./data/compiled_expression_IM767.txt", 
                           header = T, sep = '\t') |>
  dplyr::rename(Gene_name = GeneID) |>
  dplyr::distinct(Gene_name, .keep_all = TRUE)

# Isolate the numeric data (Everything except Gene_name)
numeric_data <- as.matrix(exp_complete[, -1])

# Calculate the "Mean Log Expression"
# Logic: Add 1 (pseudocount) -> Log10 -> Average across tissues
exp_complete$Mean_Log10_Exp <- rowMeans(log10(numeric_data + 1))

# Get the "Max Log Expression"
# Logic: Add 1 (pseudocount) -> Log10 -> Max across tissues
exp_complete$Max_Log10_Exp <- rowMaxs(log10(numeric_data + 1))

# Get the expression breadth
# Logic: Count the number of instances where expression is higher than a threshold
# defined by 1 CPM
CPM_thr <- 1

exp_complete$Exp_breadth <- apply(X = numeric_data, MARGIN = 1, 
                                  FUN = function(x)
                                    {
                                    sum(x > CPM_thr)
                                  })

# Geometric Mean
exp_complete$Geom_Mean_CPM <- 10^(exp_complete$Mean_Log10_Exp) - 1

# Check the result
head(exp_complete[, c("Gene_name", "Mean_Log10_Exp", "Max_Log10_Exp", 
                      "Geom_Mean_CPM", "Exp_breadth")])

# Saving the data for future usage

write.csv(exp_complete, file = "./results/Expression_Profiles_Summary.csv", 
          row.names = FALSE)

# --- Memory cleanup: expression intermediates ---
rm(numeric_data, CPM_thr)
gc()

## ============================================================================
## RESULTS 1 — CUB metrics indicate strong codon usage bias (related to expression)
##   Produces:
##     Figure 1A  RSCU bar plot per amino acid (`codon_usage_barplot.pdf`)
##     Figure 1B  Parity Rule 2 plot (`pr2_plot.pdf`)
##     Table S1   G-test heterogeneity per amino acid
##     Cited values:  ENC range 31.84–59.00, mean 53.43
##                     G-test deviation in 80.2% of genes (20213/25188)
##                     CDC range 0.0401–0.5427, mean 0.1266, median 0.1143
##                     18,722 / 22,556 genes CDC-significant at FDR<0.05
## ============================================================================

## *****************************************************************************
## 4) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results",
                          aa_group = aa_chemistry_df)

# Creation of integrated data ----

integrated_data <- dplyr::left_join(exp_complete |> dplyr::select(Gene_name, 
                                                                  Max_Log10_Exp,
                                                                  Mean_Log10_Exp,
                                                                  Exp_breadth,
                                                                  Geom_Mean_CPM), 
                                    cub_results$enc_results, 
                                    by = dplyr::join_by(Gene_name)) |>
  na.exclude()

# Add gene length (CDS length in codons and nucleotides)
codon_columns <- names(codon_usage)[names(codon_usage) != "Gene_name"]

gene_lengths <- codon_usage |>
  dplyr::mutate(
    Total_Codons = rowSums(across(all_of(codon_columns)), na.rm = TRUE),
    CDS_length_nt = Total_Codons * 3,  # nucleotides
    CDS_length_aa = Total_Codons       # amino acids (codons)
  ) |>
  dplyr::select(Gene_name, Total_Codons, CDS_length_nt, CDS_length_aa)

integrated_data <- integrated_data |>
  left_join(gene_lengths, by = "Gene_name")

# Adding GC content variables

integrated_data <- integrated_data |>
  left_join(cub_results$gc_results, by = "Gene_name") |>
  data.frame() # Strip attributes

# Memory cleanup: gene length intermediates ---
rm(gene_lengths, codon_columns)
gc()

## *****************************************************************************
## 5) CDC-based analysis ----
## _____________________________________________________________________________

# Full integration with the pipeline
integrated_data <- integrate_cdc_analysis(codon_usage, 
                                      genetic_code_dna_long, 
                                      integrated_data, 
                                      n_bootstrap = 10000,
                                      n_cores = parallel::detectCores() - 1)

# Re-plotting CDC-based neutrality plot highlighting the significant genes with CDC ----

# Merge ENC, GC3s, and CDC results
integrated_data <- integrated_data |>
  dplyr::mutate(
    CDC_significant = !is.na(p_adj) & p_adj < 0.05,
    CDC_category = dplyr::case_when(
      is.na(p_value) ~ "No CDC data",
      p_adj < 0.001 ~ "p < 0.001",
      p_adj < 0.01 ~ "p < 0.01",
      p_adj < 0.05 ~ "p < 0.05",
      TRUE ~ "Not significant"
    )
  )

# Count significant genes
n_sig <- sum(integrated_data$CDC_significant, na.rm = TRUE)
n_total <- sum(!is.na(integrated_data$p_adj))
pct_sig <- 100 * n_sig / n_total

cat(sprintf("Found %d / %d (%.1f%%) genes with significant CDC (FDR < 0.05)\n", 
            n_sig, n_total, pct_sig))
## ============================================================================
## RESULTS 2 — Codon usage bias scales with gene expression
##   Produces:
##     Figure 2  GAM prediction of CDC vs Max expression × Exp_breadth,
##               controlling for CDS length (`GAM_Interaction_Predictions_CDC.pdf`)
##     Cited value:  GAM deviance explained = 54%
##   Also creates `Expression_Group` (Top 5% / Middle 90% / Bottom 5%) which is
##   used by Sections 9, 11, and 12.
##   Section 5.5 (polymorphism pre-load) is hoisted here because the Section 6
##   GAM and the Wright MSD block (Section 4 below) both filter on Pi_mean_4fold.
## ============================================================================

## 5.5) Polymorphism data preload (Pi_mean_4fold needed by Section 6+) ----
## _____________________________________________________________________________
# Section 12 historically loaded polymorphism data and joined it into
# integrated_data. Section 6 GAMs and Section 8.3.4 msd_data filter on
# Pi_mean_4fold, so the join is hoisted here. Section 12 keeps the per-feature
# positional decomposition and downstream analyses; only the by-gene join is
# moved.

pi_data <- fread(input = "data/all_chromosomes.bygene.pi.txt")

pi_data <- pi_data |>
  dplyr::select(Chr, Gene, contains("mean"),
                contains("Sites"), contains("Pi_sum"), contains("Poly")) |>
  dplyr::mutate(Gene = paste0("MgIM767.", pi_data[['Gene']])) |>
  dplyr::rename(Gene_name = Gene)

n_pre_pi_join <- nrow(integrated_data)
integrated_data <- integrated_data |>
  dplyr::left_join(pi_data, by = "Gene_name") |>
  na.exclude()
cat(sprintf("integrated_data: %d -> %d genes after left_join(pi_data) + na.exclude() (dropped %d)\n",
            n_pre_pi_join, nrow(integrated_data), n_pre_pi_join - nrow(integrated_data)))
## 6) Modeling relationship between CDC and Expression profiles ----
## _____________________________________________________________________________

# Set of predictors we care about

predictors <- c('Max_Log10_Exp', 'Mean_Log10_Exp', 'Exp_breadth', 
                'CDS_length_nt')

corrr::correlate(integrated_data[, predictors], method = "spearman") |> shave()

# We preserve predictors with correlation less than 0.75

predictors <-  c('Max_Log10_Exp', 'Exp_breadth', 'CDS_length_nt')

# Generate Table ----
# Using bind_rows directly on the list
justification_list <- lapply(predictors, analyze_nonlinearity, 
                             data = integrated_data,
                             resp = "CDC")
justification_table <- dplyr::bind_rows(justification_list)

write.csv(justification_table, 
          "results/Linearity_Justification_Table.csv", 
          row.names = FALSE)

# Generate Visuals (Safe Loop) ----
plot_list <- list()

# Check that 0 or 1 are not between CDC values
# Enforce shrinkage to ensure compatibility with betar

integrated_data <- as.data.table(integrated_data) # Coerce to data.table
integrated_data[CDC == 1, CDC := 0.9999]
integrated_data[CDC == 0, CDC := 0.0001]

for (pred in predictors) {
    
  form_gam <- as.formula(paste0("CDC ~ s(", pred, ")"))
  model_gam <- gam(form_gam, data = integrated_data, 
                   family = betar(link = "logit"),
                   method = "REML",
                   select = T)
  
  p <- gratia::draw(model_gam, residuals = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "red", 
               alpha = 0.5) +
    labs(title = paste0("Partial Effect on CDC (logit scale): ", 
                        pred)) +
    theme_custom()
  
  plot_list[[pred]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 3, scales = "free")
ggsave("results/GAM_Partial_Effects_Gratia.pdf", combined_plot, width = 12, 
       height = 4)

# GAM final models ----

# Given the non-linearity effect of the predictors, we are going to model them
# using GAM models

# Competing models

# Model 0: Null
m_null <- gam(CDC ~ 1,
              data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
              family = betar(link = "logit"), 
              method = "REML")

# Model 1: Additive (Independent effects)
# Hypothesis: Each predictor affects CUB independently.

m_additive <- gam(CDC ~ s(Max_Log10_Exp) + s(Exp_breadth) + s(CDS_length_nt),
                  data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
                  family = betar(link = "logit"), 
                  method = "REML",
                  select = T)

# Model 2: Expression Interaction (The "Trade-off" Hypothesis)
# Hypothesis: High expression only forces strict CUB if the gene is broad.

m_interaction <- gam(CDC ~ te(Max_Log10_Exp, Exp_breadth) + s(CDS_length_nt),
                     data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
                     family = betar(link = "logit"), 
                     method = "REML",
                     select = T)

# Model 3: Complex (Full Interaction)
# Hypothesis: Length and expression interact in complex ways."
m_complex <- gam(CDC ~ te(Max_Log10_Exp, Exp_breadth, CDS_length_nt),
                 data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
                 family = betar(link = "logit"), 
                 method = "REML")

model_list <- list(Null = m_null,
                   Additive = m_additive, 
                   Interaction_Exp = m_interaction, 
                   Interaction_Com = m_complex)

# Select the best model

selection_table <- do.call(rbind, lapply(names(model_list), function(n) {
  m <- model_list[[n]]
  data.frame(Model = n,
             AIC = AIC(m),
             Deviance_Expl = summary(m)$dev.expl,
             R_sq = summary(m)$r.sq)
}))

# Selected model: m_interaction

# We hold CDS Length constant at the mean to isolate the interaction
# Visual of the predictions of the model
p_effects <- plot_predictions(m_interaction, 
                              condition = c("Max_Log10_Exp", "Exp_breadth"), 
                              newdata = datagrid(
                                CDS_length_nt = mean((integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0))$CDS_length_nt)),
                              type = "response") + 
  geom_rug(data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
           aes(x = Max_Log10_Exp), 
           sides = "b", alpha = 0.05, inherit.aes = FALSE) +
  # THEME & LABELS
  theme_custom() + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  labs(y = "Predicted CDC",
       x = "Max Expression (Log10 CPM)")

ggsave("./results/GAM_Interaction_Predictions_CDC.pdf", 
       plot = p_effects, width = 10, height = 6)
# Define expression groups: Top 5% vs Bottom 5% (extreme comparison) ----

top_5_cutoff <- quantile(integrated_data$Max_Log10_Exp, probs = 0.95)
bottom_5_cutoff <- quantile(integrated_data$Max_Log10_Exp, probs = 0.05)

integrated_data$Expression_Group <- case_when(
  integrated_data$Max_Log10_Exp >= top_5_cutoff ~ "Top 5%",
  integrated_data$Max_Log10_Exp <= bottom_5_cutoff ~ "Bottom 5%",
  TRUE ~ "Middle 90%"
)

# Confounding out-based analysis (detendred CDC) ----

# Assesing significance of expression over the detrended residuals

cat("\n=== Kruskal-Wallis Test: Detrended ENC Residuals across Groups ===\n")

kw_detrended <- kruskal.test(CDC_detrended ~ Expression_Group, 
                             data = integrated_data)

# Plotting and assessing significance using Dunn

print(kw_detrended)
if (kw_detrended$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Perform Dunn's test with FDR correction
  dunn_result_detrended <- dunn.test::dunn.test(
    x = integrated_data$CDC_detrended,
    g = integrated_data$Expression_Group,
    method = "bh",
    kw = TRUE,
    label = TRUE,
    wrap = FALSE,
    table = TRUE,
    list = FALSE,
    altp = TRUE
  )
} else {
  cat("\nNo significant difference among groups (p >= 0.05)\n")
  cat("Post-hoc tests not necessary.\n")
}

# Ploting box plot

p_boxplot_detrended <- ggplot(integrated_data, aes(x = Expression_Group, 
                                                   y = CDC_detrended, 
                                                   fill = Expression_Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(outlier.alpha = 0.3) +
  # geom_boxplot(outlier.alpha = 0.3) +
  # geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "Detrended CDC Residuals by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR",
       y = "CDC (lenght corrected)",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/Detrended_CDC_by_expression_group.pdf", 
       p_boxplot_detrended, width = 8, height = 6)

# Median and CI version

plot_data <- integrated_data |>
  dplyr::mutate(Exp_Group = factor(Expression_Group, 
                                   levels = c("Bottom 5%", "Middle 90%", "Top 5%"))) |>
  dplyr::filter(!is.na(Exp_Group))

# Define the Comparisons for the plot
my_comparisons <- list(c("Bottom 5%", "Middle 90%"), 
                       c("Middle 90%", "Top 5%"), 
                       c("Bottom 5%", "Top 5%"))

# Create the Plot
p_medians <- ggplot(plot_data, aes(x = Exp_Group, y = CDC_detrended)) +
  
  # A. Median and 95% CI (Bootstrap)
  # We use the custom function defined above
  stat_summary(fun.data = median_cl_boot, 
               geom = "errorbar", width = 0.15, linewidth = 0.8, color = "black") +
  stat_summary(fun = median, geom = "point", size = 3.5, aes(color = Exp_Group)) +
  
  # B. Reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  
  # C. Formatting
  scale_color_manual(values = c("#377EB8", "#999999", "#E41A1C")) +
  
  # Zoom in (Adjust these limits if your medians are slightly different than means)
  coord_cartesian(ylim = c(NA, 0.01)) + 
  
  labs(y = "CDC Residuals (length corrected)",
       x = NULL) +
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        axis.title.y = element_text(size = 11),
        plot.subtitle = element_text(size = 10, color = "gray30"))

# Save
ggsave("./results/CDC_detrended_Medians_CI.pdf", p_medians, width = 4, height = 3.5)

# Get CDC values for each group
top5_cdc_de <- integrated_data |> filter(Expression_Group == "Top 5%") |> pull(CDC_detrended)
middle_cdc_de <- integrated_data |> filter(Expression_Group == "Middle 90%") |> pull(CDC_detrended)
bottom5_cdc_de <- integrated_data |> filter(Expression_Group == "Bottom 5%") |> pull(CDC_detrended)

# Calculate effect sizes
if (length(top5_cdc_de) > 0 && length(middle_cdc_de) > 0) {
  d_top_middle_de <- cohens_d_calc(top5_cdc_de, middle_cdc_de)
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle_de))
}

if (length(top5_cdc_de) > 0 && length(bottom5_cdc_de) > 0) {
  d_top_bottom_de <- cohens_d_calc(top5_cdc_de, bottom5_cdc_de)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom_de))
}

if (length(middle_cdc_de) > 0 && length(bottom5_cdc_de) > 0) {
  d_middle_bottom_de <- cohens_d_calc(middle_cdc_de, bottom5_cdc_de)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom_de))
}

cat("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

# Memory cleanup: Section 5-6 plots and intermediate objects ---
# Keeping: selection_table, kw_detrended.
rm(m_null, m_additive, m_complex, m_interaction, model_list,
   justification_list, justification_table, plot_list, combined_plot,
   p_effects,
   p_boxplot_detrended, p_medians,
   plot_data, my_comparisons,
   top5_cdc_de, middle_cdc_de, bottom5_cdc_de,
   n_sig, n_total, pct_sig,
   top_5_cutoff, bottom_5_cutoff)
gc()
## ============================================================================
## RESULTS 3 — Codon Adaptation Index (CAI) discriminates highly expressed genes
##   Produces:
##     Figure S3 supporting data — CAI by expression group
##     Cited values:  CAI range 0.54–0.89, mean 0.70
## ============================================================================

## 7) Calculate Codon Adaptation Index (CAI) ----
## _____________________________________________________________________________

# Define reference set: Genes which are constitutively highly expressed
# Example: Elongation factors
reference_genes <- read.table(file = 'data/CAI_Reference_Set_Mguttatus.txt')[, 1]

message(sprintf("Using %d highly expressed genes as reference set with relevant functional annotations\n", 
            length(reference_genes)))

# Calculate CAI for all genes
cai_results <- calculate_cai(
  codon_counts = codon_usage,
  reference_genes = reference_genes,
  genetic_code = genetic_code_dna_long
)

# Extract CAI values and merge with expression data
cai_values <- cai_results$cai_values
w_table <- cai_results$w_table

# Merge CAI with expression and integrated data
n_pre_cai_join <- nrow(integrated_data)
integrated_data <- integrated_data |>
  left_join(cai_values, by = "Gene_name")
cat(sprintf("integrated_data: %d -> %d genes after left_join(cai_values) (CAI NA: %d)\n",
            n_pre_cai_join, nrow(integrated_data),
            sum(is.na(integrated_data$CAI))))
rm(n_pre_cai_join)

# CAI by Expression_Group: summary, Kruskal-Wallis test, Cohen's d, boxplot ----
cai_by_group <- integrated_data |>
  dplyr::group_by(Expression_Group) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_CAI = mean(CAI, na.rm = TRUE),
    median_CAI = median(CAI, na.rm = TRUE),
    sd_CAI = sd(CAI, na.rm = TRUE),
    mean_ENC = mean(ENC, na.rm = TRUE),
    .groups = "drop"
  )
print(cai_by_group)

kw_test <- kruskal.test(CAI ~ Expression_Group, data = integrated_data)
print(kw_test)

top_cai    <- integrated_data |> dplyr::filter(Expression_Group == "Top 5%")    |> pull(CAI)
middle_cai <- integrated_data |> dplyr::filter(Expression_Group == "Middle 90%") |> pull(CAI)
bottom_cai <- integrated_data |> dplyr::filter(Expression_Group == "Bottom 5%") |> pull(CAI)

if (length(top_cai) > 0 && length(middle_cai) > 0) {
  cat(sprintf("CAI Top 5%% vs Middle 90%%: d = %.3f\n", cohens_d_calc(top_cai, middle_cai)))
}
if (length(top_cai) > 0 && length(bottom_cai) > 0) {
  cat(sprintf("CAI Top 5%% vs Bottom 5%%: d = %.3f\n", cohens_d_calc(top_cai, bottom_cai)))
}
if (length(middle_cai) > 0 && length(bottom_cai) > 0) {
  cat(sprintf("CAI Middle 90%% vs Bottom 5%%: d = %.3f\n", cohens_d_calc(middle_cai, bottom_cai)))
}

p_cai_boxplot <- ggplot(integrated_data,
                        aes(x = Expression_Group, y = CAI, fill = Expression_Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C",
                               "Bottom 5%" = "#377EB8",
                               "Middle 90%" = "#999999")) +
  labs(title = "Codon Adaptation Index by Expression Level",
       subtitle = "Diamond = mean, box = median +/- IQR",
       y = "CAI (Codon Adaptation Index)",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")
ggsave("./results/CAI_by_expression_group.pdf", p_cai_boxplot, width = 8, height = 6)

plot_data_cai <- integrated_data |>
  dplyr::mutate(Exp_Group = factor(Expression_Group,
                                   levels = c("Bottom 5%", "Middle 90%", "Top 5%"))) |>
  dplyr::filter(!is.na(Exp_Group))

p_cai_median <- ggplot(plot_data_cai, aes(x = Exp_Group, y = CAI)) +
  stat_summary(fun.data = median_cl_boot,
               geom = "errorbar", width = 0.15, linewidth = 0.8, color = "black") +
  stat_summary(fun = median, geom = "point", size = 4, aes(color = Exp_Group)) +
  scale_color_manual(values = c("Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999",
                                "Top 5%" = "#E41A1C")) +
  labs(y = "CAI (Codon Adaptation Index)", x = NULL) +
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        panel.grid.major.x = element_blank())
ggsave("./results/CAI_by_expression_group_Median_CI.pdf", p_cai_median,
       width = 4, height = 3.5)

# Memory cleanup: Section 7 CAI intermediates ---
# Keeping: cai_results, w_table, cai_by_group, kw_test
rm(p_cai_boxplot, p_cai_median, plot_data_cai,
   top_cai, middle_cai, bottom_cai, reference_genes)
gc()
## ============================================================================
## RESULTS 4 — Candidate optimal codons; selection acts only on the elite
##   Produces:
##     Figure 4   Codon share across the expression range, by family
##                (`codon_share_expression_regimes.pdf`) — replaces the former
##                ROC-SEMPPR trajectory panel, which drew Referee 1's objection
##                that "most of the curve dynamics are beyond the extent of the
##                data points".  Every point here is an observed pooled share.
##     Figure 4b  Regime contrast per codon, with 95% CIs
##                (`preferred_codon_regime_contrast.pdf`)
##     Table 2    Candidate optimal codons by family
##                (`preferred_codon_detection_by_family.csv`)
##     Table S-x  Per-codon detail, elite-vs-bulk contrast, ROC concordance
##                (`preferred_codon_detection_per_codon.csv`,
##                 `preferred_codon_elite_vs_bulk.csv`,
##                 `preferred_codon_concordance_vs_ROC.csv`)
##     Cited values:  21/21 families resolved; 21/21 preferred codons GC-ending
##                     (16 C-ending, 5 G-ending)
##                     19/21 show a drift-to-selection regime reversal
##                     18/19 agreement with the archived ROC-SEMPPR call
##                     (the exception is Val: GTC here vs GTG under ROC)
##
##   SUPPLEMENTARY (AnaCoDa is archived upstream; retained as a cross-check):
##     Table S-y  Intron / intergenic stationary nucleotide frequencies
##                (`Mguttatus_intron_derived_dM.csv`, `Mguttatus_intergenic_derived_dM.csv`)
##     S_ROC / L_ROC translational-load axis (Section 8.3)
##   The AnaCoDa MCMC fit is run externally (see comments at lines 919–933 of
##   full_analysis.R); this script reads the saved posterior summaries.
## ============================================================================

## 8) Optimal-codon detection and selection coefficients ----
## _____________________________________________________________________________

# 8.2) Candidate optimal codons — expression-regime detection (PRIMARY) ----
#
# This replaces the former AnaCoDa call `which.min(Delta_eta)`.  The rationale
# is documented in full at the top of src/detect_preferred_codons.R; in brief:
#
#   * Codon share is not monotone in expression.  GC3-ending codons decline
#     steadily across the drift-dominated bulk of genes and reverse sharply in
#     the extreme upper tail, where Ns is large enough for selection on codon
#     usage to be effective.  A single linear expression slope averages the two
#     regimes and returns the drift answer (T-ending codons), which is a
#     statement about mutation bias, not optimality.
#   * The optimal codon is therefore identified by the REGIME CONTRAST: the
#     derivative of a smooth fit at the 99th expression percentile minus the
#     derivative at the median, tested against the joint model covariance.
#     A codon qualifies when its trajectory bends upward specifically where
#     selection becomes effective.
#   * GC12 (non-degenerate positions) controls regional base composition, so a
#     codon cannot be called preferred merely because its gene sits in a
#     GC-rich neighbourhood.  This is what separates codon preference from
#     GC-biased gene conversion.  The call is near-invariant to dropping it
#     (20/21 families unchanged), which is itself the answer to that referee
#     concern.
#
# Note this is the same quantity ROC-SEMPPR reports (which codon wins as
# phi -> infinity), estimated on the scale of the data instead of by
# extrapolating a sigmoid beyond the observed phi range.

preferred_detection <- detect_preferred_codons(
  codon_counts    = codon_usage,
  gene_meta       = integrated_data,
  genetic_code    = genetic_code_dna_long,
  expression_var  = "Mean_Log10_Exp",   # standalone trend; see project convention
  composition_var = "GC12",
  tail_quantile   = 0.99,
  bulk_quantile   = 0.50,
  verbose         = FALSE
)

preferred_codons <- preferred_detection$preferred |>
  dplyr::filter(Call == "preferred") |>
  dplyr::transmute(
    AA    = Family,
    aa    = Family,
    Codon = Preferred_Codon,
    # Ordering statistic with ROC's sign convention (lower = more preferred),
    # so merge_2_and_4_to_6_fold() collapses six-fold families unchanged.
    eta   = -Slope_contrast
  ) |>
  as.data.frame()

# Exporting preferred codons for polymorphism-based analysis ----

write.table(x = preferred_codons$Codon,
            file = './results/preferred_codons.txt',
            quote = F, row.names = F, col.names = F)

data.table::fwrite(preferred_detection$codon_table,
                   "./results/preferred_codon_detection_per_codon.csv")
data.table::fwrite(preferred_detection$preferred,
                   "./results/preferred_codon_detection_by_family.csv")
data.table::fwrite(preferred_detection$regime,
                   "./results/preferred_codon_elite_vs_bulk.csv")

# Unified object for downstream analyses (Wright MSD, CA/PCA biplots).
# Kept in the shape create_enhanced_biplot() expects.
preferred_codons_call <- preferred_codons |>
  dplyr::rename(Preferred_Codons = Codon,
                Amino_Acid = AA,
                Family = aa) |>
  dplyr::mutate(Source = "Expression_regime")

# Back-compatible alias for helpers that still refer to the old name.
preferred_codons_roc <- preferred_codons_call

message(sprintf("✓ Preferred codons (expression-regime): %d families resolved; 3rd base %s",
                nrow(preferred_codons_call),
                paste(names(table(substr(preferred_codons_call$Preferred_Codons, 3, 3))),
                      table(substr(preferred_codons_call$Preferred_Codons, 3, 3)),
                      sep = "=", collapse = " ")))

# --- PRIMARY VALIDATION: tRNA gene supply ------------------------------------
# Independent of both the detector and ROC-SEMPPR: tRNA gene content uses no
# expression data and no codon-usage data.  (ROC-SEMPPR is NOT independent —
# the run its coefficients come from was given the empirical expression, so its
# agreement is reproducibility, not corroboration.  See supplementary_anacoda.R.)
#
# Supply is the tAI weight (dos Reis et al. 2004), which applies the wobble
# penalties.  Read it by family size:
#   * 13 families are supported — 7 because the genome encodes a dedicated
#     Watson-Crick G34 reader for the called C-ending codon (and no inosine
#     route), 6 on actual C34-vs-T34 copy numbers;
#   * 8 are NOT supported: inosine-dominated families where the C-ending codon
#     is served only by I34. In 5 of them (Ala, Arg_4, Ile, Thr, Val) the
#     W(NNC)/W(NNT) ratio is pinned at 1-s(I:C)=0.720 by arithmetic, whatever
#     the copy numbers, so no data enters; the other 3 carry a single token G34
#     gene against 12-22 A34. Nothing is contradicted, but nothing is confirmed
#     either -- and this is exactly where the generic tAI s-values are least
#     reliable (s(I:C) fitted to yeast; I:C is held to be an efficient pairing).
#     The paper's "mostly C-ending" headline leans on four-fold families, so
#     this limitation has to be stated rather than averaged away.

trna_supply <- trna_supply_by_codon(
  trna_file    = "data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long
)

trna_validation <- validate_preferred_by_trna(
  preferred_detection$preferred, trna_supply, genetic_code_dna_long
)

data.table::fwrite(trna_supply,               "./results/trna_supply_by_codon.csv")
data.table::fwrite(trna_validation$by_family, "./results/trna_validation_by_family.csv")

cat(sprintf(
  paste0("[tRNA validation] %d/%d families supported (%d by anticodon inventory, ",
         "%d by copy number); %d not supported (inosine-dominated, of which %d have the ",
         "C:T ratio exactly pinned by 1-s(I:C)); %d contradicted.\n"),
  trna_validation$summary$n_supported, trna_validation$summary$n_families,
  trna_validation$summary$n_supported_inv, trna_validation$summary$n_supported_copy,
  trna_validation$summary$n_not_supported, trna_validation$summary$n_ratio_pinned,
  trna_validation$summary$n_contradicted
))
cat(sprintf("  Inosine-dominated (no tRNA support for the C-ending call): %s\n",
            paste(trna_validation$summary$inosine_dominated_families, collapse = ", ")))

plot_trna_supply_vs_preference(
  trna_supply[Family %in% preferred_detection$preferred$Family],
  preferred_detection$preferred,
  output_file = "./results/trna_supply_vs_preference.pdf"
)

preferred_bins <- codon_share_by_expression_bin(
  codon_usage, integrated_data, genetic_code_dna_long
)
data.table::fwrite(preferred_bins, "./results/codon_share_by_expression_bin.csv")

plot_codon_regimes(preferred_bins,
                   output_file = "./results/codon_share_expression_regimes.pdf")
plot_preferred_codon_slopes(preferred_detection$codon_table,
                            preferred_detection$preferred,
                            output_file = "./results/preferred_codon_regime_contrast.pdf")
## RESULTS 5 — Wright two-allele MSD framework, drift barrier, GO enrichment
##   Produces:
##     Figure 5   A) S_Wright distribution with drift threshold (S=1)
##                B) L_ROC density by drift / nearly neutral / selection group
##                C) S_eta density by group
##                (`Drift_barrier_overview.pdf`)
##     Table S2   GO enrichment for genes with S>1
##                (`Go_enrichment_selection_S_Wright.csv`,
##                 `Go_enrichment_load_ROC_eff.csv`)
##     Table S3   Top genes by translational load / S_Wright
##                (`Top_genes_strong_selection_load.csv`,
##                 `Top_genes_strong_selection_S_Wright.csv`)
##     Cited values:  18,273 / 22,355 genes below drift threshold (81.74%)
##                     Mean S = 0.60, median S = 0.55, range −∞ to 2.80
##                     4,082 genes with S > 1; 17,041 with 0 < S < 1; 1,232 with S < 0
## ============================================================================

# 8.3.4) Wright's MSD framework ----
#
# Three metrics presented alongside each other:
#   L_ROC:         per-gene translational load (phi-scaled |Δη|); used for GO/BGS.
#   ROC_eff:         per-gene codon-usage efficacy (-mean η at 4-fold sites); no phi.
#   S_Wright:      per-gene selection coefficient from Wright Q inversion; population-genetic units.
#
# Drift barrier (S_BARRIER):
#   Defined as the median S_Wright of genes at/above the expression inflection
#   — the first expression level where the GAM lower CI of Q exceeds Q_neutral.
#   This anchors the threshold to the biological inflection.
#
# Branch A: empirical U, V from introns
# Branch B: empirical U, V from low-expression neutral pool.
# S_BARRIER derivation follows per-gene S_Wright computation below.

# Per-gene preferred-base frequency at 4-fold sites
# Wright is a per-site framework. For each 4-fold AA family, the preferred
# 3rd-position nucleotide is set by the preferred codon detected in Section 8.2
# (e.g. Ala -> GCC -> C, Val -> GTC -> C).  Q_pref_base is then the per-gene
# fraction of 4-fold sites carrying that per-AA preferred base, summed over
# the 8 four-fold families.
#
# Note this framework needs only the IDENTITY of the preferred base: S_Wright
# is inverted from each gene's observed Q at 4-fold sites, not from any ROC
# parameter.  That is why the drift barrier, the selection group and the GO
# enrichment survive the move off AnaCoDa unchanged.

fourfold_families_msd <- c("Ala", "Gly", "Pro", "Thr", "Val",
                           "Leu_4", "Ser_4", "Arg_4")
fourfold_codons_msd <- names(genetic_code_dna_long)[
  genetic_code_dna_long %in% fourfold_families_msd
]
# In preferred_codons_call, `Amino_Acid` holds the family label (e.g. "Ala",
# "Leu_4") from genetic_code_dna_long.
preferred_per_AA <- preferred_codons_call |>
  dplyr::filter(Amino_Acid %in% fourfold_families_msd) |>
  dplyr::transmute(AA = Amino_Acid,
                   Preferred_Codon = Preferred_Codons,
                   Preferred_Base  = substr(Preferred_Codons, 3, 3))

# Arg_4 used to be absent here: AnaCoDa's preferred Arg codon (AGG) belongs to
# Arg_2, so the 4-fold Arg family had no ROC call and CGC was substituted by
# hand as a "defensible default".  The expression-regime detector resolves
# Arg_4 directly and returns CGC (regime contrast +0.43, FDR 5e-04), so the
# hand-patch is no longer needed — the value is now estimated, not assumed.
# The check is kept as a guard rather than a fallback.

missing_fourfold_families <- setdiff(fourfold_families_msd, preferred_per_AA$AA)
if (length(missing_fourfold_families) > 0) {
  stop(sprintf(
    "Missing preferred-base mapping for 4-fold families: %s",
    paste(missing_fourfold_families, collapse = ", ")
  ))
}

cat("[Wright MSD] Per-AA preferred 3rd-position base:\n")
print(preferred_per_AA, row.names = FALSE)

# Map every 4-fold codon to its family and to whether its 3rd-position base
# matches the per-AA preferred base.
fourfold_codon_table <- data.frame(
  Codon = fourfold_codons_msd,
  AA    = unname(genetic_code_dna_long[fourfold_codons_msd]),
  Base3 = substr(fourfold_codons_msd, 3, 3),
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(preferred_per_AA |> dplyr::select(AA, Preferred_Base),
                   by = "AA") |>
  dplyr::mutate(is_preferred = Base3 == Preferred_Base)

if (anyNA(fourfold_codon_table$Preferred_Base)) {
  stop("Preferred-base mapping is incomplete for one or more fourfold codons.")
}

# Per-gene preferred-base count over 4-fold sites.
codon_4fold_counts <- as.data.frame(
  codon_usage[, c("Gene_name", fourfold_codons_msd), with = FALSE]
)
preferred_codon_set <- fourfold_codon_table$Codon[fourfold_codon_table$is_preferred]
N_4fold_sites    <- rowSums(codon_4fold_counts[, fourfold_codons_msd])
N_preferred_base <- rowSums(codon_4fold_counts[, preferred_codon_set, drop = FALSE])

gene_Q_4fold <- data.frame(
  Gene_name        = codon_4fold_counts$Gene_name,
  N_4fold_sites    = N_4fold_sites,
  N_preferred_base = N_preferred_base,
  Q_pref_base      = ifelse(N_4fold_sites > 0, N_preferred_base / N_4fold_sites, NA_real_)
)

# Merge codon-derived Q with integrated data for MSD analysis.
# 
# IMPORTANT: Two independent denominator sources:
#   1. N_4fold_sites:  Directly from codon counts (FASTA data)
#                      Transparent, authoritative for Q calculations.
#   2. Sites_4fold:    From polymorphism data (VCF) paired with Pi_sum_4fold
#                      Used ONLY for π calculations where available.
#
# These may differ due to:
#   - Different genomic regions (VCF may have missing data at some loci)
#   - Different site definitions (frame filters, quality thresholds)
#   - Different coverage (polymorphism data has inherent missingness)
#
# Strategy:
#   - Use N_4fold_sites (codon-derived) for Q and Q_bin calculations
#   - Use Sites_4fold (VCF-derived) with Pi_sum_4fold for π_bin
#   - For per-gene π (Pi_mean_4fold), use whichever denominator is available

n_pre_msd_join <- nrow(integrated_data)
n_gene_Q_4fold <- nrow(gene_Q_4fold)
msd_data <- integrated_data |>
  dplyr::select(Gene_name, Pi_mean_4fold, Mean_Log10_Exp,
                Max_Log10_Exp, Exp_breadth, CDS_length_nt, Sites_4fold,
                Pi_sum_4fold) |>
  dplyr::inner_join(gene_Q_4fold, by = "Gene_name") |>
  dplyr::filter(!is.na(Q_pref_base), !is.na(Pi_mean_4fold),
                !is.na(Mean_Log10_Exp), !is.na(Exp_breadth),
                N_4fold_sites >= 20)
cat(sprintf("msd_data: integrated_data %d x gene_Q_4fold %d -> %d genes after inner_join + NA/Site filters\n",
            n_pre_msd_join, n_gene_Q_4fold, nrow(msd_data)))
rm(n_pre_msd_join, n_gene_Q_4fold)

# Diagnostic: quantify discrepancies and document
n_4fold_match <- sum(msd_data$N_4fold_sites == msd_data$Sites_4fold, na.rm = TRUE)
n_4fold_total <- nrow(msd_data)
pct_match <- 100 * n_4fold_match / n_4fold_total

if (pct_match < 100) {
  cat(sprintf(
    "[Wright MSD] INFO: 4-fold site count discrepancy between codon and VCF data detected.\n"
  ))
  cat(sprintf(
    "  %d of %d genes (%.1f%%) have N_4fold_sites == Sites_4fold.\n",
    n_4fold_match, n_4fold_total, pct_match
  ))
  
  # Summarize differences
  msd_data$site_count_diff <- msd_data$N_4fold_sites - msd_data$Sites_4fold
  n_discrepant <- sum(msd_data$site_count_diff != 0, na.rm = TRUE)
  
  if (n_discrepant > 0) {
    disc_subset <- msd_data |>
      dplyr::filter(site_count_diff != 0)
    
    cat(sprintf(
      "  %d genes show differences (mean: %.1f sites, max: %.0f sites, median: %.0f sites).\n",
      n_discrepant,
      mean(abs(disc_subset$site_count_diff), na.rm = TRUE),
      max(abs(disc_subset$site_count_diff), na.rm = TRUE),
      median(abs(disc_subset$site_count_diff), na.rm = TRUE)
    ))
    cat(sprintf(
      "  Mean relative difference: %.1f%% of N_4fold_sites.\n\n",
      mean(abs(disc_subset$site_count_diff) / pmax(disc_subset$N_4fold_sites, 1) * 100, na.rm = TRUE)
    ))
  }
  
  cat(sprintf(
    "  RESOLUTION: Using N_4fold_sites (codon-derived) for Q and Q_bin.\n"
  ))
  cat(sprintf(
    "  Pi calculations use Sites_4fold (VCF-derived, paired with Pi_sum_4fold).\n\n"
  ))
}

cat(sprintf(
  "[Wright MSD] %d genes with usable 4-fold per-base Q and pi (>=20 4-fold sites).\n",
  nrow(msd_data)
))

# Generate neutral mutation parameters from intronic SFS ----
# Fits the two-allele Wright model to the intronic site-frequency spectrum to
# recover alpha (4N*u toward C/G) and beta (4N*v away from C/G) for each
# nucleotide system.  Output is written to results/ so Branch A can read it;
# if either SFS file is missing the write is skipped and Branch A falls back
# to whatever file already exists on disk (or NA if none).
#
# SFS GENERATION (when sfs_introns_{G,C}.csv are missing):
#   The intronic SFS files are produced from a population VCF + GFF3 by
#     miscellanea_code/filter_vcf_for_introns.py
#   wrapped for SLURM at
#     Bash_scripts/run_sfs_introns.sh <vcf.gz> <gff3> <out_dir>
#   Below, if both SFS CSVs are missing AND a local VCF is reachable (via
#   the SFS_INTRONS_VCF env var, e.g. `export SFS_INTRONS_VCF=...`), we
#   invoke the Python script directly so the rest of main.R can proceed
#   end-to-end on a fresh dataset. On the cluster, prefer the sbatch path.

sfs_G_file <- "./data/sfs_introns_G.csv"
sfs_C_file <- "./data/sfs_introns_C.csv"

if (!file.exists(sfs_G_file) || !file.exists(sfs_C_file)) {
  sfs_vcf_env <- Sys.getenv("SFS_INTRONS_VCF", unset = "")
  gff_path    <- "./data/Mguttatusvar_IM767_887_v2.1.gene.gff3"
  py_path     <- "./miscellanea_code/filter_vcf_for_introns.py"
  if (nzchar(sfs_vcf_env) && file.exists(sfs_vcf_env) &&
      file.exists(gff_path) && file.exists(py_path)) {
    cat(sprintf("[SFS gen] Building intronic SFS from %s ...\n", sfs_vcf_env))
    stream_cmd <- if (grepl("\\.gz$", sfs_vcf_env)) "zcat" else "cat"
    cmd <- sprintf("%s %s | python3 %s --stream --gff %s",
                   stream_cmd, shQuote(sfs_vcf_env),
                   shQuote(py_path), shQuote(gff_path))
    exit_code <- system(cmd)
    if (exit_code == 0) {
      # Python writes sfs_introns_{G,C}.csv to CWD; move them under ./data/.
      for (sfs_fname in c("sfs_introns_G.csv", "sfs_introns_C.csv")) {
        if (file.exists(sfs_fname)) {
          file.rename(sfs_fname, file.path("data", sfs_fname))
        }
      }
      cat("[SFS gen] sfs_introns_{G,C}.csv written to ./data/.\n")
    } else {
      cat(sprintf("[SFS gen] python3 invocation failed (exit %d).\n", exit_code))
    }
    rm(stream_cmd, cmd, exit_code)
  } else {
    cat("[SFS gen] SFS CSVs missing and no local VCF available.\n")
    cat("  To regenerate: sbatch Bash_scripts/run_sfs_introns.sh <vcf.gz> <gff3> data/\n")
    cat("  Or set env var SFS_INTRONS_VCF=/path/to/vcf.gz before sourcing main.R.\n")
  }
  rm(sfs_vcf_env, gff_path, py_path)
}

if (file.exists(sfs_G_file) && file.exists(sfs_C_file)) {
  neutral_params_sfs <- load_and_estimate_neutral_params(sfs_G_file, sfs_C_file)

  neutral_params_df <- data.frame(
    Parameter = c("alpha_G", "beta_G", "alpha_C", "beta_C",
                  "pi_G_expected", "pi_C_expected"),
    Value = c(neutral_params_sfs$alpha_G, neutral_params_sfs$beta_G,
              neutral_params_sfs$alpha_C, neutral_params_sfs$beta_C,
              neutral_params_sfs$pi_G_expected, neutral_params_sfs$pi_C_expected),
    Description = c("4N*u for G (unpreferred->preferred)",
                    "4N*v for G (preferred->unpreferred)",
                    "4N*u for C (unpreferred->preferred)",
                    "4N*v for C (preferred->unpreferred)",
                    "Expected nucleotide diversity at G sites",
                    "Expected nucleotide diversity at C sites")
  )

  write.csv(neutral_params_df, "./results/neutral_mutation_parameters.csv",
            row.names = FALSE)
  cat("[SFS fit] neutral_mutation_parameters.csv written.\n")
  rm(neutral_params_sfs, neutral_params_df)
} else {
  cat("[SFS fit] One or both SFS files not found; skipping regeneration.\n")
  cat(sprintf("  G: %s\n  C: %s\n", sfs_G_file, sfs_C_file))
}
rm(sfs_G_file, sfs_C_file)

# Branch A: estimate U, V from SFS-derived intronic neutral mutation parameters ----
# The VCF-based two-allele approach averages Q and pi over polymorphic sites only
# (invariant sites are absent from VCF records). This inflates both estimates:
# the resulting conditional Q (~0.32) is far from the observed intronic C content
# (~0.16), and pi is inflated ~5x, producing U and V that are biologically wrong.
# The SFS fit uses the full allele-frequency spectrum at intronic sites; its
# implied Q = V/(U+V) is consistent with observed nucleotide composition and is
# the correct calibration source.
#
# Convention in neutral_mutation_parameters.csv:
#   alpha_X = 4N*u for X  (unpreferred->preferred, i.e., rate TOWARD C/G)
#   beta_X  = 4N*v for X  (preferred->unpreferred, i.e., rate AWAY FROM C/G)
#
# These are the Beta-distribution shape parameters for the DIPLOID Wright
# stationary: f(x) ∝ x^(V-1) (1-x)^(U-1), with E[x] = V/(U+V) and
# E[2x(1-x)] = 2UV/((U+V)(U+V+1)). They must be passed to wright_pi/wright_Q
# DIRECTLY. An earlier version of this block divided by 2 (treating them as
# 2N-scaled), which left Q unchanged (V/(U+V) is scale-invariant) but halved
# pi_neutral_theory — producing a ~2x mismatch against observed intronic pi.
# Validated below: with the direct assignment, pi_theory matches the observed
# unconditional intronic 2-allele pi within ~7%.

neutral_params <- tryCatch(
  read.csv("results/neutral_mutation_parameters.csv"),
  error = function(e) NULL
)

if (!is.null(neutral_params)) {
  get_np <- function(nm) {
    v <- neutral_params$Value[neutral_params$Parameter == nm]
    if (length(v)) v[1L] else NA_real_
  }

  alpha_C <- get_np("alpha_C"); beta_C <- get_np("beta_C")
  alpha_G <- get_np("alpha_G"); beta_G <- get_np("beta_G")

  V_intron <- alpha_C   # 4N*u toward preferred (C)
  U_intron <- beta_C    # 4N*v away from preferred (C)

  Q_implied_C <- V_intron / (U_intron + V_intron)
  Q_implied_G <- alpha_G / (alpha_G + beta_G)

  cat(sprintf("[Branch A] U_intron = %.6f, V_intron = %.6f, V/U = %.3f\n",
              U_intron, V_intron, V_intron / U_intron))
  cat(sprintf(
    "[Branch A] Implied Q_C = %.4f  (cross-check: observed intronic C content ~0.16)\n",
    Q_implied_C))
  cat(sprintf(
    "[Branch A - G approx check] Implied Q_G = %.4f  (C used for both systems)\n",
    Q_implied_G))

  # --- Validation: theoretical Q and pi at S=0 vs observed intronic values ----
  # Q is scale-invariant in (U, V); pi is not. The fix above (no /2) is verified
  # by comparing wright_pi(0, U, V) to the observed unconditional intronic
  # 2-allele pi (sum of polymorphic pi over total intron bp).
  intron_two_path <- "data/intron_2_allele.csv"
  feat_pi_path    <- "data/all_chromosomes.pi_per_gene_feature.txt"
  if (file.exists(intron_two_path) && file.exists(feat_pi_path)) {
    two_a <- data.table::fread(intron_two_path)
    pf_a  <- data.table::fread(feat_pi_path)
    pf_a[, Gene_norm := gsub("^MgIM767\\.", "", Gene)]
    pf_a[, Gene_norm := gsub("\\.v2\\.1$",   "", Gene_norm)]
    intron_bp <- pf_a[Feature_Type == "intron" & Degeneracy == "all",
                      .(total_intron_bp = sum(Sites)), by = Gene_norm]
    data.table::setnames(intron_bp, "Gene_norm", "Gene")
    cal <- merge(two_a, intron_bp, by = "Gene")

    # Observed unconditional intronic pi under the 2-allele system
    pi_C_obs_intron <- sum(cal$pi_2allele_C * cal$n_sites) /
                       sum(cal$total_intron_bp)
    pi_G_obs_intron <- sum(cal$pi_2allele_G * cal$n_sites) /
                       sum(cal$total_intron_bp)

    # Theoretical values from current (post-fix) U, V
    Q_theory_C  <- V_intron / (U_intron + V_intron)
    pi_theory_C <- wright_pi(0, U = U_intron, V = V_intron)
    pi_theory_G <- wright_pi(0, U = beta_G,   V = alpha_G)

    intron_C_freq_known <- 0.16   # documented intronic C content (project memory)
    intron_G_freq_known <- 0.17   # documented intronic G content (project memory)

    cat("\n[Branch A - validation] expectation vs observation at S = 0:\n")
    cat(sprintf("  Q_C : theory = %.4f  | observed intron C content ~ %.2f\n",
                Q_theory_C, intron_C_freq_known))
    cat(sprintf("  Q_G : theory = %.4f  | observed intron G content ~ %.2f\n",
                Q_implied_G, intron_G_freq_known))
    cat(sprintf("  pi_C: theory = %.6f | observed = %.6f  (ratio %.2fx)\n",
                pi_theory_C, pi_C_obs_intron, pi_C_obs_intron / pi_theory_C))
    cat(sprintf("  pi_G: theory = %.6f | observed = %.6f  (ratio %.2fx)\n",
                pi_theory_G, pi_G_obs_intron, pi_G_obs_intron / pi_theory_G))

    rm(two_a, pf_a, intron_bp, cal, pi_C_obs_intron, pi_G_obs_intron,
       Q_theory_C, pi_theory_C, pi_theory_G,
       intron_C_freq_known, intron_G_freq_known,
       intron_two_path, feat_pi_path)
  } else {
    cat("[Branch A - validation] intron_2_allele.csv or pi_per_gene_feature.txt missing; pi validation skipped.\n")
  }

  rm(neutral_params, get_np, alpha_C, beta_C, alpha_G, beta_G,
     Q_implied_C, Q_implied_G)
} else {
  U_intron <- NA_real_; V_intron <- NA_real_
  cat("[Branch A] results/neutral_mutation_parameters.csv not found; intron calibration skipped.\n")
}
# Branch B: estimate U, V from a low-expression near-neutral pool ----
# Use the bottom expression decile (rather than bottom L_ROC quartile) so
# the "neutral" sample is selected on a covariate independent of the L_ROC
# scale. Site-weight Q for noise reduction; weight pi by 4-fold site counts.

exp_q10 <- quantile(msd_data$Mean_Log10_Exp, 0.10, na.rm = TRUE)
neutral_pool <- msd_data |> dplyr::filter(Mean_Log10_Exp <= exp_q10)

Q_neutral_obs  <- with(neutral_pool,
                       sum(N_preferred_base) / sum(N_4fold_sites))
pi_neutral_obs <- with(neutral_pool,
                       sum(Pi_sum_4fold,   na.rm = TRUE) /
                       sum(Sites_4fold,    na.rm = TRUE))

cat(sprintf(
  "[Branch B] Near-neutral pool (Mean_Log10_Exp <= %.3f, n = %d): Q_obs = %.3f, pi_obs = %.4f\n",
  exp_q10, nrow(neutral_pool), Q_neutral_obs, pi_neutral_obs
))

UV_emp <- wright_solve_UV(Q_neutral_obs, pi_neutral_obs)
U_emp  <- UV_emp["U"]; V_emp <- UV_emp["V"]
cat(sprintf("[Branch B] Empirical U = %.4f, V = %.4f, V/U = %.3f\n",
            U_emp, V_emp, V_emp / U_emp))

S_grid_emp <- seq(0, 8, length.out = 200)
wright_emp <- data.frame(
  S       = S_grid_emp,
  Q       = wright_Q(S_grid_emp,  U = U_emp, V = V_emp),
  pi_site = wright_pi(S_grid_emp, U = U_emp, V = V_emp)
)

# Two-state π calibration (preferred vs unpreferred sites) ----
# Use the same neutral_pool (low-expression decile) but pull the
# two-state per-gene heterozygosity `pi_2allele` from the precomputed
# `data/Two_allele_pi.csv`.  Solve (U, V) from Q_neutral_obs and the
# site-weighted `pi_neutral_two` when feasible and keep as a fallback
# reference alongside the original empirical (U_emp, V_emp).
# File has 5 columns: Gene, n_pref_notpref, q_pref, p_notpref, pi_2allele.
# Properly select and rename — old 2-name colnames() assignment silently
# mapped column 2 (site counts) to pi_2allele, discarding q_pref entirely.
pi_data_operational <- tryCatch({
  d <- read.csv("data/pi_operational.csv")
  dplyr::mutate(d,
    Gene_name = paste0("MgIM767.", Gene)
  )
}, error = function(e) NULL)

pi_data_operational <- pi_data_operational |>
  dplyr::filter(MeanLog10_exp != 0)

if (!is.null(pi_data_operational)) {
  # Drop any columns from neutral_pool that pi_data_operational is about to
  # contribute. Otherwise, on a rerun within the same R session, msd_data
  # already carries pi_2allele (joined in below at the "Join two-allele pi"
  # block) and the inner_join here would produce pi_2allele.x / pi_2allele.y,
  # breaking the bare `pi_2allele` reference in the filter.
  op_cols <- setdiff(names(pi_data_operational), "Gene_name")
  neutral_pool_pi <- neutral_pool |>
    dplyr::select(-tidyselect::any_of(op_cols)) |>
    dplyr::inner_join(pi_data_operational, by = "Gene_name") |>
    dplyr::filter(!is.na(pi_2allele), !is.na(q_pref))
  if (nrow(neutral_pool_pi) > 0) {
    # Q_neutral from two-allele preferred frequency, weighted by two-allele site count.
    Q_neutral_two  <- with(neutral_pool_pi,
                           weighted.mean(q_pref, n_pref_notpref, na.rm = TRUE))
    # pi from two-allele heterozygosity, weighted by two-allele site count.
    pi_neutral_two <- with(neutral_pool_pi,
                           weighted.mean(pi_2allele, n_pref_notpref, na.rm = TRUE))
  } else {
    Q_neutral_two  <- NA_real_
    pi_neutral_two <- NA_real_
  }
} else {
  Q_neutral_two  <- NA_real_
  pi_neutral_two <- NA_real_
}

# Diagnostic output for two-state calibration
if (is.null(pi_data_operational)) {
  cat("[Branch B - two-state] data/Two_allele_pi.csv not found or unreadable.\n")
} else {
  n_pi <- if (exists("neutral_pool_pi")) nrow(neutral_pool_pi) else 0
  cat(sprintf("[Branch B - two-state] matched genes in neutral pool: %d\n", n_pi))
  cat(sprintf("[Branch B - two-state] Q_neutral_two = %s, pi_neutral_two = %s\n",
              ifelse(is.na(Q_neutral_two),  "NA", format(Q_neutral_two,  digits = 6)),
              ifelse(is.na(pi_neutral_two), "NA", format(pi_neutral_two, digits = 6))))
}

# Validate and solve for (U, V) using the two-state Q and π.
if (is.finite(Q_neutral_two) && is.finite(pi_neutral_two) && pi_neutral_two > 0) {
  hardy_max_two <- 2 * Q_neutral_two * (1 - Q_neutral_two)
  cat(sprintf("[Branch B - two-state] Q_neutral_two = %.6f, Hardy_max = %.6f\n",
              Q_neutral_two, hardy_max_two))
  if (pi_neutral_two < hardy_max_two) {
    UV_emp_two <- wright_solve_UV(Q_neutral_two, pi_neutral_two)
    U_emp_two <- UV_emp_two["U"]; V_emp_two <- UV_emp_two["V"]
    cat(sprintf("[Branch B - two-state π] Empirical U2 = %.6f, V2 = %.6f\n", U_emp_two, V_emp_two))
    cat("[Branch B - two-state] two-state UV calibration SUCCESS; using U_emp_two/V_emp_two for π→S inversions.\n")
  } else {
    U_emp_two <- NA_real_; V_emp_two <- NA_real_
    warning("two-state pi_neutral is outside Hardy bound; skipping two-state UV solve")
    cat(sprintf("[Branch B - two-state] pi_neutral_two = %.6f >= Hardy_max = %.6f; skipping solve.\n",
                pi_neutral_two, hardy_max_two))
  }
} else {
  U_emp_two <- NA_real_; V_emp_two <- NA_real_
}
# S_BARRIER: median S_Wright of genes at/above the Q-vs-expression inflection.
#   Genes with S_Wright_signed >= S_BARRIER are in the "selection group".
# thr_sel: L_ROC of the 50th-highest gene (top-50 load group for GO/BGS).
#
# Per-gene S_Wright: inverted from each gene's Q at 4-fold sites via
#   wright_invert_Q(Q; U_emp, V_emp).  Sign-aware version kept as
#   S_Wright_signed; genes below neutral Q are flagged is_drift = TRUE.

# Theoretical neutral-pi reference ----
# S_BARRIER: fixed at the 2N_e*s = 1 threshold (selection dominates drift) ----
# The GAM-inflection derived value (~0.54) classified genes where drift is still
# the dominant force as being under selection. S = 1 is the canonical boundary
# where selection overcomes drift. The inflection-based derivation is archived.

S_BARRIER        <- 1
S_BARRIER_source <- "fixed_4Ns_gt_1"

pi_neutral_theory <- wright_pi(0,
  U = if (is.finite(U_intron)) U_intron else U_emp,
  V = if (is.finite(V_intron)) V_intron else V_emp
)
S_BARRIER_advisor <- 0.1   # JK's Mathematica reference; kept for comparison.

cat(sprintf(
  "\n[Wright MSD] pi_neutral_theory = %.5f, pi_neutral_obs = %.5f\n",
  pi_neutral_theory, pi_neutral_obs
))

# Per-gene S_Wright ----
# Compute the SIGNED inversion first (full-information diagnostic), then
# floor at zero for the operational column used downstream.
# Use two-state calibration if available, otherwise fall back to original calibration.

# This is the entry point of the U and V parameters. Several options are available
# to explore. For compatibility, we are going to adopt intron-based if present. 
# The second source in hierarchy are the parameters derived from the two-allele
# system. And finally, the approximated parameters based on the regular pi.

U_gene_calib <- if (is.finite(U_intron) && is.finite(V_intron)) U_intron else
                if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) U_emp_two else U_emp
V_gene_calib <- if (is.finite(U_intron) && is.finite(V_intron)) V_intron else
                if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) V_emp_two else V_emp
cat(sprintf("[UV calib] U_gene_calib = %.6f, V_gene_calib = %.6f  [source: %s]\n",
            U_gene_calib, V_gene_calib,
            if (is.finite(U_intron) && is.finite(V_intron)) "Branch A (intron)" else
            if (exists("U_emp_two") && is.finite(U_emp_two)) "Branch B (two-state)" else
            "Branch B (raw empirical)"))

msd_data$S_Wright_signed <- vapply(msd_data$Q_pref_base, function(q) {
  tryCatch(wright_invert_Q(q, U = U_gene_calib, V = V_gene_calib),
           error = function(e) NA_real_)
}, numeric(1))
msd_data$is_drift     <- !is.na(msd_data$S_Wright_signed) &
                          msd_data$S_Wright_signed < S_BARRIER
msd_data$S_Wright_raw <- pmax(msd_data$S_Wright_signed, 0)   # operational

cat(sprintf(
  "\n[Wright MSD] Per-gene S_Wright_signed: min = %.3f, q05 = %.3f, median = %.3f, q95 = %.3f, max = %.3f, NA = %d / %d\n",
  min(msd_data$S_Wright_signed, na.rm = TRUE),
  quantile(msd_data$S_Wright_signed, 0.05, na.rm = TRUE),
  median(msd_data$S_Wright_signed, na.rm = TRUE),
  quantile(msd_data$S_Wright_signed, 0.95, na.rm = TRUE),
  max(msd_data$S_Wright_signed, na.rm = TRUE),
  sum(is.na(msd_data$S_Wright_signed)), nrow(msd_data)
))
cat(sprintf(
  "[Wright MSD] %d / %d genes (%.1f%%) flagged is_drift (Q < Q_neutral; floored to S_Wright_raw = 0).\n",
  sum(msd_data$is_drift, na.rm = TRUE), nrow(msd_data),
  100 * mean(msd_data$is_drift, na.rm = TRUE)
))

# Join two-allele pi into msd_data ----
# pi_2allele (per-gene heterozygosity under the two-allele model) was read
# during the two-state UV calibration above (Branch B). Join here so it is
# available for downstream Wright comparisons alongside S_Wright_signed.
if (!is.null(pi_data_operational)) {
  # Defensive: drop any pre-existing pi_2allele on msd_data so reruns in the
  # same R session do not produce .x/.y suffixes from the left_join below.
  msd_data <- msd_data |>
    dplyr::select(-tidyselect::any_of("pi_2allele")) |>
    dplyr::left_join(pi_data_operational |> dplyr::select(Gene_name, pi_2allele),
                     by = "Gene_name")
}


# Binned S_Wright table ----
# 30 site-weighted ntile bins used for the pi-consistency validation and the
# diversity-hump figure (bin_sw). Use two-state calibration if available.
# The former ROC_eff-binned companion table (bin_roc) moved to
# supplementary_anacoda.R along with the rest of the AnaCoDa material; the
# pi-consistency check now runs on the S_Wright bins, which is the more apt
# test anyway — it validates Wright predictions against Wright-binned data.

U_bin_calib <- if (is.finite(U_intron) && is.finite(V_intron)) U_intron else
               if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) U_emp_two else U_emp
V_bin_calib <- if (is.finite(U_intron) && is.finite(V_intron)) V_intron else
               if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) V_emp_two else V_emp


# Also build S_Wright-binned table for the two-panel diversity hump
# Use the SIGNED S_Wright here so negative values (mutation-dominated
# genes) are preserved in the binning for the diversity-hump diagnostic.
bin_sw <- msd_data |>
  dplyr::filter(!is.na(S_Wright_signed), !is.na(Q_pref_base)) |>
  dplyr::arrange(S_Wright_signed) |>
  dplyr::mutate(Ssw_bin = ntile(S_Wright_signed, 30)) |>
  dplyr::group_by(Ssw_bin) |>
  dplyr::summarize(
    n_genes       = dplyr::n(),
    mean_S_Wright = mean(S_Wright_signed),
    sites_total   = sum(N_4fold_sites),
    Q_bin         = sum(N_preferred_base) / sum(N_4fold_sites),
    pi_bin        = sum(Pi_sum_4fold, na.rm = TRUE) /
                    sum(Sites_4fold,  na.rm = TRUE),
    pi_se = sd(sum(Pi_sum_4fold, na.rm = TRUE) /
               sum(Sites_4fold,  na.rm = TRUE)) / sqrt(n()),
    .groups = "drop"
  )
bin_sw$pi_se <- sqrt(bin_sw$pi_bin * (1 - bin_sw$pi_bin / 2) /
                     pmax(bin_sw$sites_total, 1))

# S_Wright inverted from each bin's pooled Q, for the pi-consistency check
# below. (mean_S_Wright above is the mean of per-gene values; S_Wright_bin is
# the bin-level inversion, which is what wright_pi() must be evaluated at.)
bin_sw$S_Wright_bin <- vapply(bin_sw$Q_bin, function(q) {
  tryCatch(wright_invert_Q(q, U = U_bin_calib, V = V_bin_calib),
           error = function(e) NA_real_)
}, numeric(1))

cat(sprintf(
  "    S_BARRIER = %.4f; genes with S_Wright_raw >= S_BARRIER: %d\n",
  S_BARRIER,
  sum(msd_data$S_Wright_raw >= S_BARRIER, na.rm = TRUE)
))

# Bin-level pi consistency: parameter-free Wright validation ----
# pi_pred_wright is wright_pi evaluated at S_Wright_bin (already inverted from
# Q_bin). chi^2 = sum((pi_bin - pi_pred_wright)^2 / pi_se^2). No fitted
# parameter, so df = number of bins. The selection-regime subset isolates the
# regime where the Wright model is identifiable.
bin_sw$pi_pred_wright <- wright_pi(bin_sw$S_Wright_bin,
                                   U = U_bin_calib, V = V_bin_calib)
bin_sw$pi_residual    <- bin_sw$pi_bin - bin_sw$pi_pred_wright

chi2_pi_terms <- (bin_sw$pi_residual)^2 /
                 pmax(bin_sw$pi_se, .Machine$double.eps)^2
chi2_pi_stat  <- sum(chi2_pi_terms, na.rm = TRUE)
chi2_pi_df    <- sum(!is.na(chi2_pi_terms))
chi2_pi_p     <- pchisq(chi2_pi_stat, df = chi2_pi_df, lower.tail = FALSE)
cat(sprintf(
  "[Validation] Bin-level pi consistency (S_Wright bins): chi^2 = %.2f / df = %d -> p = %.3g\n",
  chi2_pi_stat, chi2_pi_df, chi2_pi_p
))

sel_bins <- bin_sw |>
  dplyr::filter(!is.na(S_Wright_bin), S_Wright_bin >= S_BARRIER)
chi2_pi_sel_stat <- if (nrow(sel_bins) > 0) {
  sum((sel_bins$pi_residual)^2 /
      pmax(sel_bins$pi_se, .Machine$double.eps)^2, na.rm = TRUE)
} else NA_real_
chi2_pi_sel_df <- nrow(sel_bins)
chi2_pi_sel_p  <- if (chi2_pi_sel_df > 0 && is.finite(chi2_pi_sel_stat)) {
  pchisq(chi2_pi_sel_stat, df = chi2_pi_sel_df, lower.tail = FALSE)
} else NA_real_
cat(sprintf(
  "[Validation] Selection-regime subset (S_Wright_bin >= %.4f): chi^2 = %s / df = %d -> p = %s\n",
  S_BARRIER,
  if (is.finite(chi2_pi_sel_stat)) sprintf("%.2f", chi2_pi_sel_stat) else "NA",
  chi2_pi_sel_df,
  if (is.finite(chi2_pi_sel_p))    sprintf("%.3g", chi2_pi_sel_p)    else "NA"
))


# 8.3.5) Drift-barrier overview ----
#
# S_Wright_signed histogram, filled by selection / nearly-neutral / drift group.
# The vertical dashed line marks S_BARRIER.
#
# Logic: the barrier is derived from Q-inflection (expression level where
# selection becomes detectable). Genes at/above inflection are classified as
# "selection" if S_Wright_raw >= S_BARRIER.
#
# The former Panels B and C — L_ROC and S_ROC densities split by this same
# S_Wright classification — are the AnaCoDa cross-metric validation and now
# live in supplementary_anacoda.R, which reads the handoff object written at
# the end of this section and reproduces the full three-panel figure.
# `SW_group` is exported there so the classification is identical.

plot_barrier <- msd_data |>
  dplyr::filter(is.finite(S_Wright_signed)) |>
  dplyr::mutate(
    SW_group = dplyr::if_else(S_Wright_signed >= S_BARRIER, 
                              "Selection", 
                              ifelse(S_Wright_signed < S_BARRIER & S_Wright_signed >= 0,
                                     "Nearly neutral", "Drift"))
  )

n_sel_barrier <- sum(plot_barrier$SW_group == "Selection")
n_nearly_neutral <- sum(plot_barrier$SW_group == "Nearly neutral")
n_drift_barrier <- sum(plot_barrier$SW_group == "Drift")

barrier_colors <- c("Selection" = "#E41A1C", "Nearly neutral" = "gray", 
                    "Drift" = "#377EB8")

# Panel A: S_Wright histogram coloured by group
p_sw_dist <- ggplot(plot_barrier, aes(x = S_Wright_signed, fill = SW_group)) +
  # Removed position="dodge" and color="white" to create a continuous histogram
  geom_histogram(
    binwidth = 0.05, 
    boundary = 0, 
    position = "stack"
  ) + 
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


# cairo_pdf is not always available on headless compute nodes, and this call
# sits BEFORE the handoff write below — a missing device would otherwise kill
# the run before supplementary_anacoda.R has anything to read. Fall back to the
# default pdf device, which changes the backend only, not the figure.
.dev_pdf <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave("./results/Drift_barrier_distribution.pdf",
       p_sw_dist, width = 8, height = 5, device = .dev_pdf)

write.csv(wright_emp,
          "./results/Wright_curve_empirical.csv",        row.names = FALSE)
write.csv(
  msd_data |> dplyr::select(Gene_name, Q_pref_base,
                            S_Wright_signed, S_Wright_raw, is_drift,
                            pi_2allele,
                            Mean_Log10_Exp, Max_Log10_Exp, N_4fold_sites),
  "./results/Wright_per_gene_S_Wright.csv", row.names = FALSE
)

# Genes clearing the drift barrier. This is the AnaCoDa-free selection group
# used by Section 12's isolation analyses (it replaces the former
# `L_ROC > thr_sel` top-50-by-load set, which moved to supplementary_anacoda.R).
# The two are close in size, but this one is defined by efficacy of selection
# rather than by the load being paid.
selection_gene_set <- msd_data |>
  dplyr::filter(!is.na(S_Wright_signed), S_Wright_signed >= S_BARRIER) |>
  dplyr::pull(Gene_name)
cat(sprintf("[Selection group] %d genes with S_Wright_signed >= %.4f\n",
            length(selection_gene_set), S_BARRIER))

# --- Handoff for supplementary_anacoda.R -------------------------------------
# The AnaCoDa cross-metric validation (S_ROC / L_ROC densities against this
# S_Wright classification, and their per-gene correlation) needs exactly the
# objects below.  Written here so the supplementary script never has to re-run
# any part of the CUB pipeline.
saveRDS(
  list(
    msd_data        = msd_data,
    plot_barrier    = plot_barrier,   # carries SW_group, the shared classification
    bin_sw          = bin_sw,
    S_BARRIER       = S_BARRIER,
    U_bin_calib     = U_bin_calib,
    V_bin_calib     = V_bin_calib,
    custom_bag      = integrated_data |> dplyr::pull(Gene_name),
    integrated_data = integrated_data,
    codon_usage     = codon_usage,
    exp_complete    = exp_complete,
    genetic_code    = genetic_code_dna_long,
    preferred       = preferred_detection$preferred,
    barrier_colors  = c("Selection" = "#E41A1C", "Nearly neutral" = "gray",
                        "Drift" = "#377EB8")
  ),
  "./results/cub_handoff_for_anacoda.rds"
)

rm(plot_barrier, n_sel_barrier, n_drift_barrier,
   barrier_colors, p_sw_dist)
write.csv(
  data.frame(
    criterion                    = "fixed_2Ns_gt_1",
    S_BARRIER                    = S_BARRIER,
    S_BARRIER_advisor            = S_BARRIER_advisor,
    U_empirical                  = U_emp,
    V_empirical                  = V_emp,
    Q_neutral_obs                = Q_neutral_obs,
    pi_neutral_obs               = pi_neutral_obs,
    pi_neutral_theory            = pi_neutral_theory,
    n_above_S_BARRIER            = sum(msd_data$S_Wright_signed >= S_BARRIER,
                                       na.rm = TRUE),
    n_drift_genes                = sum(msd_data$is_drift, na.rm = TRUE),
    frac_drift_genes             = mean(msd_data$is_drift, na.rm = TRUE),
    chi2_pi_stat                 = chi2_pi_stat,
    chi2_pi_df                   = chi2_pi_df,
    chi2_pi_p                    = chi2_pi_p,
    chi2_pi_sel_stat             = if (is.finite(chi2_pi_sel_stat)) chi2_pi_sel_stat else NA_real_,
    chi2_pi_sel_df               = chi2_pi_sel_df,
    chi2_pi_sel_p                = if (is.finite(chi2_pi_sel_p))    chi2_pi_sel_p    else NA_real_
  ),
  "./results/Wright_threshold_adopted.csv", row.names = FALSE
)

# Memory cleanup -- keep: U_emp, V_emp, Q_neutral_obs,
# pi_neutral_obs, pi_neutral_theory, S_BARRIER, S_BARRIER_advisor,
# msd_data, bin_sw, integrated_data.
rm(codon_4fold_counts, N_4fold_sites, N_preferred_base, gene_Q_4fold,
   preferred_codon_set, fourfold_codon_table, preferred_per_AA,
   S_grid_emp,
   neutral_pool, neutral_pool_pi, pi_data_operational,
   Q_neutral_two, pi_neutral_two,
   chi2_pi_terms, sel_bins)
gc()

# 8.4) GO-enrichment for the selection group ----
#
#   S_Wright_signed >= S_BARRIER -> "selection" group (drift-barrier genes)
#
# The former group (a) — top 50 by L_ROC, the "load-paying" set — was defined
# from an AnaCoDa quantity and moved to supplementary_anacoda.R.  The two
# groups are close in size (top-50 by load vs the ~47-49 genes clearing the
# drift barrier), so the S_Wright group is the natural AnaCoDa-free
# replacement, and it is the metric this paper actually argues about:
# S_Wright measures efficacy of selection, L_ROC measures the load being paid.

custom_bag <- integrated_data |> dplyr::pull(Gene_name)


# (b) Selection group: S_Wright >= S_BARRIER ------------------------------
subset_selection <- msd_data |>
  dplyr::filter(!is.na(S_Wright_raw), S_Wright_raw >= S_BARRIER) |>
  dplyr::pull(Gene_name)

GO_results_selection <- gost(query = subset_selection,
                             organism = "gp__q7VP_EAck_dZk",
                             multi_query = FALSE, significant = TRUE,
                             correction_method = "fdr",
                             domain_scope = "custom", custom_bg = custom_bag,
                             user_threshold = 0.05)
write.csv(x = GO_results_selection$result |> dplyr::select(-parents),
          file = "./results/Go_enrichment_selection_S_Wright.csv",
          quote = TRUE, row.names = FALSE)
cat(sprintf("[GO] Selection group (S_Wright >= %.4f): n = %d genes\n",
            S_BARRIER, length(subset_selection)))

# Backwards-compatible aliases now point at the S_Wright selection group
# (they previously pointed at the L_ROC load group, which moved out).
subset_strongly_shaped_by_s <- subset_selection
GO_results <- GO_results_selection
write.csv(x = GO_results$result |> dplyr::select(-parents),
          file = "./results/Go_enrichment.csv",
          quote = TRUE, row.names = FALSE)

# 8.4b) GO enrichment dot-plot visualisation ----
#
# Shows ALL significant terms after removing overly generic ones
# (term_size > 500). Plot height scales automatically with term count.
# Gene ratio (precision) drives the x-axis; dot size = overlap count;
# colour = -log10(FDR p-value).

.go_dotplot <- function(go_result, title, max_term_size = 500) {
  empty_plot <- function(sub) {
    ggplot() +
      labs(title = title, subtitle = sub) +
      theme_void()
  }

  if (is.null(go_result) || nrow(go_result$result) == 0)
    return(list(plot = empty_plot("No significant GO terms"), nterms = 0L))

  df <- go_result$result |>
    dplyr::filter(term_size <= max_term_size) |>
    dplyr::arrange(p_value) |>
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      gene_ratio  = precision,
      term_label  = sapply(
        term_name,
        function(x) paste(strwrap(x, width = 45), collapse = "\n")
      )
    )

  if (nrow(df) == 0)
    return(list(plot = empty_plot("No specific terms (all filtered as generic)"), nterms = 0L))

  p <- ggplot(df, aes(x = gene_ratio,
                      y = reorder(term_label, gene_ratio),
                      size = intersection_size,
                      colour = neg_log10_p)) +
    geom_point() +
    scale_colour_viridis_c(
      name   = expression(-log[10](italic(p))),
      option = "plasma",
      begin  = 0.2, end = 0.95
    ) +
    scale_size_continuous(name = "Genes", range = c(2, 9)) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.05, 0.2))
    ) +
    labs(title = title, x = "Gene ratio", y = NULL) +
    theme_custom() +
    theme(
      axis.text.y       = element_text(size = 8, lineheight = 1.15),
      plot.title        = element_text(size = 9, face = "bold"),
      legend.key.height = unit(0.5, "cm")
    )

  list(plot = p, nterms = nrow(df))
}

.go_height <- function(n) max(5, 2 + n * 0.32)

res_sel <- .go_dotplot(
  GO_results_selection,
  sprintf("Population-genetic selection  (S_Wright >= %.1f,  n = %d)",
          S_BARRIER, length(subset_selection))
)

ggsave("./results/GO_dotplot_selection.pdf",
       res_sel$plot,  width = 8, height = .go_height(res_sel$nterms))

cat(sprintf("[GO plots] Selection: %d terms shown\n", res_sel$nterms))

rm(res_sel, .go_dotplot, .go_height)
gc()

# 8.5) Top genes by L_ROC (load) and by S_Wright (selection) -----------------------------------

detailed_annotation_full <- read.delim(
  "data/Mguttatusvar_IM767_887_v2.1.annotation_info.txt",
  header = TRUE, sep = "\t", comment.char = "", quote = "",
  fill = TRUE, na.strings = ""
) |>
  dplyr::select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline) |>
  dplyr::distinct()

# The former "top 50 by L_ROC (load-paying)" table moved to
# supplementary_anacoda.R along with L_ROC itself.

# Selection group: S_Wright >= S_BARRIER
top_selection <- msd_data |>
  dplyr::filter(!is.na(S_Wright_raw), S_Wright_raw >= S_BARRIER) |>
  dplyr::arrange(desc(S_Wright_raw)) |>
  dplyr::select(Gene_name, S_Wright_raw, Mean_Log10_Exp) |>
  dplyr::left_join(detailed_annotation_full,
                   by = c("Gene_name" = "locusName"))
write.csv(top_selection,
          "./results/Top_genes_strong_selection_S_Wright.csv",
          quote = TRUE, row.names = FALSE)
cat(sprintf("[Top genes] S_Wright_raw >= %.4f: %d genes (Q-inflection-derived selection group)\n",
            S_BARRIER, nrow(top_selection)))
## ============================================================================
## RESULTS 6 — M. guttatus preferred codons vs other plants
##   Produces:
##     Figure 3  Cross-species preferred codon comparison
##               (`plant_codon_preference_comparison_colored.pdf`,
##                `plant_preferred_codons_comparison.csv`)
## ============================================================================

## 9) Comparing preferred codon of Mimulus guttatus to other plants ----

# Use w_table from CAI analysis (already calculated preferred codons)
cat("Using optimal codons from corrected reference set...\n")

# Get preferred codons (those with relative_adaptiveness == 1.0)
preferred_codons_comparative <- preferred_codons |>
  dplyr::mutate(Codon_RNA = gsub("T", "U", Codon)) |>
  dplyr::select(Amino_Acid = AA, aa, Codon_RNA, eta)

# Collapse amino acids with six codons back into six, based on relative adaptiveness

preferred_codons_comparative <- preferred_codons_comparative |>
  dplyr::mutate(AA_root = sapply(preferred_codons_comparative$Amino_Acid, 
                                 function(x) 
                                 {
                                   unlist(strsplit(x, "_"))[1]
                                 }))

# `eta` is carried through from Section 8.2 (negated regime contrast, so that
# lower = more preferred, matching the sign convention merge_2_and_4_to_6_fold()
# expects).  It no longer comes from the AnaCoDa posterior.

preferred_codons_mg <- merge_2_and_4_to_6_fold(
  preferred_codons_comparative,
  "AA_root"
)

# Add M. guttatus to the global plant comparison table
mg_prefs <- preferred_codons_mg |>
  dplyr::select(Amino_Acid, Mimulus_guttatus = Codon_RNA)

plant_codons_extended <- model_plants_PC |>
  left_join(mg_prefs, by = "Amino_Acid") |>
  na.omit()

# Reorder columns
plant_codons_extended <- plant_codons_extended |>
  dplyr::select(Group, Amino_Acid, Arabidopsis_thaliana, Populus_trichocarpa, 
                Mimulus_guttatus, Physcomitrella_patens, Synonymous_Codons)

# Save extended table
write.csv(plant_codons_extended, "./results/plant_preferred_codons_comparison.csv", 
          row.names = FALSE, quote = FALSE)
# Create a data frame for the plot
plot_data <- data.frame()

# Species order: Arabidopsis, Populus, Physcomitrella, then Mimulus
species_order <- c("Arabidopsis_thaliana", "Populus_trichocarpa", 
                   "Physcomitrella_patens", "Mimulus_guttatus")
species_labels <- c("A. thaliana", "P. trichocarpa", "P. patens", "M. guttatus")

for (aa in sort(unique(plant_codons_extended$Amino_Acid))) {
  # Determine chemistry group
  aa_group <- "Other"
  for (grp in names(aa_chemistry)) {
    if (aa %in% aa_chemistry[[grp]]) {
      aa_group <- gsub("_", " ", grp)  # Convert underscores to spaces here
      break
    }
  }
  
  aa_data <- plant_codons_extended |> dplyr::filter(Amino_Acid == aa)
  
  # Get preferred codons for each species
  codons_list <- list()
  for (sp in species_order) {
    if (sp %in% colnames(plant_codons_extended)) {
      codon_str <- aa_data[[sp]][1]
      if (!is.na(codon_str) && codon_str != "") {
        codons_list[[sp]] <- unique(unlist(strsplit(codon_str, "/")))
      } else {
        codons_list[[sp]] <- character(0)
      }
    }
  }
  
  # For each species, add their preferred codons
  for (i in 1:length(species_order)) {
    sp <- species_order[i]
    sp_label <- species_labels[i]
    
    if (length(codons_list[[sp]]) > 0) {
      codon_text <- paste(codons_list[[sp]], collapse = "/")
      
      # Determine color for M. guttatus column
      if (sp == "Mimulus_guttatus") {
        # Check which species M. guttatus shares with
        mg_codons <- codons_list[["Mimulus_guttatus"]]
        at_codons <- codons_list[["Arabidopsis_thaliana"]]
        pt_codons <- codons_list[["Populus_trichocarpa"]]
        pp_codons <- codons_list[["Physcomitrella_patens"]]
        
        shares_with <- c()
        if (length(intersect(mg_codons, at_codons)) > 0) shares_with <- c(shares_with, "Arabidopsis")
        if (length(intersect(mg_codons, pt_codons)) > 0) shares_with <- c(shares_with, "Populus")
        if (length(intersect(mg_codons, pp_codons)) > 0) shares_with <- c(shares_with, "Physcomitrella")
        
        # Assign color based on sharing pattern
        if (length(shares_with) == 0) {
          codon_color <- "Unique"
        } else if (length(shares_with) == 3) {
          codon_color <- "All_three"
        } else if (length(shares_with) == 2) {
          codon_color <- "Two_species"
        } else {
          # Shares with only one species
          if ("Arabidopsis" %in% shares_with) {
            codon_color <- "Only_Arabidopsis"
          } else if ("Populus" %in% shares_with) {
            codon_color <- "Only_Populus"
          } else {
            codon_color <- "Only_Physcomitrella"
          }
        }
      } else {
        # For other species, use their own color
        codon_color <- sp_label
      }
      
      plot_data <- rbind(plot_data,
                         data.frame(
                           Amino_Acid = aa,
                           Chemistry = aa_group,  # Already converted above
                           Species = sp_label,
                           Codon = codon_text,
                           Color_Category = codon_color,
                           stringsAsFactors = FALSE
                         ))
    }
  }
}

# Set factor levels for proper ordering
plot_data$Species <- factor(plot_data$Species, levels = species_labels)
plot_data$Chemistry <- factor(plot_data$Chemistry, 
                              levels = c("Nonpolar Aliphatic", "Aromatic", 
                                         "Polar Uncharged", "Positively Charged", 
                                         "Negatively Charged", "Other"))

# Define colors
color_palette <- c(
  "A. thaliana" = "#E41A1C",           # Red for Arabidopsis
  "P. trichocarpa" = "#377EB8",        # Blue for Populus
  "P. patens" = "#4DAF4A",             # Green for Physcomitrella
  "Only_Arabidopsis" = "#E41A1C",      # Red - shares only with Arabidopsis
  "Only_Populus" = "#377EB8",          # Blue - shares only with Populus
  "Only_Physcomitrella" = "#4DAF4A",   # Green - shares only with Physcomitrella
  "Two_species" = "#FF7F00",           # Orange - shares with two species
  "All_three" = "#984EA3",             # Purple - shares with all three
  "Unique" = "#999999"                 # Gray - unique to M. guttatus
)

# Create the plot
p_comparison <- ggplot(plot_data, aes(x = Species, y = Amino_Acid, label = Codon)) +
  geom_tile(aes(fill = Color_Category), color = "white", size = 1, alpha = 0.3) +
  geom_text(size = 3, fontface = "bold") +
  scale_fill_manual(values = color_palette,
                    labels = c("A. thaliana" = "A. thaliana",
                               "P. trichocarpa" = "P. trichocarpa",
                               "P. patens" = "P. patens",
                               "Only_Arabidopsis" = "M.g. shares with Arabidopsis only",
                               "Only_Populus" = "M.g. shares with Populus only",
                               "Only_Physcomitrella" = "M.g. shares with Physcomitrella only",
                               "Two_species" = "M.g. shares with two species",
                               "All_three" = "M.g. shares with all three",
                               "Unique" = "M.g. unique preference"),
                    name = "") +
  facet_grid(Chemistry ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Preferred Codon Usage Across Plant Species",
       subtitle = "M. guttatus (rightmost column) colored by sharing pattern with other species",
       x = "", y = "") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
        axis.text.y = element_text(size = 10),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 11),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        panel.grid = element_blank())

ggsave("./results/plant_codon_preference_comparison_colored.pdf", p_comparison, 
       width = 12, height = 16)
## ============================================================================
## RESULTS 7 — tRNA / codon-anticodon correspondence
##   Produces:
##     Figure S5 supporting data (aa-tRNA correspondence)
##     Cited value:  Spearman r = 0.761 (p = 3.83 × 10⁻⁴) between
##                   amino acid frequencies and tRNA gene copy number
## ============================================================================

## 11) tRNA abundance correlation analysis ----
## _____________________________________________________________________________

# NOTE: The genome-wide analysis (Analysis 1) correlates raw codon frequencies
# with tRNA supply. Because this genome has AT-rich mutational bias, the most
# frequent codons genome-wide are AT-ending (due to mutation, not selection).
# The within-family analysis corrects for this by examining proportions within
# each amino acid family. The top-expression tier analysis (Analysis 2) further
# isolates the selection signal by focusing on genes under strongest selection.

# Analysis 1: Genome-wide (baseline, with proper wobble rules)

cat("\n=== Analysis 1: Genome-wide tRNA-Codon Correlation (Baseline) ===\n")

tRNA_copynumber_results <- tRNA_codon_correlation(
  codon_counts = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_copynumber",
  test_method = "spearman",
  mode = "by.copy.number",
  wobble_mode = "conservative"  # eukaryotic rules: A34→I34 modification
)

aa_trna_check <- check_aa_frequency_vs_tRNA_supply(
  codon_usage = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/aa_trna_sanity_check"
)
## ============================================================================
## RESULTS 8 — Selection on codon usage elevates synonymous diversity
##   The McVean & Charlesworth expectation: where selection on codon usage is
##   effective it opposes mutation pressure at 4-fold sites, holding variants at
##   intermediate frequency and RAISING synonymous diversity rather than
##   depressing it.  This result never depended on position within genes.
##   Produces:
##     Figure 7A  Synonymous & non-synonymous π vs expression rank
##                (`pi_4fold_by_expression_rank.pdf`)
##     Figure 7B  Preferred codon frequency vs expression group
##                (`Frequency_preferred_by_expression_group_Median_CI.pdf`)
##     Figure 6B  Preferred-codon-frequency landscape vs expression × gene length
##                (`Preferred_freq_contour_exp_x_length_broad.pdf`,
##                 `Preferred_freq_contour_exp_x_length_narrow.pdf`)
##
##   MOVED OUT: the former Figure 7C (intron vs exon π by distance from gene
##   start) and the first-300 bp decomposition are within-gene positional
##   analyses and now live in paper2_linked_selection.R.
## ============================================================================

## 12) Polymorphism data integration ----
## _____________________________________________________________________________
# By-gene polymorphism (Pi_mean_4fold etc.) was preloaded in Section 5.5 because
# Section 6 GAMs and Section 8.3.4 msd_data depend on those columns.

# Memory cleanup: polymorphism raw data (now joined into integrated_data) ---
rm(pi_data)


# 12.1) Expression-ranked 4-fold π analysis (Kelly replication) ----
# Bin genes into groups of ~1000 ranked by Mean_Log10_Exp, calculate
# weighted mean 4-fold nucleotide diversity within each bin.

cat(sprintf("\n=== Section 12.1: Expression-ranked 4-fold pi | integrated_data N = %d (with Mean_Log10_Exp: %d, with Pi_sum_4fold: %d) ===\n",
            nrow(integrated_data),
            sum(!is.na(integrated_data$Mean_Log10_Exp)),
            sum(!is.na(integrated_data$Pi_sum_4fold))))

bin_size <- 1000

# Check if mutation-type columns are available (requires extended calculate_pi.py)
mutation_types <- c("AC", "AG", "AT", "CG", "CT", "GT")
has_mutation_types <- all(paste0("Pi_sum_4fold_", mutation_types) %in% 
                            names(integrated_data))

# Rank genes by Mean_Log10_Exp and create bins
pi_by_expression <- integrated_data |>
  dplyr::filter(!(Gene_name %in% selection_gene_set)) |>
  dplyr::arrange(Mean_Log10_Exp) |>
  dplyr::mutate(
    Rank = dplyr::row_number(),
    Exp_Bin = ceiling(Rank / bin_size)
  ) |>
  dplyr::group_by(Exp_Bin) |>
  dplyr::summarize(
    n_genes = n(),
    mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
    # Weighted mean π at 4-fold sites: total π_sum / total sites
    total_pi_sum_4fold = sum(Pi_sum_4fold, na.rm = TRUE),
    total_sites_4fold = sum(Sites_4fold, na.rm = TRUE),
    weighted_pi_4fold = total_pi_sum_4fold / total_sites_4fold,
    # Also compute individual-gene SD for error bars
    sd_pi_4fold = sd(Pi_mean_4fold, na.rm = TRUE),
    se_pi_4fold = sd_pi_4fold / sqrt(n()),
    .groups = "drop"
  )

# Selection group lives in the bin immediately after the last expression bin,
# so the index stays correct under any upstream filter changes.
sel_bin_id <- max(pi_by_expression$Exp_Bin) + 1L

sel_cat <- integrated_data |>
  dplyr::filter(Gene_name %in% selection_gene_set) |>
  dplyr::summarize(
    Exp_Bin = sel_bin_id,
    n_genes = n(),
    mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
    # Weighted mean π at 4-fold sites: total π_sum / total sites
    total_pi_sum_4fold = sum(Pi_sum_4fold, na.rm = TRUE),
    total_sites_4fold = sum(Sites_4fold, na.rm = TRUE),
    weighted_pi_4fold = total_pi_sum_4fold / total_sites_4fold,
    # Also compute individual-gene SD for error bars
    sd_pi_4fold = sd(Pi_mean_4fold, na.rm = TRUE),
    se_pi_4fold = sd_pi_4fold / sqrt(n()),
    .groups = "drop"
  )

pi_by_expression <- pi_by_expression |>
  rbind(sel_cat) # Final bin holds the S_Wright selection group

cat("\n=== 4-fold π by Expression Rank (groups of ~1000 genes) ===\n")
print(pi_by_expression)

# Plot: replicating advisor's graph (interval plot with individual SDs)
p_pi_by_expression <- ggplot(pi_by_expression, 
                             aes(x = Exp_Bin, y = weighted_pi_4fold)) +
  geom_point(size = 3, color = "#377EB8") +
  geom_errorbar(aes(ymin = weighted_pi_4fold - se_pi_4fold,
                    ymax = weighted_pi_4fold + se_pi_4fold),
                width = 0.3, color = "#377EB8") +
  labs(
    title = expression(paste("4-fold Nucleotide Diversity (", pi, 
                             ") by Expression Level")),
    subtitle = "Genes ranked by Mean Log10 Expression, binned in groups of ~1000",
    x = "Expression level category",
    y = expression(paste("nuc_diversity (4 fold)"))
  ) +
  scale_x_continuous(breaks = seq_len(max(pi_by_expression$Exp_Bin))) +
  theme_custom()

ggsave("./results/pi_4fold_by_expression_rank.pdf", 
       p_pi_by_expression, width = 10, height = 6)

cat("✓ Saved: ./results/pi_4fold_by_expression_rank.pdf\n")

# Breakdown by segregating base pair type at 4-fold sites ----
if (has_mutation_types) {
  
  cat("\n=== 4-fold π by Mutation Type and Expression Rank ===\n")
  
  # Calculate per-mutation-type pi component within each expression bin
  # Component = sum(Pi_sum_type) / sum(Sites_4fold) → additive decomposition
  pi_by_mutation <- integrated_data |>
    dplyr::filter(!(Gene_name %in% selection_gene_set)) |>
    dplyr::arrange(Mean_Log10_Exp) |>
    dplyr::mutate(
      Rank = dplyr::row_number(),
      Exp_Bin = ceiling(Rank / bin_size)
    ) |>
    dplyr::group_by(Exp_Bin) |>
    dplyr::summarize(
      n_genes = n(),
      mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
      total_sites_4fold = sum(Sites_4fold, na.rm = TRUE),
      pi_AC = sum(Pi_sum_4fold_AC, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AG = sum(Pi_sum_4fold_AG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AT = sum(Pi_sum_4fold_AT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CG = sum(Pi_sum_4fold_CG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CT = sum(Pi_sum_4fold_CT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_GT = sum(Pi_sum_4fold_GT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      .groups = "drop"
    )
  
  pi_sel_group <- integrated_data |>
    dplyr::filter(Gene_name %in% selection_gene_set) |>
    dplyr::summarize(
      Exp_Bin = sel_bin_id,
      n_genes = n(),
      mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
      total_sites_4fold = sum(Sites_4fold, na.rm = TRUE),
      pi_AC = sum(Pi_sum_4fold_AC, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AG = sum(Pi_sum_4fold_AG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AT = sum(Pi_sum_4fold_AT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CG = sum(Pi_sum_4fold_CG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CT = sum(Pi_sum_4fold_CT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_GT = sum(Pi_sum_4fold_GT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      .groups = "drop"
    )
  
  pi_by_mutation <- pi_by_mutation |>
    rbind(pi_sel_group)
  
  # Verify additive decomposition
  pi_check <- pi_by_expression |>
    dplyr::left_join(
      pi_by_mutation |> 
        dplyr::mutate(sum_components = pi_AC + pi_AG + pi_AT + pi_CG + pi_CT + pi_GT) |>
        dplyr::select(Exp_Bin, sum_components), 
      by = "Exp_Bin"
    )
  cat(sprintf("Additive check (max |total - sum_components|): %.2e\n",
              max(abs(pi_check$weighted_pi_4fold - pi_check$sum_components), na.rm = TRUE)))
  
  # Pivot to long format for plotting
  pi_mutation_long <- pi_by_mutation |>
    tidyr::pivot_longer(
      cols = starts_with("pi_"),
      names_to = "Mutation_Type",
      values_to = "Pi_component",
      names_prefix = "pi_"
    )
  
  # Plot 1: All mutation types overlaid
  p_pi_by_mutation <- ggplot(pi_mutation_long, 
                             aes(x = Exp_Bin, y = Pi_component, 
                                 color = Mutation_Type)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    scale_color_brewer(palette = "Set2", name = "Segregating\nBases") +
    labs(
      title = expression(paste("4-fold ", pi, " by Segregating Base Pair and Expression")),
      subtitle = "Additive components of 4-fold diversity by mutation type",
      x = "Expression level category",
      y = expression(paste(pi, " component (4-fold)"))
    ) +
    scale_x_continuous(breaks = seq_len(max(pi_by_mutation$Exp_Bin))) +
    theme_custom() +
    theme(legend.position = "right")
  
  ggsave("./results/pi_4fold_by_mutation_type_and_expression.pdf", 
         p_pi_by_mutation, width = 12, height = 6)
  
  cat("✓ Saved: ./results/pi_4fold_by_mutation_type_and_expression.pdf\n")
  
} else {
  cat("\nNote: Mutation-type columns not found in pi data.\n")
  cat("Re-run calculate_pi.py (extended version) to generate per-mutation-type output.\n")
  cat("Required columns: Pi_sum_4fold_AC, Pi_sum_4fold_AG, ..., Pi_sum_4fold_GT\n")
}

# Memory cleanup
rm(p_pi_by_expression, bin_size, mutation_types)

# 12.2) Tracking frequency of preferred allele as a function of expression ----

cat(sprintf("\n=== Section 12.2: Preferred-allele frequency vs expression | integrated_data N = %d ===\n",
            nrow(integrated_data)))

preferred_data <- read.delim("./data/all_chromosomes.codon_frequencies_preferred.txt",
                             stringsAsFactors = FALSE) |>
  dplyr::mutate(Gene = paste0("MgIM767.", Gene))

preferred_data <- preferred_data |>
  dplyr::select(Gene, Preferred_Freq) |>
  dplyr::rename(Gene_name = Gene) |>
  dplyr::group_by(Gene_name) |>
  summarize(
    Mean_preferred_freq = mean(Preferred_Freq)
  ) |>
  ungroup()

integrated_data <- integrated_data |>
  left_join(preferred_data) |>
  na.exclude()

# Memory cleanup: preferred frequency data (now joined into integrated_data) ---
rm(preferred_data)

# GAM wrappers

predictors_p <- c("Max_Log10_Exp", "Exp_breadth", "Total_Codons")

preferred_nonlinearity_results <- analyze_nonlinearity_suite(
  resp = "Mean_preferred_freq",
  predictors = predictors_p,
  data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0),
  family = betar(link = "logit")
)

preferred_models <- fit_codon_gam_suite(
  data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0),
  model_list = alist(
    Null        = ~ 1,
    Additive    = ~ s(Max_Log10_Exp) + s(Exp_breadth) + s(Total_Codons),
    Interaction = ~ te(Max_Log10_Exp, Exp_breadth) + s(Total_Codons),
    Complex     = ~ te(Max_Log10_Exp, Exp_breadth, Total_Codons)
  ),
  response_var = "Mean_preferred_freq",
  family = betar(link = "logit") 
)

preferred_selection <- get_model_selection_stats(preferred_models)
preferred_selection_winner <- preferred_selection |> 
  dplyr::filter(AIC == min(AIC))

run_posteriori_gam_analysis(model = preferred_models[["Complex"]], 
                            data = integrated_data |> 
                              dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0),
                            focal_pred = "Max_Log10_Exp", 
                            interact_pred = "Exp_breadth",
                            third_pred = NULL,
                            response_name = "Mean_preferred_freq",
                            prefix = "Preferred")

summary(preferred_models[["Complex"]])

# Assessing significance of expression over the detrended residuals

confounder_model_gam <- gam(Mean_preferred_freq ~ s(CDS_length_nt),
                            data = integrated_data,
                            family = betar(link = "logit"))

integrated_data$Mean_preferred_freq_detrended <- residuals(confounder_model_gam, 
                                           type = "response")

cat("\n=== Kruskal-Wallis Test: Frequency of preferred codons across Groups ===\n")

kw_preferred_freq <- kruskal.test(Mean_preferred_freq_detrended ~ Expression_Group, 
                                  data = integrated_data)

# Plotting and assessing significance using Dunn

print(kw_preferred_freq)
if (kw_preferred_freq$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Perform Dunn's test with FDR correction
  dunn_result_detrended <- dunn.test::dunn.test(
    x = integrated_data$Mean_preferred_freq_detrended,
    g = integrated_data$Expression_Group,
    method = "bh",
    kw = TRUE,
    label = TRUE,
    wrap = FALSE,
    table = TRUE,
    list = FALSE,
    altp = TRUE
  )
} else {
  cat("\nNo significant difference among groups (p >= 0.05)\n")
  cat("Post-hoc tests not necessary.\n")
}

# Extract groups for effect size calculations
top5_preferred <- integrated_data |>
  dplyr::filter(Expression_Group == "Top 5%") |>
  dplyr::pull(Mean_preferred_freq_detrended)

middle_preferred <- integrated_data |>
  dplyr::filter(Expression_Group == "Middle 90%") |>
  dplyr::pull(Mean_preferred_freq_detrended)

bottom5_preferred <- integrated_data |>
  dplyr::filter(Expression_Group == "Bottom 5%") |>
  dplyr::pull(Mean_preferred_freq_detrended)

# Calculate effect sizes
if (length(top5_preferred) > 0 && length(middle_preferred) > 0) {
  d_top_middle_preferred <- cohens_d_calc(top5_preferred, middle_preferred)
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle_preferred))
}

if (length(top5_preferred) > 0 && length(bottom5_preferred) > 0) {
  d_top_bottom_preferred <- cohens_d_calc(top5_preferred, bottom5_preferred)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom_preferred))
}

if (length(middle_preferred) > 0 && length(bottom5_preferred) > 0) {
  d_middle_bottom_preferred <- cohens_d_calc(middle_preferred, bottom5_preferred)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom_preferred))
}

# Median and CI

plot_data_pref <- integrated_data |>
  dplyr::mutate(Exp_Group = factor(Expression_Group, 
                                   levels = c("Bottom 5%", "Middle 90%", "Top 5%"))) |>
  dplyr::filter(!is.na(Exp_Group))

p_preferred_median <- ggplot(plot_data_pref, aes(x = Exp_Group, y = Mean_preferred_freq_detrended)) +
  
  # A. Median and 95% Bootstrap CI
  stat_summary(fun.data = median_cl_boot, 
               geom = "errorbar", width = 0.15, size = 0.8, color = "black") +
  stat_summary(fun = median, geom = "point", size = 4, aes(color = Exp_Group)) +
  
  # B. Formatting
  scale_color_manual(values = c("Bottom 5%" = "#377EB8", 
                                "Middle 90%" = "#999999", 
                                "Top 5%" = "#E41A1C")) +
  
  labs(y = "Median Frequency of Preferred Codons (length corrected)",
       x = NULL) + # Remove X label as groups are self-explanatory
  
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        panel.grid.major.x = element_blank())

# Save
ggsave("./results/Frequency_preferred_by_expression_group_Median_CI.pdf", 
       p_preferred_median, width = 5, height = 6)

# Execute for 16% model
p_surface_pref <- plot_selection_surface(
  model = preferred_models[["Complex"]], 
  data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0),
  response_name = "Mean_preferred_freq"
)

# 12.3) Contour plot: preferred codon frequency ~ expression x gene length ----
# Shows joint effect of expression and gene length on the frequency of
# ROC-preferred codons across ALL amino acid families (not only 4-fold).
# Mean_preferred_freq is the per-gene average across all codon positions.
# Also produces separate contour surfaces for C-ending and G-ending
# preferred codons, which reduces noise from mixing two distinct nucleotide
# biases that selection must act against.

cat(sprintf("\n=== Section 12.3: Contour, preferred codon freq ~ Expression x Length | integrated_data N = %d (with Mean_preferred_freq: %d) ===\n",
            nrow(integrated_data),
            sum(!is.na(integrated_data$Mean_preferred_freq))))

# --- 12.4a: Overall (all preferred codons pooled) ---

contour_data <- integrated_data |>
  dplyr::filter(!is.na(Mean_preferred_freq),
                !is.na(Max_Log10_Exp),
                Total_Codons > 0) |>
  dplyr::mutate(log10_length = log10(Total_Codons))
cat(sprintf("contour_data: integrated_data %d -> %d genes after NA filter on (Mean_preferred_freq, Max_Log10_Exp, Total_Codons > 0)\n",
            nrow(integrated_data), nrow(contour_data)))

contour_gam <- mgcv::gam(
  Mean_preferred_freq ~ te(Max_Log10_Exp, Exp_breadth, log10_length, k = c(10, 10)),
  data = contour_data,
  family = betar(link = "logit")
)
cat("GAM surface R-sq(adj) [all preferred]: ", summary(contour_gam)$r.sq, "\n")

pred_grid_pref_broad <- expand.grid(
  Max_Log10_Exp = seq(min(contour_data$Max_Log10_Exp, na.rm = TRUE),
                      max(contour_data$Max_Log10_Exp, na.rm = TRUE),
                      length.out = 200),
  log10_length  = seq(min(contour_data$log10_length, na.rm = TRUE),
                      max(contour_data$log10_length, na.rm = TRUE),
                      length.out = 200),
  Exp_breadth = 33) # Holding breadth constant (broadly expressed genes)

pred_grid_pref_narrow <- expand.grid(
  Max_Log10_Exp = seq(min(contour_data$Max_Log10_Exp, na.rm = TRUE),
                      max(contour_data$Max_Log10_Exp, na.rm = TRUE),
                      length.out = 200),
  log10_length  = seq(min(contour_data$log10_length, na.rm = TRUE),
                      max(contour_data$log10_length, na.rm = TRUE),
                      length.out = 200),
  Exp_breadth = 1) # Holding breadth constant (narrowly expressed genes)

# Shared (x, y) grid for plotting; the two predictions only differ in Exp_breadth
pred_grid_pref <- pred_grid_pref_broad[, c("Max_Log10_Exp", "log10_length")]
pred_grid_pref$Predicted_broad  <- predict(contour_gam, newdata = pred_grid_pref_broad,
                                           type = "response")
pred_grid_pref$Predicted_narrow <- predict(contour_gam, newdata = pred_grid_pref_narrow,
                                           type = "response")

# Broad
p_pref_contour <- ggplot(pred_grid_pref,
                         aes(x = Max_Log10_Exp, y = log10_length)) +
  geom_raster(aes(fill = Predicted_broad), interpolate = TRUE) +
  geom_contour(aes(z = Predicted_broad), colour = "grey30",
               linewidth = 0.4, bins = 12) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Freq.\npreferred"
  ) +
  labs(
    title = "Frequency of ROC-Preferred Codons (All Amino Acids)",
    subtitle = "GAM-predicted surface across expression and gene length",
    x = expression(log[10](Max~Expression~CPM)),
    y = expression(log[10](Gene~Length~"(codons)"))
  ) +
  theme_custom() +
  theme(legend.position = "right",
        panel.grid = element_blank())

ggsave("./results/Preferred_freq_contour_exp_x_length_broad.pdf",
       p_pref_contour, width = 9, height = 7)
cat("Saved: ./results/Preferred_freq_contour_exp_x_length_broad.pdf\n")

# Narrow

p_pref_contour <- ggplot(pred_grid_pref,
                         aes(x = Max_Log10_Exp, y = log10_length)) +
  geom_raster(aes(fill = Predicted_narrow), interpolate = TRUE) +
  geom_contour(aes(z = Predicted_narrow), colour = "grey30",
               linewidth = 0.4, bins = 12) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Freq.\npreferred"
  ) +
  labs(
    title = "Frequency of ROC-Preferred Codons (All Amino Acids)",
    subtitle = "GAM-predicted surface across expression and gene length",
    x = expression(log[10](Max~Expression~CPM)),
    y = expression(log[10](Gene~Length~"(codons)"))
  ) +
  theme_custom() +
  theme(legend.position = "right",
        panel.grid = element_blank())

ggsave("./results/Preferred_freq_contour_exp_x_length_narrow.pdf",
       p_pref_contour, width = 9, height = 7)
cat("Saved: ./results/Preferred_freq_contour_exp_x_length_narrow.pdf\n")
## ============================================================================
## RESULTS 9 — Diversity hump: GC-segregating sites carry excess heterozygosity
##   Produces:
##     Figure 7D  Excess heterozygosity at GC-segregating 4-fold sites
##                vs matched intronic controls (`Hump_Hypothesis_Confirmation.pdf`)
##     Cited values:  Paired t-test t₁₃ = 39.76, P < 2.9 × 10⁻¹⁵
##                     Mean Δπ ≈ 0.0107 (≈1% extra heterozygosity at GC sites)
## ============================================================================

## *****************************************************************************
## 14) Diversity across different genomic compartment ----
## _____________________________________________________________________________

pi_compartment <- read.table(file = "./data/all_chromosomes.pi_by_compartment.txt",
                             header = T)

cat(sprintf("\n=== Section 14: Diversity across genomic compartments | pi_compartment N = %d rows (compartments x nucleotide categories; not gene-keyed) ===\n",
            nrow(pi_compartment)))

# HUMP EFFECT TEST ----
# 1. Aggregate Data by "Selection Potential"
# We group C, G, and CG as "GC_Segregating" (Where selection acts)
# We keep AT as "AT_Only" (Where selection is absent/invisible)
# IMPORTANT: exclude Nuc_Category == "all" to avoid double-counting
# (the "all" row is the total that already contains C+G+AT+CG)

hump_test_data <- pi_compartment |>
  dplyr::filter(Compartment %in% c("nonfirst_exon_4fold", "intron"),
                Nuc_Category != "all") |>
  dplyr::mutate(Site_Type = ifelse(Nuc_Category == "AT", "AT_Only", "GC_Segregating")) |>
  dplyr::group_by(Compartment, Site_Type) |>
  summarise(
    Total_Pi = sum(Pi_sum),
    Total_Sites = sum(Sites),
    Weighted_Pi = Total_Pi / Total_Sites,
    .groups = "drop"
  )

print(hump_test_data)

ggplot(hump_test_data, aes(x = Compartment, y = Weighted_Pi, fill = Site_Type)) +
  geom_col(position = "dodge") +
  labs(title = "Testing the McVean Hump Hypothesis",
       subtitle = "Does opposing selection boost diversity at GC sites?",
       y = "Weighted Nucleotide Diversity (Pi)",
       x = "Genomic Compartment") +
  scale_fill_manual(values = c("AT_Only" = "gray70", "GC_Segregating" = "firebrick")) +
  theme_custom()

ggsave("./results/diversity_hump_test_by_compartment.pdf", width = 7, height = 5)

# 2) Formal testing

# Build per-chromosome data for proper mean ± CI visualization
hump_per_chrom <- pi_compartment |>
  dplyr::filter(Compartment %in% c("nonfirst_exon_4fold", "intron"),
                Nuc_Category != "all") |>
  dplyr::mutate(Site_Type = ifelse(Nuc_Category == "AT", "AT_Only", "GC_Segregating")) |>
  dplyr::group_by(Chromosome, Compartment, Site_Type) |>
  dplyr::summarize(
    Pi_weighted = sum(Pi_sum) / sum(Sites),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Compartment_label = factor(
      ifelse(Compartment == "intron", "Intron", "Exon (4-fold)"),
      levels = c("Intron", "Exon (4-fold)")
    )
  )

# Summary: mean ± 95% CI across chromosomes
hump_summary <- hump_per_chrom |>
  dplyr::group_by(Compartment_label, Site_Type) |>
  dplyr::summarize(
    Mean_Pi = mean(Pi_weighted),
    SD_Pi = sd(Pi_weighted),
    n = dplyr::n(),
    SE_Pi = SD_Pi / sqrt(n),
    CI_lo = Mean_Pi - qt(0.975, n - 1) * SE_Pi,
    CI_hi = Mean_Pi + qt(0.975, n - 1) * SE_Pi,
    .groups = "drop"
  )

# Plot: Mean ± 95% CI with individual chromosome points and connecting lines
p_hump <- ggplot(hump_summary,
                 aes(x = Compartment_label, y = Mean_Pi,
                     color = Site_Type, group = Site_Type)) +
  # Connecting lines between compartments
  geom_line(linewidth = 1.2) +
  # 95% CI error bars
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi),
                width = 0.08, linewidth = 0.8) +
  # Mean points (larger, filled)
  geom_point(size = 4) +
  # Individual chromosome values (smaller, semi-transparent)
  geom_point(data = hump_per_chrom,
             aes(x = Compartment_label, y = Pi_weighted,
                 color = Site_Type, group = Site_Type),
             position = position_dodge(width = 0.15),
             size = 1.5, alpha = 0.4, shape = 16) +
  # Formatting
  scale_color_manual(
    values = c("AT_Only" = "#999999", "GC_Segregating" = "#C0392B"),
    labels = c("AT_Only" = "AT-only sites",
               "GC_Segregating" = "GC-segregating sites (C + G + CG)")
  ) +
  labs(
    title = "Evidence for Weak Selection (The 'Hump' Effect)",
    subtitle = "Mean ± 95% CI across 14 chromosomes | Points = individual chromosomes",
    y = expression("Per-site Nucleotide Diversity (" * pi * ")"),
    x = NULL,
    color = NULL
  ) +
  theme_custom() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 10))

ggsave("./results/Hump_Hypothesis_Confirmation.pdf", p_hump, width = 8, height = 6)

# Statistical test

# 1. Prepare Data with Aggregation step
# IMPORTANT: exclude Nuc_Category == "all" before grouping to avoid
# contaminating GC_Segregating with the total row (which already includes AT)
paired_test_data <- pi_compartment |>
  dplyr::filter(Compartment %in% c("nonfirst_exon_4fold", "intron"),
                Nuc_Category != "all") |>
  
  # Create the new categories (now "all" is excluded, so only C/G/CG -> GC_Segregating)
  dplyr::mutate(Site_Type = ifelse(Nuc_Category == "AT", "AT_Only", "GC_Segregating")) |>
  
  # Aggregate the multiple GC rows (C, G, CG) into one value per Chromosome
  dplyr::group_by(Chromosome, Compartment, Site_Type) |>
  dplyr::summarise(
    # Recalculate weighted mean: Sum of Pi / Sum of Sites
    Pi_mean = sum(Pi_sum) / sum(Sites), 
    .groups = "drop"
  ) |>
  
  # Now pivot (guaranteed to be unique now)
  tidyr::pivot_wider(names_from = Compartment, values_from = Pi_mean) |>
  
  # Calculate the boost
  dplyr::mutate(Diversity_Boost = nonfirst_exon_4fold - intron) |>
  stats::na.omit()

# 2. Run the Paired T-test
# Compare if the boost in GC sites is larger than the boost in AT sites
boost_comparison <- paired_test_data |>
  dplyr::select(Chromosome, Site_Type, Diversity_Boost) |>
  tidyr::pivot_wider(names_from = Site_Type, values_from = Diversity_Boost)

t_test_result <- stats::t.test(boost_comparison$GC_Segregating, 
                               boost_comparison$AT_Only, 
                               paired = TRUE, 
                               alternative = "greater")

print("=== Paired Test: Does Selection Boost GC Diversity More than AT? ===")
print(t_test_result)

# Memory cleanup: Section 14 plot objects and raw data ---
# Keeping: t_test_result, boost_comparison, hump_test_data
rm(pi_compartment, hump_per_chrom, hump_summary, p_hump,
   paired_test_data)
gc()
## 16) GAM models for codon-based analysis ----
## _____________________________________________________________________________


# 16.1: Load per-codon data and filter to 4-fold degenerate sites ----

cat("Loading codon frequency data...\n")
codon_data_raw <- fread("data/all_chromosomes.codon_frequencies_preferred.txt",
                        showProgress = FALSE)

# Identify 4-fold degenerate amino acid families by AA code + codon prefix:
#   Simple 4-fold: Ala (A), Gly (G), Pro (P), Thr (T), Val (V)
#   Split-family 4-fold: Leu_4 (L,CT*), Ser_4 (S,TC*), Arg_4 (R,CG*)
fourfold_simple <- c("A", "G", "P", "T", "V")
codon_data_raw[, codon_prefix := substr(Ref_Codon, 1, 2)]

codon_4fold <- codon_data_raw[
  (AA %in% fourfold_simple) |
  (AA == "L" & codon_prefix == "CT") |
  (AA == "S" & codon_prefix == "TC") |
  (AA == "R" & codon_prefix == "CG")
]
rm(codon_data_raw)
gc()

# 16.2: Determine polymorphism status at each 4-fold site ----
# A site is polymorphic if >1 codon (differing at 3rd position) has freq > 0.
# Strategy: remove zero-freq entries from the Frequencies string;
# if a semicolon remains, >= 2 codons are segregating -> polymorphic.

cat("Classifying poly / mono at each site...\n")

freq_clean <- codon_4fold$Frequencies
# Remove zero-frequency entries: ";XXX:0.000" (middle or end positions)
freq_clean <- gsub(";[ACGT]{3}:0\\.000", "", freq_clean)
# Remove zero-frequency entries: "XXX:0.000;" (beginning position)
freq_clean <- gsub("[ACGT]{3}:0\\.000;", "", freq_clean)
# Remove zero-frequency entries: "XXX:0.000" (standalone / last remaining)
freq_clean <- gsub("[ACGT]{3}:0\\.000",  "", freq_clean)

codon_4fold[, is_poly := as.integer(grepl(";", freq_clean, fixed = TRUE))]

cat(sprintf("  Polymorphic: %s / %s (%.2f%%)\n",
            format(sum(codon_4fold$is_poly), big.mark = ","),
            format(nrow(codon_4fold), big.mark = ","),
            100 * mean(codon_4fold$is_poly)))

# 16.3: Merge expression data and compute normalized predictors ----

# Amino acid full names for reporting
aa_name_map <- c(A = "Ala", G = "Gly", L = "Leu", P = "Pro",
                 R = "Arg", S = "Ser", T = "Thr", V = "Val")
codon_4fold[, AA_name := aa_name_map[AA]]

# Gene name matching — codon data uses e.g. "01G000100",
# integrated_data uses "MgIM767.01G000100"
codon_4fold[, Gene_clean := paste0("MgIM767.", Gene)]

gene_info_16 <- integrated_data[, c("Gene_name", 
                                    "Mean_Log10_Exp", "Total_Codons")]

codon_4fold <- merge(codon_4fold, gene_info_16,
                     by.x = "Gene_clean", by.y = "Gene_name",
                     all.x = FALSE)

# Within-gene position (dist_norm / dist_z) was constructed here but never
# entered any model downstream; it moved to paper2_linked_selection.R with the
# rest of the positional analyses.

codon_4fold[, exp_norm := Mean_Log10_Exp]

# Store expression range for prediction grids (contour plots)
exp_range_16 <- range(codon_4fold$exp_norm, na.rm = TRUE)

# Standardised predictors (z-scored) for comparable coefficient magnitudes
# Raw predictors are retained for contour-plot prediction grids
codon_4fold[, exp_z  := as.numeric(scale(exp_norm))]

cat(sprintf("  After expression merge: %s sites in %s genes\n",
            format(nrow(codon_4fold), big.mark = ","),
            format(length(unique(codon_4fold$Gene_clean)), big.mark = ",")))
cat(sprintf("  exp_norm: mean=%.3f, sd=%.3f\n",
            mean(codon_4fold$exp_norm),  sd(codon_4fold$exp_norm)))

# 3rd-position base of the preferred codon at each site
# (downstream gene-level GAMs in 16.11 split sites by C-ending vs G-ending).
codon_4fold[, pref_base3 := substr(Preferred_Codon, 3, 3)]

# 16.4: Compute per-site π at the 3rd codon position ----
## For 4-fold degenerate sites, only the 3rd position is degenerate.
## π = n/(n-1) * (1 - Σ p_i²)  where p_i are nucleotide frequencies
## at the 3rd position, n = total number of alleles sampled.
## We extract the 3rd nucleotide from each codon in the Frequencies string,
## aggregate frequencies by nucleotide, then compute heterozygosity.

cat("\n\n=== Computing per-site π at 3rd codon position ===\n")

# Parse the Frequencies column into long format for efficient computation
codon_4fold[, site_id := .I]

# Fast vectorised parsing: split Frequencies by ";", then by ":"
freq_entries <- strsplit(codon_4fold$Frequencies, ";", fixed = TRUE)
freq_long_16 <- data.table(
  site_id = rep(codon_4fold$site_id, lengths(freq_entries)),
  entry   = unlist(freq_entries)
)
freq_long_16[, codon := sub(":.*", "", entry)]
freq_long_16[, freq  := as.numeric(sub(".*:", "", entry))]
freq_long_16[, nuc3  := substr(codon, 3, 3)]

# Aggregate frequency by 3rd-position nucleotide within each site
nuc_agg <- freq_long_16[, .(nuc_freq = sum(freq)), by = .(site_id, nuc3)]

# Determine sample size from Codon_Variants (sum of counts)
variant_entries <- strsplit(codon_4fold$Codon_Variants, ";", fixed = TRUE)
var_long_16 <- data.table(
  site_id = rep(codon_4fold$site_id, lengths(variant_entries)),
  entry   = unlist(variant_entries)
)
var_long_16[, count := as.integer(sub(".*:", "", entry))]
site_n <- var_long_16[, .(n_alleles = sum(count)), by = site_id]

# Per-site π = n/(n-1) * (1 - Σ p_i²)
site_pi_16 <- nuc_agg[, .(sum_p2 = sum(nuc_freq^2)), by = site_id]
site_pi_16 <- merge(site_pi_16, site_n, by = "site_id")
site_pi_16[, pi_site := ifelse(n_alleles > 1,
                               (n_alleles / (n_alleles - 1)) * (1 - sum_p2),
                               0)]

codon_4fold <- merge(codon_4fold, site_pi_16[, .(site_id, pi_site)],
                     by = "site_id", all.x = TRUE)

# Clean up large intermediates
rm(freq_entries, freq_long_16, nuc_agg, variant_entries, var_long_16,
   site_n, site_pi_16)
gc()
# 16.11) Gene-level preferred codon frequency GAMs ----
## Site-level GAMs give low R^2 because each codon position is dominated by
## stochastic noise.  Gene-level aggregation averages this noise out and yields
## a surface comparable to the gene-level CDC ~ expression + length models.
## We aggregate separately for C-ending and G-ending ROC-preferred sites.

cat("\n=== Section 16.11: Gene-Level Preferred Codon Frequency GAMs ===\n")

# Aggregate per gene: mean Preferred_Freq, split by ending
gene_pf_CG <- codon_4fold[, .(
  mean_pf   = mean(Preferred_Freq, na.rm = TRUE),
  n_sites   = .N,
  exp_norm  = Mean_Log10_Exp[1],
  log_len   = log10(Total_Codons[1])
), by = .(Gene_clean, pref_base3)]

gene_pf_C <- gene_pf_CG[pref_base3 == "C" & n_sites >= 10]
gene_pf_G <- gene_pf_CG[pref_base3 == "G" & n_sites >= 10]

cat(sprintf("  Genes with >= 10 C-ending 4-fold sites: %s\n",
            format(nrow(gene_pf_C), big.mark = ",")))
cat(sprintf("  Genes with >= 10 G-ending 4-fold sites: %s\n",
            format(nrow(gene_pf_G), big.mark = ",")))

# GAM: mean preferred freq ~ s(expression) + s(log gene length)
gam_gene_C <- gam(
  mean_pf ~ s(exp_norm, k = 10) + s(log_len, k = 10),
  data = gene_pf_C, method = "REML"
)
gam_gene_G <- gam(
  mean_pf ~ s(exp_norm, k = 10) + s(log_len, k = 10),
  data = gene_pf_G, method = "REML"
)

cat("\n--- Gene-level GAM: C-ending preferred codons ---\n")
print(summary(gam_gene_C))
cat("\n--- Gene-level GAM: G-ending preferred codons ---\n")
print(summary(gam_gene_G))

# Partial-effect plots: expression smooth from each model
exp_grid_gene <- data.frame(
  exp_norm = seq(min(gene_pf_C$exp_norm), max(gene_pf_C$exp_norm),
                 length.out = 200),
  log_len  = median(gene_pf_C$log_len)
)

pred_gene_C <- predict(gam_gene_C, newdata = exp_grid_gene, se.fit = TRUE)
pred_gene_G <- predict(gam_gene_G, newdata = exp_grid_gene, se.fit = TRUE)

gene_trend_df <- data.frame(
  exp_norm = rep(exp_grid_gene$exp_norm, 2),
  fit      = c(pred_gene_C$fit, pred_gene_G$fit),
  se       = c(pred_gene_C$se.fit, pred_gene_G$se.fit),
  Group    = rep(c("C-ending preferred", "G-ending preferred"),
                 each = nrow(exp_grid_gene))
)
gene_trend_df$lo <- gene_trend_df$fit - 1.96 * gene_trend_df$se
gene_trend_df$hi <- gene_trend_df$fit + 1.96 * gene_trend_df$se

p_gene_cg <- ggplot(gene_trend_df,
                     aes(x = exp_norm, y = fit, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("C-ending preferred" = "#2171B5",
                                 "G-ending preferred" = "#238B45")) +
  scale_fill_manual(values = c("C-ending preferred" = "#2171B5",
                                "G-ending preferred" = "#238B45")) +
  labs(
    title = "Gene-Level GAM: Mean Preferred Codon Frequency vs Expression",
    subtitle = paste0(
      "s(expression) + s(log10 gene length). Evaluated at median gene length.\n",
      sprintf("C-ending R2=%.3f, G-ending R2=%.3f",
              summary(gam_gene_C)$r.sq, summary(gam_gene_G)$r.sq)
    ),
    x = expression(log[10](Expression)),
    y = "Gene-Mean Preferred Codon Frequency",
    color = NULL, fill = NULL
  ) +
  theme_custom() +
  theme(legend.position = c(0.75, 0.85))

ggsave("./results/gene_level_pref_freq_C_vs_G.pdf",
       p_gene_cg, width = 8, height = 6)
cat("  Saved: ./results/gene_level_pref_freq_C_vs_G.pdf\n")

cat("\n=== Section 16.11: Gene-Level Preferred Codon Frequency GAMs (2D Contours) ===\n")

# Fit the GAMs using a 2D tensor product (Expression x Gene Length)
gam_gene_C_2d <- gam(
  mean_pf ~ te(exp_norm, log_len, k = c(8, 8)),
  data = gene_pf_C, method = "REML"
)
gam_gene_G_2d <- gam(
  mean_pf ~ te(exp_norm, log_len, k = c(8, 8)),
  data = gene_pf_G, method = "REML"
)

cat("\n--- 2D GAM: C-ending preferred codons ---\n")
print(summary(gam_gene_C_2d))
cat("\n--- 2D GAM: G-ending preferred codons ---\n")
print(summary(gam_gene_G_2d))

# Create a high-resolution 2D prediction grid
exp_range <- range(c(gene_pf_C$exp_norm, gene_pf_G$exp_norm), na.rm = TRUE)
len_range <- range(c(gene_pf_C$log_len, gene_pf_G$log_len), na.rm = TRUE)

pred_grid_global <- expand.grid(
  exp_norm = seq(exp_range[1], exp_range[2], length.out = 150),
  log_len  = seq(len_range[1], len_range[2], length.out = 150)
)

# Predict the surfaces for both C and G
pred_grid_global$pf_C <- as.numeric(predict(gam_gene_C_2d, newdata = pred_grid_global))
pred_grid_global$pf_G <- as.numeric(predict(gam_gene_G_2d, newdata = pred_grid_global))

# 4. Pivot longer so we can facet the plot in ggplot2
contour_global_long <- tidyr::pivot_longer(
  pred_grid_global,
  cols = c(pf_C, pf_G),
  names_to = "Ending",
  values_to = "Predicted_Pref_Freq"
) |>
  dplyr::mutate(
    Ending = ifelse(Ending == "pf_C", "C-ending preferred", "G-ending preferred")
  )

# Build the Contour Plot
p_global_contour <- ggplot(contour_global_long, aes(x = exp_norm, y = log_len)) +
  geom_raster(aes(fill = Predicted_Pref_Freq), interpolate = TRUE) +
  geom_contour(aes(z = Predicted_Pref_Freq), colour = "grey30", linewidth = 0.4, bins = 12) +
  facet_wrap(~ Ending) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Predicted\nPref Freq"
  ) +
  labs(
    title = "Global Gene-Level Preferred Codon Frequency",
    subtitle = "Interaction between Expression and Gene Length. Note the distinct topologies.",
    x = expression(log[10](Expression)),
    y = expression(log[10]("Gene length in codons"))
  ) +
  theme_custom() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# Save the plot
ggsave("./results/global_gene_level_pref_freq_contour.pdf", 
       p_global_contour, width = 12, height = 6)
cat("  Saved: ./results/global_gene_level_pref_freq_contour.pdf\n")
## ============================================================================
## END OF PAPER-REPLICATION PIPELINE
##
## To reproduce additional analyses (model selection, diagnostic plots,
## alternative parameterizations, etc.) run `full_analysis.R`, which is a
## verbatim copy of the historical pipeline.
## ============================================================================
