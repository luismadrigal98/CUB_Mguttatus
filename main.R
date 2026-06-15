##' @title CUB in Mimulus guttatus — Paper-replication pipeline
##'
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' @date   2026-05-18
##' ____________________________________________________________________________
##'
##' This script contains ONLY the code required to reproduce the figures,
##' tables, and quantitative claims in:
##'
##'     Madrigal-Roca, L. J. & Kelly, J. K. (2026).
##'     The dynamics of "silent" variation in Mimulus guttatus:
##'     Codon usage bias and linked selection.
##'
##' Section headers below map directly onto the paper's Results subsections.
##' For each section, the produced figures, tables, and cited values are listed
##' in the header comment.  Exploratory analyses, alternative parameterizations,
##' diagnostic plots, model-selection runs, and validation checks that did not
##' make it into the final manuscript live in `full_analysis.R` (verbatim copy
##' of the historical pipeline) — please run that script if you need the full
##' computational record.
##'
##' Pipeline order is preserved from full_analysis.R.  Inter-section
##' dependencies (e.g., Pi_mean_4fold pre-load in Section 5.5, preferred_codons
##' from Section 8 used by Sections 9, 11, 12, 14, 16) are unchanged.
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
                        'AnaCoDa', 'rtracklayer', 'tidyverse',
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
## RESULTS 4 — ROC-SEMPPR identifies preferred codons; selection acts only on the elite
##   Produces:
##     Table 1    Intron / intergenic stationary nucleotide frequencies
##                (`Mguttatus_intron_derived_dM.csv`, `Mguttatus_intergenic_derived_dM.csv`)
##     Figure 4   ROC-SEMPPR codon trajectories vs phi (dM-fixed-with-phi-introns)
##                (`ROC_codon_trajectories.pdf`)
##     Cited values:  74% of preferred codons C-ending (14/19)
##                     ρ = -0.329 (p = 0.011) between ROC cost Δη and
##                     expression-driven codon frequency change across 59 codons
##                     16/19 amino acid families show ROC-predicted codon
##                     increase with expression (14 statistically significant)
##   The full AnaCoDa MCMC fit is run externally (see comments at lines 919–933
##   of full_analysis.R); this script reads the saved posterior summaries.
## ============================================================================

## 8) AnaCoDa-based analysis ----
## _____________________________________________________________________________

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

# 8.2) Getting the preferred codon from the best model (dM-fixed-with_phi) ----
# We chose intron-based models because introns are more appropriate neutral
# baselines given the AT strand bias we see in coding sequences

# Using chain 1 results (independent chains are indistinguishable)

eta_data <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv")

preferred_codons <- sapply(unique(eta_data$AA), function(x) {
  AA <- x
  aa_subset <- eta_data[eta_data$AA == AA, ]
  preferred <- aa_subset[which.min(aa_subset$Mean), "Codon"]
  preferred
})

preferred_codons <- data.frame(AA = sapply(preferred_codons, 
                                           function(x) genetic_code_dna_long[[x]]),
  aa = names(preferred_codons),
  Codon = preferred_codons)

# Exporting preferred codons for polymorphism-based analysis ----

write.table(x = preferred_codons$Codon, 
            file = './results/preferred_codons.txt', 
            quote = F, row.names = F, col.names = F)

# Create unified preferred_codons_roc object for downstream analyses (CA/PCA biplots)
# This format is compatible with create_enhanced_biplot() function
preferred_codons_roc <- preferred_codons |>
  dplyr::rename(Preferred_Codons = Codon,
                Amino_Acid = AA,
                Family = aa) |>
  dplyr::mutate(Source = "ROC_SEMPPR")

message(sprintf("✓ Preferred codons from ROC model: %d amino acids\n", nrow(preferred_codons_roc) - 1)) # Serine is split in 2
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
# 3rd-position nucleotide is set by the AnaCoDa-determined preferred codon
# (e.g. Ala -> GCC -> C, Val -> GTG -> G).  Q_pref_base is then the per-gene
# fraction of 4-fold sites carrying that per-AA preferred base, summed over
# the 8 four-fold families.

fourfold_families_msd <- c("Ala", "Gly", "Pro", "Thr", "Val",
                           "Leu_4", "Ser_4", "Arg_4")
fourfold_codons_msd <- names(genetic_code_dna_long)[
  genetic_code_dna_long %in% fourfold_families_msd
]
# NOTE: in preferred_codons_roc, `Amino_Acid` holds the family label (e.g.
# "Ala", "Leu_4") from genetic_code_dna_long, while `Family` holds AnaCoDa's
# 1-letter AA code. We filter on Amino_Acid to match fourfold_families_msd.
preferred_per_AA <- preferred_codons_roc |>
  dplyr::filter(Amino_Acid %in% fourfold_families_msd) |>
  dplyr::transmute(AA = Amino_Acid,
                   Preferred_Codon = Preferred_Codons,
                   Preferred_Base  = substr(Preferred_Codons, 3, 3))

# Arg_4 is missing from the ROC output (the AnaCoDa preferred Arg codon AGG
# belongs to Arg_2). In Mimulus the dominant Arg_4 codon at high expression
# is CGC; its 3rd-position base C aligns with the GC-preference shared by
# all other 4-fold families. Using C as the Arg_4 preferred base is a
# defensible default — flagged explicitly so the manuscript can address it.
if (!"Arg_4" %in% preferred_per_AA$AA) {
  preferred_per_AA <- dplyr::bind_rows(
    preferred_per_AA,
    data.frame(AA = "Arg_4", Preferred_Codon = "CGC", Preferred_Base = "C")
  )
  cat("[Wright MSD] Arg_4 absent from ROC output; using CGC (preferred base = C) as default.\n")
}

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
  dplyr::select(Gene_name, ROC_eff, ROC_eff_4, L_ROC, Pi_mean_4fold, Mean_Log10_Exp,
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


# Binned ROC_eff and S_Wright tables ----
# 30 site-weighted ntile bins used for the pi-consistency validation (bin_roc)
# and the diversity-hump figure (bin_sw). Use two-state calibration if available.

U_bin_calib <- if (is.finite(U_intron) && is.finite(V_intron)) U_intron else
               if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) U_emp_two else U_emp
V_bin_calib <- if (is.finite(U_intron) && is.finite(V_intron)) V_intron else
               if (exists("U_emp_two") && is.finite(U_emp_two) && is.finite(V_emp_two)) V_emp_two else V_emp

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

# thr_sel: L_ROC threshold = load of the 50th-highest L_ROC gene ----
# "Outlier load group" for GO / BGS isolation = top 50 by translational load.

thr_sel <- sort(integrated_data$L_ROC, decreasing = TRUE, na.last = NA)[50]
if (!is.finite(thr_sel)) stop("thr_sel could not be determined: fewer than 50 finite L_ROC values.")
attr(thr_sel, "criterion") <- "top50_L_ROC"
attr(thr_sel, "U_empirical") <- U_gene_calib
attr(thr_sel, "V_empirical") <- V_gene_calib

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
bin_roc$pi_pred_wright <- wright_pi(bin_roc$S_Wright_bin,
                                    U = U_bin_calib, V = V_bin_calib)
bin_roc$pi_residual    <- bin_roc$pi_bin - bin_roc$pi_pred_wright

chi2_pi_terms <- (bin_roc$pi_residual)^2 /
                 pmax(bin_roc$pi_se, .Machine$double.eps)^2
chi2_pi_stat  <- sum(chi2_pi_terms, na.rm = TRUE)
chi2_pi_df    <- sum(!is.na(chi2_pi_terms))
chi2_pi_p     <- pchisq(chi2_pi_stat, df = chi2_pi_df, lower.tail = FALSE)
cat(sprintf(
  "[Validation] Bin-level pi consistency (ROC_eff bins): chi^2 = %.2f / df = %d -> p = %.3g\n",
  chi2_pi_stat, chi2_pi_df, chi2_pi_p
))

sel_bins <- bin_roc |>
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

# 8.3.5) Three-panel drift-barrier overview ----
#
# Panel A: S_Wright_signed histogram, filled by selection/drift group.
#   The vertical dashed line marks S_BARRIER (GAM-inflection threshold).
# Panel B: L_ROC density split by the same S_Wright classification.
# Panel C: ROC_eff density split by the same S_Wright classification.
#
# Logic: the barrier is derived from Q-inflection (expression level where
# selection becomes detectable). Genes at/above inflection are classified as
# "selection" if S_Wright_raw >= S_BARRIER (median of selected genes).
# Overlaying on L_ROC and ROC_eff provides cross-metric validation.

plot_barrier <- msd_data |>
  dplyr::filter(is.finite(S_Wright_signed), is.finite(L_ROC), is.finite(ROC_eff)) |>
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

rm(plot_barrier, n_sel_barrier, n_drift_barrier,
   barrier_colors, p_sw_dist, p_lroc_split, p_roc_eff_split, p_barrier_overview)

write.csv(wright_emp,
          "./results/Wright_curve_empirical.csv",        row.names = FALSE)
write.csv(
  msd_data |> dplyr::select(Gene_name, ROC_eff, ROC_eff_4, L_ROC, Q_pref_base,
                            S_Wright_signed, S_Wright_raw, is_drift,
                            pi_2allele,
                            Mean_Log10_Exp, Max_Log10_Exp, N_4fold_sites),
  "./results/Wright_per_gene_ROC_eff_S_Wright.csv", row.names = FALSE
)
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
    thr_sel                      = as.numeric(thr_sel),
    n_above_thr_sel              = sum(integrated_data$L_ROC > as.numeric(thr_sel),
                                       na.rm = TRUE),
    n_above_S_BARRIER            = sum(msd_data$S_Wright_signed >= S_BARRIER,
                                       na.rm = TRUE),
    n_drift_genes                = sum(msd_data$is_drift, na.rm = TRUE),
    frac_drift_genes             = mean(msd_data$is_drift, na.rm = TRUE),
    chi2_pi_stat                 = chi2_pi_stat,
    chi2_pi_df                   = chi2_pi_df,
    chi2_pi_p                    = chi2_pi_p,
    chi2_pi_sel_stat             = if (is.finite(chi2_pi_sel_stat)) chi2_pi_sel_stat else NA_real_,
    chi2_pi_sel_df               = chi2_pi_sel_df,
    chi2_pi_sel_p                = if (is.finite(chi2_pi_sel_p))    chi2_pi_sel_p    else NA_real_,
    cor_ROC_eff_4_S_Wright_spearman = cor_roc_eff_4_spearman,
    cor_ROC_eff_4_S_Wright_pearson  = cor_roc_eff_4_pearson
  ),
  "./results/Wright_threshold_adopted.csv", row.names = FALSE
)

# Memory cleanup -- keep: thr_sel, U_emp, V_emp, Q_neutral_obs,
# pi_neutral_obs, pi_neutral_theory, S_BARRIER, S_BARRIER_advisor,
# msd_data, bin_roc, bin_sw, integrated_data.
rm(codon_4fold_counts, N_4fold_sites, N_preferred_base, gene_Q_4fold,
   preferred_codon_set, fourfold_codon_table, preferred_per_AA,
   S_grid_emp,
   neutral_pool, neutral_pool_pi, pi_data_operational,
   Q_neutral_two, pi_neutral_two,
   chi2_pi_terms, sel_bins, per_gene_pool)
gc()

# 8.4) GO-enrichment for two selection groups ----
#
#   (a) Top 50 by L_ROC  -> "load-paying" group (Rubisco/photosynthesis enrichment)
#   (b) S_Wright_signed >= S_BARRIER -> "selection" group (drift-barrier genes)

custom_bag <- integrated_data |> dplyr::pull(Gene_name)

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

# Backwards-compatible alias
subset_strongly_shaped_by_s <- subset_load_paying
GO_results <- GO_results_load
write.csv(x = GO_results$result |> dplyr::select(-parents),
          file = "./results/Go_enrichment.csv",
          quote = TRUE, row.names = FALSE)

# 8.4b) GO enrichment dot-plot visualisation ----
#
# Filters out overly generic terms (term_size > 500) and shows the top N
# most significant remaining terms. Gene ratio (precision) drives the x-axis;
# dot size = overlap count; colour = -log10(FDR p-value).

.go_dotplot <- function(go_result, title, max_terms = 20, max_term_size = 500) {
  if (is.null(go_result) || nrow(go_result$result) == 0) {
    return(
      ggplot() +
        labs(title = title, subtitle = "No significant GO terms after filtering") +
        theme_void()
    )
  }

  df <- go_result$result |>
    dplyr::filter(term_size <= max_term_size) |>
    dplyr::arrange(p_value) |>
    dplyr::slice_head(n = max_terms) |>
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      gene_ratio  = precision,
      term_label  = sapply(
        term_name,
        function(x) paste(strwrap(x, width = 45), collapse = "\n")
      )
    )

  if (nrow(df) == 0) {
    return(
      ggplot() +
        labs(title = title, subtitle = "No specific GO terms (all filtered as generic)") +
        theme_void()
    )
  }

  ggplot(df, aes(x = gene_ratio,
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
}

p_go_load <- .go_dotplot(
  GO_results_load,
  sprintf("High translational load  (top 50 by L_ROC,  n = %d)", length(subset_load_paying))
)

p_go_sel <- .go_dotplot(
  GO_results_selection,
  sprintf("Population-genetic selection  (S_Wright ≥ %.1f,  n = %d)",
          S_BARRIER, length(subset_selection))
)

ggsave("./results/GO_dotplot_load.pdf",      p_go_load, width = 8, height = 6)
ggsave("./results/GO_dotplot_selection.pdf", p_go_sel,  width = 8, height = 6)

p_go_combined <- (p_go_load | p_go_sel) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("./results/GO_dotplot_combined.pdf", p_go_combined, width = 16, height = 7)

rm(p_go_load, p_go_sel, p_go_combined, .go_dotplot)
gc()

# 8.5) Top genes by L_ROC (load) and by S_Wright (selection) -----------------------------------

detailed_annotation_full <- read.delim(
  "data/Mguttatusvar_IM767_887_v2.1.annotation_info.txt",
  header = TRUE, sep = "\t", comment.char = "", quote = "",
  fill = TRUE, na.strings = ""
) |>
  dplyr::select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline) |>
  dplyr::distinct()

# Top 50 by L_ROC (load-paying)
top_L_ROC <- integrated_data |>
  dplyr::arrange(desc(L_ROC)) |>
  dplyr::select(Gene_name, L_ROC)

top_L_ROC <- top_L_ROC[1:50,]

top_L_ROC <- top_L_ROC |>
  left_join(detailed_annotation_full, by = join_by("Gene_name" == "locusName"))

write.csv(top_L_ROC,
          "./results/Top_genes_strong_selection_load.csv",
          quote = TRUE, row.names = FALSE)

# Selection group: S_Wright >= S_BARRIER
top_selection <- msd_data |>
  dplyr::filter(!is.na(S_Wright_raw), S_Wright_raw >= S_BARRIER) |>
  dplyr::arrange(desc(S_Wright_raw)) |>
  dplyr::select(Gene_name, S_Wright_raw, ROC_eff_4, L_ROC, Mean_Log10_Exp) |>
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
  dplyr::select(Amino_Acid = AA, aa, Codon_RNA)

# Collapse amino acids with six codons back into six, based on relative adaptiveness

preferred_codons_comparative <- preferred_codons_comparative |>
  dplyr::mutate(AA_root = sapply(preferred_codons_comparative$Amino_Acid, 
                                 function(x) 
                                 {
                                   unlist(strsplit(x, "_"))[1]
                                 }))

preferred_codons_comparative <- preferred_codons_comparative |>
  left_join(eta_data |> dplyr::select(AA, Mean) |> dplyr::rename(aa = AA, eta = Mean))

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
## RESULTS 8 — Polymorphism reveals linked selection and codon-usage selection
##   Produces:
##     Figure 7A  Synonymous & non-synonymous π vs expression rank
##                (`pi_4fold_by_expression_rank.pdf`)
##     Figure 7B  Preferred codon frequency vs expression group
##                (`Frequency_preferred_by_expression_group_Median_CI.pdf`)
##     Figure 7C  Intron vs exon (4-fold) π by distance from gene start
##                (`pi_by_gene_distance.pdf`)
##     Figure 6B  Preferred-codon-frequency landscape vs expression × gene length
##                (`Preferred_freq_contour_exp_x_length_broad.pdf`,
##                 `Preferred_freq_contour_exp_x_length_narrow.pdf`)
## ============================================================================

## 12) Polymorphism data integration ----
## _____________________________________________________________________________
# By-gene polymorphism (Pi_mean_4fold etc.) was preloaded in Section 5.5 because
# Section 6 GAMs and Section 8.3.4 msd_data depend on those columns. This
# section continues with the per-feature positional and decomposed analyses.

# Memory cleanup: polymorphism raw data (now joined into integrated_data) ---
rm(pi_data)

# 12.0) Positional decomposition of pi: first 300 bp vs after 300 bp ----
# Compute per-gene 4-fold and 0-fold pi split at the 300 bp CDS boundary.
# This captures Hill-Robertson interference patterns: linked selection near
# the 5' end (translational ramp) vs the gene body.

pi_feature <- data.table::fread("data/all_chromosomes.pi_per_gene_feature.txt")

# Harmonize gene names: feature file has "MgIM767.01G000100.v2.1",
# integrated_data has "MgIM767.01G000100"
pi_feature[, Gene := sub("\\.v[0-9.]+$", "", Gene)]

cat(sprintf("\n=== Section 12.0: Positional pi decomposition | %d unique genes in feature file (no AnaCoDa filter applied) ===\n",
            length(unique(pi_feature$Gene))))

# Per-exon CDS sizes (from "all" degeneracy) and degeneracy-specific pi
exon_all <- pi_feature[Feature_Type == "exon" & Degeneracy == "all",
                       .(Gene, Feature_Num, Exon_Sites = Sites)]
exon_4fold <- pi_feature[Feature_Type == "exon" & Degeneracy == "4-fold",
                         .(Gene, Feature_Num, Sites_4f = Sites, Pi_sum_4f = Pi_sum)]
exon_0fold <- pi_feature[Feature_Type == "exon" & Degeneracy == "0-fold",
                         .(Gene, Feature_Num, Sites_0f = Sites, Pi_sum_0f = Pi_sum)]

exon_data <- merge(exon_all, exon_4fold, by = c("Gene", "Feature_Num"), all.x = TRUE)
exon_data <- merge(exon_data, exon_0fold, by = c("Gene", "Feature_Num"), all.x = TRUE)
exon_data[is.na(Sites_4f), c("Sites_4f", "Pi_sum_4f") := 0]
exon_data[is.na(Sites_0f), c("Sites_0f", "Pi_sum_0f") := 0]

# Order exons and compute cumulative CDS position
data.table::setorder(exon_data, Gene, Feature_Num)
exon_data[, cum_end := cumsum(Exon_Sites), by = Gene]
exon_data[, cum_start := cum_end - Exon_Sites + 1L, by = Gene]

# Fraction of each exon falling within the first 300 bp of CDS
bp_cutoff <- 300
exon_data[, frac_first := data.table::fifelse(
  cum_end <= bp_cutoff, 1.0,
  data.table::fifelse(cum_start > bp_cutoff, 0.0,
                      (bp_cutoff - cum_start + 1) / Exon_Sites)
)]

# Aggregate per gene: first 300 bp vs after 300 bp
pi_300bp <- exon_data[, .(
  Sites_4fold_first300  = sum(Sites_4f * frac_first),
  Pi_sum_4fold_first300 = sum(Pi_sum_4f * frac_first),
  Sites_4fold_after300  = sum(Sites_4f * (1 - frac_first)),
  Pi_sum_4fold_after300 = sum(Pi_sum_4f * (1 - frac_first)),
  Sites_0fold_first300  = sum(Sites_0f * frac_first),
  Pi_sum_0fold_first300 = sum(Pi_sum_0f * frac_first),
  Sites_0fold_after300  = sum(Sites_0f * (1 - frac_first)),
  Pi_sum_0fold_after300 = sum(Pi_sum_0f * (1 - frac_first))
), by = .(Gene_name = Gene)]

# Merge into integrated_data
integrated_data <- integrated_data |>
  dplyr::left_join(pi_300bp, by = "Gene_name")

cat(sprintf("300 bp decomposition: %d genes matched\n",
            sum(!is.na(integrated_data$Sites_4fold_first300))))

# 12.0b) π by distance from gene start — exon (4-fold) vs intron ----
# Assigns each feature a cumulative genomic position by interleaving exons and
# introns in Feature_Num order: exon k → genomic rank 2k-1, intron k → 2k.
# Feature bp length is approximated by "all"-degeneracy Sites (surveyed sites).
# Exons use 4-fold π; introns use all-site π (the only class they carry).

cat(sprintf("\n=== Section 12.0b: pi vs distance from gene start | %d unique genes (pi_feature; no AnaCoDa filter) ===\n",
            length(unique(pi_feature$Gene))))

feat_pos <- pi_feature[Degeneracy == "all",
                        .(Gene, Feature_Type, Feature_Num, feat_bp = Sites)]

feat_pos[Feature_Type == "exon",   genomic_order := 2L * Feature_Num - 1L]
feat_pos[Feature_Type == "intron", genomic_order := 2L * Feature_Num]
data.table::setorder(feat_pos, Gene, genomic_order)

# Cumulative distance from gene start (approximation via surveyed-site lengths)
feat_pos[, cum_end   := cumsum(feat_bp),           by = Gene]
feat_pos[, cum_start := cum_end - feat_bp,          by = Gene]
feat_pos[, midpoint  := (cum_start + cum_end) / 2L, by = Gene]

# 1 kb windows up to 10 kb, last bin 10–20 kb
dist_breaks <- c(0, seq(1000, 10000, 1000), 20000)
dist_labels <- c(as.character(seq(1000, 10000, 1000)), "10k-20k")

feat_pos[, dist_bin := cut(midpoint, breaks = dist_breaks, labels = dist_labels,
                            right = TRUE, include.lowest = FALSE)]

# π values at the appropriate degeneracy class per feature type
pi_exon_4f <- pi_feature[Feature_Type == "exon"   & Degeneracy == "4-fold",
                          .(Gene, Feature_Type, Feature_Num, Pi_sum, n_pi = Sites)]
pi_int_all <- pi_feature[Feature_Type == "intron" & Degeneracy == "all",
                          .(Gene, Feature_Type, Feature_Num, Pi_sum, n_pi = Sites)]
pi_vals_dist <- data.table::rbindlist(list(pi_exon_4f, pi_int_all))

feat_pi_dist <- merge(
  pi_vals_dist,
  feat_pos[, .(Gene, Feature_Type, Feature_Num, dist_bin)],
  by = c("Gene", "Feature_Type", "Feature_Num"), all.x = TRUE
)

# Weighted mean π per feature type × distance window
pi_by_dist <- feat_pi_dist[n_pi > 0 & !is.na(dist_bin), .(
  n_features = .N,
  Pi_mean    = sum(Pi_sum) / sum(n_pi)
), by = .(Feature_Type, dist_bin)]

pi_by_dist[, dist_bin := factor(dist_bin, levels = dist_labels)]
data.table::setorder(pi_by_dist, Feature_Type, dist_bin)

plot_pi_by_dist <- ggplot(
    pi_by_dist,
    aes(x = dist_bin, y = Pi_mean, fill = Feature_Type)
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c("exon" = "#1F4E79", "intron" = "#E07830"),
    labels = c("exon" = "Exon (4-fold)", "intron" = "Intron")
  ) +
  labs(
    x = "Distance from gene start (bp)",
    y = "Nucleotide Diversity",
    fill = NULL
  ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./results/pi_by_gene_distance.pdf", plot_pi_by_dist, width = 10, height = 6)

rm(feat_pos, pi_exon_4f, pi_int_all, pi_vals_dist, feat_pi_dist, pi_by_dist,
   plot_pi_by_dist, dist_breaks, dist_labels)

rm(pi_feature, exon_all, exon_4fold, exon_0fold, exon_data, pi_300bp)
gc()

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
  dplyr::filter(L_ROC < thr_sel) |>
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
  dplyr::filter(L_ROC > thr_sel) |>
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
  rbind(sel_cat) # Final bin holds genes with L_ROC > thr_sel (Wright-calibrated)

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
    dplyr::filter(L_ROC < thr_sel) |>
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
    dplyr::filter(L_ROC > thr_sel) |>
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
## ============================================================================
## RESULTS 10 — Translational ramp: preferred codons cluster near the 5' end
##   Produces:
##     Figure 6A  Preferred codon frequency vs codon position by expression tier
##                (`ramp_by_expression_polymorphism.pdf`)
##     Also runs the gene-level GAM that produces the supplementary contour
##     `global_gene_level_pref_freq_contour.pdf` (Fig S10 / Fig 6B alt).
## ============================================================================

# ******************************************************************************
# 15) Testing the translational ramp hypothesis ----
# ______________________________________________________________________________

# 15.1) Polymorphism-based ramp models ----
# Per-codon preferred-allele frequencies in the first 200 codons; aggregate into
# 5-codon windows; fit beta-regression bam() with random intercept by gene.
# We subsample 3,000 genes (seed = 1998) — the gene RE makes the full set
# prohibitively slow without changing the population-level smooth.

poly_data <- fread(
  "data/all_chromosomes.codon_frequencies_preferred.txt",
  select = c("Gene", "Codon_Pos", "Preferred_Freq", "Non_Preferred_Freq"),
  showProgress = FALSE
)

poly_data <- poly_data |>
  dplyr::mutate(Gene_clean = paste0("MgIM767.", Gene)) |>
  dplyr::rename(Position = Codon_Pos) |>
  dplyr::filter(Position <= 200)

poly_with_exp <- poly_data |>
  dplyr::left_join(
    integrated_data |>
      dplyr::select(Gene_name, Max_Log10_Exp, Exp_breadth),
    by = c("Gene_clean" = "Gene_name")
  ) |>
  dplyr::filter(!is.na(Max_Log10_Exp), !is.na(Exp_breadth)) |>
  dplyr::mutate(
    Exp_Z = as.numeric(scale(Max_Log10_Exp)),
    Breadth_Z = as.numeric(scale(Exp_breadth)),
    Gene_clean = factor(Gene_clean)
  )

cat(sprintf("Loaded %d codon positions from %d genes\n",
            nrow(poly_with_exp), length(unique(poly_with_exp$Gene_clean))))

n_subsample_genes <- 3000
all_genes_15 <- unique(poly_with_exp$Gene_clean)
set.seed(1998)
sampled_genes_15 <- sample(all_genes_15,
                           size = min(n_subsample_genes, length(all_genes_15)))
poly_with_exp <- poly_with_exp |>
  dplyr::filter(Gene_clean %in% sampled_genes_15) |>
  dplyr::mutate(Gene_clean = droplevels(Gene_clean))

cat(sprintf("After subsampling: %d codon positions from %d genes (seed = 1998)\n",
            nrow(poly_with_exp), length(unique(poly_with_exp$Gene_clean))))

window_size <- 5
poly_agg <- poly_with_exp |>
  dplyr::mutate(Window = ceiling(Position / window_size)) |>
  dplyr::group_by(Gene_clean, Window, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    Position_mid = mean(Position),
    Preferred_Freq_mean = mean(Preferred_Freq, na.rm = TRUE),
    n_codons = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::filter(Position_mid <= 200) |>
  dplyr::mutate(
    Preferred_Freq_beta = dplyr::case_when(
      Preferred_Freq_mean <= 0.001 ~ 0.001,
      Preferred_Freq_mean >= 0.999 ~ 0.999,
      TRUE ~ Preferred_Freq_mean
    )
  )

if (any(poly_agg$Preferred_Freq_beta <= 0 | poly_agg$Preferred_Freq_beta >= 1)) {
  stop("Beta regression requires values strictly between 0 and 1")
}

fit_ramp_poly <- bam(
  Preferred_Freq_beta ~
    s(Position_mid, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  data = poly_agg, family = betar(),
  method = "fREML", discrete = TRUE, nthreads = 1
)

fit_ramp_int_poly <- bam(
  Preferred_Freq_beta ~
    s(Position_mid, k = 10, bs = "tp") +
    s(Position_mid, by = Exp_Z, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  data = poly_agg, family = betar(),
  method = "fREML", discrete = TRUE, nthreads = 1
)

# Plot 1: Global ramp shape from polymorphism data
pred_positions <- data.frame(
  Position_mid = seq(5, 200, by = 2),
  Exp_Z = 0,
  Breadth_Z = 0,
  Gene_clean = poly_agg$Gene_clean[1]
)

pred_ramp <- predict(fit_ramp_poly, newdata = pred_positions,
                     type = "response", se.fit = TRUE,
                     exclude = "s(Gene_clean)",
                     unconditional = TRUE)

pred_positions$fit   <- pred_ramp$fit
pred_positions$lower <- pred_ramp$fit - 1.96 * pred_ramp$se.fit
pred_positions$upper <- pred_ramp$fit + 1.96 * pred_ramp$se.fit

plot_ramp_poly <- ggplot(pred_positions, aes(x = Position_mid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = fit), color = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Translational Ramp: Population Polymorphism Data",
    subtitle = "Frequency of preferred codons across positions",
    x = "Codon Position",
    y = "Mean Preferred Codon Frequency",
    caption = "Ribbon = +/-1.96 SE | Data from population genomics"
  ) +
  theme_bw(base_size = 12)

ggsave("./results/translational_ramp_polymorphism.pdf",
       plot_ramp_poly, width = 10, height = 6)

# Plot 2: Ramp by expression level
pred_grid_exp <- expand.grid(
  Position_mid = seq(5, 200, by = 5),
  Exp_Z = c(-1.5, 0, 1.5),
  Breadth_Z = 0
) |>
  dplyr::mutate(
    Gene_clean = poly_agg$Gene_clean[1],
    Exp_level = factor(
      Exp_Z,
      levels = c(-1.5, 0, 1.5),
      labels = c("Low Expression", "Medium Expression", "High Expression")
    )
  )

pred_exp <- predict(fit_ramp_int_poly, newdata = pred_grid_exp,
                    type = "response", se.fit = TRUE,
                    exclude = "s(Gene_clean)")

pred_grid_exp$fit   <- pred_exp$fit
pred_grid_exp$lower <- pred_exp$fit - 1.96 * pred_exp$se.fit
pred_grid_exp$upper <- pred_exp$fit + 1.96 * pred_exp$se.fit

plot_ramp_exp_poly <- ggplot(pred_grid_exp,
                             aes(x = Position_mid, color = Exp_level, fill = Exp_level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Ramp Shape by Expression: Population Data",
    subtitle = "Does selection for preferred codons vary by expression level?",
    x = "Codon Position",
    y = "Mean Preferred Codon Frequency",
    color = NULL, fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("./results/ramp_by_expression_polymorphism.pdf",
       plot_ramp_exp_poly, width = 10, height = 6)

# Memory cleanup: Section 15 large intermediates ---
rm(poly_data, poly_with_exp, poly_agg,
   pred_positions, pred_ramp, plot_ramp_poly,
   pred_grid_exp, pred_exp, plot_ramp_exp_poly,
   sampled_genes_15, all_genes_15, n_subsample_genes, window_size)
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

# Normalised distance: Codon_Pos / (Total_Codons - 1) -> [0, 1]
codon_4fold[, dist_norm := Codon_Pos / pmax(Total_Codons - 1L, 1L)]
codon_4fold[dist_norm > 1, dist_norm := 1]                      # safety cap

codon_4fold[, exp_norm := Mean_Log10_Exp]

# Store expression range for prediction grids (contour plots)
exp_range_16 <- range(codon_4fold$exp_norm, na.rm = TRUE)

# Standardised predictors (z-scored) for comparable coefficient magnitudes
# Raw predictors are retained for contour-plot prediction grids
codon_4fold[, dist_z := as.numeric(scale(dist_norm))]
codon_4fold[, exp_z  := as.numeric(scale(exp_norm))]

cat(sprintf("  After expression merge: %s sites in %s genes\n",
            format(nrow(codon_4fold), big.mark = ","),
            format(length(unique(codon_4fold$Gene_clean)), big.mark = ",")))
cat(sprintf("  dist_norm: mean=%.3f, sd=%.3f  |  exp_norm: mean=%.3f, sd=%.3f\n",
            mean(codon_4fold$dist_norm), sd(codon_4fold$dist_norm),
            mean(codon_4fold$exp_norm),  sd(codon_4fold$exp_norm)))

# 3rd-position base of the ROC-preferred codon at each site
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
