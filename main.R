##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' 
##' @date 02/24/2026
##' ____________________________________________________________________________

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
                        'betareg', 'brms', 'cmdstanr')

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

# Full integration with your pipeline
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

# Calculate expected ENC under mutation-drift equilibrium (Wright 1990)
# ENC_expected = 2 + GC3s + 29/(GC3s^2 + (1-GC3s)^2)
gc3s_range <- seq(0, 1, by = 0.01)
enc_expected <- 2 + gc3s_range + 29 / (gc3s_range^2 + (1 - gc3s_range)^2)

expected_curve <- data.frame(
  GC3s = gc3s_range,
  ENC_expected = enc_expected
)

# Create enhanced ENC plot
p_enc_cdc <- ggplot() +
  # Background: non-significant genes
  geom_point(data = integrated_data |> filter(!CDC_significant | is.na(CDC_significant)),
             aes(x = GC3s, y = ENC), 
             color = "gray70", alpha = 0.3, size = 0.8) +
  # Foreground: CDC-significant genes
  geom_point(data = integrated_data |> filter(CDC_significant),
             aes(x = GC3s, y = ENC, color = CDC_category), 
             size = 2, alpha = 0.7) +
  # Expected neutrality curve (Wright 1990)
  geom_line(data = expected_curve, 
            aes(x = GC3s, y = ENC_expected), 
            color = "black", linewidth = 1.2, linetype = "solid") +
  scale_color_manual(
    values = c("p < 0.001" = "#d73027",     # Dark red
               "p < 0.01" = "#fc8d59",      # Orange
               "p < 0.05" = "#fee08b"),     # Yellow
    name = "CDC Significance"
  ) +
  labs(
    title = "ENC Plot with CDC-Significant Genes Highlighted",
    subtitle = sprintf("Black curve = expected ENC under mutation-drift equilibrium\nGenes below curve = selection for codon bias\n%d genes (%.1f%%) have significant CDC", 
                       n_sig, pct_sig),
    x = "GC3s (GC content at synonymous 3rd codon positions)",
    y = "ENC (Effective Number of Codons)"
  ) +
  ylim(20, 61) +
  xlim(0, 1) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave("./results/ENC_plot_CDC_highlighted.pdf", p_enc_cdc, 
       width = 11, height = 7)

message("Enhanced ENC plot saved: ./results/ENC_plot_CDC_highlighted.pdf\n\n")

# Analyze CDC-significant genes: are they below the curve (under selection)?
cat("=== Position Analysis: CDC-Significant Genes Relative to Neutrality Curve ===\n")

# Calculate deviation from expected ENC
integrated_data <- integrated_data |>
  mutate(
    ENC_expected = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2),
    ENC_deviation = ENC - ENC_expected,
    Below_curve = ENC_deviation < 0
  )

# Compare CDC-significant vs non-significant genes
cdc_position_summary <- integrated_data |>
  dplyr::filter(!is.na(CDC_significant)) |>
  dplyr::group_by(CDC_significant) |>
  summarize(
    n = n(),
    mean_ENC = mean(ENC, na.rm = TRUE),
    mean_ENC_expected = mean(ENC_expected, na.rm = TRUE),
    mean_deviation = mean(ENC_deviation, na.rm = TRUE),
    pct_below_curve = 100 * sum(Below_curve, na.rm = TRUE) / n(),
    mean_CDC = mean(CDC, na.rm = TRUE)
  )

print(cdc_position_summary)

# Statistical tests
if (n_sig > 0) {
  cat("\n=== Statistical Comparisons ===\n")
  
  # Test if CDC-significant genes have different ENC deviation
  wilcox_enc <- wilcox.test(
    ENC_deviation ~ CDC_significant,
    data = integrated_data |> filter(!is.na(CDC_significant))
  )
  cat(sprintf("ENC deviation (CDC-sig vs non-sig): W = %.0f, p = %.2e\n", 
              wilcox_enc$statistic, wilcox_enc$p.value))
  
  # Test if more CDC-significant genes are below the curve
  below_curve_table <- table(
    integrated_data |> dplyr::filter(!is.na(CDC_significant)) |> dplyr::select(CDC_significant, Below_curve)
  )
  fisher_test <- fisher.test(below_curve_table)
  cat(sprintf("Position relative to curve (Fisher test): OR = %.2f, p = %.2e\n", 
              fisher_test$estimate, fisher_test$p.value))
}

## *****************************************************************************
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

integrated_data <- integrated_data |>
  dplyr::left_join(pi_data, by = "Gene_name") |>
  na.exclude()

rm(pi_data)

## *****************************************************************************
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

# VISUALIZING THE RATE OF CHANGE (The Slope)

p_slopes <- plot_slopes(m_interaction, 
                        variables = "Max_Log10_Exp",
                        condition = c("Max_Log10_Exp", "Exp_breadth"), 
                        newdata = datagrid(
                          CDS_length_nt = mean((integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0))$CDS_length_nt)),
                        type = "response") +
  # RED LINE = Zero Slope (Where the relationship flattens/saturates)
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_rug(data = integrated_data |> dplyr::filter(Exp_breadth > 0,                                           Pi_mean_4fold > 0), 
           aes(x = Max_Log10_Exp), 
           sides = "b", alpha = 0.05, inherit.aes = FALSE) +
  
  # THEME & LABELS
  theme_custom() +
  labs(y = "Slope (Change in CDC / Change in Exp)",
       x = "Max Expression (Log10 CPM)")

ggsave("./results/GAM_Interaction_Slopes.pdf", p_slopes, width = 10, height = 6)

# FORMAL HYPOTHESIS TESTING
# Question: Is the rate of bias accumulation significantly different?

# Calculate average slopes for specific breadths (e.g., Narrow vs Broad)
slopes <- avg_slopes(
  m_interaction,
  variables = "Max_Log10_Exp",
  by = "Exp_breadth", 
  newdata = datagrid(Exp_breadth = c(1, 25, 29))
)

# Plotting detrended ENC against expression

confounder_model_gam <- gam(CDC ~ s(CDS_length_nt),
                            data = integrated_data,
                            family = betar(link = "logit"))

integrated_data$CDC_detrended <- residuals(confounder_model_gam, 
                                           type = "response")

p_detrended <- ggplot(integrated_data, aes(x = Max_Log10_Exp, 
                                           y = CDC_detrended)) +
  # Use density to handle the 22k points
  geom_pointdensity(adjust = 0.1) + 
  scale_color_viridis_c() +
  
  # This allows the J-shape (Selection threshold) to appear
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
  
  # Add a horizontal line at 0 (No deviation from length expectation)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  labs(
    title = "Isolated Effect of Expression on CDC",
    subtitle = "Residuals after removing the effect of CDS Length",
    y = "Deviation from Expected CDC (given Length)",
    x = "Max Log10 Expression"
  ) +
  theme_custom()

ggsave("./results/CDC_detrended_vs_expression_curved.pdf", p_detrended, 
       width = 8, height = 6)

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
# Keeping: selection_table, test_broad_vs_narrow, test_medium_vs_narrow,
#          kw_detrended, cdc_position_summary, slopes, contrasts
rm(m_null, m_additive, m_complex, m_interaction, model_list,
   justification_list, justification_table, plot_list, combined_plot,
   p_effects, p_slopes,
   m_noise, confounder_model_gam, p_detrended,
   p_boxplot_detrended, p_medians, p_constraint, p_enc_cdc,
   expected_curve, gc3s_range, enc_expected,
   plot_data, my_comparisons,
   top5_cdc_de, middle_cdc_de, bottom5_cdc_de,
   n_sig, n_total, pct_sig,
   top_5_cutoff, bottom_5_cutoff)
gc()

## *****************************************************************************
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
integrated_data <- integrated_data |>
  left_join(cai_values, by = "Gene_name")

cat("\n=== CAI vs Expression Level ===\n")
# Compare CAI across expression groups
cai_by_group <- integrated_data |>
  dplyr::group_by(Expression_Group) |>
  summarise(
    n = n(),
    mean_CAI = mean(CAI, na.rm = TRUE),
    median_CAI = median(CAI, na.rm = TRUE),
    sd_CAI = sd(CAI, na.rm = TRUE),
    mean_ENC = mean(ENC, na.rm = TRUE)
  )

print(cai_by_group)

# Statistical tests for three groups
message("\n=== Kruskal-Wallis Test: CAI across All Three Groups ===\n")
message("H0: All three groups have the same median CAI\n")
kw_test <- kruskal.test(CAI ~ Expression_Group, data = integrated_data)
print(kw_test)

if (kw_test$p.value < 0.05) {
  message("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  message("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Perform Dunn's test with FDR correction
  dunn_result <- dunn.test::dunn.test(
    x = integrated_data$CAI,
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
  message("\nNo significant difference among groups (p >= 0.05)\n")
  message("Post-hoc tests not necessary.\n")
}

# Additional pairwise effect sizes
message("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

# Get CAI values for each group
top_cai <- integrated_data |> dplyr::filter(Expression_Group == "Top 5%") |> pull(CAI)
middle_cai <- integrated_data |> dplyr::filter(Expression_Group == "Middle 90%") |> pull(CAI)
bottom_cai <- integrated_data |> dplyr::filter(Expression_Group == "Bottom 5%") |> pull(CAI)

# If groups don't exist with new names, try old names
if (length(top_cai) == 0) {
  top_cai <- integrated_data |> dplyr::filter(Expression_Group == "Top 5% Expressed") |> pull(CAI)
}
if (length(bottom_cai) == 0) {
  bottom_cai <- integrated_data |> dplyr::filter(Expression_Group == "Bottom 95%") |> pull(CAI)
}

# Calculate effect sizes
if (length(top_cai) > 0 && length(middle_cai) > 0) {
  d_top_middle <- cohens_d_calc(top_cai, middle_cai)
  message(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle))
}

if (length(top_cai) > 0 && length(bottom_cai) > 0) {
  d_top_bottom <- cohens_d_calc(top_cai, bottom_cai)
  message(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom))
}

if (length(middle_cai) > 0 && length(bottom_cai) > 0) {
  d_middle_bottom <- cohens_d_calc(middle_cai, bottom_cai)
  message(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom))
}

message("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

# Plot CAI by expression group
p_cai_boxplot <- ggplot(integrated_data, aes(x = Expression_Group, y = CAI, fill = Expression_Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(outlier.alpha = 0.3)  +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "Codon Adaptation Index by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR",
       y = "CAI (Codon Adaptation Index)",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/CAI_by_expression_group.pdf", p_cai_boxplot, width = 8, height = 6)

# Median + CI

plot_data_cai <- integrated_data |>
  dplyr::mutate(Exp_Group = factor(Expression_Group, 
                                   levels = c("Bottom 5%", "Middle 90%", "Top 5%"))) |>
  dplyr::filter(!is.na(Exp_Group))

# 3. Create the Plot
p_cai_median <- ggplot(plot_data_cai, aes(x = Exp_Group, y = CAI)) +
  
  # A. Median and 95% CI
  stat_summary(fun.data = median_cl_boot, 
               geom = "errorbar", width = 0.15, linewidth = 0.8, color = "black") +
  stat_summary(fun = median, geom = "point", size = 4, aes(color = Exp_Group)) +
  
  # B. Formatting
  # Using your specific colors
  scale_color_manual(values = c("Bottom 5%" = "#377EB8", 
                                "Middle 90%" = "#999999", 
                                "Top 5%" = "#E41A1C")) +
  
  labs(y = "CAI (Codon Adaptation Index)",
       x = NULL) + # Remove X label as groups are self-explanatory
  
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        panel.grid.major.x = element_blank())

# Save
ggsave("./results/CAI_by_expression_group_Median_CI.pdf", p_cai_median, 
       width = 4, height = 3.5)

# Memory cleanup: Section 7 CAI intermediates ---
# Keeping: cai_results, w_table, cai_by_group, kw_test, comparison_table, wilcox_test
rm(p_cai_boxplot, p_cai_median, plot_data_cai,
   top_cai, middle_cai, bottom_cai, reference_genes)
gc()

## *****************************************************************************
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

# Optional: Additional analysis on the intermediate data
# For example, cluster genomic windows by mutational spectrum
if (!is.null(dM_results$intermediates$introns$nuc_filtered) &
    !is.null(dM_results$intermediates$intergenic$nuc_filtered)) {
  
  # Prepare window data for clustering (optional advanced analysis)
  window_data_introns <- dM_results$intermediates$introns$nuc_filtered |>
    dplyr::select(window_idx, pi_A, pi_C, pi_G, pi_T)
  
  window_data_intergenic <- dM_results$intermediates$intergenic$nuc_filtered |>
    dplyr::select(window_idx, pi_A, pi_C, pi_G, pi_T)
  
  # PCA summary
  pca_introns <- prcomp(
    x = as.matrix(window_data_introns[, c("pi_A", "pi_C", "pi_G", "pi_T")]),
    center = TRUE,
    scale. = TRUE
  )
  
  pca_intergenic <- prcomp(
    x = as.matrix(window_data_intergenic[, c("pi_A", "pi_C", "pi_G", "pi_T")]),
    center = TRUE,
    scale. = TRUE
  )
  
  cat("\nPCA of nucleotide composition (introns):\n")
  print(summary(pca_introns))
  
  cat("\nPCA of nucleotide composition (intergenic):\n")
  print(summary(pca_intergenic))
  
  # GMM clustering (optional - may find no evidence for multiple clusters)
  clusters_localM_introns <- make_clusters(data = window_data_introns[, -1], G = 1:10)
  clusters_localM_intergenic <- make_clusters(data = window_data_intergenic[, -1], G = 1:10)
}

# Memory cleanup: dM estimation intermediates ---
# Keeping: dM_results (mutation rate estimates)
if (exists("window_data_introns")) rm(window_data_introns, window_data_intergenic,
                                       pca_introns, pca_intergenic,
                                       clusters_localM_introns, clusters_localM_intergenic)
gc()

# No evidence for more than one cluster (genome-wide mutational pressure)

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

# Memory cleanup: convergence diagnostics, phi comparisons, trace objects ---
# Keeping: dM_fixed_with_phi_conv, dM_fixed_intergenic, dM_fixed_with_phi_intergenic,
#          phi_hat_dM_fixed_intergenic, phi_hat_dM_fixed_with_phi_intergenic, exp_complete
rm(phi_dM_fixed_with_phi, phi_dM_fixed_intergenic,
   phi_dM_fixed_with_phi_intergenic,
   parameters_objects, plot_data, p1, p2, p3, acf_val, acf_df,
   codon_index, run_dirs)
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

# 8.3) Extracting selection estimates from the best model (dM-fixed-with_phi) ----

genome <- initializeGenomeObject(file = 'data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa',
                                 match.expression.by.id = TRUE,
                                 observed.expression.file = 'data/compiled_expression_IM767.txt') # Warnings are expected if genes are missing from expression file

parameter_object <- loadParameterObject(file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/R_objects/parameter.Rda")

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
rm(genome, parameter_object, selection_coeff,
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
selection_genes <- final_analysis_data |> dplyr::filter(Mean.log10.Phi >= 1.5)
cor_selection <- cor.test(selection_genes$Mean.log10.Phi, 
                          selection_genes$Intrinsic_Inefficiency, method = "spearman",
                          exact = F)

# Memory cleanup: section 8.3.1 plot objects ---
# Keeping: final_analysis_data, cor_load, cor_eff, cor_selection
rm(p_load, p_optim, selection_genes)
gc()

# 8.3.2) Analyzing the correlation between total selective pressure and CAI and CDC ----

integrated_data <- integrated_data |>
  dplyr::left_join(selection_metrics, by = "Gene_name") |>
  dplyr::filter(!is.na(L_ROC))

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

# 8.3.3) Final visualization ----

# 1. Data Preparation
plot_data <- integrated_data |>
  dplyr::mutate(
    # Expression metric
    Log_Phi = Mean.log10.Phi,
    
    # Calculate raw Phi to isolate Intrinsic Inefficiency
    # Using 10^Mean.log10.Phi approximates the scale
    Phi_raw = 10^Mean.log10.Phi,
    Intrinsic_Ineff = L_ROC / Phi_raw,
    
    # Log of L_ROC
    Log_L_ROC = log10(L_ROC + 1e-6), # Add small constant to avoid log(0)

    # Log Transform Length
    Log_Length = log10(Total_Codons)
  ) |>
  dplyr::filter(!is.na(ENC), !is.na(Total_Codons), !is.na(CAI), !is.na(Log_Phi))

# --- Print Stats for the Manuscript ---
cat("\n=== Selection Efficacy Stats (For Manuscript) ===\n")

# Define common color scale limits
phi_range <- range(plot_data$Log_Phi, na.rm = TRUE)

# Panel A: ROC_eff distribution (signed efficacy at 4-fold sites; linear scale)
y_max_anno <- 10000

p1 <- ggplot(plot_data |> dplyr::filter(!is.na(ROC_eff)), aes(x = ROC_eff)) +
  geom_histogram(bins = 100, fill = "#69b3a2", color = "white", linewidth = 0.05) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000),
                     labels = scales::comma_format(accuracy = 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off", ylim = c(0, NA)) +
  labs(x = expression(S[ROC] ~ "(signed AnaCoDa efficacy, 4-fold sites)"),
       y = "Gene Count (Log1p Scale)") +
  theme_custom() +
  theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))

# Panel B: L_ROC distribution
p2 <- ggplot(plot_data, aes(x = L_ROC)) +
  geom_histogram(bins = 100, fill = "#69b3a2", color = "white", linewidth = 0.05) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000), 
                     labels = scales::comma_format(accuracy = 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off", ylim = c(0, NA)) + 
  labs(x = expression(Log[10](L[ROC])), y = "Gene Count (Log1p Scale)") +
  theme_custom()

# Combine all 4 panels using patchwork
combined_plot <- p1 | p2 + plot_annotation(tag_levels = 'A')

# Saved with a larger dimension to accommodate the 2x2 grid nicely
ggsave("results/Selection_Landscape_Final.pdf", combined_plot, width = 14, height = 8)

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

msd_data <- integrated_data |>
  dplyr::select(Gene_name, ROC_eff, ROC_eff_4, L_ROC, Pi_mean_4fold, Mean_Log10_Exp,
                Max_Log10_Exp, Exp_breadth, CDS_length_nt, Sites_4fold,
                Pi_sum_4fold) |>
  dplyr::inner_join(gene_Q_4fold, by = "Gene_name") |>
  dplyr::filter(!is.na(Q_pref_base), !is.na(Pi_mean_4fold),
                !is.na(Mean_Log10_Exp), !is.na(Exp_breadth),
                N_4fold_sites >= 20)

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

sfs_G_file <- "./data/sfs_introns_G.csv"
sfs_C_file <- "./data/sfs_introns_C.csv"

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
  neutral_pool_pi <- neutral_pool |>
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

# Visualizations based on this data (q a two-allele system) ----

plot_data <- pi_data_operational |>
  # Remove rows with missing data
  dplyr::filter(!is.na(MeanLog10_exp), !is.na(q_pref)) |>
  
  dplyr::mutate(
    exp_bin = as.numeric(as.character(
      cut(MeanLog10_exp, 
          # Extend breaks to 2.4, then catch the rest with Inf
          breaks = c(seq(0, 2.4, by = 0.2), Inf), 
          # Shift labels to start at 0 (lower bounds)
          labels = seq(0, 2.4, by = 0.2),
          include.lowest = TRUE,
          right = FALSE) 
    ))
  ) |>
  
  dplyr::filter(!is.na(exp_bin)) |>
  
  dplyr::group_by(exp_bin) |>
  dplyr::summarise(
    size = n(),
    mean_q_pref = mean(q_pref),
    ci = 1.96 * (sd(q_pref) / sqrt(n())),
    .groups = 'drop'
  )

q_by_exp_plot <- ggplot(plot_data, aes(x = exp_bin, y = mean_q_pref)) +
  # Add error bars first so the points sit on top of them
  geom_errorbar(aes(ymin = mean_q_pref - ci, ymax = mean_q_pref + ci),
                width = 0.05, color = "#004b87") +
  
  # Add the points
  geom_point(color = "#004b87", size = 3) +
  
  # Format axes and labels to match the image
  scale_x_continuous(breaks = seq(0.0, 2.4, by = 0.2)) +
  labs(x = "exp bin", y = "q_pref") +
  
  # Apply your custom theme
  theme_custom()

ggsave("./results/Q_by_exp.pdf", plot = q_by_exp_plot, width = 8,
       height = 6)

# Pi 4-fold as a function of expression

pi_4_fold_by_exp <- pi_data_operational |>
  # Remove rows with missing data
  dplyr::filter(!is.na(MeanLog10_exp), !is.na(q_pref)) |>
  
  dplyr::mutate(
    exp_bin = as.numeric(as.character(
      cut(MeanLog10_exp, 
          # Extend breaks to 2.4, then catch the rest with Inf
          breaks = c(seq(0, 2.4, by = 0.2), Inf), 
          # Shift labels to start at 0 (lower bounds)
          labels = seq(0, 2.4, by = 0.2),
          include.lowest = TRUE,
          right = FALSE) 
    ))
  ) |>
  
  dplyr::filter(!is.na(exp_bin)) |>
  
  dplyr::group_by(exp_bin) |>
  dplyr::summarise(
    size = n(),
    mean_pi_4_fold = mean(pi_4fold),
    ci = 1.96 * (sd(pi_4fold) / sqrt(n())),
    .groups = 'drop'
  )

p_pi_by_expression <- ggplot(pi_4_fold_by_exp, 
                             aes(x = exp_bin, 
                                 y = mean_pi_4_fold)) +
  geom_point(size = 3, color = "#377EB8") +
  geom_errorbar(aes(ymin = mean_pi_4_fold - ci,
                    ymax = mean_pi_4_fold + ci),
                width = 0.1, color = "#377EB8") +
  labs(
    title = expression(paste("Mean 4-fold Nucleotide Diversity (", pi, 
                             ") by Expression Level")),
    x = "Expression level category",
    y = expression(paste("nuc_diversity (4 fold)"))
  ) +
  scale_x_continuous(breaks = seq(0.0, 2.4, by = 0.2)) +
  theme_custom()

ggsave("./results/pi_4fold_by_expression_asym.pdf", 
       p_pi_by_expression, width = 10, height = 6)

# Including non-synonimous pi values

pi_non_syn <- read.csv('data/NonSynPi.csv') |>
  dplyr::mutate(Gene_name = paste0("MgIM767.", Gene)) |>
  dplyr::rename(pi_non_syn = pi)

pi_data_operational <- pi_data_operational |>
  left_join(pi_non_syn |> dplyr::select(Gene_name, pi_non_syn))

pi_non_syn_by_exp <- pi_data_operational |>
  # Remove rows with missing data
  dplyr::filter(!is.na(MeanLog10_exp), !is.na(q_pref), !is.na(pi_non_syn)) |>
  
  dplyr::mutate(
    exp_bin = as.numeric(as.character(
      cut(MeanLog10_exp, 
          # Extend breaks to 2.4, then catch the rest with Inf
          breaks = c(seq(0, 2.4, by = 0.2), Inf), 
          # Shift labels to start at 0 (lower bounds)
          labels = seq(0, 2.4, by = 0.2),
          include.lowest = TRUE,
          right = FALSE) 
    ))
  ) |>
  
  dplyr::filter(!is.na(exp_bin)) |>
  
  dplyr::group_by(exp_bin) |>
  dplyr::summarise(
    size = n(),
    mean_pi_non_syn = mean(pi_non_syn),
    ci = 1.96 * (sd(pi_non_syn) / sqrt(n())),
    .groups = 'drop'
  )

p_pi_by_expression <- ggplot(pi_non_syn_by_exp, 
                             aes(x = exp_bin, 
                                 y = mean_pi_non_syn)) +
  geom_point(size = 3, color = "#377EB8") +
  geom_errorbar(aes(ymin = mean_pi_non_syn - ci,
                    ymax = mean_pi_non_syn + ci),
                width = 0.1, color = "#377EB8") +
  labs(
    title = expression(paste("Mean Non_synonimous Nucleotide Diversity (", pi, 
                             ") by Expression Level")),
    x = "Expression level category",
    y = expression(paste("nuc_diversity (non synonimous)"))
  ) +
  scale_x_continuous(breaks = seq(0.0, 2.4, by = 0.2)) +
  theme_custom()

ggsave("./results/pi_non_syn_by_expression_asym.pdf", 
       p_pi_by_expression, width = 10, height = 6)

# Combined pi plot ----

# 1. Create a unified summary dataframe
combined_pi_summary <- pi_data_operational |>
  dplyr::filter(!is.na(MeanLog10_exp), !is.na(pi_4fold), !is.na(pi_non_syn)) |>
  dplyr::mutate(
    exp_bin = as.numeric(as.character(
      cut(MeanLog10_exp, 
          breaks = c(seq(0, 2.4, by = 0.2), Inf), 
          labels = seq(0, 2.4, by = 0.2),
          include.lowest = TRUE, right = FALSE) 
    ))
  ) |>
  dplyr::filter(!is.na(exp_bin)) |>
  dplyr::group_by(exp_bin) |>
  dplyr::summarise(
    size = n(),
    mean_4fold = mean(pi_4fold),
    ci_4fold = 1.96 * (sd(pi_4fold) / sqrt(n())),
    mean_nonsyn = mean(pi_non_syn),
    ci_nonsyn = 1.96 * (sd(pi_non_syn) / sqrt(n())),
    .groups = 'drop'
  )

# Pivot data for faceting
long_pi_summary <- combined_pi_summary |>
  pivot_longer(
    cols = c(mean_4fold, mean_nonsyn),
    names_to = "metric",
    values_to = "mean_val"
  ) |>
  dplyr::mutate(
    # Match the correct CI to the correct metric
    ci = ifelse(metric == "mean_4fold", ci_4fold, ci_nonsyn),
    # Create clean labels for the graph
    metric_label = ifelse(metric == "mean_4fold", "4-fold Diversity", "Non-synonymous Diversity")
  )

# The Faceted Plot
p_faceted <- ggplot(long_pi_summary, aes(x = exp_bin, y = mean_val, color = metric_label)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_val - ci, ymax = mean_val + ci), width = 0.1) +
  facet_wrap(~ metric_label, ncol = 1, scales = "free_y") + # 'free_y' prevents compression
  scale_color_manual(values = c("4-fold Diversity" = "#377EB8", "Non-synonymous Diversity" = "#E41A1C")) +
  scale_x_continuous(breaks = seq(0.0, 2.4, by = 0.2)) +
  labs(
    title = expression(paste("Mean Nucleotide Diversity (", pi, ") by Expression Level")),
    x = "Expression level category",
    y = expression(paste("Nucleotide Diversity (", pi, ")"))
  ) +
  theme_custom() +
  theme(legend.position = "none") # Facet labels make the legend redundant

ggsave("./results/pi_combined_faceted.pdf", p_faceted, width = 10, height = 8)

# Operational thresholds for the two selection groups ----

# 1. Calculate dynamic transformation parameters based on your data ranges
min_4f <- min(combined_pi_summary$mean_4fold, na.rm = TRUE)
max_4f <- max(combined_pi_summary$mean_4fold, na.rm = TRUE)
min_ns <- min(combined_pi_summary$mean_nonsyn, na.rm = TRUE)
max_ns <- max(combined_pi_summary$mean_nonsyn, na.rm = TRUE)

# Calculate the slope (b_scale) and intercept (a_shift)
b_scale <- (max_4f - min_4f) / (max_ns - min_ns)
a_shift <- min_4f - b_scale * min_ns

# 2. Plot the aligned dual-axis graph
p_dual_aligned <- ggplot(combined_pi_summary, aes(x = exp_bin)) +
  
  # --- PRIMARY METRIC (4-fold) ---
  geom_point(aes(y = mean_4fold, color = "4-fold"), size = 3) +
  geom_errorbar(aes(y = mean_4fold, 
                    ymin = mean_4fold - ci_4fold, 
                    ymax = mean_4fold + ci_4fold, 
                    color = "4-fold"), width = 0.1) +
  
  # --- SECONDARY METRIC (Non-synonymous) ---
  # Apply both the shift and the scale so the trends map perfectly onto the primary axis space
  geom_point(aes(y = a_shift + b_scale * mean_nonsyn, color = "Non-synonymous"), size = 3) +
  geom_errorbar(aes(y = a_shift + b_scale * mean_nonsyn, 
                    ymin = a_shift + b_scale * (mean_nonsyn - ci_nonsyn), 
                    ymax = a_shift + b_scale * (mean_nonsyn + ci_nonsyn), 
                    color = "Non-synonymous"), width = 0.1) +
  
  # --- AXIS SCALING ---
  scale_y_continuous(
    name = expression(paste("4-fold Diversity (", pi, ")")),
    # Back-transform the secondary axis labels so they show the true non-synonymous values
    sec.axis = sec_axis(~ (. - a_shift) / b_scale, 
                        name = expression(paste("Non-synonymous Diversity (", pi, ")")))
  ) +
  scale_color_manual(name = "Metric", 
                     values = c("4-fold" = "#377EB8", "Non-synonymous" = "#E41A1C")) +
  scale_x_continuous(breaks = seq(0.0, 2.4, by = 0.2)) +
  labs(
    title = expression(paste("Mean Nucleotide Diversity (", pi, ") by Expression Level")),
    x = "Expression level category"
  ) +
  theme_custom() +
  theme(legend.position = "top")

ggsave("./results/pi_combined_dual_axis_aligned.pdf", p_dual_aligned, width = 10, height = 6)

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
S_BARRIER_source <- "fixed_2Ns_gt_1"

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
  msd_data <- msd_data |>
    dplyr::left_join(pi_data_operational |> dplyr::select(Gene_name, pi_2allele),
                     by = "Gene_name")
}


# Figure: Q vs expression with GAM fit ----
gam_q_exp <- mgcv::gam(
  Q_pref_base ~ s(Mean_Log10_Exp, k = 10),
  data    = msd_data |> dplyr::filter(!is.na(Q_pref_base), !is.na(Mean_Log10_Exp)),
  weights = N_4fold_sites
)
exp_grid_full <- seq(min(msd_data$Mean_Log10_Exp, na.rm = TRUE),
                     max(msd_data$Mean_Log10_Exp, na.rm = TRUE),
                     length.out = 300)
gam_full      <- predict(gam_q_exp,
                          newdata = data.frame(Mean_Log10_Exp = exp_grid_full),
                          type = "response", se.fit = TRUE)
gam_ribbon_df <- data.frame(
  x   = exp_grid_full,
  fit = gam_full$fit,
  lo  = gam_full$fit - 1.96 * gam_full$se.fit,
  hi  = gam_full$fit + 1.96 * gam_full$se.fit
)

p_q_expression <- ggplot(msd_data |> dplyr::filter(!is.na(Q_pref_base)),
                          aes(x = Mean_Log10_Exp, y = Q_pref_base)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "n") +
  geom_ribbon(data = gam_ribbon_df, aes(x = x, ymin = lo, ymax = hi),
              inherit.aes = FALSE, fill = "#E41A1C", alpha = 0.20) +
  geom_line(data = gam_ribbon_df, aes(x = x, y = fit),
            inherit.aes = FALSE, color = "#E41A1C", linewidth = 0.9) +
  geom_hline(yintercept = Q_neutral_obs, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Preferred-base frequency (Q) vs expression",
    subtitle = sprintf("Q_neutral = %.3f | S_BARRIER = %.4f (2N_e*s = 1 threshold)",
                       Q_neutral_obs, S_BARRIER),
    x = expression(Log[10](expression)),
    y = "Q (preferred-base freq., 4-fold sites)"
  ) +
  theme_custom()
ggsave("./results/Wright_Q_vs_expression_diagnostic.pdf",
       p_q_expression, width = 7, height = 5, device = cairo_pdf)

rm(gam_q_exp, exp_grid_full, gam_full, gam_ribbon_df, p_q_expression)

# Three-metric correlation summary ----
cat("\n=== Spearman correlations among selection metrics ===\n")
metric_df <- msd_data |>
  dplyr::select(ROC_eff, L_ROC, S_Wright_signed, Mean_Log10_Exp) |>
  dplyr::filter(dplyr::if_all(dplyr::everything(), is.finite))
cat(sprintf("  n = %d genes with all four metrics finite\n", nrow(metric_df)))
print(round(cor(metric_df, method = "spearman"), 3))
rm(metric_df)

# --- Binned ROC_eff and S_Wright tables ----------------------------------------
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

# Bin-level pi consistency test (parameter-free Wright validation)

# Parameter-free validation: does Wright's Q -> pi map predict pi_bin
# from the bin-level S_Wright_bin (already inverted from Q_bin)?
#
# Replacement is a parameter-free check: does Wright's Q -> pi map predict
# pi_bin from the bin-level S_Wright_bin (already inverted from Q_bin)?
#
#   pi_pred_wright = wright_pi(S_Wright_bin; U_emp, V_emp)
#   pi_residual    = pi_bin - pi_pred_wright
#   chi2_pi        = sum( pi_residual^2 / pi_se^2 )
#
# No fitted parameter, so df = n_bins.  Drift bins (negative S_Wright_bin)
# reflect residual selection in the empirical neutral pool that nudges
# (U_emp, V_emp) -- they surface as inflated pi_residual at low ROC_eff bins.
# We therefore also report the chi^2 restricted to the selection regime
# (S_Wright_bin >= S_BARRIER), the regime where the model is identifiable.

# Bin-level pi consistency using ROC_eff bins (parameter-free Wright validation)
bin_roc$pi_pred_wright <- wright_pi(bin_roc$S_Wright_bin,
                                    U = U_bin_calib, V = V_bin_calib)
bin_roc$pi_residual    <- bin_roc$pi_bin - bin_roc$pi_pred_wright

chi2_pi_terms <- (bin_roc$pi_residual)^2 /
                 pmax(bin_roc$pi_se, .Machine$double.eps)^2
chi2_pi_stat  <- sum(chi2_pi_terms, na.rm = TRUE)
chi2_pi_df    <- sum(!is.na(chi2_pi_terms))
chi2_pi_p     <- pchisq(chi2_pi_stat, df = chi2_pi_df, lower.tail = FALSE)
cat(sprintf(
  "\n[Validation] Bin-level pi consistency (ROC_eff bins): chi^2 = %.2f / df = %d -> p = %.3g\n",
  chi2_pi_stat, chi2_pi_df, chi2_pi_p
))

sel_bins <- bin_roc |>
  dplyr::filter(!is.na(S_Wright_bin), S_Wright_bin >= S_BARRIER)
chi2_pi_sel_stat <- if (nrow(sel_bins) > 0) {
  sum((sel_bins$pi_residual)^2 /
      pmax(sel_bins$pi_se, .Machine$double.eps)^2,
      na.rm = TRUE)
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

# Bin-level pi consistency figure (ROC_eff bins)
p_pi_validation <- ggplot(bin_roc, aes(x = mean_ROC_eff, y = pi_bin)) +
  geom_errorbar(aes(ymin = pi_bin - 2 * pi_se,
                    ymax = pi_bin + 2 * pi_se),
                width = 0, color = "#377EB8", alpha = 0.6) +
  geom_point(aes(size = sites_total),
             color = "#377EB8", alpha = 0.85) +
  geom_line(aes(y = pi_pred_wright),
            color = "#E41A1C", linewidth = 0.8) +
  geom_hline(yintercept = pi_neutral_obs,
             linetype = "dotted", color = "grey40") +
  scale_size_continuous(name = "4-fold sites / bin",
                        range = c(2, 6),
                        labels = scales::comma_format()) +
  labs(
    title    = "Bin-level pi consistency: observed pi vs Wright pi(S_Wright_bin)",
    subtitle = sprintf(
      "All-bin chi^2 = %.1f / df %d (p = %.2g) | sel-only chi^2 = %s / df %d (p = %s) | red: Wright pi(S_Wright_bin)",
      chi2_pi_stat, chi2_pi_df, chi2_pi_p,
      if (is.finite(chi2_pi_sel_stat)) sprintf("%.1f", chi2_pi_sel_stat) else "NA",
      chi2_pi_sel_df,
      if (is.finite(chi2_pi_sel_p))    sprintf("%.2g", chi2_pi_sel_p)    else "NA"
    ),
    x = expression(mean(S[ROC]) ~ "(per bin, 30 site-weighted bins)"),
    y = expression(pi[4*fold] ~ "(site-weighted)")
  ) +
  theme_custom()
ggsave("./results/Wright_pi_vs_S_Wright_binned.pdf",
       p_pi_validation, width = 7, height = 5, device = cairo_pdf)

# Per-gene ROC_eff_4 vs S_Wright sanity-check scatter ----
# CRITICAL: Filter to genes where BOTH ROC_eff_4 and S_Wright_raw are non-zero
# to avoid noise contamination from genes with no detectable selection signal.
# This avoids correlating noise with noise, preserving signal-to-noise in the
# subset of genes actually exhibiting codon usage optimization.

per_gene_pool <- msd_data |>
  dplyr::filter(!is.na(S_Wright_signed), !is.na(ROC_eff_4), 
                N_4fold_sites >= 50,        # Non-zero S_Wright (positive selection)
                ROC_eff_4 != 0)            # Non-zero ROC_eff (detectable codon preference)

cor_roc_eff_4_spearman <- cor(per_gene_pool$ROC_eff_4, per_gene_pool$S_Wright_signed,
                          method = "spearman")
cor_roc_eff_4_pearson  <- cor(per_gene_pool$ROC_eff_4, per_gene_pool$S_Wright_signed,
                          method = "pearson")

cat(sprintf(
  "[Validation] cor(ROC_eff_4, S_Wright_raw): Spearman = %+.3f, Pearson = %+.3f (n = %d)\n",
  cor_roc_eff_4_spearman, cor_roc_eff_4_pearson, nrow(per_gene_pool)
))
cat(sprintf("Filtered to: S_Wright_raw > 0 & ROC_eff_4 != 0 & N_4fold_sites >= 50\n"))

p_ROC_eff_4_vs_SWright <- ggplot(per_gene_pool,
                              aes(x = ROC_eff_4, y = S_Wright_signed)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Gene\nCount") +
  geom_smooth(method = "lm", color = "#E41A1C", linewidth = 0.9, se = TRUE) +
  geom_hline(yintercept = S_BARRIER, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0,         linetype = "dotted", color = "grey60") +
  labs(
    title = expression(paste(S[ROC][",4"], " vs ", S[Wright], "  (non-zero genes only)")),
    subtitle = sprintf(
      "Spearman rho = %+.3f (n = %d, S_Wright & ROC_eff) | S_BARRIER = %.4f",
      cor_roc_eff_4_spearman, nrow(per_gene_pool), S_BARRIER
    ),
    x = expression(S[ROC] ~ "(AnaCoDa efficacy, 4-fold)"),
    y = expression(S[Wright] ~ "(per-gene, non-negative)")
  ) +
  theme_custom()
ggsave("./results/Wright_ROC_eff_4_vs_S_Wright.pdf",
       p_ROC_eff_4_vs_SWright, width = 7, height = 5, device = cairo_pdf)

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

# 8.3.6) Diversity hump: pi vs S_Wright ----
#
# Wright's pi(S) is non-monotone: plateau at S=0, small rise (the "hump"),
# then collapse as selection drives preferred alleles toward fixation.
# 30 site-weighted bins of S_Wright with theory overlay (red) and S_BARRIER (dashed).

# Relative diversity: pi_bin / pi_neutral_obs, echoing McVean & Charlesworth (1999)
# who plot pi/pi(S=0) vs N_e*s.  Our calibrated U/V = (1-Q_neutral)/Q_neutral
# ≈ 2.9, placing Mimulus in the μ₁₀/μ₀₁ ≈ 3 regime (dotted curve in M&C Fig 1),
# where weak selection transiently raises diversity above neutral before the
# preferred allele is fixed at high S.

mu_ratio <- (1 - Q_neutral_obs) / Q_neutral_obs   # U/V = μ_away / μ_toward

# Use signed S_Wright for the diversity-hump but display it on a
# relative 0->1 axis for visual comparison while preserving negative values
# in the binning.  Compute a relative S for plotting: S_rel = (mean_S - min)/(max-min).
bin_sw_rel <- bin_sw |>
  dplyr::mutate(
    pi_moe = 1.96 * pi_se,
    pi_rel = pi_bin / pi_neutral_obs,
    pi_rel_lo = (pi_bin - 2 * pi_se) / pi_neutral_obs,
    pi_rel_hi = (pi_bin + 2 * pi_se) / pi_neutral_obs
  )

# Compute relative scaling based on the observed binned S range
min_S <- min(bin_sw_rel$mean_S_Wright, na.rm = TRUE)
max_S <- max(bin_sw_rel$mean_S_Wright, na.rm = TRUE)
if (!is.finite(min_S) || !is.finite(max_S) || max_S == min_S) {
  stop("Cannot compute relative S scale for diversity hump: invalid S range.")
}
bin_sw_rel <- bin_sw_rel |>
  dplyr::mutate(
    S_rel = (mean_S_Wright - min_S) / (max_S - min_S)
  )

# Theory curve evaluated on the signed S range then mapped to relative x
S_seq_hump <- seq(min_S, max_S, length.out = 500)
theory_sw  <- data.frame(
  x_rel  = (S_seq_hump - min_S) / (max_S - min_S),
  pi_rel = wright_pi(S_seq_hump, U = U_bin_calib, V = V_bin_calib) / pi_neutral_obs
)

# Compute S_BARRIER on the same relative scale (clamped to [0,1])
S_BARRIER_rel <- (S_BARRIER - min_S) / (max_S - min_S)
S_BARRIER_rel <- pmin(pmax(S_BARRIER_rel, 0), 1)

p_hump_sw <- ggplot(bin_sw_rel, aes(x = mean_S_Wright, y = pi_bin)) +
  geom_vline(xintercept = S_BARRIER,
             linetype = "dashed", color = "#377EB8") +
  geom_errorbar(aes(x = mean_S_Wright, ymin = pi_bin - pi_moe, ymax = pi_bin + pi_moe),
                width = 0.1, color = "#377EB8", alpha = 0.6) +
  geom_point(aes(x = mean_S_Wright, size = sites_total),
             color = "#377EB8", alpha = 0.9) +
  scale_size_continuous(name = "4-fold sites / bin",
                        range = c(2, 6),
                        labels = scales::comma_format()) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    x = "Mean(S[Wright])",
    y = expression(pi[4*fold])
  ) +
  theme_custom()

ggsave("./results/Diversity_hump_S_Wright.pdf",
       p_hump_sw, width = 7, height = 5.5, device = cairo_pdf)

write.csv(bin_roc |> dplyr::select(ROC_eff_bin, n_genes, mean_ROC_eff, sites_total,
                                   Q_bin, S_Wright_bin, pi_bin, pi_se,
                                   pi_pred_wright, pi_residual),
          "./results/Diversity_hump_pi_vs_ROC_eff_binned.csv",
          row.names = FALSE)
write.csv(bin_sw_rel |> dplyr::select(Ssw_bin, n_genes, mean_S_Wright, sites_total,
                                      Q_bin, pi_bin, pi_se, pi_rel),
          "./results/Diversity_hump_pi_vs_S_Wright_binned.csv",
          row.names = FALSE)

rm(p_hump_sw, theory_sw, bin_sw_rel, S_seq_hump, mu_ratio)
gc()

# Persist artefacts.
write.csv(wright_fid,
          "./results/Wright_curve_advisor_fiducial.csv", row.names = FALSE)
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
   p_advisor_Q, p_advisor_pi,
   S_grid_fid, S_grid_emp,
   neutral_pool, neutral_pool_pi, pi_data_operational,
   Q_neutral_two, pi_neutral_two, hardy_max_two,
   chi2_pi_terms, sel_bins,
   per_gene_pool, p_pi_validation, p_ROC_eff_4_vs_SWright)
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

# Compatibility alias for downstream sections expecting `detailed_annotation`
detailed_annotation <- top_ROC_eff |>
  dplyr::transmute(locusName = Gene_name, Best.hit.arabi.name,
                   Best.hit.arabi.defline)

# 8.6) Goodness-of-fit test based on Anacoda predictions ----

gof_results <- run_gof_analysis(
  mutation_file  = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Mutation.csv",
  selection_file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv",
  phi_file       = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/gene_expression.txt",
  codon_counts_long = codon_freq_long,
  test           = "chisq",
  min_aa_total   = 5,
  output_prefix  = "./results/ROC_model_goodness_of_fit"
)

# 8.7) Independent multinomial validation of ROC preferred codons ----

# Validates ROC preferred codons via two complementary approaches:
#
#   1) Targeted binomial GLM: For each AA, tests whether the ROC-predicted
#      preferred codon's usage increases significantly with expression.
#      A positive expression coefficient confirms the ROC prediction.
#
#   2) Full multinomial models: Fits per-AA multinomial regressions and
#      checks whether the ROC preferred codon shows positive Delta_P,
#      plus dEta vs Delta_P correlation across all codons.
#
# These tests do NOT assume the ROC mutational background derived from
# introns — they ask whether the ROC model's predictions are borne out
# by the expression-frequency relationship in the data.

message("\n=== Independent Multinomial Validation of ROC Preferred Codons ===\n")

# Prepare expression and length data
expr_for_multinom <- exp_complete |>
  dplyr::mutate(Exp_log10 = Max_Log10_Exp) |>
  dplyr::select(Gene_name, Exp_log10) |>
  dplyr::rename(Gene = Gene_name)

length_for_multinom <- integrated_data |>
  dplyr::select(Gene_name, CDS_length_nt) |>
  dplyr::rename(Gene = Gene_name) |>
  as.data.frame()

# Load CSP parameters for quantitative dEta correlation
csp_for_validation <- load_csp_parameters(
  mutation_file  = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Mutation.csv",
  selection_file = "./results/MCMC_results/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv"
)

# Run validation
multinomial_validation <- run_multinomial_validation(
  codon_counts_long = codon_freq_long,
  expression_df     = expr_for_multinom,
  length_df         = length_for_multinom,
  roc_preferred_df  = preferred_codons_roc,
  csp_df            = csp_for_validation,
  min_aa_total      = 5,
  output_prefix     = "./results/ROC_multinomial_validation"
)

# Memory cleanup: multinomial validation intermediates ---
# Keeping: multinomial_validation
rm(expr_for_multinom, length_for_multinom, csp_for_validation)
gc()

# Memory cleanup: large codon frequency data ---
# Keeping: gof_results, multinomial_validation
rm(codon_freq_long)
gc()

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

# Get all sense codons from global genetic_code_dna_long (excluding stops)
all_codons_rna <- gsub("T", "U", names(genetic_code_dna_long)[!genetic_code_dna_long %in% c("STOP", "Trp", "Met")])

# Initialize matrix
species <- c("Arabidopsis_thaliana", "Populus_trichocarpa", "Mimulus_guttatus", "Physcomitrella_patens")
codon_matrix <- matrix(0, nrow = length(species), ncol = length(all_codons_rna))
rownames(codon_matrix) <- species
colnames(codon_matrix) <- all_codons_rna

# Fill matrix
for (sp_idx in 1:length(species)) {
  sp_name <- species[sp_idx]
  if (sp_name == "Mimulus_guttatus") {
    preferred <- preferred_codons_mg$Codon_RNA
  } else {
    preferred <- plant_codons_extended[[sp_name]]
  }
  
  # Mark preferred codons as 1
  for (codon in preferred) {
    # Handle multiple codons separated by /
    codons_split <- unlist(strsplit(codon, "/"))
    for (c in codons_split) {
      if (c %in% all_codons_rna) {
        codon_matrix[sp_name, c] <- 1
      }
    }
  }
}

# Build similarity matrix
n_species <- length(species)
similarity_matrix <- matrix(0, nrow = n_species, ncol = n_species)
rownames(similarity_matrix) <- species
colnames(similarity_matrix) <- species

for (i in 1:n_species) {
  for (j in 1:n_species) {
    similarity_matrix[i, j] <- jaccard_similarity(codon_matrix[i, ], codon_matrix[j, ])
  }
}

# Convert to distance matrix
distance_matrix <- as.dist(1 - similarity_matrix)

# Hierarchical clustering
hc <- hclust(distance_matrix, method = "average")

# Save dendrogram
pdf("./results/plant_codon_preference_dendrogram.pdf", width = 10, height = 7)
par(mar = c(5, 4, 4, 2))
plot(hc, main = "Plant Species Clustering by Codon Preference Similarity",
     xlab = "Species", ylab = "Distance (1 - Jaccard Similarity)",
     sub = paste("Based on preferred codon usage in", length(all_codons_rna), "sense codons"),
     cex.main = 1.3)
dev.off()

# Create unrooted phylogram using ape package
tree <- as.phylo(hc)

# Save unrooted tree
pdf("./results/plant_codon_preference_unrooted.pdf", width = 10, height = 10)
par(mar = c(1, 1, 3, 1))
plot(tree, type = "unrooted", main = "Unrooted Tree: Codon Preference Similarity",
     cex = 1.2, lab4ut = "axial", edge.width = 2)
dev.off()

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

# Memory cleanup: Section 9 comparative analysis intermediates ---
# Keeping: eta_data, similarity_matrix, hc, tree
rm(preferred_codons_comparative, preferred_codons_mg, mg_prefs,
   plant_codons_extended, codon_matrix,
   distance_matrix, plot_data, p_comparison,
   all_codons_rna, species)
gc()

## *****************************************************************************
## 10) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________

# 10.1) CA Analysis -

codon_usage_m <- as.matrix(codon_usage[, -1])
rownames(codon_usage_m) <- codon_usage[[1]]
colnames(codon_usage_m) <- names(codon_usage)[-1]

codon_usage_CA <- CA(X = codon_usage_m, graph = FALSE)

# Extract CA coordinates and merge with expression data
codon_usage_CA_coord <- as.data.frame(codon_usage_CA$row$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", rownames(codon_usage_CA$row$coord)))

names(codon_usage_CA_coord)[names(codon_usage_CA_coord) %in% c("Dim 1", "Dim 2", "Dim 3")] <- 
  c("Dim.1", "Dim.2", "Dim.3")

codon_usage_CA_coord <- integrated_data |>
  dplyr::left_join(codon_usage_CA_coord, by = "Gene_name") |>
  dplyr::mutate(Expression_Group = as.character(Expression_Group))

# Filter to extreme expression groups for MANOVA
codon_usage_CA_coord_extremes <- codon_usage_CA_coord |>
  dplyr::filter(Expression_Group %in% c("Top 5%", "Bottom 5%"))

# Prepare gene data for biplot
gene_data_ca <- codon_usage_CA_coord |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\nCA Analysis: Preferred vs Non-preferred Codons ---\n")
p_ca <- create_preference_biplot(
  ordination_result = codon_usage_CA,
  gene_data = gene_data_ca,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 3),
  arrow_scale = 1.0,
  title = NULL,
  subtitle = NULL,
  output_file = "./results/CA_preference_biplot.pdf"
)

# Analyze loading direction
ca_loading_test <- analyze_codon_loading_direction(
  ordination_result = codon_usage_CA,
  preferred_codons = preferred_codons_roc,
  dim = 1
)

# MANOVA test for CA dimension separation
cat("\n=== MANOVA Test: CA Dimensions by Expression Group ===\n")
ca_manova <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ Expression_Group, 
                    data = codon_usage_CA_coord_extremes)
print(summary(ca_manova))

# Univariate Wilcoxon tests with FDR correction
cat("\n=== Univariate Tests for CA Dimensions by Expression (FDR corrected) ===\n")
ca_wilcox_results <- data.frame(
  Dimension = character(),
  W = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (dim_name in c("Dim.1", "Dim.2", "Dim.3")) {
  if (dim_name %in% names(codon_usage_CA_coord_extremes)) {
    wtest <- wilcox.test(as.formula(paste(dim_name, "~ Expression_Group")), 
                         data = codon_usage_CA_coord_extremes)
    ca_wilcox_results <- rbind(ca_wilcox_results, 
                                data.frame(Dimension = dim_name,
                                          W = wtest$statistic,
                                          p_value = wtest$p.value))
  }
}

# Apply FDR correction
ca_wilcox_results$p_adj <- p.adjust(ca_wilcox_results$p_value, method = "fdr")
ca_wilcox_results$Significance <- ifelse(ca_wilcox_results$p_adj < 0.05, "***", "")

for (i in seq_len(nrow(ca_wilcox_results))) {
  cat(sprintf("%s: W = %.2f, p = %.4f, p_adj = %.4f %s\n", 
              ca_wilcox_results$Dimension[i], 
              ca_wilcox_results$W[i], 
              ca_wilcox_results$p_value[i],
              ca_wilcox_results$p_adj[i],
              ca_wilcox_results$Significance[i]))
}

# 10.2) PCA Analysis ----

rscu_m <- as.matrix(cub_results$rscu_results[, -1])
rownames(rscu_m) <- cub_results$rscu_results[[1]]
colnames(rscu_m) <- names(cub_results$rscu_results)[-1]

rscu_PCA <- PCA(rscu_m, graph = FALSE)

# Extract PCA coordinates and merge with expression data
rscu_PCA_coord <- as.data.frame(rscu_PCA$ind$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", rownames(rscu_PCA$ind$coord)))

rscu_PCA_coord <- integrated_data |>
  dplyr::left_join(rscu_PCA_coord, by = "Gene_name") |>
  dplyr::mutate(Expression_Group = as.character(Expression_Group))

# Filter to extreme expression groups for MANOVA
rscu_PCA_coord_extremes <- rscu_PCA_coord |>
  dplyr::filter(Expression_Group %in% c("Top 5%", "Bottom 5%"))

# Prepare gene data for biplot
gene_data_pca <- rscu_PCA_coord |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\nPCA Analysis: Preferred vs Non-preferred Codons ---\n")
p_pca <- create_preference_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 2),
  arrow_scale = 1.5,
  title = NULL,
  subtitle = NULL,
  output_file = "./results/PCA_preference_biplot.pdf"
)

# Analyze loading direction
pca_loading_test <- analyze_codon_loading_direction(
  ordination_result = rscu_PCA,
  preferred_codons = preferred_codons_roc,
  dim = 1
)

# MANOVA test for PCA dimension separation
cat("\n=== MANOVA Test: PCA Dimensions by Expression Group ===\n")
pca_manova <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ Expression_Group, 
                     data = rscu_PCA_coord_extremes)
print(summary(pca_manova))

# Univariate Wilcoxon tests with FDR correction
cat("\n=== Univariate Tests for PCA Dimensions by Expression (FDR corrected) ===\n")
pca_wilcox_results <- data.frame(
  Dimension = character(),
  W = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (dim_name in c("Dim.1", "Dim.2", "Dim.3")) {
  wtest <- wilcox.test(as.formula(paste(dim_name, "~ Expression_Group")), 
                       data = rscu_PCA_coord_extremes)
  pca_wilcox_results <- rbind(pca_wilcox_results, 
                               data.frame(Dimension = dim_name,
                                         W = wtest$statistic,
                                         p_value = wtest$p.value))
}

# Apply FDR correction
pca_wilcox_results$p_adj <- p.adjust(pca_wilcox_results$p_value, method = "fdr")
pca_wilcox_results$Significance <- ifelse(pca_wilcox_results$p_adj < 0.05, "***", "")

for (i in seq_len(nrow(pca_wilcox_results))) {
  cat(sprintf("%s: W = %.2f, p = %.4f, p_adj = %.4f %s\n", 
              pca_wilcox_results$Dimension[i], 
              pca_wilcox_results$W[i], 
              pca_wilcox_results$p_value[i],
              pca_wilcox_results$p_adj[i],
              pca_wilcox_results$Significance[i]))
}

# 10.3) Selection (L_ROC) Based Analysis ----

# Test whether genes under strong vs weak selection show distinct codon patterns

if (exists("selection_metrics") && "L_ROC" %in% names(selection_metrics)) {
  
  cat("\nCA/PCA Analysis by Selection Load ---\n")
  
  # Create S_load-based groups
  s_quantiles <- quantile(selection_metrics$L_ROC, 
                          probs = c(0.05, 0.95), na.rm = TRUE)
  
  selection_groups <- selection_metrics |>
    dplyr::mutate(
      S_Group = dplyr::case_when(
        L_ROC >= s_quantiles[2] ~ "High Selection (Top 5%)",
        L_ROC <= s_quantiles[1] ~ "Low Selection (Bottom 5%)",
        TRUE ~ "Intermediate"
      )
    ) |>
    dplyr::select(Gene_name, L_ROC, S_Group)
  
  cat("\nSelection Load Group Distribution:\n")
  print(table(selection_groups$S_Group))
  
  # CA analysis by selection load
  codon_usage_CA_coord_S <- codon_usage_CA_coord |>
    dplyr::left_join(selection_groups, by = "Gene_name") |>
    dplyr::filter(S_Group %in% c("High Selection (Top 5%)", "Low Selection (Bottom 5%)"))
  
  gene_data_ca_S <- codon_usage_CA_coord_S |>
    dplyr::select(Gene_name, expression_group = S_Group)
  
  p_ca_S <- create_preference_biplot(
    ordination_result = codon_usage_CA,
    gene_data = gene_data_ca_S,
    preferred_codons = preferred_codons_roc,
    dims = c(1, 2),
    arrow_scale = 1.0,
    title = "CA Biplot: By Selection Load",
    subtitle = "High vs Low S_load genes",
    output_file = "./results/CA_preference_biplot_S_load.pdf"
  )
  
  # PCA analysis by selection load
  rscu_PCA_coord_S <- rscu_PCA_coord |>
    dplyr::left_join(selection_groups, by = "Gene_name") |>
    dplyr::filter(S_Group %in% c("High Selection (Top 5%)", "Low Selection (Bottom 5%)"))
  
  gene_data_pca_S <- rscu_PCA_coord_S |>
    dplyr::select(Gene_name, expression_group = S_Group)
  
  p_pca_S <- create_preference_biplot(
    ordination_result = rscu_PCA,
    gene_data = gene_data_pca_S,
    preferred_codons = preferred_codons_roc,
    dims = c(1, 2),
    arrow_scale = 1.5,
    title = "PCA Biplot: By Selection Load",
    subtitle = "High vs Low S_load genes",
    output_file = "./results/PCA_preference_biplot_S_load.pdf"
  )
  
  # MANOVA tests for S_load grouping
  cat("\n=== MANOVA Test: CA Dimensions by Selection Load ===\n")
  ca_manova_S <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ S_Group, 
                        data = codon_usage_CA_coord_S)
  print(summary(ca_manova_S))
  
  cat("\n=== MANOVA Test: PCA Dimensions by Selection Load ===\n")
  pca_manova_S <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ S_Group, 
                         data = rscu_PCA_coord_S)
  print(summary(pca_manova_S))
  
  # Univariate Wilcoxon tests with FDR correction for CA by S_load
  cat("\n=== Univariate Tests for CA Dimensions by S_load (FDR corrected) ===\n")
  ca_wilcox_S_results <- data.frame(
    Dimension = character(),
    W = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (dim_name in c("Dim.1", "Dim.2", "Dim.3")) {
    if (dim_name %in% names(codon_usage_CA_coord_S)) {
      wtest <- wilcox.test(as.formula(paste(dim_name, "~ S_Group")), 
                           data = codon_usage_CA_coord_S)
      ca_wilcox_S_results <- rbind(ca_wilcox_S_results, 
                                    data.frame(Dimension = dim_name,
                                              W = wtest$statistic,
                                              p_value = wtest$p.value))
    }
  }
  
  ca_wilcox_S_results$p_adj <- p.adjust(ca_wilcox_S_results$p_value, method = "fdr")
  ca_wilcox_S_results$Significance <- ifelse(ca_wilcox_S_results$p_adj < 0.05, "***", "")
  
  for (i in seq_len(nrow(ca_wilcox_S_results))) {
    cat(sprintf("%s: W = %.2f, p = %.4f, p_adj = %.4f %s\n", 
                ca_wilcox_S_results$Dimension[i], 
                ca_wilcox_S_results$W[i], 
                ca_wilcox_S_results$p_value[i],
                ca_wilcox_S_results$p_adj[i],
                ca_wilcox_S_results$Significance[i]))
  }
  
  # Univariate Wilcoxon tests with FDR correction for PCA by S_load
  cat("\n=== Univariate Tests for PCA Dimensions by S_load (FDR corrected) ===\n")
  pca_wilcox_S_results <- data.frame(
    Dimension = character(),
    W = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (dim_name in c("Dim.1", "Dim.2", "Dim.3")) {
    wtest <- wilcox.test(as.formula(paste(dim_name, "~ S_Group")), 
                         data = rscu_PCA_coord_S)
    pca_wilcox_S_results <- rbind(pca_wilcox_S_results, 
                                   data.frame(Dimension = dim_name,
                                             W = wtest$statistic,
                                             p_value = wtest$p.value))
  }
  
  pca_wilcox_S_results$p_adj <- p.adjust(pca_wilcox_S_results$p_value, method = "fdr")
  pca_wilcox_S_results$Significance <- ifelse(pca_wilcox_S_results$p_adj < 0.05, "***", "")
  
  for (i in seq_len(nrow(pca_wilcox_S_results))) {
    cat(sprintf("%s: W = %.2f, p = %.4f, p_adj = %.4f %s\n", 
                pca_wilcox_S_results$Dimension[i], 
                pca_wilcox_S_results$W[i], 
                pca_wilcox_S_results$p_value[i],
                pca_wilcox_S_results$p_adj[i],
                pca_wilcox_S_results$Significance[i]))
  }
  
  cat("\n✓ Selection load-based analysis complete\n")
  
} else {
  cat("\n⚠ selection_coeff_intensity not found - skipping S_load-based analysis\n")
}

cat("\n✓ CA/PCA ordination analysis complete\n")
cat("  Key outputs: CA_preference_biplot.pdf, PCA_preference_biplot.pdf\n\n")

# Memory cleanup: Section 10 large matrices and plot objects ---
# Keeping: cub_results, codon_usage_CA, rscu_PCA, ca_loading_test, pca_loading_test,
#          ca_manova, pca_manova, ca_wilcox_results, pca_wilcox_results
rm(codon_usage_m, codon_usage_CA_coord,
   rscu_m, rscu_PCA_coord,
   gene_data_ca, gene_data_pca, p_ca, p_pca)
if (exists("selection_groups")) rm(selection_groups, codon_usage_CA_coord_S,
                                   rscu_PCA_coord_S, gene_data_ca_S, gene_data_pca_S,
                                   p_ca_S, p_pca_S)
gc()

## *****************************************************************************
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

# Analysis 2: Top expression tier (isolating selection signal)

cat("\n=== Analysis 2: Top 5% Expressed Genes - tRNA Correlation ===\n")
cat("Rationale: If selection shapes codon usage to match tRNA supply,\n")
cat("the signal should be STRONGEST in highly expressed genes.\n\n")

# Subset codon usage to top 5% expressed genes
top_expressed_genes <- integrated_data |>
  dplyr::filter(Expression_Group == "Top 5%") |>
  dplyr::pull(Gene_name)

codon_usage_top <- codon_usage |>
  dplyr::filter(Gene_name %in% top_expressed_genes)

cat(sprintf("Top expression tier: %d genes (5%% cutoff = %.2f log10 CPM)\n",
            nrow(codon_usage_top), top_5_cutoff))

tRNA_top_results <- tRNA_codon_correlation(
  codon_counts = codon_usage_top,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_top_expressed",
  test_method = "spearman",
  mode = "by.copy.number",
  wobble_mode = "conservative",
  is_genome_wide = FALSE
)

# Background comparison: rest of genes
rest_expressed_genes <- integrated_data |>
  dplyr::filter(Expression_Group != "Top 5%") |>
  dplyr::pull(Gene_name)

codon_usage_rest <- codon_usage |>
  dplyr::filter(Gene_name %in% rest_expressed_genes)

tRNA_rest_results <- tRNA_codon_correlation(
  codon_counts = codon_usage_rest,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_rest_genes",
  test_method = "spearman",
  mode = "by.copy.number",
  wobble_mode = "conservative",
  is_genome_wide = FALSE
)

# Comparison: Top 5% vs Rest within-family correlations

cat("\n=== Expression Tier Comparison: Within-Family tRNA-Codon Correlations ===\n")
cat("Testing whether top-expressed genes show stronger tRNA co-adaptation.\n\n")

top_cors <- sapply(tRNA_top_results$correlation_results$per_amino_acid,
                   function(x) x$estimate)
rest_cors <- sapply(tRNA_rest_results$correlation_results$per_amino_acid,
                    function(x) x$estimate)

common_aas <- intersect(names(top_cors), names(rest_cors))

if (length(common_aas) >= 3) {
  tier_comparison <- data.frame(
    AA = common_aas,
    Top5_r = top_cors[common_aas],
    Rest_r = rest_cors[common_aas],
    Delta_r = top_cors[common_aas] - rest_cors[common_aas],
    row.names = NULL
  )
  tier_comparison <- tier_comparison[order(-tier_comparison$Delta_r), ]

  cat("Within-family correlation comparison (sorted by Delta):\n")
  print(tier_comparison, row.names = FALSE)

  # Paired Wilcoxon: are top-gene correlations systematically stronger?
  paired_test <- wilcox.test(tier_comparison$Top5_r, tier_comparison$Rest_r,
                             paired = TRUE, alternative = "greater")
  cat(sprintf("\nPaired Wilcoxon test (Top 5%% > Rest): V = %.0f, p = %.4f\n",
              paired_test$statistic, paired_test$p.value))
  cat(sprintf("Mean Delta_r (top - rest): %.3f\n", mean(tier_comparison$Delta_r, na.rm = TRUE)))

  # Visualization: comparison plot
  p_tier_compare <- ggplot(tier_comparison, aes(x = Rest_r, y = Top5_r)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(aes(label = AA), size = 3, max.overlaps = 20) +
    labs(
      title = "tRNA-Codon Correlation: Top 5% vs Rest (Within-Family)",
      subtitle = sprintf(
        "Paired Wilcoxon p = %.4f | Mean Delta_r = %.3f\nAbove diagonal = stronger correlation in top-expressed genes",
        paired_test$p.value, mean(tier_comparison$Delta_r, na.rm = TRUE)
      ),
      x = "Within-family Spearman r (Rest 95%)",
      y = "Within-family Spearman r (Top 5%)"
    ) +
    theme_custom() +
    coord_equal()

  ggsave("./results/tRNA_correlation_tier_comparison.pdf",
         p_tier_compare, width = 8, height = 8)

  # Save comparison data
  write.csv(tier_comparison,
            "./results/tRNA_correlation_tier_comparison.csv",
            row.names = FALSE)
  cat("Saved: ./results/tRNA_correlation_tier_comparison.pdf\n")
  cat("Saved: ./results/tRNA_correlation_tier_comparison.csv\n")
} else {
  cat("Not enough common amino acid families for tier comparison.\n")
}

# Translational Accuracy Hypothesis

# Load tRNA data for the pairing analysis
tRNA_data <- fread("./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt")

# Get codon supply from the copy number analysis (already calculated)
codon_supply <- tRNA_copynumber_results$codon_supply

# Create complete codon status table using ROC-derived preferred codons
# First, get all sense codons from genetic code
all_sense_codons <- data.frame(
  Codon = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
            "TAT", "TAC", "TGT", "TGC", "TGG",
            "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
            "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
            "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
            "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
            "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"),
  Amino_Acid = c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser",
                 "Tyr", "Tyr", "Cys", "Cys", "Trp",
                 "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro",
                 "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
                 "Ile", "Ile", "Ile", "Met", "Thr", "Thr", "Thr", "Thr",
                 "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg",
                 "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala",
                 "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly")
)

# Mark preferred codons from ROC model (those in preferred_codons_roc)
roc_preferred <- preferred_codons_roc$Preferred_Codons
all_sense_codons$Status <- ifelse(all_sense_codons$Codon %in% roc_preferred, 
                                   "Preferred", "Non-Preferred")

# Convert to RNA format for pairing function
roc_codon_status <- all_sense_codons
roc_codon_status$Codon <- gsub("T", "U", roc_codon_status$Codon)

n_preferred <- sum(roc_codon_status$Status == "Preferred")
n_total <- nrow(roc_codon_status)

cat("Using ROC model (AnaCoDa) to classify all codons:\n")
cat("  Preferred:", n_preferred, "codons (lowest eta = highest fitness)\n")
cat("  Non-Preferred:", n_total - n_preferred, "codons\n")
cat("  Total:", n_total, "sense codons\n\n")

# Run the translational accuracy test
pairing_analysis <- classify_codon_anticodon_pairing(
  tRNA_data = tRNA_data,
  codon_supply = codon_supply,
  preferred_codons = roc_codon_status[, c("Codon", "Amino_Acid", "Status")],
  output_dir = "./results/tRNA_analysis_pairing",
  save_results = TRUE
)

# Memory cleanup: Section 11 tRNA intermediate objects ---
# Keeping: tRNA_copynumber_results, tRNA_top_results, tRNA_rest_results, 
#          aa_trna_check, pairing_analysis, preferred_codons_roc, tier_comparison
rm(tRNA_data, codon_supply,
   all_sense_codons, roc_codon_status, roc_preferred,
   n_preferred, n_total, top_expressed_genes, rest_expressed_genes,
   codon_usage_top, codon_usage_rest)
if (exists("top_cors")) rm(top_cors, rest_cors, common_aas)
if (exists("p_tier_compare")) rm(p_tier_compare, paired_test)
gc()

## *****************************************************************************
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

cat("\n=== Contour Plot: Preferred Codon Frequency ~ Expression x Gene Length ===\n")

# --- 12.4a: Overall (all preferred codons pooled) ---

contour_data <- integrated_data |>
  dplyr::filter(!is.na(Mean_preferred_freq),
                !is.na(Max_Log10_Exp),
                Total_Codons > 0) |>
  dplyr::mutate(log10_length = log10(Total_Codons))

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

pred_grid_pref$Predicted_broad <- predict(contour_gam, newdata = pred_grid_pref_broad,
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

# 12.3b: Split by C-ending vs G-ending preferred codons ----
# Re-read the per-position file to compute separate per-gene frequencies

preferred_raw <- read.delim("./data/all_chromosomes.codon_frequencies_preferred.txt",
                            stringsAsFactors = FALSE) |>
  dplyr::mutate(Gene_name = paste0("MgIM767.", Gene))

# Classify preferred codons by their 3rd nucleotide
preferred_raw$Ending <- substr(preferred_raw$Preferred_Codon,
                               nchar(preferred_raw$Preferred_Codon),
                               nchar(preferred_raw$Preferred_Codon))

# Skip stop codons and methionine/tryptophan (single-codon families)
preferred_raw <- preferred_raw |>
  dplyr::filter(!AA %in% c("*", "Stop", "M", "W"))

# Per-gene average preferred freq split by ending nucleotide
pref_by_ending <- preferred_raw |>
  dplyr::group_by(Gene_name, Ending) |>
  dplyr::summarize(
    Mean_pref_freq = mean(Preferred_Freq, na.rm = TRUE),
    n_positions = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = Ending,
    values_from = c(Mean_pref_freq, n_positions),
    names_sep = "_"
  )

contour_split <- integrated_data |>
  dplyr::left_join(pref_by_ending, by = "Gene_name") |>
  dplyr::filter(!is.na(Max_Log10_Exp), Total_Codons > 0) |>
  dplyr::mutate(log10_length = log10(Total_Codons))

# Build contour for each ending type
for (ending in c("C", "G")) {
  col_name <- paste0("Mean_pref_freq_", ending)
  if (!col_name %in% names(contour_split)) next

  sub <- contour_split |>
    dplyr::filter(!is.na(.data[[col_name]]),
                  .data[[col_name]] > 0,
                  .data[[col_name]] < 1)

  if (nrow(sub) < 200) {
    message(sprintf("  %s-ending: skipped (%d genes)", ending, nrow(sub)))
    next
  }

  gam_split <- mgcv::gam(
    as.formula(paste0(col_name, " ~ te(Max_Log10_Exp, Exp_breadth, log10_length, k = c(10, 10))")),
    data = sub,
    family = betar(link = "logit")
  )
  cat(sprintf("GAM R-sq(adj) [%s-ending]: %.3f\n", ending, summary(gam_split)$r.sq))

  grid_split <- expand.grid(
    Max_Log10_Exp = seq(min(sub$Max_Log10_Exp, na.rm = TRUE),
                        max(sub$Max_Log10_Exp, na.rm = TRUE),
                        length.out = 200),
    log10_length  = seq(min(sub$log10_length, na.rm = TRUE),
                        max(sub$log10_length, na.rm = TRUE),
                        length.out = 200),
    Exp_breadth = 33 # Holding this constant (broad genes)
  )
  grid_split$Predicted <- predict(gam_split, newdata = grid_split, type = "response")

  p_split <- ggplot(grid_split, aes(x = Max_Log10_Exp, y = log10_length)) +
    geom_raster(aes(fill = Predicted), interpolate = TRUE) +
    geom_contour(aes(z = Predicted), colour = "grey30",
                 linewidth = 0.4, bins = 12) +
    scale_fill_gradientn(
      colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                  "#FEE8C8", "#FDBB84", "#E34A33"),
      name = paste0("Freq.\n", ending, "-end")
    ) +
    labs(
      title = sprintf("Freq. of ROC-Preferred %s-Ending Codons (All AA Families)", ending),
      subtitle = "GAM-predicted surface across expression and gene length",
      x = expression(log[10](Max~Expression~CPM)),
      y = expression(log[10](Gene~Length~"(codons)"))
    ) +
    theme_custom() +
    theme(legend.position = "right", panel.grid = element_blank())

  outfile <- sprintf("./results/Preferred_freq_contour_%s_ending.pdf", ending)
  ggsave(outfile, p_split, width = 9, height = 7)
  cat(sprintf("Saved: %s\n", outfile))
}

rm(preferred_raw, pref_by_ending, contour_split)

# 12.4) Background selection test: nonsynonymous pi vs expression ----
# Under background selection, genes under stronger purifying selection
# (= more highly expressed) should have lower pi at nonsynonymous (0-fold)
# sites. This is independent of the codon usage bias signal.
#
# We decompose 4-fold pi into three components:
#   (a) Overall 4-fold pi (linked selection signal)
#   (b) 4-fold pi at the first 300 bp of CDS (ramp / Hill-Robertson zone)
#   (c) 4-fold pi after the first 300 bp (gene body)
# The "selection group" (L_ROC > thr_sel) is isolated as an independent point.

cat("\n=== Background Selection: Nonsynonymous Pi vs Expression ===\n")

bgs_data <- integrated_data |>
  dplyr::filter(Sites_0fold >= 10, Pi_mean_0fold >= 0,
                !is.na(Mean_Log10_Exp)) |>
  dplyr::mutate(log10_length = log10(CDS_length_nt))

# Spearman correlation: pi_0fold vs expression
cor_0fold <- cor.test(bgs_data$Mean_Log10_Exp, bgs_data$Pi_mean_0fold,
                      method = "spearman", exact = FALSE)
cat(sprintf("Spearman rho (Pi_0fold ~ Expression): %.4f (p = %.2e, n = %d)\n",
            cor_0fold$estimate, cor_0fold$p.value, nrow(bgs_data)))

# --- Binned analysis: expression-ranked bins (excluding selection group) ---
bgs_nonsel <- bgs_data |> dplyr::filter(L_ROC < thr_sel)

bgs_binned <- bgs_nonsel |>
  dplyr::arrange(Mean_Log10_Exp) |>
  dplyr::mutate(Rank = dplyr::row_number(),
                Exp_Bin = ceiling(Rank / 1000)) |>
  dplyr::group_by(Exp_Bin) |>
  dplyr::summarize(
    n_genes           = dplyr::n(),
    mean_expression   = mean(Mean_Log10_Exp, na.rm = TRUE),
    # Overall 4-fold and 0-fold
    pi_4fold          = sum(Pi_sum_4fold, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
    pi_0fold          = sum(Pi_sum_0fold, na.rm = TRUE) / sum(Sites_0fold, na.rm = TRUE),
    # 4-fold decomposed by position
    pi_4fold_first300 = sum(Pi_sum_4fold_first300, na.rm = TRUE) /
                        sum(Sites_4fold_first300, na.rm = TRUE),
    pi_4fold_after300 = sum(Pi_sum_4fold_after300, na.rm = TRUE) /
                        sum(Sites_4fold_after300, na.rm = TRUE),
    .groups = "drop"
  )

# --- Selection group (L_ROC > thr_sel) as isolated point ---
bgs_sel <- bgs_data |> dplyr::filter(L_ROC > thr_sel)

sel_point <- bgs_sel |>
  dplyr::summarize(
    n_genes           = dplyr::n(),
    mean_expression   = mean(Mean_Log10_Exp, na.rm = TRUE),
    pi_4fold          = sum(Pi_sum_4fold, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
    pi_0fold          = sum(Pi_sum_0fold, na.rm = TRUE) / sum(Sites_0fold, na.rm = TRUE),
    pi_4fold_first300 = sum(Pi_sum_4fold_first300, na.rm = TRUE) /
                        sum(Sites_4fold_first300, na.rm = TRUE),
    pi_4fold_after300 = sum(Pi_sum_4fold_after300, na.rm = TRUE) /
                        sum(Sites_4fold_after300, na.rm = TRUE)
  )

cat(sprintf("Selection group (L_ROC > thr_sel = %.3f): n = %d, mean_exp = %.2f\n",
            as.numeric(thr_sel),
            sel_point$n_genes, sel_point$mean_expression))

# Compute pi_0fold / pi_4fold ratio
bgs_binned$Pi_ratio_0_4 <- bgs_binned$pi_0fold / bgs_binned$pi_4fold

# --- Pivot to long format for the combined plot ---
bgs_long <- bgs_binned |>
  tidyr::pivot_longer(
    cols = c(pi_4fold, pi_0fold, pi_4fold_first300, pi_4fold_after300),
    names_to = "Site_Class", values_to = "Pi"
  ) |>
  dplyr::mutate(
    Site_Class = dplyr::recode(Site_Class,
      "pi_0fold"          = "0-fold (nonsynonymous)",
      "pi_4fold"          = "4-fold (overall)",
      "pi_4fold_first300" = "4-fold (first 300 bp)",
      "pi_4fold_after300" = "4-fold (after 300 bp)"
    ),
    Source = "Expression bins"
  )

sel_long <- sel_point |>
  tidyr::pivot_longer(
    cols = c(pi_4fold, pi_0fold, pi_4fold_first300, pi_4fold_after300),
    names_to = "Site_Class", values_to = "Pi"
  ) |>
  dplyr::mutate(
    Site_Class = dplyr::recode(Site_Class,
      "pi_0fold"          = "0-fold (nonsynonymous)",
      "pi_4fold"          = "4-fold (overall)",
      "pi_4fold_first300" = "4-fold (first 300 bp)",
      "pi_4fold_after300" = "4-fold (after 300 bp)"
    ),
    Source = "Selection group (S > thr_sel)"
  )

bgs_colors <- c("0-fold (nonsynonymous)" = "#E41A1C",
                "4-fold (overall)"        = "#377EB8",
                "4-fold (first 300 bp)"   = "#984EA3",
                "4-fold (after 300 bp)"   = "#4DAF4A")

p_bgs <- ggplot(bgs_long,
                aes(x = mean_expression, y = Pi, color = Site_Class)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.6) +
  # Selection group as distinct larger points
  geom_point(data = sel_long,
             aes(x = mean_expression, y = Pi, color = Site_Class),
             shape = 18, size = 5, stroke = 1.2) +
  scale_color_manual(values = bgs_colors) +
  labs(
    title = expression("Linked Selection: " * pi * " at 0-fold and 4-fold Sites vs Expression"),
    subtitle = "4-fold decomposed: first 300 bp (ramp/HRI zone) vs gene body | Diamonds = selection group (S > thr_sel)",
    x = expression(Mean~log[10](Expression)),
    y = expression("Weighted " * pi),
    color = NULL
  ) +
  theme_custom() +
  theme(legend.position = "top")

ggsave("./results/Pi_0fold_4fold_vs_expression.pdf", p_bgs, width = 11, height = 7)
cat("Saved: ./results/Pi_0fold_4fold_vs_expression.pdf\n")

# --- 12.4b) C-frequency contribution to pi increase in selection group ---
# The rise in 4-fold pi for the selection category is expected if selection
# increases C frequency at 4-fold sites (ROC-preferred codons are GC-ending),
# creating more C/non-C heterozygosity. We show this via mutation-type
# decomposition: C-involving pairs (AC, CG, CT) vs non-C pairs (AG, AT, GT).

if (has_mutation_types) {
  # Mutation-type pi components for the binned data
  bgs_mut_binned <- bgs_nonsel |>
    dplyr::arrange(Mean_Log10_Exp) |>
    dplyr::mutate(Rank = dplyr::row_number(),
                  Exp_Bin = ceiling(Rank / 1000)) |>
    dplyr::group_by(Exp_Bin) |>
    dplyr::summarize(
      mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
      pi_AC = sum(Pi_sum_4fold_AC, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AG = sum(Pi_sum_4fold_AG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AT = sum(Pi_sum_4fold_AT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CG = sum(Pi_sum_4fold_CG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CT = sum(Pi_sum_4fold_CT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_GT = sum(Pi_sum_4fold_GT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(Source = "Expression bins")

  # Selection group mutation-type decomposition
  sel_mut <- bgs_sel |>
    dplyr::summarize(
      mean_expression = mean(Mean_Log10_Exp, na.rm = TRUE),
      pi_AC = sum(Pi_sum_4fold_AC, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AG = sum(Pi_sum_4fold_AG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_AT = sum(Pi_sum_4fold_AT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CG = sum(Pi_sum_4fold_CG, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_CT = sum(Pi_sum_4fold_CT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE),
      pi_GT = sum(Pi_sum_4fold_GT, na.rm = TRUE) / sum(Sites_4fold, na.rm = TRUE)
    ) |>
    dplyr::mutate(Source = "Selection group (S > thr_sel)")

  bgs_mut_all <- dplyr::bind_rows(bgs_mut_binned, sel_mut)

  # Tag C-involving vs non-C pairs
  bgs_mut_long <- bgs_mut_all |>
    tidyr::pivot_longer(
      cols = starts_with("pi_"),
      names_to = "Mutation_Type", values_to = "Pi_component",
      names_prefix = "pi_"
    ) |>
    dplyr::mutate(
      Involves_C = ifelse(grepl("C", Mutation_Type), "C-involving", "Non-C"),
      Is_Selection = (Source == "Selection group (S > thr_sel)")
    )

  # Stacked area by C-involvement
  bgs_c_summary <- bgs_mut_long |>
    dplyr::group_by(mean_expression, Source, Involves_C) |>
    dplyr::summarize(Pi_total = sum(Pi_component), .groups = "drop")

  # Split for bins vs selection group
  bgs_c_bins <- bgs_c_summary |> dplyr::filter(Source == "Expression bins")
  bgs_c_sel  <- bgs_c_summary |> dplyr::filter(Source == "Selection group (S > thr_sel)")

  p_c_contrib <- ggplot(bgs_c_bins,
                        aes(x = mean_expression, y = Pi_total, fill = Involves_C)) +
    geom_area(alpha = 0.7, position = "stack") +
    # Selection group as stacked bars
    geom_col(data = bgs_c_sel,
             aes(x = mean_expression, y = Pi_total, fill = Involves_C),
             width = 0.08, alpha = 0.9, position = "stack") +
    scale_fill_manual(values = c("C-involving" = "#377EB8", "Non-C" = "#BDBDBD")) +
    labs(
      title = expression("Decomposition of 4-fold " * pi *
                          ": C-involving vs Non-C Segregating Pairs"),
      subtitle = "Bar at right = selection group (S > thr_sel); C-involvement drives pi increase under selection",
      x = expression(Mean~log[10](Expression)),
      y = expression(pi * " component (4-fold)"),
      fill = NULL
    ) +
    theme_custom() +
    theme(legend.position = "top")

  ggsave("./results/Pi_4fold_C_contribution_vs_expression.pdf",
         p_c_contrib, width = 10, height = 6)
  cat("Saved: ./results/Pi_4fold_C_contribution_vs_expression.pdf\n")

  # Per mutation-type line plot with selection group highlighted
  p_mut_sel <- ggplot(bgs_mut_long |> dplyr::filter(!Is_Selection),
                      aes(x = mean_expression, y = Pi_component,
                          color = Mutation_Type)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 0.6) +
    # Selection group as diamonds
    geom_point(data = bgs_mut_long |> dplyr::filter(Is_Selection),
               aes(x = mean_expression, y = Pi_component, color = Mutation_Type),
               shape = 18, size = 4) +
    scale_color_brewer(palette = "Set2", name = "Segregating\nBases") +
    labs(
      title = expression("4-fold " * pi * " by Mutation Type: Expression Bins + Selection Group"),
      subtitle = "Diamonds = selection group (S > thr_sel); note elevated AC, CG, CT (C-involving pairs)",
      x = expression(Mean~log[10](Expression)),
      y = expression(pi * " component (4-fold)")
    ) +
    theme_custom() +
    theme(legend.position = "right")

  ggsave("./results/Pi_4fold_mutation_type_with_selection_group.pdf",
         p_mut_sel, width = 11, height = 6)
  cat("Saved: ./results/Pi_4fold_mutation_type_with_selection_group.pdf\n")

  rm(bgs_mut_binned, sel_mut, bgs_mut_all, bgs_mut_long,
     bgs_c_summary, bgs_c_bins, bgs_c_sel, p_c_contrib, p_mut_sel)
}

# GAM contour: pi at 0-fold sites ~ expression x gene length
bgs_contour_data <- bgs_data |>
  dplyr::filter(Pi_mean_0fold > 0)

bgs_gam <- mgcv::gam(
  Pi_mean_0fold ~ te(Max_Log10_Exp, Exp_breadth, log10_length,
                     k = c(10, 10)),
  data = bgs_contour_data,
  family = Gamma(link = "log")
)
cat(sprintf("GAM R-sq(adj) [Pi_0fold surface]: %.3f\n", summary(bgs_gam)$r.sq))

bgs_grid <- expand.grid(
  Max_Log10_Exp = seq(min(bgs_contour_data$Max_Log10_Exp, na.rm = TRUE),
                      max(bgs_contour_data$Max_Log10_Exp, na.rm = TRUE),
                      length.out = 200),
  log10_length  = seq(min(bgs_contour_data$log10_length, na.rm = TRUE),
                      max(bgs_contour_data$log10_length, na.rm = TRUE),
                      length.out = 200),
  Exp_breadth = 33 # Fixed for broadly expressed genes (bigger category)
)

bgs_grid$Predicted <- predict(bgs_gam, newdata = bgs_grid, type = "response")

p_bgs_contour <- ggplot(bgs_grid, aes(x = Max_Log10_Exp, y = log10_length)) +
  geom_raster(aes(fill = Predicted), interpolate = TRUE) +
  geom_contour(aes(z = Predicted), colour = "grey30",
               linewidth = 0.4, bins = 12) +
  scale_fill_gradientn(
    colours = rev(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                    "#FEE8C8", "#FDBB84", "#E34A33")),
    name = expression(pi["0-fold"])
  ) +
  labs(
    title = expression("Nonsynonymous Diversity (" * pi["0-fold"] * ") Surface"),
    subtitle = "Background selection: pi decreases with expression (purifying selection)",
    x = expression(log[10](Max~Expression~CPM)),
    y = expression(log[10](Gene~Length~"(nt)"))
  ) +
  theme_custom() +
  theme(legend.position = "right", panel.grid = element_blank())

ggsave("./results/Pi_0fold_contour_exp_x_length.pdf",
       p_bgs_contour, width = 9, height = 7)
cat("Saved: ./results/Pi_0fold_contour_exp_x_length.pdf\n")

# Memory cleanup: Section 12 plot objects and temporary subsets ---
# Keeping: pi_nonlinearity_results, pi_models, pi_selection,
#          pi_nonlinearity_results_s, pi_models_s, pi_selection_s,
#          preferred_nonlinearity_results, preferred_models, preferred_selection,
#          confounding_results, kw_preferred_freq
rm(pi_selection_winner, pi_secection_winner_s,
   preferred_selection_winner, confounder_model_gam,
   predictors, predictors_s, predictors_p,
   top5_preferred, middle_preferred, bottom5_preferred,
   p_boxplot_preferred, p_preferred_median, p_surface_pref, plot_data_pref,
   contour_data, contour_gam, pred_grid_pref, p_pref_contour,
   bgs_data, bgs_binned, bgs_long, bgs_contour_data, bgs_nonsel, bgs_sel,
   bgs_gam, bgs_grid, p_bgs, p_ratio, p_bgs_contour, cor_0fold,
   sel_point, sel_long, bgs_colors, has_mutation_types)
gc()

## *****************************************************************************
## 13) Nucleotide composition at 4-fold degenerate sites ----
## _____________________________________________________________________________
#
# This analysis provides a SFS-independent test for selection on codon usage.
# If selection favors GC-ending codons (as indicated by ROC model), then GC
# content at 4-fold degenerate sites should INCREASE with expression level,
# opposing the AT-rich mutational bias.
#
# 4-fold degenerate sites: all four nucleotides at the 3rd position are
# synonymous (Ala, Gly, Pro, Thr, Val, Leu_4, Ser_4, Arg_4).

cat("\n=== Nucleotide Composition at 4-Fold Degenerate Sites ===\n")

# Identify 4-fold degenerate amino acid families
fourfold_families <- c("Ala", "Gly", "Pro", "Thr", "Val", "Leu_4", "Ser_4", "Arg_4")
fourfold_codons <- names(genetic_code_dna_long)[genetic_code_dna_long %in% fourfold_families]

cat(sprintf("4-fold degenerate families: %s\n", paste(fourfold_families, collapse = ", ")))
cat(sprintf("Codons included: %d\n", length(fourfold_codons)))

# Calculate per-gene nucleotide composition at 4-fold sites
fourfold_cols <- intersect(fourfold_codons, names(codon_usage))
fourfold_usage <- codon_usage[, c("Gene_name", fourfold_cols), with = FALSE]

fourfold_long <- fourfold_usage |>
  tidyr::pivot_longer(-Gene_name, names_to = "Codon", values_to = "Count") |>
  dplyr::mutate(Third_base = substr(Codon, 3, 3))

fourfold_per_gene <- fourfold_long |>
  dplyr::group_by(Gene_name) |>
  dplyr::summarise(
    Total_4fold = sum(Count),
    N_A = sum(Count[Third_base == "A"]),
    N_T = sum(Count[Third_base == "T"]),
    N_C = sum(Count[Third_base == "C"]),
    N_G = sum(Count[Third_base == "G"]),
    GC3_4fold = (N_G + N_C) / Total_4fold,
    Freq_A = N_A / Total_4fold,
    Freq_T = N_T / Total_4fold,
    Freq_C = N_C / Total_4fold,
    Freq_G = N_G / Total_4fold,
    .groups = "drop"
  )

# Merge with expression data
fourfold_with_expr <- integrated_data |>
  dplyr::select(Gene_name, Max_Log10_Exp, Mean_Log10_Exp,
                Expression_Group,
                Geom_Mean_CPM, CDS_length_nt, L_ROC) |>
  dplyr::inner_join(fourfold_per_gene, by = "Gene_name") |>
  dplyr::filter(Total_4fold >= 20)  # Minimum for reliable composition estimates

cat(sprintf("Genes with >= 20 4-fold codons: %d / %d (%.1f%%)\n",
            nrow(fourfold_with_expr), nrow(integrated_data),
            100 * nrow(fourfold_with_expr) / nrow(integrated_data)))

# Summary by expression group
fourfold_summary <- fourfold_with_expr |>
  dplyr::group_by(Expression_Group) |>
  dplyr::summarise(
    N = dplyr::n(),
    Mean_GC3 = mean(GC3_4fold, na.rm = TRUE),
    SD_GC3 = sd(GC3_4fold, na.rm = TRUE),
    Median_GC3 = median(GC3_4fold, na.rm = TRUE),
    Mean_A = mean(Freq_A, na.rm = TRUE),
    Mean_T = mean(Freq_T, na.rm = TRUE),
    Mean_C = mean(Freq_C, na.rm = TRUE),
    Mean_G = mean(Freq_G, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nGC3 at 4-fold degenerate sites by expression group:\n")
print(as.data.frame(fourfold_summary))

# Per-nucleotide composition vs expression
nuc_long <- fourfold_with_expr |>
  tidyr::pivot_longer(
    cols = c(Freq_A, Freq_T, Freq_C, Freq_G),
    names_to = "Nucleotide",
    values_to = "Frequency"
  ) |>
  dplyr::mutate(
    Nucleotide = gsub("Freq_", "", Nucleotide),
    Bias_direction = ifelse(Nucleotide %in% c("G", "C"),
                            "GC (ROC-preferred)", "AT (mutation-favored)")
  )

p_nuc_4fold <- ggplot(nuc_long, aes(x = Mean_Log10_Exp, y = Frequency,
                                     color = Nucleotide)) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = c("A" = "#E41A1C", "T" = "#FF7F00",
                                "C" = "#377EB8", "G" = "#4DAF4A")) +
  geom_vline(xintercept = 1.8, colour="black", linetype = "longdash") +
  labs(
    title = "Nucleotide Composition at 4-fold Sites vs Expression",
    subtitle = "Opposing trends (GC up, AT down) indicate selection opposing mutational bias",
    x = "Mean Expression (Log10 CPM)",
    y = "Nucleotide Frequency at 4-fold Sites"
  ) +
  theme_custom() +
  theme(legend.position = "right")

ggsave("./results/nucleotide_comp_4fold_vs_expression.pdf",
       p_nuc_4fold, width = 10, height = 7)

# Save 4-fold composition data
write.csv(fourfold_with_expr,
          "./results/fourfold_degenerate_composition.csv",
          row.names = FALSE)

# Create the binary group based on the inflection point
fourfold_with_expr <- fourfold_with_expr |>
  dplyr::mutate(
    Inflection_Group = ifelse(Mean_Log10_Exp >= 1.8, "High (>= 2.0)", "Low/Med (< 2.0)")
  )

#Check the summary statistics of C frequency across these two groups
inflection_summary <- fourfold_with_expr |>
  dplyr::group_by(Inflection_Group) |>
  dplyr::summarise(
    N = dplyr::n(),
    Mean_Freq_C = mean(Freq_C, na.rm = TRUE),
    Median_Freq_C = median(Freq_C, na.rm = TRUE),
    SD_Freq_C = sd(Freq_C, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Summary of Freq_C by 2.8 Inflection Point ===\n")
print(as.data.frame(inflection_summary))

# Perform the Wilcoxon rank-sum test (equivalent to Kruskal-Wallis for 2 groups)
wilcox_c_4fold <- wilcox.test(Freq_C ~ Inflection_Group, 
                              data = fourfold_with_expr)

cat(sprintf("\nWilcoxon test (Freq_C by Inflection): W = %.2f, p = %.2e\n",
            wilcox_c_4fold$statistic, wilcox_c_4fold$p.value))

# Cleanup Section 13b intermediates
rm(fourfold_long, fourfold_per_gene, fourfold_summary,
   fourfold_cols, fourfold_codons, fourfold_families,
   nuc_long, kw_gc3_4fold, m_gc3_4fold, p_gc3_4fold, p_nuc_4fold)
# Keep fourfold_with_expr for potential downstream use
gc()

## *****************************************************************************
## 14) Diversity across different genomic compartment ----
## _____________________________________________________________________________

pi_compartment <- read.table(file = "./data/all_chromosomes.pi_by_compartment.txt",
                             header = T)

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
# Keeping: overall_pi_stats, t_test_result, boost_comparison, hump_test_data
rm(pi_compartment, compartment_order, pi_compart,
   hump_per_chrom, hump_summary, p_hump,
   paired_test_data, plot_data, additivity_check)
gc()

# ******************************************************************************
# 15) Testing the translational ramp hypothesis ----
# ______________________________________________________________________________

# 15.1) Reference based analysis (historic) ----
binary_preferred <- codons_to_preferred_state_bernoulli(trans, 
                                                        as.character(preferred_codons$Codon))

binary_preferred <- binary_preferred |>
  left_join(integrated_data |> dplyr::select(Gene_name, Max_Log10_Exp, 
                                             Exp_breadth), 
            by = join_by("GeneID" == "Gene_name")) |>
  na.exclude()

clean_data <- binary_preferred |>
  dplyr::filter(!is.na(Max_Log10_Exp), !is.na(Exp_breadth)) %>%
  # Focus on the Translational Ramp region (first 200 codons)
  dplyr::filter(Position <= 400) |> 
  dplyr::mutate(
    # Standardize predictors (Mean=0, SD=1) for faster MCMC
    Exp_Z = as.numeric(scale(Max_Log10_Exp)),
    Breadth_Z = as.numeric(scale(Exp_breadth)),
    GeneID = factor(GeneID),
    Is_Preferred = as.integer(Is_Preferred)
  )

# Subsampling genes to make MCMC manageable

set.seed(1998)
target_genes <- sample(unique(clean_data$GeneID), 3000)
model_data <- clean_data |> dplyr::filter(GeneID %in% target_genes)

model_data$GeneID <- droplevels(model_data$GeneID)

# Aggregate data
window_size <- 5
model_data_agg <- model_data |>
  dplyr::mutate(Window = ceiling(Position / window_size)) |>
  dplyr::group_by(GeneID, Window, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    Position_mid = mean(Position),
    n_preferred = sum(Is_Preferred),
    n_total = n(),
    .groups = "drop"
  ) |>
  dplyr::filter(n_total > 0, Position_mid <= 200)

# Model 1: Null (no ramp) ----

fit_null_ml <- bam(
  cbind(n_preferred, n_total - n_preferred) ~ 
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(GeneID, bs = "re"),  # Random effect
  
  data = model_data_agg,
  family = binomial(link = "logit"),
  method = "fREML",  # Fast REML
  discrete = TRUE,   # Faster for large datasets
  nthreads = 1      # Parallel
)

# Model 2: Global Ramp (Position Effect Only) ----

fit_ramp_ml <- bam(
  cbind(n_preferred, n_total - n_preferred) ~ 
    s(Position_mid, k = 10, bs = "tp") +  # The ramp
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(GeneID, bs = "re"),
  
  data = model_data_agg,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

# Model 3: Ramp x Expression interaction ----

fit_ramp_int_ml <- bam(
  cbind(n_preferred, n_total - n_preferred) ~ 
    s(Position_mid, k = 10, bs = "tp") +
    s(Position_mid, by = Exp_Z, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(GeneID, bs = "re"),
  
  data = model_data_agg,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

# Model comparison (Likelihood ratio test) ----

# Test 1: Does position matter? (Ramp vs Null)
anova(fit_null_ml, fit_ramp_ml, test = "Chisq")

# Test 2: Does ramp vary by expression?
anova(fit_ramp_ml, fit_ramp_int_ml, test = "Chisq")

# AIC comparison
AIC(fit_null_ml, fit_ramp_ml, fit_ramp_int_ml)

# Visualization

pred_positions <- data.frame(
  Position_mid = seq(5, 200, by = 2),
  Exp_Z = 0,
  Breadth_Z = 0,
  GeneID = model_data_agg$GeneID[1],  # Reference gene for RE
  n_total = 5
)

pred_ramp <- predict(fit_ramp_ml, newdata = pred_positions, 
                     type = "link", se.fit = TRUE, exclude = "s(GeneID)")

pred_positions$fit <- plogis(pred_ramp$fit)
pred_positions$lower <- plogis(pred_ramp$fit - 1.96 * pred_ramp$se.fit)
pred_positions$upper <- plogis(pred_ramp$fit + 1.96 * pred_ramp$se.fit)

plot_ramp_ml <- ggplot(pred_positions, aes(x = Position_mid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = fit), color = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Translational Ramp: ML Estimate",
    subtitle = "Does preferred codon usage increase with position?",
    x = "Codon Position",
    y = "P(Preferred Codon)",
    caption = "Ribbon = ±1.96 SE (approx 95% CI)"
  ) +
  theme_custom()

ggsave("./results/translational_ramp_ml.pdf", plot_ramp_ml, width = 10, 
       height = 6)

pred_grid_exp <- expand.grid(
  Position_mid = seq(5, 200, by = 5),
  Exp_Z = c(-1.5, 0, 1.5),
  Breadth_Z = 0,
  n_total = 5
) |>
  dplyr::mutate(
    GeneID = model_data_agg$GeneID[1],
    Exp_level = factor(
      Exp_Z,
      levels = c(-1.5, 0, 1.5),
      labels = c("Low", "Medium", "High")
    )
  )

pred_exp <- predict(fit_ramp_int_ml, newdata = pred_grid_exp, 
                    type = "link", se.fit = TRUE, 
                    exclude = "s(GeneID)",
                    unconditional = TRUE)

pred_grid_exp$fit <- plogis(pred_exp$fit)
pred_grid_exp$lower <- plogis(pred_exp$fit - 1.96 * pred_exp$se.fit)
pred_grid_exp$upper <- plogis(pred_exp$fit + 1.96 * pred_exp$se.fit)

plot_ramp_exp_ml <- ggplot(pred_grid_exp, 
                           aes(x = Position_mid, color = Exp_level, fill = Exp_level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Ramp Shape by Expression Level (ML)",
    x = "Codon Position",
    y = "P(Preferred Codon)",
    color = "Expression",
    fill = "Expression"
  ) +
  theme_custom() +
  theme(legend.position = "top")

ggsave("./results/ramp_by_expression_ml.pdf", plot_ramp_exp_ml, width = 10, 
       height = 6)

# POST-HOC: Contrast Test (Positions 1-50 vs 51-200) ----

contrast_data <- model_data_agg |>
  dplyr::mutate(Region = ifelse(Position_mid <= 50, "Ramp", "Body")) |>
  dplyr::group_by(GeneID, Region, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    n_preferred = sum(n_preferred),
    n_total = sum(n_total),
    .groups = "drop"
  )

fit_contrast_ml <- bam(
  cbind(n_preferred, n_total - n_preferred) ~ 
    Region + Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(GeneID, bs = "re"),
  
  data = contrast_data,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

summary(fit_contrast_ml)

# Extract coefficient for Region
region_coef <- coef(fit_contrast_ml)["RegionRamp"]
region_se <- summary(fit_contrast_ml)$se[names(coef(fit_contrast_ml)) == "RegionRamp"]

cat("\n=== Ramp vs Body Contrast ===\n")
cat(sprintf("Log-odds difference (Ramp - Body): %.3f ± %.3f\n", 
            region_coef, region_se))
cat(sprintf("Z-score: %.2f\n", region_coef / region_se))
cat(sprintf("P-value: %.2e\n", 2 * pnorm(-abs(region_coef / region_se))))

# Interpret on probability scale
baseline <- coef(fit_contrast_ml)["(Intercept)"]
prob_body <- plogis(baseline)
prob_ramp <- plogis(baseline + region_coef)

cat(sprintf("\nP(Preferred | Body): %.3f\n", prob_body))
cat(sprintf("P(Preferred | Ramp): %.3f\n", prob_ramp))
cat(sprintf("Difference: %.3f (%.1f%% change)\n", 
            prob_ramp - prob_body, 
            100 * (prob_ramp - prob_body) / prob_body))

# Memory cleanup: Section 15.1 large binary matrix and plot objects ---
# Keeping: fit_null_ml, fit_ramp_ml, fit_ramp_int_ml, fit_contrast_ml
rm(binary_preferred, clean_data, model_data, model_data_agg, target_genes,
   contrast_data, pred_positions, pred_ramp,
   plot_ramp_ml, pred_grid_exp, pred_exp, plot_ramp_exp_ml,
   region_coef, region_se, baseline, prob_body, prob_ramp)
gc()

# 15.2) Polymorphism based (contemporaneous) ----

poly_data <- fread(
  "data/all_chromosomes.codon_frequencies_preferred.txt",
  select = c("Gene", "Codon_Pos", "Preferred_Freq", "Non_Preferred_Freq"),
  showProgress = TRUE
)

poly_data <- poly_data |>
  dplyr::mutate(
    Gene_clean = paste0("MgIM767.", Gene)  # Add prefix
  ) |>
  dplyr::rename(Position = Codon_Pos) |>
  dplyr::filter(Position <= 200)  # Focus on first 200 codons

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
            nrow(poly_with_exp), 
            length(unique(poly_with_exp$Gene_clean))))

# Subsample genes to make GAM with gene random effects tractable ----
# The full dataset (~20K genes) makes s(Gene_clean, bs="re") extremely slow.
# We randomly sample 3,000 genes (seed = 1998 for reproducibility) —
# sufficient for stable estimation of population-level fixed + smooth effects.

n_subsample_genes <- 3000
all_genes_15 <- unique(poly_with_exp$Gene_clean)
cat(sprintf("Subsampling %d / %d genes for GAM fitting (seed = 1998)...\n",
            n_subsample_genes, length(all_genes_15)))

set.seed(1998)
sampled_genes_15 <- sample(all_genes_15, size = min(n_subsample_genes, length(all_genes_15)))
poly_with_exp <- poly_with_exp |>
  dplyr::filter(Gene_clean %in% sampled_genes_15) |>
  dplyr::mutate(Gene_clean = droplevels(Gene_clean))

cat(sprintf("After subsampling: %d codon positions from %d genes\n",
            nrow(poly_with_exp),
            length(unique(poly_with_exp$Gene_clean))))

# Aggregate using same window_size

poly_agg <- poly_with_exp |>
  dplyr::mutate(Window = ceiling(Position / window_size)) |>
  dplyr::group_by(Gene_clean, Window, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    Position_mid = mean(Position),
    # Weighted average of preferred frequency
    Preferred_Freq_mean = mean(Preferred_Freq, na.rm = TRUE),
    n_codons = n(),
    .groups = "drop"
  ) |>
  dplyr::filter(Position_mid <= 200) |>
  dplyr::mutate(
    # Simple boundary adjustment - just avoid exact 0 and 1
    Preferred_Freq_beta = case_when(
      Preferred_Freq_mean <= 0.001 ~ 0.001,
      Preferred_Freq_mean >= 0.999 ~ 0.999,
      TRUE ~ Preferred_Freq_mean
    )
  )

# Verify no exact 0s or 1s
if(any(poly_agg$Preferred_Freq_beta <= 0 | poly_agg$Preferred_Freq_beta >= 1)) {
  stop("ERROR: Beta regression requires values strictly between 0 and 1")
}

# MODEL 1: Null (no position effect) ----
fit_null_poly <- bam(
  Preferred_Freq_beta ~ 
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  
  data = poly_agg,
  family = betar(),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

# MODEL 2: Global ramp ----
fit_ramp_poly <- bam(
  Preferred_Freq_beta ~ 
    s(Position_mid, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  
  data = poly_agg,
  family = betar(),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

# MODEL 3: Ramp × Expression interaction ----
fit_ramp_int_poly <- bam(
  Preferred_Freq_beta ~ 
    s(Position_mid, k = 10, bs = "tp") +
    s(Position_mid, by = Exp_Z, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  
  data = poly_agg,
  family = betar(),
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

# Model comparison ----

cat("\n=== Model Comparison: Polymorphism Data ===\n")
anova(fit_null_poly, fit_ramp_poly, test = "Chisq")
anova(fit_ramp_poly, fit_ramp_int_poly, test = "Chisq")

aic_comparison <- AIC(fit_null_poly, fit_ramp_poly, fit_ramp_int_poly)
print(aic_comparison)

# Visualizations ----

# PLOT 1: Ramp shape from polymorphism data
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

pred_positions$fit <- pred_ramp$fit
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
    caption = "Ribbon = ±1.96 SE | Data from population genomics"
  ) +
  theme_bw(base_size = 12)

ggsave("./results/translational_ramp_polymorphism.pdf", 
       plot_ramp_poly, width = 10, height = 6)

# PLOT 2: Ramp by expression level
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

pred_grid_exp$fit <- pred_exp$fit
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
    color = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("./results/ramp_by_expression_polymorphism.pdf", 
       plot_ramp_exp_poly, width = 10, height = 6)

# Contrast the ramp vs the body

contrast_data <- poly_agg |>
  dplyr::mutate(Region = ifelse(Position_mid <= 50, "Ramp", "Body")) |>
  dplyr::group_by(Gene_clean, Region, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    Preferred_Freq_mean = mean(Preferred_Freq_mean, na.rm = TRUE),
    .groups = "drop"
  )

fit_contrast_poly <- bam(
  Preferred_Freq_mean ~ 
    Region + Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  
  data = contrast_data,
  family = betar(),
  method = "fREML",
  discrete = TRUE,
  nthreads = 4
)

summary(fit_contrast_poly)

# Extract statistics
region_coef <- coef(fit_contrast_poly)["RegionRamp"]
region_se <- summary(fit_contrast_poly)$se[names(coef(fit_contrast_poly)) == "RegionRamp"]

cat("\n=== Polymorphism Data: Ramp vs Body Contrast ===\n")
cat(sprintf("Difference (Ramp - Body): %.4f ± %.4f\n", region_coef, region_se))
cat(sprintf("t-statistic: %.2f\n", region_coef / region_se))
cat(sprintf("P-value: %.2e\n", 2 * pt(-abs(region_coef / region_se), 
                                      df = fit_contrast_poly$df.residual)))

baseline <- coef(fit_contrast_poly)["(Intercept)"]
cat(sprintf("\nMean Preferred Freq (Body): %.4f\n", baseline))
cat(sprintf("Mean Preferred Freq (Ramp): %.4f\n", baseline + region_coef))
cat(sprintf("Relative difference: %.1f%%\n", 
            100 * region_coef / baseline))

# Final memory cleanup: Section 15.2 large intermediates ---
# Keeping: trans, codon_usage, selection_metrics,
#          fit_null_poly, fit_ramp_poly, fit_ramp_int_poly, fit_contrast_poly
rm(poly_data, poly_with_exp, poly_agg, window_size,
   contrast_data,
   pred_positions, pred_ramp, plot_ramp_poly,
   pred_grid_exp, pred_exp, plot_ramp_exp_poly,
   region_coef, region_se, baseline,
   all_genes_15, sampled_genes_15, n_subsample_genes, n_threads_15)
gc()

## *****************************************************************************
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

# 16.5: C-ending vs G-ending preferred codon decomposition ----
# Classify each 4-fold site by 3rd base of its preferred codon
codon_4fold[, pref_base3 := substr(Preferred_Codon, 3, 3)]

n_C_sites <- sum(codon_4fold$pref_base3 == "C")
n_G_sites <- sum(codon_4fold$pref_base3 == "G")

cat(sprintf("  4-fold sites with C-ending preferred codon: %s\n",
            format(n_C_sites, big.mark = ",")))
cat(sprintf("  4-fold sites with G-ending preferred codon: %s\n",
            format(n_G_sites, big.mark = ",")))

# 16.5a: GAM fits for C-ending and G-ending ----
# GAM with tensor product smooth of position x expression
# Using te() to capture the interaction and individual non-linear effects

cat("\nFitting GAMs: Preferred_Freq ~ te(dist_norm, exp_norm)\n")

# Subsample for GAM fitting if datasets are very large (>500k)
# GAM with te() on millions of rows is memory-intensive
max_gam_sites <- 500000L

set.seed(1998)
if (n_C_sites > max_gam_sites) {
  idx_C <- sample.int(n_C_sites, max_gam_sites)
  gam_data_C <- codon_4fold[pref_base3 == "C"][idx_C]
  cat(sprintf("  C-ending: subsampled %s sites for GAM\n",
              format(max_gam_sites, big.mark = ",")))
} else {
  gam_data_C <- codon_4fold[pref_base3 == "C"]
}

if (n_G_sites > max_gam_sites) {
  idx_G <- sample.int(n_G_sites, max_gam_sites)
  gam_data_G <- codon_4fold[pref_base3 == "G"][idx_G]
  cat(sprintf("  G-ending: subsampled %s sites for GAM\n",
              format(max_gam_sites, big.mark = ",")))
} else {
  gam_data_G <- codon_4fold[pref_base3 == "G"]
}

gam_pf_C <- bam(
  Preferred_Freq ~ te(dist_norm, exp_norm, k = c(8, 8)),
  data = gam_data_C,
  method = "fREML",
  discrete = TRUE,
  nthreads = 4
)

gam_pf_G <- bam(
  Preferred_Freq ~ te(dist_norm, exp_norm, k = c(8, 8)),
  data = gam_data_G,
  method = "fREML",
  discrete = TRUE,
  nthreads = 4
)

cat("\n--- GAM Summary: C-ending preferred codons ---\n")
print(summary(gam_pf_C))
cat("\n--- GAM Summary: G-ending preferred codons ---\n")
print(summary(gam_pf_G))

# 16.5b: GAM contour plots for C-ending vs G-ending ----

pred_grid_cg <- expand.grid(
  dist_norm = seq(0, 1, length.out = 200),
  exp_norm  = seq(exp_range_16[1], exp_range_16[2], length.out = 200)
)

pred_grid_cg$pf_C <- predict(gam_pf_C, newdata = pred_grid_cg)
pred_grid_cg$pf_G <- predict(gam_pf_G, newdata = pred_grid_cg)

# C-ending contour
p_pf_C_contour <- ggplot(pred_grid_cg, aes(x = dist_norm, y = exp_norm)) +
  geom_raster(aes(fill = pf_C), interpolate = TRUE) +
  geom_contour(aes(z = pf_C), colour = "grey30", linewidth = 0.4, bins = 12) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Pref Freq"
  ) +
  labs(
    title = "GAM: Preferred Codon Frequency — C-ending Preferred Codons",
    subtitle = "CUB selection signal: expected to increase at high expression",
    x = "Normalised distance from gene start",
    y = expression(log[10](Expression))
  ) +
  theme_custom() +
  theme(legend.position = "right", panel.grid = element_blank())

# G-ending contour
p_pf_G_contour <- ggplot(pred_grid_cg, aes(x = dist_norm, y = exp_norm)) +
  geom_raster(aes(fill = pf_G), interpolate = TRUE) +
  geom_contour(aes(z = pf_G), colour = "grey30", linewidth = 0.4, bins = 12) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Pref Freq"
  ) +
  labs(
    title = "GAM: Preferred Codon Frequency — G-ending Preferred Codons",
    subtitle = "No selection support: expected to decrease (mutational decay to AT)",
    x = "Normalised distance from gene start",
    y = expression(log[10](Expression))
  ) +
  theme_custom() +
  theme(legend.position = "right", panel.grid = element_blank())

ggsave("./results/gam_4fold_pref_freq_C_ending_contour.pdf",
       p_pf_C_contour, width = 8, height = 7)
ggsave("./results/gam_4fold_pref_freq_G_ending_contour.pdf",
       p_pf_G_contour, width = 8, height = 7)

cat("  Saved: ./results/gam_4fold_pref_freq_C_ending_contour.pdf\n")
cat("  Saved: ./results/gam_4fold_pref_freq_G_ending_contour.pdf\n")

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

# 16.12) Per-amino-acid analysis of G-ending preferred codons ----
## Only 4 amino acids have G-ending preferred codons:
##   Glu (GAG), Lys (AAG), Arg_2 (AGG), Val (GTG)
## At these sites the ROC-preferred codon ends in G, but G is competing
## against C within the same amino acid family.  The aggregate "G-ending
## decline" (Section 16.10) could be confounded by C winning the within-
## family competition everywhere.  If genuine CUB acts on G at these four
## amino acids, Preferred_Freq should still increase with expression.
##
## DATA ISSUE — SPLIT FAMILIES:
## The ROC model assigns one preferred codon per amino acid *letter*.
## For Arg: ROC preferred = AGG (2-fold family, AGR).
## But 4-fold Arg sites are CGN — AGG can NEVER appear there.
## => Preferred_Freq = 0 at ALL CGN-Arg sites, and pref_base3 = "G"
##    from substr("AGG", 3, 3).
## This makes Arg a structural artefact: its "G-ending preferred" signal
## is uninformative (always zero).  We exclude Arg from the per-AA
## decomposition and flag it in the output.
##
## The same logic applies to split Leu (CTC preferred, 4-fold = CTN -> OK)
## and split Ser (TCC preferred, 4-fold = TCN -> OK).

cat("\n=== Section 16.12: Per-Amino-Acid Analysis of Preferred Codons ===\n")

# Identify which amino acid families have G-ending preferred codons
g_preferred_aas <- preferred_codons$aa[substr(preferred_codons$Codon, 3, 3) == "G"]
cat(sprintf("  Amino acids with G-ending preferred codon: %s\n",
            paste(g_preferred_aas, collapse = ", ")))
cat(sprintf("  Codons: %s\n",
            paste(preferred_codons$Codon[substr(preferred_codons$Codon, 3, 3) == "G"],
                  collapse = ", ")))

# Flag and exclude Arg: the 4-fold family (CGN) cannot contain the
# ROC-preferred codon (AGG), so Preferred_Freq = 0 everywhere.
n_arg_before <- nrow(codon_4fold[AA_name == "Arg"])
cat(sprintf("\n  WARNING: Excluding Arg (CGN) — Preferred_Codon=AGG cannot appear\n"))
cat(sprintf("  at CGN sites (%s sites removed).\n",
            format(n_arg_before, big.mark = ",")))
codon_4fold_noarg <- codon_4fold[AA_name != "Arg"]

# Recompute aggregate C/G counts without Arg
n_C_noarg <- sum(codon_4fold_noarg$pref_base3 == "C")
n_G_noarg <- sum(codon_4fold_noarg$pref_base3 == "G")
cat(sprintf("  After exclusion: C-ending=%s, G-ending=%s (change: -%s G sites)\n",
            format(n_C_noarg, big.mark = ","),
            format(n_G_noarg, big.mark = ","),
            format(n_arg_before - (n_G_noarg - n_G_sites + n_arg_before),
                   big.mark = ",")))

# ---- 16.12a: Gene-level aggregation per amino acid ----
# For each amino acid family, aggregate Preferred_Freq per gene and model vs expression
# Using codon_4fold_noarg (Arg excluded — see note above)

# Build per-AA, per-gene mean preferred freq
gene_pf_by_aa <- codon_4fold_noarg[, .(
  mean_pf   = mean(Preferred_Freq, na.rm = TRUE),
  n_sites   = .N,
  exp_norm  = Mean_Log10_Exp[1],
  log_len   = log10(Total_Codons[1]),
  pref_end  = pref_base3[1]        # same for all sites in this AA
), by = .(Gene_clean, AA_name)]

# Keep genes with >= 5 sites per amino acid
gene_pf_by_aa <- gene_pf_by_aa[n_sites >= 5]

# Spearman correlation: mean preferred freq ~ expression, per amino acid
aa_corr_table <- gene_pf_by_aa[, {
  ct <- cor.test(mean_pf, exp_norm, method = "spearman", exact = FALSE)
  .(rho = ct$estimate, p = ct$p.value, n_genes = .N,
    n_sites_total = sum(n_sites), pref_end = pref_end[1])
}, by = AA_name]

setorder(aa_corr_table, pref_end, -rho)

cat("\nSpearman correlation: gene-mean Preferred_Freq ~ expression, per amino acid\n")
cat("(genes with >= 5 sites per AA family)\n\n")
print(aa_corr_table, digits = 4)

write.csv(aa_corr_table,
          "./results/pref_freq_by_aa_expression_correlations.csv",
          row.names = FALSE)
cat("\nSaved: ./results/pref_freq_by_aa_expression_correlations.csv\n")

# ---- 16.12b: GAM trends for the 4 G-ending preferred AAs ----

cat("\n--- GAM trends for G-ending preferred amino acids ---\n")

# Filter to amino acids with G-ending preferred codons
# Map the preferred_codons$aa labels to the AA_name used in codon_4fold
g_aa_names <- preferred_codons$AA[substr(preferred_codons$Codon, 3, 3) == "G"]
# These are full names like "Val", "Glu" etc.

# Also get the equivalent names used in AA_name column of codon_4fold
# AA_name was created by aa_name_map which only maps A,G,L,P,R,S,T,V
# For the 4 G-ending preferred: Val (V), Glu (not 4-fold), Lys (not 4-fold),
# Arg_2 (not 4-fold -- it's AGG, 2-fold)
# Wait -- the 4-fold data only includes 4-fold degenerate families!
# Let's check which G-ending preferred AAs are actually in codon_4fold

g_end_in_data <- codon_4fold[pref_base3 == "G", unique(AA_name)]
cat(sprintf("  G-ending preferred AAs in 4-fold data: %s\n",
            paste(g_end_in_data, collapse = ", ")))

# Gene-level data for G-ending preferred AAs
gene_pf_g_aas <- gene_pf_by_aa[pref_end == "G"]

if (nrow(gene_pf_g_aas) > 0 && length(unique(gene_pf_g_aas$AA_name)) > 0) {
  
  # GAM per AA for G-ending
  gam_results_g <- list()
  for (aa in unique(gene_pf_g_aas$AA_name)) {
    aa_data <- gene_pf_g_aas[AA_name == aa]
    if (nrow(aa_data) >= 50) {
      gam_aa <- gam(mean_pf ~ s(exp_norm, k = 8), data = aa_data, method = "REML")
      gam_results_g[[aa]] <- gam_aa
      cat(sprintf("  %s: n_genes=%d, R2=%.4f, deviance_expl=%.1f%%\n",
                  aa, nrow(aa_data),
                  summary(gam_aa)$r.sq,
                  summary(gam_aa)$dev.expl * 100))
    }
  }
  
  # Prediction curves for each G-ending AA
  if (length(gam_results_g) > 0) {
    exp_pred_range <- range(gene_pf_g_aas$exp_norm, na.rm = TRUE)
    exp_pred_seq <- seq(exp_pred_range[1], exp_pred_range[2], length.out = 200)
    
    g_trend_list <- lapply(names(gam_results_g), function(aa) {
      nd <- data.frame(exp_norm = exp_pred_seq)
      pp <- predict(gam_results_g[[aa]], newdata = nd, se.fit = TRUE)
      pref_codon <- preferred_codons$Codon[
        preferred_codons$AA == aa | preferred_codons$aa == aa]
      if (length(pref_codon) == 0) pref_codon <- "?"
      data.frame(
        exp_norm = exp_pred_seq,
        fit = pp$fit, se = pp$se.fit,
        AA = aa,
        label = sprintf("%s (%s)", aa, pref_codon[1])
      )
    })
    g_trend_df <- do.call(rbind, g_trend_list)
    g_trend_df$lo <- g_trend_df$fit - 1.96 * g_trend_df$se
    g_trend_df$hi <- g_trend_df$fit + 1.96 * g_trend_df$se
    
    p_g_aa_trends <- ggplot(g_trend_df,
                             aes(x = exp_norm, y = fit, color = label, fill = label)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
      geom_line(linewidth = 1.2) +
      labs(
        title = "G-ending Preferred Codons: Gene-Mean Preferred Freq vs Expression",
        subtitle = paste(
          "If CUB supports G at these amino acids, preferred freq should increase.",
          "\nDecline -> C within-family competitor is winning / mutational decay."
        ),
        x = expression(log[10](Expression)),
        y = "Gene-Mean Preferred Codon Frequency",
        color = "Amino acid (preferred codon)",
        fill  = "Amino acid (preferred codon)"
      ) +
      theme_custom() +
      theme(legend.position = "bottom")
    
    ggsave("./results/G_ending_preferred_per_AA_trends.pdf",
           p_g_aa_trends, width = 9, height = 6)
    cat("  Saved: ./results/G_ending_preferred_per_AA_trends.pdf\n")
  }
} else {
  cat("  No G-ending preferred amino acids found in 4-fold data.\n")
  cat("  (This is expected if only Val has G-ending preferred among 4-fold families,\n")
  cat("   and Glu/Lys/Arg_2 are 2-fold.)\n")
}

# ---- 16.12c: Comparison panel: all AAs, colored by preferred ending ----

cat("\n--- All amino acid families: preferred freq trend by ending type ---\n")

# GAM per AA for ALL amino acids in 4-fold data
gam_results_all_aa <- list()
for (aa in unique(gene_pf_by_aa$AA_name)) {
  aa_data <- gene_pf_by_aa[AA_name == aa]
  if (nrow(aa_data) >= 50) {
    gam_aa <- gam(mean_pf ~ s(exp_norm, k = 8), data = aa_data, method = "REML")
    gam_results_all_aa[[aa]] <- list(
      model = gam_aa,
      pref_end = aa_data$pref_end[1],
      n_genes = nrow(aa_data)
    )
  }
}

# Build prediction dataframe
if (length(gam_results_all_aa) > 0) {
  exp_pred_range_all <- range(gene_pf_by_aa$exp_norm, na.rm = TRUE)
  exp_pred_seq_all <- seq(exp_pred_range_all[1], exp_pred_range_all[2],
                          length.out = 200)
  
  all_aa_trend_list <- lapply(names(gam_results_all_aa), function(aa) {
    obj <- gam_results_all_aa[[aa]]
    nd <- data.frame(exp_norm = exp_pred_seq_all)
    pp <- predict(obj$model, newdata = nd, se.fit = TRUE)
    data.frame(
      exp_norm = exp_pred_seq_all,
      fit      = pp$fit,
      se       = pp$se.fit,
      AA       = aa,
      pref_end = obj$pref_end,
      n_genes  = obj$n_genes
    )
  })
  all_aa_trend_df <- do.call(rbind, all_aa_trend_list)
  all_aa_trend_df$lo <- all_aa_trend_df$fit - 1.96 * all_aa_trend_df$se
  all_aa_trend_df$hi <- all_aa_trend_df$fit + 1.96 * all_aa_trend_df$se
  
  # Color by preferred ending: blue = C-ending, green/orange = G-ending
  all_aa_trend_df$Ending <- ifelse(all_aa_trend_df$pref_end == "C",
                                    "C-ending preferred", "G-ending preferred")
  
  p_all_aa <- ggplot(all_aa_trend_df,
                      aes(x = exp_norm, y = fit, color = AA, linetype = Ending)) +
    geom_line(linewidth = 1.0) +
    scale_linetype_manual(values = c("C-ending preferred" = "solid",
                                      "G-ending preferred" = "dashed")) +
    facet_wrap(~ Ending, ncol = 2, scales = "free_y") +
    labs(
      title = "Gene-Mean Preferred Codon Freq vs Expression: All 4-fold AA Families",
      subtitle = "Solid = C-ending preferred, Dashed = G-ending preferred",
      x = expression(log[10](Expression)),
      y = "Gene-Mean Preferred Codon Frequency",
      color = "Amino acid"
    ) +
    guides(linetype = "none") +
    theme_custom() +
    theme(legend.position = "bottom")
  
  ggsave("./results/pref_freq_all_AA_by_ending.pdf",
         p_all_aa, width = 12, height = 6)
  cat("  Saved: ./results/pref_freq_all_AA_by_ending.pdf\n")
  
  # Print summary table
  aa_summary <- do.call(rbind, lapply(names(gam_results_all_aa), function(aa) {
    obj <- gam_results_all_aa[[aa]]
    ss <- summary(obj$model)
    # Low and high expression predictions
    nd_lo <- data.frame(exp_norm = quantile(exp_pred_seq_all, 0.05))
    nd_hi <- data.frame(exp_norm = quantile(exp_pred_seq_all, 0.95))
    delta <- predict(obj$model, newdata = nd_hi) -
             predict(obj$model, newdata = nd_lo)
    data.frame(
      AA       = aa,
      Pref_End = obj$pref_end,
      N_genes  = obj$n_genes,
      R2       = round(ss$r.sq, 4),
      Dev_Expl = round(ss$dev.expl * 100, 1),
      Delta_lo_hi = round(delta, 4),
      Direction = ifelse(delta > 0, "UP", "DOWN")
    )
  }))
  setorder(as.data.table(aa_summary), Pref_End, -Delta_lo_hi)
  cat("\nSummary: GAM-predicted change (low -> high expression) per amino acid:\n\n")
  print(aa_summary, row.names = FALSE)
  
  write.csv(aa_summary,
            "./results/pref_freq_per_AA_expression_direction.csv",
            row.names = FALSE)
  cat("\nSaved: ./results/pref_freq_per_AA_expression_direction.csv\n")
}

# ---- 16.12d: GAM contour plots per amino acid ----
# Contour surface: Preferred_Freq ~ f(expression, gene length) for each AA
# These give a richer picture than the 1-D trend curves

cat("\n--- GAM contour plots per amino acid (expression x gene length) ---\n")

# Build gene-level data with both predictors per AA
gene_pf_for_contour <- codon_4fold_noarg[, .(
  mean_pf   = mean(Preferred_Freq, na.rm = TRUE),
  n_sites   = .N,
  exp_norm  = Mean_Log10_Exp[1],
  log_len   = log10(Total_Codons[1]),
  pref_end  = pref_base3[1]
), by = .(Gene_clean, AA_name)]
gene_pf_for_contour <- gene_pf_for_contour[n_sites >= 5]

# Prediction grid (shared across AAs)
exp_range_aa <- range(gene_pf_for_contour$exp_norm, na.rm = TRUE)
len_range_aa <- range(gene_pf_for_contour$log_len, na.rm = TRUE)

pred_grid_aa <- expand.grid(
  exp_norm = seq(exp_range_aa[1], exp_range_aa[2], length.out = 100),
  log_len  = seq(len_range_aa[1], len_range_aa[2], length.out = 100)
)

# Fit te() GAM per AA and collect predictions
aa_list <- sort(unique(gene_pf_for_contour$AA_name))
contour_preds <- list()

for (aa in aa_list) {
  aa_data <- gene_pf_for_contour[AA_name == aa]
  if (nrow(aa_data) < 200) next
  
  gam_aa_2d <- gam(
    mean_pf ~ te(exp_norm, log_len, k = c(6, 6)),
    data = aa_data, method = "REML"
  )
  
  ss <- summary(gam_aa_2d)
  cat(sprintf("  %s: n=%d, R2=%.4f, dev_expl=%.1f%%\n",
              aa, nrow(aa_data), ss$r.sq, ss$dev.expl * 100))
  
  preds <- predict(gam_aa_2d, newdata = pred_grid_aa)
  contour_preds[[aa]] <- data.frame(
    pred_grid_aa,
    pf_pred  = as.numeric(preds),
    AA       = aa,
    pref_end = aa_data$pref_end[1]
  )
}

contour_all <- do.call(rbind, contour_preds)
contour_all$Ending <- ifelse(contour_all$pref_end == "C",
                              "C-ending preferred", "G-ending preferred")

# Faceted contour plot: one panel per amino acid
p_aa_contours <- ggplot(contour_all, aes(x = exp_norm, y = log_len)) +
  geom_raster(aes(fill = pf_pred), interpolate = TRUE) +
  geom_contour(aes(z = pf_pred), colour = "grey30", linewidth = 0.3, bins = 10) +
  facet_wrap(~ AA, ncol = 4, scales = "free") +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
                "#FEE8C8", "#FDBB84", "#E34A33"),
    name = "Predicted\nPref Freq"
  ) +
  labs(
    title = "GAM: Gene-Mean Preferred Codon Frequency per Amino Acid",
    subtitle = paste(
      "te(expression, log10 gene length). Arg excluded (split-family artefact).",
      "\nWarm colours = higher preferred codon frequency."
    ),
    x = expression(log[10](Expression)),
    y = expression(log[10](Gene~length~"in"~codons))
  ) +
  theme_custom() +
  theme(legend.position = "right", panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("./results/gam_pref_freq_per_AA_contour.pdf",
       plot = p_aa_contours, width = 14, height = 10)
cat("  Saved: ./results/gam_pref_freq_per_AA_contour.pdf\n")

# Also make individual high-resolution contour PDFs for each AA
for (aa in unique(contour_all$AA)) {
  # 1. Subset data
  aa_df <- contour_all[contour_all$AA == aa, ]
  ending_label <- aa_df$Ending[1]
  
  # 2. Build plot
  p_single <- ggplot(aa_df, aes(x = exp_norm, y = log_len)) +
    geom_raster(aes(fill = pf_pred), interpolate = TRUE) +
    geom_contour(aes(z = pf_pred), colour = "grey30", linewidth = 0.4, bins = 12) +
    scale_fill_gradientn(
      colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#FEE8C8", "#FDBB84", "#E34A33"),
      name = "Pref Freq"
    ) +
    labs(
      title = sprintf("GAM: Preferred Codon Frequency — %s (%s)", aa, ending_label),
      x = expression(log[10](Expression)),
      # Use quotes around "in" to prevent the reserved keyword error
      y = expression(log[10](Gene~length~"in"~codons)) 
    ) +
    theme_custom() +
    theme(legend.position = "right", panel.grid = element_blank())
  
  # 3. Save plot
  ggsave(sprintf("./results/gam_pref_freq_%s_contour.pdf", tolower(aa)), 
         plot = p_single, width = 8, height = 7)
}

cat(sprintf("  Saved %d individual AA contour plots\n", 
            length(unique(contour_all$AA))))

# Section 16 cleanup ----
rm(codon_4fold, codon_4fold_noarg, gene_info_16, freq_clean, fourfold_simple,
   aa_name_map,
   model_all_4fold, pred_grid_16, p_contour_filled, p_contour_lines,
   logistic_table, logistic_table_z, result_all, result_all_z,
   results_aa, results_aa_z, aa_levels_16, exp_range_16,
   model_pi_all, pred_grid_pi, p_pi_contour, pi_table, result_pi_all, results_pi_aa,
   model_pf_all, pred_grid_pf, p_pf_contour, pf_table, result_pf_all, results_pf_aa,
   gam_pf_C, gam_pf_G, gam_data_C, gam_data_G,
   pred_grid_cg, p_pf_C_contour, p_pf_G_contour,
   trend_df, trend_long, p_cg_trends,
   mid_df, p_cg_mid, pred_C_mid, pred_G_mid,
   pf_cg_table_z, result_pf_C_z, result_pf_G_z,
   n_C_sites, n_G_sites, n_C_noarg, n_G_noarg, n_arg_before, exp_seq,
   gene_pf_CG, gene_pf_C, gene_pf_G, gam_gene_C, gam_gene_G,
   exp_grid_gene, pred_gene_C, pred_gene_G, gene_trend_df, p_gene_cg,
   gene_pf_by_aa, aa_corr_table, gene_pf_g_aas,
   all_aa_trend_df, aa_summary, p_all_aa, p_g_aa_trends,
   gene_pf_for_contour, pred_grid_aa, contour_preds, contour_all,
   p_aa_contours, exp_range_aa, len_range_aa, aa_list)
gc()

message("=== Analysis complete. Memory cleaned. ===")