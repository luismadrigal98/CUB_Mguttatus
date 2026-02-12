##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
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
                        'betareg', 'brms')

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
  na.omit() |>
  distinct(Gene_name, .keep_all = TRUE)

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
    CDC_significant = !is.na(p_value) & p_value < 0.05,
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

write.csv(justification_table, "results/Linearity_Justification_Table.csv", 
          row.names = FALSE)

# Generate Visuals (Safe Loop) ----
plot_list <- list()

for (pred in predictors) {
    
  form_gam <- as.formula(paste0("CDC ~ s(", pred, ")"))
  model_gam <- gam(form_gam, data = integrated_data, 
                   family = betar(link = "logit"),
                   method = "REML",
                   select = T)
  
  p <- gratia::draw(model_gam, residuals = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", 
               alpha = 0.5) +
    labs(title = paste0("Partial Effect on CDC (logit scale): ", pred)) +
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
              data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
              family = betar(link = "logit"), 
              method = "REML")

# Model 1: Additive (Independent effects)
# Hypothesis: Each predictor affects CUB independently.

m_additive <- gam(CDC ~ s(Max_Log10_Exp) + s(Exp_breadth) + s(CDS_length_nt),
                  data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
                  family = betar(link = "logit"), 
                  method = "REML",
                  select = T)

# Model 2: Expression Interaction (The "Trade-off" Hypothesis)
# Hypothesis: High expression only forces strict CUB if the gene is broad.

m_interaction <- gam(CDC ~ te(Max_Log10_Exp, Exp_breadth) + s(CDS_length_nt),
                     data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
                     family = betar(link = "logit"), 
                     method = "REML",
                     select = T)

# Model 3: Complex (Full Interaction)
# Hypothesis: Length and expression interact in complex ways."
m_complex <- gam(CDC ~ te(Max_Log10_Exp, Exp_breadth, CDS_length_nt),
                 data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
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
                                CDS_length_nt = mean(clean_data$CDS_length_nt)),
                              type = "response") + 
  geom_rug(data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
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
                          CDS_length_nt = mean(clean_data$CDS_length_nt)),
                        type = "response") +
  # RED LINE = Zero Slope (Where the relationship flattens/saturates)
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_rug(data = integrated_data |> dplyr::filter(Exp_breadth > 0), 
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

# Pairwise Tests (Comparing the slopes)
# We test if the relationship between Exp and CDC is stronger in one group.
message("Test: Difference in Slopes (Broad vs Narrow) ---\n")
test_broad_vs_narrow <- hypotheses(slopes, hypothesis = "b3 - b1 = 0")
print(test_broad_vs_narrow)

message("Test: Difference in Slopes (Medium vs Narrow) ---\n")
test_medium_vs_narrow <- hypotheses(slopes, hypothesis = "b2 - b1 = 0")
print(test_medium_vs_narrow)

# VARIANCE ANALYSIS (Evolutionary Constraint)
# Question: Does breadth constrain the variation in codon bias?

# Extract Residuals (The "Noise" or deviation from the expected CDC)
# If residuals are high, other factors (mutation, drift) are overpowering expression.
integrated_data$bias_deviation <- abs(residuals(m_interaction, type = "response"))

# Visualize the Constraint
p_constraint <- ggplot(integrated_data, aes(x = Exp_breadth, y = bias_deviation)) +
  geom_point(alpha = 0.1, color = "gray60") +
  
  # Fit a smooth trend line to the noise
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              color = "firebrick", fill = "firebrick") +
  
  theme_custom() +
  labs(y = "Absolute Residuals (Deviation from Expected CDC)",
       x = "Expression Breadth (Number of Tissues)")

ggsave("./results/Evolutionary_Constraint_Plot.pdf", 
       p_constraint, width = 8, height = 6)

# Statistical Test of Constraint
# Gamma regression is ideal for modeling variance (always positive)
m_noise <- gam(bias_deviation ~ s(Exp_breadth), 
               data = clean_data, 
               family = Gamma(link = "log"))

print(summary(m_noise))

# Compare predictions at extreme expression levels, holding length constant
contrasts <- avg_comparisons(
  m_interaction,
  variables = "Max_Log10_Exp",
  newdata = datagrid(
    Max_Log10_Exp = c(
      quantile(integrated_data$Max_Log10_Exp, 0.05),  # Bottom 5%
      quantile(integrated_data$Max_Log10_Exp, 0.95)   # Top 5%
    ),
    CDS_length_nt = mean(integrated_data$CDS_length_nt),
    Exp_breadth = c(1, 29)
  )
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
        axis.title.y = element_text(size = 11),,
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

# 7.1) Making bar plot per aminoacid to show differences in use ----

ggplot(data = w_table |> dplyr::filter(amino_acid != "Met" & amino_acid != "Trp" &
                                         amino_acid != "STOP"), 
       mapping = aes(x = reorder(codon, relative_adaptiveness), 
                                     y = relative_adaptiveness)) +
  geom_segment(aes(xend = codon, y = 0, yend = relative_adaptiveness), 
               color = "gray70", linewidth = 1) +
  geom_point(aes(color = relative_adaptiveness), size = 4) +
  scale_color_gradient2(
    low = "#d73027", 
    mid = "#fee090", 
    high = "#1a9850", 
    midpoint = 0.5,
    name = "w"
  ) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(~amino_acid, scales = "free", ncol = 4) +
  coord_flip() +
  theme_custom() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(y = "Relative Adaptiveness (w)", x = "Codon")

ggsave("./results/codons_relative_adaptiveness.pdf", 
       width = 12, height = 10)

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

# 7.2) Compare absolute codon frequencies: Top 5% vs Rest ----
# This shows that raw frequencies differ, but not all differences are due to 
# selection some codons are frequent simply because their amino acids are frequent

# Get gene lists
top5_genes <- integrated_data |>
  filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

rest_genes <- integrated_data |>
  filter(Expression_Group %in% c("Middle 90%", "Bottom 5%")) |>
  pull(Gene_name)

# Calculate absolute frequencies (sum of codon counts)
codon_cols <- setdiff(names(codon_usage), "Gene_name")

freq_top5 <- codon_usage |>
  dplyr::filter(Gene_name %in% top5_genes) |>
  dplyr::select(all_of(codon_cols)) |>
  summarise(across(everything(), sum, na.rm = TRUE)) |>
  pivot_longer(everything(), names_to = "Codon", values_to = "Count_Top5")

freq_rest <- codon_usage |>
  dplyr::filter(Gene_name %in% rest_genes) |>
  dplyr::select(all_of(codon_cols)) |>
  summarise(across(everything(), sum, na.rm = TRUE)) |>
  pivot_longer(everything(), names_to = "Codon", values_to = "Count_Rest")

# Combine and calculate proportions
freq_comparison <- freq_top5 |>
  left_join(freq_rest, by = "Codon") |>
  dplyr::mutate(
    Total_Top5 = sum(Count_Top5),
    Total_Rest = sum(Count_Rest),
    Freq_Top5 = Count_Top5 / Total_Top5,
    Freq_Rest = Count_Rest / Total_Rest,
    Freq_Diff = Freq_Top5 - Freq_Rest
  )

# Add amino acid information
freq_comparison <- freq_comparison |>
  mutate(AA = genetic_code_dna_long[Codon]) |>
  filter(!is.na(AA) & AA != "STOP")

# Convert to long format for plotting
freq_long <- freq_comparison |>
  dplyr::select(Codon, AA, Freq_Top5, Freq_Rest) |>
  pivot_longer(cols = c(Freq_Top5, Freq_Rest), 
               names_to = "Group", 
               values_to = "Frequency") |>
  dplyr::mutate(Group = dplyr::recode(Group, 
                               "Freq_Top5" = "Top 5%",
                               "Freq_Rest" = "Rest")) |>
  dplyr::filter(AA != "Trp" & AA != "Met")

# Create lollipop plot
p_freq_comparison <- ggplot(freq_long, 
                            aes(x = reorder(Codon, Frequency), 
                                y = Frequency, 
                                color = Group)) +
  geom_line(aes(group = Codon), color = "gray80", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", "Rest" = "#377EB8"),
                     name = "Gene Group") +
  facet_wrap(~AA, scales = "free", ncol = 4) +
  coord_flip() +
  theme_custom() +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top"
  ) +
  labs(
    y = "Frequency",
    x = "Codon",
    title = "Raw Codon Frequencies: Top 5% vs Rest",
    subtitle = "Not all differences reflect selection - some codons are frequent due to amino acid composition"
  )

ggsave("./results/codon_frequency_top5_vs_rest.pdf", 
       p_freq_comparison, width = 12, height = 14)

# Summary statistics
cat(sprintf("\nTop 5%% genes: %d genes, %d total codons\n", 
            length(top5_genes), freq_comparison$Total_Top5[1]))
cat(sprintf("Rest genes: %d genes, %d total codons\n", 
            length(rest_genes), freq_comparison$Total_Rest[1]))

# 7.2.1) Enriched codons usage ----

# Get enriched codons (w = 1.0 from CAI)
enriched_codons_vec <- w_table |>
  dplyr::filter(relative_adaptiveness == 1.0) |>
  dplyr::pull(codon)

# Merge codon usage with expression groups
codon_usage_with_groups <- codon_usage |>
  dplyr::left_join(integrated_data |> dplyr::select(Gene_name, Expression_Group), 
                   by = "Gene_name")

# Filter to top 5% and rest (bottom 95%)
top5_genes <- codon_usage_with_groups |> dplyr::filter(Expression_Group == "Top 5%")
rest_genes <- codon_usage_with_groups |> dplyr::filter(Expression_Group != "Top 5%")


# Calculate for both groups
cat("Calculating enrichment of codon per amino acid...\n")

enriched_aa <- count_preferred_by_aa(top5_genes, enriched_codons_vec, 
                                     genetic_code_dna_long)
enriched_aa$Group <- "Selected (Top 5%)"

rest_aa <- count_preferred_by_aa(rest_genes, enriched_codons_vec, 
                                 genetic_code_dna_long)
rest_aa$Group <- "Rest (Bottom 95%)"

# Combine for comparison
comparison_table <- enriched_aa |>
  dplyr::select(Amino_Acid, N_synonymous, Preferred_codons, 
                Selected_count = Preferred_count, 
                Selected_prop = Prop_preferred) |>
  dplyr::left_join(
    rest_aa |> dplyr::select(Amino_Acid, 
                             Rest_count = Preferred_count,
                             Rest_prop = Prop_preferred),
    by = "Amino_Acid"
  ) |>
  dplyr::mutate(
    Difference = Selected_prop - Rest_prop,
    Fold_enrichment = Selected_prop / Rest_prop
  ) |>
  dplyr::arrange(dplyr::desc(Difference))

# Save table
write.csv(comparison_table, "./results/enriched_codon_usage.csv",
          row.names = FALSE)

# Wilcoxon test
wilcox_test <- wilcox.test(comparison_table$Selected_prop, 
                           comparison_table$Rest_prop,
                           paired = TRUE)

# Create visualization
p_comparison <- ggplot(comparison_table, 
                       aes(x = reorder(Amino_Acid, -Difference))) +
  geom_segment(aes(xend = Amino_Acid, y = Rest_prop, yend = Selected_prop),
               color = "gray70", linewidth = 1) +
  geom_point(aes(y = Selected_prop, color = "Top 5%"), 
             size = 3, shape = 16) +
  geom_point(aes(y = Rest_prop, color = "Rest 95%"), 
             size = 3, shape = 16) +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                "Rest 95%" = "#377EB8"),
                     name = "") +
  labs(x = "Amino Acid (ordered by difference)",
       y = "Proportion of Enriched Codons") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("./results/preferred_codon_usage_comparison.pdf", p_comparison,
       width = 10, height = 6)

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

# 1. Filter for complete cases (Intersection of Leaf and Bud)
# We strictly remove genes with 0 counts in either tissue
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

# Load validation functions
source("./src/roc_model_validation.R")

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
  mutation_file = "./results/MCMC_results/results_dM_fixed_with_phi/run_1/Parameter_est/Cluster_1_Mutation.csv",
  selection_file = "./results/MCMC_results/results_dM_fixed_with_phi/run_1/Parameter_est/Cluster_1_Selection.csv",
  codon_freq_df = codon_freq_long,
  expression_df = expr_data,
  output_file = "./results/ROC_codon_trajectories.pdf",
  n_bins = 10
)

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

# Convergence achieved.

# Checking correlation with empirical values

phi_hat_dM_fixed_intergenic <- read.csv(file = "results/MCMC_results/results_dM_fixed_intergenic/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed_intergenic <- exp_complete |>
  left_join(phi_hat_dM_fixed_intergenic, by = join_by("Gene_name" == "GeneID")) |>
  na.exclude()

cor.test(phi_dM_fixed_intergenic$Mean.log10.Phi, 
         phi_dM_fixed_intergenic$Mean_Log10_Exp)

# Inappropriate fix to empirical data

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

# Convergence: TRUE (intergenic-based mutation bias is not adequate)

# Checking the correlation between phi and empirical values

phi_hat_dM_fixed_with_phi_intergenic <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi_intergenic_final/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed_with_phi_intergenic <- exp_complete |>
  left_join(phi_hat_dM_fixed_with_phi_intergenic, by = join_by("Gene_name" == "GeneID")) |>
  na.exclude()

cor.test(phi_dM_fixed_with_phi_intergenic$Mean.log10.Phi, 
         phi_dM_fixed_with_phi_intergenic$Mean_Log10_Exp)

# Exploring and overlapping traces to rule out bad mixing ----

parameters_objects <- list(run1 = loadParameterObject(file = paste(run_dirs[1], "R_objects/parameter.Rda", sep = '/')),
                           run2 = loadParameterObject(file = paste(run_dirs[2], "R_objects/parameter.Rda", sep = '/')),
                           run3 = loadParameterObject(file = paste(run_dirs[3], "R_objects/parameter.Rda", sep = '/')))

parameters_objects <- list(run1 = listRDA(paste(run_dirs[1], "R_objects/parameter.Rda", sep = '/')),
                           run2 = listRDA(paste(run_dirs[2], "R_objects/parameter.Rda", sep = '/')),
                           run3 = listRDA(paste(run_dirs[3], "R_objects/parameter.Rda", sep = '/')))

# 1. Select one of the high deviants codons
codon_index <- 36  

# 2. Create the Master Data Frame
plot_data <- extract_traces(parameters_objects, codon_index)
plot_data <- subset(plot_data, Iteration > max(plot_data$Iteration) * 0.5)

# 3. PLOT 1: The Trace Overlay (The "Are we mixing?" Plot)
p1 <- ggplot(plot_data, aes(x = Iteration, y = Value, color = Run)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  theme_custom() +
  labs(title = paste("Trace Overlay: Codon Index", codon_index),
       subtitle = "Flat lines at different levels = Good Mixing + Multimodality",
       y = "Selection Cost (Delta Eta)",
       x = "MCMC Sample") +
  theme(legend.position = "top")

# 4. PLOT 2: Density Overlay (The "Distinct Realities" Plot)
p2 <- ggplot(plot_data, aes(x = Value, fill = Run)) +
  geom_density(alpha = 0.5) +
  theme_custom() +
  labs(title = paste("Posterior Density: Codon Index", codon_index),
       subtitle = "Non-overlapping peaks = Model Misspecification / Broken Energy Landscape",
       x = "Selection Cost (Delta Eta)",
       y = "Density") +
  theme(legend.position = "top")

# 5. PLOT 3: Autocorrelation (The "Stickiness" Check)
# We calculate ACF for just Run 1 to prove it's not "sticky"
acf_val <- acf(subset(plot_data, Run == "run1")$Value, plot = FALSE, lag.max = 40)
acf_df <- data.frame(Lag = acf_val$lag, ACF = acf_val$acf)

p3 <- ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -0.05, linetype = "dashed", color = "red") +
  theme_custom() +
  labs(title = "Autocorrelation (Run 1)",
       subtitle = "Rapid drop to 0 = Efficient Sampling (Not Bad Mixing)",
       y = "Autocorrelation")

# Display Plots
ggsave(filename = paste0("./results/ROC_dM_fixed_with_phi_intergenic_codon_", codon_index, "_trace_overlay.pdf"),
       plot = p1, width = 10, height = 4)
ggsave(filename = paste0("./results/ROC_dM_fixed_with_phi_intergenic_codon_", codon_index, "_density_overlay.pdf"),
       plot = p2, width = 8, height = 4)
ggsave(filename = paste0("./results/ROC_dM_fixed_with_phi_intergenic_codon_", codon_index, "_acf.pdf"),
       plot = p3, width = 6, height = 4)

# 8.2) Getting the preferred codon from the best model (dM-fixed-with_phi) ----

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

message(sprintf("✓ Preferred codons from ROC model: %d amino acids\n", nrow(preferred_codons_roc)))

# 8.3) Extracting selection estimates from the best model (dM-fixed-with_phi) ----

genome <- initializeGenomeObject(file = 'data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa',
                                 match.expression.by.id = TRUE,
                                 observed.expression.file = 'data/observed_expression_multitissue.csv') 

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

# Total selection intensity
total_selection_intensity <- rowSums(counts_aligned * abs(sel_aligned), na.rm = TRUE)

# S_ROC: The average selection coefficient acting on a codon in this gene
S_ROC <- total_selection_intensity / n_synonymous_codons

# L_ROC: Total Load (The total fitness cost of this gene's sequence)
L_ROC <- rowSums(counts_aligned * load_matrix_per_codon, na.rm = TRUE)

# Lprime_ROC: Average Load per site (Length-corrected)
Lprime_ROC <- L_ROC / n_synonymous_codons

selection_metrics <- data.frame(
  Gene_name = common_genes,
  
  # Use this for Drift Barrier Plots (The "Hump" x-axis)
  # Theoretical Threshold: If S_avg < 1, Drift Dominates. If S_avg > 1, Selection Dominates.
  S_ROC = S_ROC, 
  
  # Use these for assessing maladaptation
  L_ROC = L_ROC,
  Lprime_ROC = Lprime_ROC,
  
  n_codons = n_synonymous_codons,
  
  row.names = common_genes
)

selection_metrics <- selection_metrics |>
  left_join(phi_hat_dM_fixed_with_phi |> dplyr::select(GeneID, Mean.log10.Phi, MeanPhi),
            by = join_by(Gene_name == GeneID))

# 8.3.1) Relationship between L' and phi ----

final_analysis_data <- selection_metrics |>
  dplyr::filter(S_ROC > 0) |>
  dplyr::mutate(
    Intrinsic_Inefficiency = S_ROC / MeanPhi
  )

p_load <- ggplot(final_analysis_data, aes(x = Mean.log10.Phi, 
                                          y = Lprime_ROC)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log10", 
                       name = "Gene Count") +
  # Trend line
  geom_smooth(method = "gam", color = "cyan", size = 1.2, se = TRUE) +
  labs(
    x = expression(bold(Log[10]("Expression" ~ (Phi)))),
    y = expression(bold("Realized Load" ~ (L[prime])))
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

cor_load <- cor.test(final_analysis_data$Mean.log10.Phi, final_analysis_data$Lprime_ROC, method = "spearman")
cor_eff <- cor.test(final_analysis_data$Mean.log10.Phi, final_analysis_data$Intrinsic_Inefficiency, method = "spearman")
# Only in tail
selection_genes <- final_analysis_data |> dplyr::filter(Mean.log10.Phi >= 1.5)
cor_selection <- cor.test(selection_genes$Mean.log10.Phi, 
                          selection_genes$Intrinsic_Inefficiency, method = "spearman",
                          exact = F)

# 8.3.2) Analyzing the correlation between total selective pressure and CAI and CDC ----

integrated_data <- integrated_data |>
  left_join(selection_metrics) |>
  na.exclude()

# Correlation between selection metric and CUB metrics

cor_S_and_bias <- corrr::correlate(x = as.matrix(integrated_data[, c(29, 6, 16, 28)]),
                                   method = "spearman")

cor.test(integrated_data$S_ROC, integrated_data$CAI)

# 8.3.3) Final visualization ----

plot_data <- integrated_data |>
  dplyr::mutate(
    # Log Transform Selection Load (Load = Total Cost per Gene)
    # Note: If you want Intensity (per codon), swap S_load for Selection_Intensity
    Log_S_ROC = log10(S_ROC + 0.01), 
    
    # Log Transform Expression (using the new clear name)
    Log_Phi = Mean.log10.Phi, 
    
    # Log Transform Length
    Log_Length = log10(Total_Codons)
  ) |>
  dplyr::filter(!is.na(ENC), !is.na(Total_Codons), !is.na(CAI))

# 4. Visualization Setup

# Define common color scale limits
phi_range <- range(plot_data$Log_Phi, na.rm = TRUE)

# Fit linear model for Panel C annotation (Load vs Length)
tail_data <- plot_data |> 
  dplyr::filter(Log_S_ROC > -0.5)

# Fit LM specifically on the tail
lm_tail <- lm(GC3s ~ Log_S_ROC, data = tail_data)
lm_tail_eq <- sprintf("Tail: y = %.2f + %.2fx\nR² = %.3f, p = %.4f",
                      coef(lm_tail)[1], 
                      coef(lm_tail)[2], 
                      summary(lm_tail)$r.squared,
                      summary(lm_tail)$coefficients["Log_S_ROC","Pr(>|t|)"])

# Panel A: Selection Load Distribution
drift_thresh <- log10(1)   
y_max_anno <- 10000

p1 <- ggplot(plot_data, aes(x = Log_S_ROC)) +
  
  # --- Background Shading ---
  annotate("rect", xmin = -Inf, xmax = drift_thresh, 
           ymin = 0, ymax = Inf, fill = "gray95", alpha = 0.8) +
  annotate("rect", xmin = drift_thresh, xmax = Inf, 
           ymin = 0, ymax = Inf, fill = "#ffe5e5", alpha = 0.5) +
  
  # --- Histogram ---
  geom_histogram(bins = 100, fill = "#69b3a2", color = "white", linewidth = 0.05) +
  
  # --- Vertical Threshold Lines ---
  geom_vline(xintercept = drift_thresh, linetype = "dotted", color = "red") +
  
  # --- Vertical Text Annotations (No Ne*s text) ---
  # Drift Label (Left side, Vertical)
  annotate("text", x = log10(0.05), y = y_max_anno/12, 
           label = "Drift Dominated", 
           color = "gray50", fontface = "bold", size = 4, 
           angle = 0) +
  
  # Strong Selection Label (Right side, Vertical)
  annotate("text", x = log10(2), y = y_max_anno/12, 
           label = "Selection", 
           color = "red", fontface = "bold", size = 4, 
           angle = 0) +
  
  # --- Custom "Separated" Rug ---
  # We draw segments explicitly at y = -0.5 (below axis) to create the gap
  geom_segment(aes(x = Log_S_ROC, xend = Log_S_ROC, 
                   y = -0.05, yend = -0.2), # Adjust these values for tick length/position
               alpha = 0.3, color = "darkgreen") +
  
  # --- Scales & Coordinates ---
  scale_y_continuous(
    trans = "log1p", 
    breaks = c(0, 10, 100, 1000, 10000), 
    labels = comma_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  
  # This allows drawing below the axis (where we put the rug)
  coord_cartesian(clip = "off", ylim = c(0, NA)) + 
  
  labs(x = expression(Log[10](S[ROC])), 
       y = "Gene Count (Log1p Scale)") +
  
  theme_custom() +
  # Add margin at bottom to ensure the new rug doesn't get cut off
  theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))

# Panel B: CAI vs S_ROC
p2 <- ggplot(plot_data, aes(x = Log_S_ROC, y = CAI)) +
  geom_point(aes(color = Log_Phi), alpha = 0.6, size = 1) +
  scale_color_viridis_c(option = "plasma", name = expression(Log[10](Phi[geom])), 
                        limits = phi_range, direction = 1) +
  geom_smooth(color = "black") +
  labs(x = expression(Log[10](S[ROC])), y = "CAI") +
  theme_custom()

# Combine
combined_plot <- (p1 | p2) + plot_annotation(tag_levels = 'A')
ggsave("results/Selection_Landscape_Final.pdf", 
       combined_plot, width = 12, height = 6)

# 8.4) GO-enrichment analysis of genes with a massive selection load ----

thr_sel <- 1

subset_strongly_shaped_by_s <- integrated_data |>
  dplyr::filter(S_ROC > thr_sel) |>
  dplyr::pull(Gene_name)

custom_bag <- integrated_data |> dplyr::pull(Gene_name)

GO_results <- gost(query = subset_strongly_shaped_by_s,
                   organism = 'gp__q7VP_EAck_dZk',
                   multi_query = F,
                   significant = T,
                   correction_method = 'fdr',
                   domain_scope = "custom",
                   custom_bg = custom_bag,
                   user_threshold = 0.05)
  
# Export results

write.csv(x = GO_results$result |> dplyr::select(-parents), 
          file = "./results/Go_enrichment.csv", quote = T, 
          row.names = F)

# 8.5) Getting top 20 genes in terms of S_load ----

subset_strongly_shaped_by_s <- integrated_data |>
  dplyr::filter(S_ROC > thr_sel) |>
  dplyr::arrange(desc(S_ROC)) |>
  dplyr::slice(1:20) |>
  dplyr::pull(Gene_name)

detailed_annotation <- read.delim(
  "data/Mguttatusvar_IM767_887_v2.1.annotation_info.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "",  # Don't treat # as comment
  quote = "",         # No quote characters (deflines may have quotes)
  fill = TRUE,        # Handle rows with varying numbers of fields
  na.strings = ""     # Treat empty strings as NA
) |>
  dplyr::select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline) |>
  dplyr::filter(locusName %in% subset_strongly_shaped_by_s) |>
  dplyr::distinct()

# Export information about top 10 genes

write.csv(x = detailed_annotation, 
          file = "./results/Top20_genes_strong_selection_load.csv", 
          quote = T, row.names = F)

# 8.6) Goodness-of-fit test based on Anacoda predictions ----

gof_results <- run_gof_analysis(
  mutation_file  = "./results/MCMC_results_backup/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Mutation.csv",
  selection_file = "./results/MCMC_results_backup/results_dM_fixed_with_phi_final/run_1/Parameter_est/Cluster_1_Selection.csv",
  phi_file       = "./results/MCMC_results_backup/results_dM_fixed_with_phi_final/run_1/Parameter_est/gene_expression.txt",
  codon_counts_long = codon_freq_long,
  test           = "chisq",
  min_aa_total   = 5,
  output_prefix  = "./results/ROC_model_goodness_of_fit"
)

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

## *****************************************************************************
## 10) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________

# 10.1) CA Analysis ---- 

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

# Prepare gene data for biplot
gene_data_ca <- codon_usage_CA_coord |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\n--- CA Analysis: Preferred vs Non-preferred Codons ---\n")
p_ca <- create_preference_biplot(
  ordination_result = codon_usage_CA,
  gene_data = gene_data_ca,
  preferred_codons = preferred_codons_roc,
  dims = c(2, 3),
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

# Prepare gene data for biplot
gene_data_pca <- rscu_PCA_coord |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\n--- PCA Analysis: Preferred vs Non-preferred Codons ---\n")
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

# 10.3) Selection (S_ROC) Based Analysis ----

# Test whether genes under strong vs weak selection show distinct codon patterns

if (exists("selection_metrics") && "S_ROC" %in% names(selection_metrics)) {
  
  cat("\n--- CA/PCA Analysis by Selection Load ---\n")
  
  # Create S_load-based groups
  s_quantiles <- quantile(selection_metrics$S_ROC, 
                          probs = c(0.05, 0.95), na.rm = TRUE)
  
  selection_groups <- selection_metrics |>
    dplyr::mutate(
      S_Group = dplyr::case_when(
        S_ROC >= s_quantiles[2] ~ "High Selection (Top 5%)",
        S_ROC <= s_quantiles[1] ~ "Low Selection (Bottom 5%)",
        TRUE ~ "Intermediate"
      )
    ) |>
    dplyr::select(Gene_name, S_ROC, S_Group)
  
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

## *****************************************************************************
## 11) tRNA abundance correlation analysis ----
## _____________________________________________________________________________

# Analysis 1: By tRNA gene copy number (traditional approach)

tRNA_copynumber_results <- tRNA_codon_correlation(
  codon_counts = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_copynumber",
  test_method = "spearman",
  mode = "by.copy.number"
)

aa_trna_check <- check_aa_frequency_vs_tRNA_supply(
  codon_usage = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/aa_trna_sanity_check"
)

# Analysis 2: Translational Accuracy Hypothesis ----

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

## *****************************************************************************
## 12) Polymorphism data integration ----
## _____________________________________________________________________________

# Pi ----

pi_data <- fread(input = "data/all_chromosomes.bygene.pi.txt")

# Homogenizing gene names to match the previous convention

pi_data <- pi_data |>
  dplyr::select(Chr, Gene, contains("Tajima"), contains("mean")) |>
  dplyr::mutate(Gene = paste0("MgIM767.", pi_data[['Gene']])) |>
  dplyr::rename(Gene_name = Gene)

# Join polymorphism data to integrated_data
integrated_data <- integrated_data |>
  dplyr::left_join(pi_data, by = "Gene_name") |>
  na.exclude()

# Is the relationship between pi and predictors of interest linear?

predictors <- c('Max_Log10_Exp', 'Exp_breadth', 'CDS_length_nt')

# Analysis based on expression
pi_nonlinearity_results <- analyze_nonlinearity_suite(
  resp = "Pi_mean_4fold",
  predictors = predictors,
  data = integrated_data |> dplyr::filter(Exp_breadth > 0),
  family = Gamma(link = "log")
)

pi_models <- fit_codon_gam_suite(
  data = integrated_data |> dplyr::filter(Exp_breadth > 0),
  response_var = "Pi_mean_4fold",
  family = Gamma(link = "log") 
)

pi_selection <- get_model_selection_stats(pi_models)
pi_selection_winner <- pi_selection |> dplyr::filter(AIC == min(AIC))

run_posteriori_gam_analysis(model = pi_models[["Complex"]], 
                            data = integrated_data |> dplyr::filter(Exp_breadth > 0),
                            response_name = "Pi_mean_4fold")

# Analysis based on selection estimates

# Geting the mutational bias (fraction of GC3 that cannot be explained by S_ROC)

confounding_results <- 
  check_confounding_vif(integrated_data, focal_pred = "S_ROC", 
                        comp_pred = "GC")
print(confounding_results$plot)

predictors_s <- c("S_ROC", "Total_Codons", "Exp_breadth")

pi_nonlinearity_results_s <- analyze_nonlinearity_suite(
  resp = "Pi_mean_4fold",
  predictors = predictors_s,
  data = integrated_data,
  family = Gamma(link = "log")
)

pi_models_s <- fit_codon_gam_suite(
  data = integrated_data,
  model_list = alist(
    Null        = ~ 1,
    Additive    = ~ s(S_ROC) + s(Total_Codons) + s(Exp_breadth),
    S_Length    = ~ te(S_ROC, Total_Codons) + s(Exp_breadth),
    S_Breadth   = ~ te(S_ROC, Exp_breadth) + s(Total_Codons),
    Complex     = ~ te(S_ROC, Total_Codons, Exp_breadth)
  ),
  response_var = "Pi_mean_4fold",
  family = Gamma(link = "log") 
)

pi_selection_s <- get_model_selection_stats(pi_models_s)
pi_secection_winner_s <- pi_selection_s |> dplyr::filter(AIC == min(AIC))

run_posteriori_gam_analysis(model = pi_models_s[["S_Length"]], 
                            data = integrated_data,
                            focal_pred = "S_ROC", 
                            interact_pred = "Exp_breadth",
                            third_pred = NULL,
                            response_name = "Pi_mean_4fold",
                            prefix = "SelectionPi")

# 12.1) Tracking frequency of preferred allele as a function of expression ----

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

# GAM wrappers

predictors_p <- c("Max_Log10_Exp", "Exp_breadth", "Total_Codons")

preferred_nonlinearity_results <- analyze_nonlinearity_suite(
  resp = "Mean_preferred_freq",
  predictors = predictors_p,
  data = integrated_data |> dplyr::filter(Exp_breadth > 0),
  family = betar(link = "logit")
)

preferred_models <- fit_codon_gam_suite(
  data = integrated_data |> dplyr::filter(Exp_breadth > 0),
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
                              dplyr::filter(Exp_breadth > 0),
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
if (kw_detrended$p.value < 0.05) {
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

# Ploting box plot

p_boxplot_preferred <- ggplot(integrated_data, aes(x = Expression_Group, 
                                                   y = Mean_preferred_freq_detrended, 
                                                   fill = Expression_Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                               "Bottom 5%" = "#377EB8",
                               "Middle 90%" = "#999999")) +
  labs(y = "Mean Frequency of Preferred Codon",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/Frequency_preferred_by_expression_group.pdf", 
       p_boxplot_preferred, width = 8, height = 6)

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

# Execute for your 16% model
p_surface_pref <- plot_selection_surface(
  model = preferred_models[["Complex"]], 
  data = integrated_data |> dplyr::filter(Exp_breadth > 0),
  response_name = "Mean_preferred_freq"
)

## *****************************************************************************
## 13) Intronic Polymorphism-Based Selection Validation ----
## _____________________________________________________________________________

# STEP 1: Check for pre-computed intronic SFS files ----
sfs_G_file <- "./data/sfs_introns_G.csv"
sfs_C_file <- "./data/sfs_introns_C.csv"

# STEP 2: Estimate neutral mutation parameters ----
neutral_params <- load_and_estimate_neutral_params(sfs_G_file, sfs_C_file)
str(neutral_params)

# Save neutral parameter estimates
neutral_params_df <- data.frame(
  Parameter = c("alpha_G", "beta_G", "alpha_C", "beta_C",
                "pi_G_expected", "pi_C_expected"),
  Value = c(neutral_params$alpha_G, neutral_params$beta_G,
            neutral_params$alpha_C, neutral_params$beta_C,
            neutral_params$pi_G_expected, neutral_params$pi_C_expected),
  Description = c("4N*u for G (unpreferred->preferred)", 
                  "4N*v for G (preferred->unpreferred)",
                  "4N*u for C (unpreferred->preferred)",
                  "4N*v for C (preferred->unpreferred)",
                  "Expected nucleotide diversity at G sites",
                  "Expected nucleotide diversity at C sites")
)

write.csv(neutral_params_df, 
          "./results/neutral_mutation_parameters.csv",
          row.names = FALSE)

cat("✓ Neutral parameters saved: ./results/neutral_mutation_parameters.csv\n\n")

# STEP 3: Calculate expected SFS for C and G ----

# To ensure a common sample size, we get the min value of n between G and C
sfs_C <- read.csv(sfs_C_file) |> dplyr:: filter(n > 90)
sfs_G <- read.csv(sfs_G_file) |> dplyr:: filter(n > 90)

target_n <- 90

obs_sfs_G <- project_sfs(sfs_G, target_n)
obs_sfs_C <- project_sfs(sfs_C, target_n)

observed_list <- list(G = obs_sfs_G, C = obs_sfs_C)

# Generate null expectations

expected_sfs <- generate_expected_counts(
  neutral_param = neutral_params, 
  observed_sfs_list = observed_list, 
  target_n = target_n
)

sfs_contrast <- data.frame(Num_seq = seq(0, target_n, 1),
                           Expected_C = expected_sfs$C,
                           Expected_G = expected_sfs$G) |>
  tidyr::pivot_longer(cols = c('Expected_C', 'Expected_G'), 
                      names_to = "Metric", values_to = "Counts")

sfs_expected_counts <- sfs_contrast |>
  dplyr::group_by(Metric) |>
  dplyr::summarise(
    Expected_p = sum(Num_seq * Counts) / sum(Counts)
  )
  
# Visualize null expectations

ggplot(sfs_contrast, aes(x = Num_seq, y = Counts, fill = Metric)) + # Changed 'color' to 'fill'
  geom_col(position = "dodge", width = 0.8) + # 'dodge' places them side-by-side; 'width' adds spacing
  labs(title = "Expected SFS for Intronic C and G Sites",
       x = "Number of Sequences with Target Allele",
       y = "Expected Count") +
  scale_fill_manual(values = c("Expected_C" = "blue", 
                               "Expected_G" = "red"),
                    labels = c("C Sites", "G Sites"),
                    name = "Site Type") +
  scale_y_log10() +
  theme_custom()

ggsave("./results/diversity_modeling/Expected_SFS_C_G.pdf",
       width = 6, height = 4)

# STEP 4: Calculate the observed SFS for C and G at third positions ----

# --- 1. Define Gene Sets ---
target_n <- 90  # Define the projection size (must match your intron analysis)

# Split the by quantiles in 20 intervals (to match with average S_roc per quantile and
# determine scaling factor)

cutoffs_exp <- quantile(integrated_data$Geom_Mean_CPM, seq(0, 1, by = 0.05),
                    na.rm = TRUE)

names_quantiles <- names(cutoffs_exp)
interval_names <- paste0(names_quantiles[-length(names_quantiles)], "-", 
                         names_quantiles[-1])

empirical_SFS <- future_lapply(X = 1:(length(cutoffs_exp) - 1), 
                           FUN = function(idx)
{
  low_thr <- ifelse(idx == 1, cutoffs_exp[idx] - 0.001, 
                    cutoffs_exp[idx]) # This ensure that the lowest value is included in subset
  high_thr <- cutoffs_exp[idx + 1]
  
  # Subset the original data as a function of the stablished thresholds
  sub_data_genes <- integrated_data |>
    dplyr::filter(Geom_Mean_CPM > low_thr & Geom_Mean_CPM <= high_thr) |>
    pull(Gene_name)
  
  # Extreact the SFS for the subset of genes
  SFS_sub <- process_gene_set_sfs(sub_data_genes, 
                                  interval_names[idx], 
                                  target_n)
  
  freq_bins <- 0:target_n
  
  obs_sfs_G <- SFS_sub$obs_sfs_G
  obs_sfs_C <- SFS_sub$obs_sfs_C
  
  # Build a data frame with all the information required for plotting 
  # NOTE: Order must match - G data with G label, C data with C label
  sfs_observed_df <- data.frame(
    Num_seq = rep(freq_bins, 2),
    Counts = c(obs_sfs_G, obs_sfs_C),
    Metric = rep(c(paste0("Observed_", interval_names[idx], "_G"), 
                   paste0("Observed_", interval_names[idx], "_C")), 
                 each = length(freq_bins))
  )
  
  # Generate the null expectations
  expected_sfs <- generate_expected_counts(
    neutral_param = neutral_params,
    observed_sfs_list = list(G = obs_sfs_G, C = obs_sfs_C), 
    target_n = target_n
  )
  
  sfs_neutral <- data.frame(
    Num_seq = rep(freq_bins, 2),
    Counts = c(expected_sfs$G, expected_sfs$C),
    Metric = rep(c(paste0("Expected_", interval_names[idx], "_G"), 
                   paste0("Expected_", interval_names[idx], "_C")), 
                 each = length(freq_bins))
  )
  
  result <- list(Observed_df = sfs_observed_df,
                 Neutral_df = sfs_neutral)
  
  return(result)
})

names(empirical_SFS) <- interval_names

# --- 7. Estimate Gamma for G and C in all groups ---
cat("\n=== Estimating Gamma (Selection Coefficient) ===\n")

gamma_estimates <- future_lapply(X = 1:length(empirical_SFS), 
                                 FUN = function(idx)
{
   # Extract the sub-list (each sub-list is an interval)
   sublist <- empirical_SFS[[idx]]
   
   # Extract relevant counts
   sfs_C <- sublist$Observed_df[grepl(pattern = "_C", 
                                     x = sublist$Observed_df$Metric), "Counts"]
   sfs_G <- sublist$Observed_df[grepl(pattern = "_G", 
                                      x = sublist$Observed_df$Metric), "Counts"]
   
   # Estimate Gamma
   Gamma_C <- estimate_gamma(sfs_C, neutral_params, "C", target_n)
   Gamma_G <- estimate_gamma(sfs_G, neutral_params, "G", target_n)
   
   result <- list(Gamma_C = Gamma_C,
                  Gamma_G = Gamma_G)
   
   return(result)
})

names(gamma_estimates) <- names(empirical_SFS)

# Extracting gamma estimates per category and nucleotide (with p-values)

gamma_summary_df <- data.frame(Category = names(gamma_estimates))

gamma_summary_df[['Gamma_C']] <- sapply(gamma_summary_df[['Category']],
                                      function(c)
{
  gamma_estimates[[c]][['Gamma_C']][['gamma']]
})

gamma_summary_df[['Gamma_C_pval']] <- sapply(gamma_summary_df[['Category']],
                                             function(c)
{
  gamma_estimates[[c]][['Gamma_C']][['p_value']]
})

gamma_summary_df[['Gamma_G']] <- sapply(gamma_summary_df[['Category']],
                                        function(c)
                                        {
                                          gamma_estimates[[c]][['Gamma_G']][['gamma']]
                                        })

gamma_summary_df[['Gamma_G_pval']] <- sapply(gamma_summary_df[['Category']],
                                             function(c)
                                             {
                                               gamma_estimates[[c]][['Gamma_G']][['p_value']]
                                             })

gamma_summary_df[['Gamma_avg']] <- (gamma_summary_df[['Gamma_C']] +
                                        gamma_summary_df[['Gamma_G']]) / 2

# Flag significant departures from neutrality
gamma_summary_df[['Sig_C']] <- gamma_summary_df$Gamma_C_pval < 0.05
gamma_summary_df[['Sig_G']] <- gamma_summary_df$Gamma_G_pval < 0.05

cat("\n=== Gamma Estimates by Expression Quantile ===\n")
print(gamma_summary_df[, c("Category", "Gamma_C", "Gamma_C_pval", "Gamma_G", "Gamma_G_pval", "Gamma_avg")])

cutoffs_exp <- quantile(integrated_data$Max_Log10_Exp, seq(0, 1, by = 0.05),
                        na.rm = TRUE)

names_quantiles <- names(cutoffs_exp)
interval_names <- paste0(names_quantiles[-length(names_quantiles)], "-", 
                         names_quantiles[-1])

# Compute both mean and SE for S_ROC
selection_summary_full <- lapply(1:(length(cutoffs_exp) - 1), function(idx)
{
  low_thr <- ifelse(idx == 1, cutoffs_exp[idx] - 0.001, 
                    cutoffs_exp[idx])
  high_thr <- cutoffs_exp[idx + 1]
  
  sub_data_genes <- integrated_data |>
    dplyr::filter(Max_Log10_Exp > low_thr & Max_Log10_Exp <= high_thr)
  
  list(
    mean_S = mean(sub_data_genes$S_ROC, na.rm = TRUE),
    se_S = sd(sub_data_genes$S_ROC, na.rm = TRUE) / sqrt(sum(!is.na(sub_data_genes$S_ROC))),
    n_genes = nrow(sub_data_genes),
    mean_exp = mean(sub_data_genes$Max_Log10_Exp, na.rm = TRUE)
  )
})

selection_summary <- sapply(selection_summary_full, `[[`, "mean_S")

cat("✓ Step 4 Complete: Data parsed and ready for plotting.\n")

scaling_df <- gamma_summary_df
scaling_df$S_ROC <- selection_summary
scaling_df$S_ROC_SE <- sapply(selection_summary_full, `[[`, "se_S")
scaling_df$N_genes <- sapply(selection_summary_full, `[[`, "n_genes")
scaling_df$Mean_Exp <- sapply(selection_summary_full, `[[`, "mean_exp")

# Add expression quantile midpoint for plotting
scaling_df$Exp_Quantile <- seq(2.5, 97.5, by = 5)

# Save the full scaling dataframe
write.csv(scaling_df, "./results/gamma_vs_Sroc_scaling_data.csv", row.names = FALSE)
cat("✓ Scaling data saved: ./results/gamma_vs_Sroc_scaling_data.csv\n")

# Calculate pi for each expression quantile using estimated gamma values
cat("\n=== Calculating Expected Pi for C and G Ending Codons ===\n")

scaling_df$Pi_C <- sapply(seq_len(nrow(scaling_df)), function(i) {
  pi_estimator_eq(
    alpha = neutral_params$alpha_C,
    beta = neutral_params$beta_C,
    gamma = scaling_df$Gamma_C[i]
  )
})

scaling_df$Pi_G <- sapply(seq_len(nrow(scaling_df)), function(i) {
  pi_estimator_eq(
    alpha = neutral_params$alpha_G,
    beta = neutral_params$beta_G,
    gamma = scaling_df$Gamma_G[i]
  )
})

scaling_df$Pi_avg <- (scaling_df$Pi_C + scaling_df$Pi_G) / 2

# Also calculate neutral (gamma=0) pi for reference
pi_neutral_C <- pi_estimator_eq(neutral_params$alpha_C, neutral_params$beta_C, 0)
pi_neutral_G <- pi_estimator_eq(neutral_params$alpha_G, neutral_params$beta_G, 0)

cat(sprintf("Neutral Pi (gamma=0): C = %.6f, G = %.6f\n", pi_neutral_C, pi_neutral_G))
cat(sprintf("Pi range across quantiles: C = [%.6f, %.6f], G = [%.6f, %.6f]\n",
            min(scaling_df$Pi_C), max(scaling_df$Pi_C),
            min(scaling_df$Pi_G), max(scaling_df$Pi_G)))

# Update saved scaling data
write.csv(scaling_df, "./results/gamma_vs_Sroc_scaling_data.csv", row.names = FALSE)
cat("✓ Updated scaling data saved with Pi estimates\n")

# STEP 5: Plot observed vs expected SFS for extreme expression groups ----
cat("\n=== STEP 5: Observed vs Expected SFS (Bottom 5% vs Top 5% Expression) ===\n")

# Extract SFS data for bottom 5% (first interval) and top 5% (last interval)
bottom5_interval <- "0%-5%"
top5_interval <- "95%-100%"

# Get the gamma values from scaling_df
gamma_bottom5_C <- scaling_df$Gamma_C[1]
gamma_bottom5_G <- scaling_df$Gamma_G[1]
gamma_top5_C <- scaling_df$Gamma_C[nrow(scaling_df)]
gamma_top5_G <- scaling_df$Gamma_G[nrow(scaling_df)]

# Build combined data frame for plotting
build_sfs_plot_data <- function(empirical_SFS, interval_name, expression_label, 
                                 gamma_C, gamma_G, neutral_params, target_n) {
  
  sfs_data <- empirical_SFS[[interval_name]]
  
  # Extract observed and neutral counts
  obs_df <- sfs_data$Observed_df
  neutral_df <- sfs_data$Neutral_df
  
  # The Metric column format is "Observed_{interval}_C" or "Observed_{interval}_G"
  # Use pattern matching to find C and G
  obs_C <- obs_df$Counts[grepl("_C$", obs_df$Metric)]
  obs_G <- obs_df$Counts[grepl("_G$", obs_df$Metric)]
  neutral_C <- neutral_df$Counts[grepl("_C$", neutral_df$Metric)]
  neutral_G <- neutral_df$Counts[grepl("_G$", neutral_df$Metric)]
  
  # Combine into plot-ready format
  freq_bins <- 0:target_n
  
  # For C sites
  plot_df_C <- data.frame(
    Frequency = freq_bins,
    Observed = obs_C,
    Neutral = neutral_C,
    Nucleotide = "C",
    Expression_Group = expression_label,
    Gamma = gamma_C
  )
  
  # For G sites
  plot_df_G <- data.frame(
    Frequency = freq_bins,
    Observed = obs_G,
    Neutral = neutral_G,
    Nucleotide = "G",
    Expression_Group = expression_label,
    Gamma = gamma_G
  )
  
  rbind(plot_df_C, plot_df_G)
}

# Build data for both expression groups
sfs_plot_bottom5 <- build_sfs_plot_data(
  empirical_SFS, bottom5_interval, "Bottom 5%",
  gamma_bottom5_C, gamma_bottom5_G, neutral_params, target_n
)

sfs_plot_top5 <- build_sfs_plot_data(
  empirical_SFS, top5_interval, "Top 5%",
  gamma_top5_C, gamma_top5_G, neutral_params, target_n
)

sfs_plot_combined <- rbind(sfs_plot_bottom5, sfs_plot_top5)

# Factor ordering
sfs_plot_combined$Expression_Group <- factor(
  sfs_plot_combined$Expression_Group,
  levels = c("Bottom 5%", "Top 5%")
)

# Create gamma annotation labels (using "gamma" instead of Unicode for PDF compatibility)
gamma_labels <- sfs_plot_combined |>
  dplyr::group_by(Nucleotide, Expression_Group) |>
  dplyr::summarize(
    Gamma = unique(Gamma),
    label = sprintf("gamma = %.2f", Gamma),
    .groups = "drop"
  )

# Verify that observed and neutral totals match (they should be scaled to same total)
cat("\n=== Verifying Observed vs Neutral Scaling ===\n")

# Full totals (including endpoints - should match)
verify_full <- sfs_plot_combined |>
  dplyr::group_by(Nucleotide, Expression_Group) |>
  dplyr::summarize(
    Total_Observed = sum(Observed),
    Total_Neutral = sum(Neutral),
    Ratio_Full = sum(Observed) / sum(Neutral),
    .groups = "drop"
  )
cat("\nFull SFS (including fixed sites at 0 and n) - should have Ratio ~ 1.0:\n")
print(verify_full)

# Polymorphic only (excluding endpoints - will differ due to SFS shape)
verify_poly <- sfs_plot_combined |>
  dplyr::filter(Frequency > 0, Frequency < target_n) |>
  dplyr::group_by(Nucleotide, Expression_Group) |>
  dplyr::summarize(
    Poly_Observed = sum(Observed),
    Poly_Neutral = sum(Neutral),
    Ratio_Poly = sum(Observed) / sum(Neutral),
    .groups = "drop"
  )
cat("\nPolymorphic sites only (0 < freq < n) - ratio > 1 indicates selection signature:\n")
print(verify_poly)

# Pivot for easier plotting (observed as bars, neutral as line)
sfs_long <- sfs_plot_combined |>
  tidyr::pivot_longer(
    cols = c(Observed, Neutral),
    names_to = "Type",
    values_to = "Count"
  )

# Keep ALL frequency bins including monomorphic sites (0 and n)
# This ensures observed and neutral totals match
sfs_all <- sfs_long

# Also create polymorphic-only subset for ratio analysis
sfs_poly <- sfs_long |>
  dplyr::filter(Frequency > 0, Frequency < target_n)

# Create the main SFS comparison plot (INCLUDING monomorphic sites)
p_sfs_comparison <- ggplot(sfs_all |> dplyr::filter(Type == "Observed"),
                           aes(x = Frequency, y = Count + 1)) +  # +1 to handle zeros on log scale
  # Observed as bars
  geom_col(aes(fill = Expression_Group), alpha = 0.7, position = "identity") +
  # Neutral expectation as blue line
  geom_line(data = sfs_all |> dplyr::filter(Type == "Neutral"),
            aes(x = Frequency, y = Count + 1),
            color = "blue", linewidth = 1.2, linetype = "solid") +
  geom_point(data = sfs_all |> dplyr::filter(Type == "Neutral"),
             aes(x = Frequency, y = Count + 1),
             color = "blue", size = 1.5) +
  # Add gamma annotations
  geom_text(data = gamma_labels,
            aes(x = target_n * 0.5, y = Inf, label = label),
            vjust = 2, hjust = 0.5, size = 4, fontface = "bold") +
  # Facet by nucleotide and expression group
  facet_grid(Nucleotide ~ Expression_Group, scales = "free_y") +
  # Styling
  scale_fill_manual(values = c("Bottom 5%" = "#377EB8", "Top 5%" = "#E41A1C")) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Site Frequency Spectrum: Observed vs Neutral Expectation (Full SFS)",
    subtitle = "Blue line = Neutral (gamma=0); Bars = Observed; Includes monomorphic sites (freq 0 and n)",
    x = "Derived Allele Count",
    y = "Number of Sites (log scale, +1)",
    fill = "Expression\nGroup"
  ) +
  theme_custom() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("./results/diversity_modeling/SFS_observed_vs_neutral_by_expression_FULL.pdf",
       plot = p_sfs_comparison, width = 10, height = 8)

cat("✓ Full SFS plot saved: ./results/diversity_modeling/SFS_observed_vs_neutral_by_expression_FULL.pdf\n")

# Also create polymorphic-only plot for comparison
p_sfs_poly <- ggplot(sfs_poly |> dplyr::filter(Type == "Observed"),
                     aes(x = Frequency, y = Count)) +
  geom_col(aes(fill = Expression_Group), alpha = 0.7, position = "identity") +
  geom_line(data = sfs_poly |> dplyr::filter(Type == "Neutral"),
            aes(x = Frequency, y = Count),
            color = "blue", linewidth = 1.2, linetype = "solid") +
  geom_point(data = sfs_poly |> dplyr::filter(Type == "Neutral"),
             aes(x = Frequency, y = Count),
             color = "blue", size = 1.5) +
  geom_text(data = gamma_labels,
            aes(x = target_n * 0.5, y = Inf, label = label),
            vjust = 2, hjust = 0.5, size = 4, fontface = "bold") +
  facet_grid(Nucleotide ~ Expression_Group, scales = "free_y") +
  scale_fill_manual(values = c("Bottom 5%" = "#377EB8", "Top 5%" = "#E41A1C")) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Site Frequency Spectrum: Polymorphic Sites Only (0 < freq < n)",
    subtitle = "Blue line = Neutral; Bars = Observed; Ratio > 1 indicates selection signature",
    x = "Derived Allele Count",
    y = "Number of Sites (log scale)",
    fill = "Expression\nGroup"
  ) +
  theme_custom() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("./results/diversity_modeling/SFS_observed_vs_neutral_by_expression_POLY.pdf",
       plot = p_sfs_poly, width = 10, height = 8)

cat("✓ SFS comparison plot saved: ./results/diversity_modeling/SFS_observed_vs_neutral_by_expression.pdf\n")

# Create a summary table of the differences
sfs_summary <- sfs_plot_combined |>
  dplyr::filter(Frequency > 0, Frequency < target_n) |>
  dplyr::group_by(Nucleotide, Expression_Group) |>
  dplyr::summarize(
    Total_Observed = sum(Observed),
    Total_Neutral = sum(Neutral),
    Obs_Neutral_Ratio = sum(Observed) / sum(Neutral),
    Mean_Freq_Observed = weighted.mean(Frequency, Observed),
    Mean_Freq_Neutral = weighted.mean(Frequency, Neutral),
    Gamma = unique(Gamma),
    .groups = "drop"
  )

cat("\n=== SFS Summary Statistics ===\n")
print(sfs_summary)

# STEP 6: Mutual Information Analysis (S_ROC vs Gamma) ----
cat("\n=== STEP 6: Mutual Information Analysis ===\n")
cat("Quantifying non-linear association between S_ROC and Gamma\n\n")

# Data Cleaning & Type Conversion ---
# Ensure variables are pure numeric vectors (not lists) and handle Infinite values
scaling_df$S_ROC <- as.numeric(unlist(scaling_df$S_ROC))
scaling_df$Gamma_avg <- as.numeric(unlist(scaling_df$Gamma_avg))

# Convert Infinite values (possible from log(0)) to NA
is.na(scaling_df$S_ROC) <- !is.finite(scaling_df$S_ROC)
is.na(scaling_df$Gamma_avg) <- !is.finite(scaling_df$Gamma_avg)

compute_mi_analysis <- function(x, y, n_bins = 10) {
  #' Compute mutual information between two continuous variables
  #' Uses discretization with equal-frequency binning
  
  # Remove NA values
  valid_idx <- complete.cases(x, y)
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  n <- length(x)
  
  # Equal-frequency binning (quantile-based)
  x_breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1)))
  y_breaks <- unique(quantile(y, probs = seq(0, 1, length.out = n_bins + 1)))
  
  # If data is too sparse for requested bins, return NA or handle gracefully
  if(length(x_breaks) < n_bins + 1 || length(y_breaks) < n_bins + 1) {
    return(list(MI=NA, NMI=NA, H_X=NA, H_Y=NA, U_Y_given_X=NA, U_X_given_Y=NA))
  }
  
  x_bins <- cut(x, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  y_bins <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  
  # Compute joint and marginal probabilities
  joint_table <- table(x_bins, y_bins)
  p_xy <- joint_table / sum(joint_table)
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  
  # Compute mutual information
  mi <- 0
  for (i in seq_along(p_x)) {
    for (j in seq_along(p_y)) {
      if (p_xy[i, j] > 0 && p_x[i] > 0 && p_y[j] > 0) {
        mi <- mi + p_xy[i, j] * log2(p_xy[i, j] / (p_x[i] * p_y[j]))
      }
    }
  }
  
  # Compute entropy for normalization
  H_x <- -sum(p_x[p_x > 0] * log2(p_x[p_x > 0]))
  H_y <- -sum(p_y[p_y > 0] * log2(p_y[p_y > 0]))
  
  # Normalized MI (ranges 0-1)
  nmi <- ifelse(H_x * H_y > 0, mi / sqrt(H_x * H_y), 0)
  
  return(list(
    MI = mi,
    NMI = nmi,
    H_X = H_x,
    H_Y = H_y,
    U_Y_given_X = ifelse(H_y > 0, mi / H_y, 0),
    U_X_given_Y = ifelse(H_x > 0, mi / H_x, 0),
    n_samples = n,
    n_bins = n_bins
  ))
}

# We calculate for different bins, but note that for N=20 points, high bins are unreliable
mi_results <- lapply(c(3, 4, 5, 8), function(nb) {
  res <- compute_mi_analysis(scaling_df$S_ROC, scaling_df$Gamma_avg, n_bins = nb)
  res$n_bins <- nb
  res
})

cat("=== Mutual Information Results (Sensitivity) ===\n")
for (res in mi_results) {
  if(!is.na(res$MI)) {
    cat(sprintf("  Bins=%2d: MI=%.4f bits, NMI=%.4f, U(Gamma|S)=%.3f\n", 
                res$n_bins, res$MI, res$NMI, res$U_Y_given_X))
  }
}

# For 20 data points (quantiles), 5 bins is the statistical maximum 
# (approx 4 points per bin). Anything higher causes overfitting (NMI -> 1.0).
optimal_bins <- 5

mi_optimal <- compute_mi_analysis(scaling_df$S_ROC, scaling_df$Gamma_avg, 
                                  n_bins = optimal_bins)

cat(sprintf("\nOptimal Binning (Robust N=%d):\n", optimal_bins))
cat(sprintf("  Mutual Information: %.4f bits\n", mi_optimal$MI))
cat(sprintf("  Normalized MI:      %.4f (0=independent, 1=perfect dependence)\n", mi_optimal$NMI))
cat(sprintf("  U(Gamma|S_ROC):     %.3f (fraction of Gamma entropy explained by S_ROC)\n", mi_optimal$U_Y_given_X))

# Permutation test for significance
cat("=== Permutation Test for MI Significance ===\n")
n_permutations <- 1000
observed_mi <- mi_optimal$MI

permuted_mi <- replicate(n_permutations, {
  y_shuffled <- sample(scaling_df$Gamma_avg)
  # Suppress warnings during permutation (sparse bins in random shuffles are common)
  suppressWarnings(compute_mi_analysis(scaling_df$S_ROC, y_shuffled, n_bins = optimal_bins)$MI)
})

# Handle potential NAs in permutations
permuted_mi <- permuted_mi[!is.na(permuted_mi)]

mi_p_value <- mean(permuted_mi >= observed_mi)
mi_zscore <- (observed_mi - mean(permuted_mi)) / sd(permuted_mi)

cat(sprintf("Observed MI: %.4f bits\n", observed_mi))
cat(sprintf("Permuted MI: %.4f ± %.4f bits (mean ± SD)\n", mean(permuted_mi), sd(permuted_mi)))
cat(sprintf("Z-score:     %.2f\n", mi_zscore))
cat(sprintf("P-value:     %.4f (from %d valid permutations)\n", mi_p_value, length(permuted_mi)))

if (mi_p_value < 0.05) {
  cat("RESULT: Significant non-random association between S_ROC and Gamma\n\n")
} else {
  cat("RESULT: Association not significantly different from random\n\n")
}

# Compare linear vs non-linear information
# Use use="complete.obs" to handle any NAs introduced by the infinite check
spearman_cor <- cor(scaling_df$S_ROC, scaling_df$Gamma_avg, method = "spearman", use = "complete.obs")
pearson_cor <- cor(scaling_df$S_ROC, scaling_df$Gamma_avg, method = "pearson", use = "complete.obs")

# Create visualization of the MI analysis
p_mi_heatmap <- ggplot(scaling_df, aes(x = S_ROC, y = Gamma_avg)) +
  stat_bin2d(bins = optimal_bins, aes(fill = after_stat(density))) +
  scale_fill_viridis_c(name = "Density", option = "plasma") +
  geom_smooth(method = "lm", color = "white", linetype = "dashed", se = FALSE) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(
    title = "Joint Distribution: S_ROC vs Gamma",
    subtitle = sprintf("MI=%.3f bits, NMI=%.3f, Pearson R²=%.3f", 
                       mi_optimal$MI, mi_optimal$NMI, pearson_cor^2),
    x = "Translational Selection (AnaCoDa S_ROC)",
    y = "Polymorphism-Based Gamma (4Nes)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("./results/MI_heatmap_Sroc_vs_Gamma.pdf", plot = p_mi_heatmap, width = 8, height = 6)

# Permutation distribution plot
permutation_df <- data.frame(MI = permuted_mi)

p_permutation <- ggplot(permutation_df, aes(x = MI)) +
  geom_histogram(bins = 30, fill = "gray70", color = "gray40", alpha = 0.7) +
  geom_vline(xintercept = observed_mi, color = "red", linewidth = 1.5, linetype = "solid") +
  annotate("text", x = observed_mi, y = Inf, label = sprintf("Observed\n(p=%.3f)", mi_p_value),
           color = "red", hjust = -0.1, vjust = 1.5, fontface = "bold") +
  labs(
    title = "Permutation Test: Mutual Information Significance",
    x = "Mutual Information (bits)",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("./results/MI_permutation_test.pdf", plot = p_permutation, width = 8, height = 5)

# Save summary statistics
mi_summary <- data.frame(
  Metric = c("Mutual_Information_bits", "Normalized_MI", "Entropy_S_ROC", "Entropy_Gamma",
             "U_Gamma_given_S", "U_S_given_Gamma", "Pearson_R", "Spearman_rho",
             "MI_pvalue", "MI_zscore", "Nonlinearity_index"),
  Value = c(mi_optimal$MI, mi_optimal$NMI, mi_optimal$H_X, mi_optimal$H_Y,
            mi_optimal$U_Y_given_X, mi_optimal$U_X_given_Y, pearson_cor, spearman_cor,
            mi_p_value, mi_zscore, mi_optimal$NMI - pearson_cor^2),
  Description = c("Information shared between S_ROC and Gamma",
                  "MI normalized by geometric mean of entropies (0-1 scale)",
                  "Entropy of S_ROC distribution",
                  "Entropy of Gamma distribution",
                  "Fraction of Gamma uncertainty reduced by knowing S_ROC",
                  "Fraction of S_ROC uncertainty reduced by knowing Gamma",
                  "Linear correlation coefficient",
                  "Monotonic correlation coefficient",
                  "P-value from permutation test",
                  "Z-score relative to null distribution",
                  "NMI - R² (indicates non-linear component)")
)

write.csv(mi_summary, "./results/MI_analysis_summary.csv", row.names = FALSE)

## *****************************************************************************
## 14) Diversity across different genomic compartment ----
## _____________________________________________________________________________

pi_compartment <- read.table(file = "./data/all_chromosomes.pi_by_compartment.txt",
                             header = T)

# 1. Sanity Check: Calculate Overall Weighted Pi per Compartment
# We filter for "all" nucleotides to avoid double counting C/G/AT breakdowns
overall_pi_stats <- pi_compartment %>%
  dplyr::filter(Nuc_Category == "all") %>%
  group_by(Compartment) %>%
  summarise(
    Total_Sites = sum(Sites, na.rm = TRUE),
    Total_Pi_Sum = sum(Pi_sum, na.rm = TRUE),
    # Weighted Average Pi = (Sum of all Pi differences) / (Total Sites)
    Weighted_Mean_Pi = Total_Pi_Sum / Total_Sites,
    # Standard Error of the mean (optional, treat chromosomes as replicates)
    SE_Pi = sd(Pi_mean) / sqrt(n()) 
  )

print("=== Overall Weighted Pi by Compartment ===")
print(overall_pi_stats)

# Visual of pi per compartment in each chromosome

# 1. Define a Logical Biological Order
# From "Most Neutral" to "Most Constrained"
compartment_order <- c(
  "intergenic", 
  "intergenic_upstream_10kb", 
  "intergenic_upstream_2kb", 
  "intron", 
  "nonfirst_exon_4fold", 
  "first_exon_4fold", 
  "exon_all"
)

# 2. Prepare the Data
plot_data <- pi_compartment |>
  # Remove the 'all' aggregate so we can see the components clearly
  dplyr::filter(Nuc_Category %in% c("AT", "C", "G")) |> 
  dplyr::mutate(
    # Apply the logical order
    Compartment = factor(Compartment, levels = compartment_order),
    # Ensure Nucleotides are ordered for consistency (AT = Baseline)
    Nuc_Category = factor(Nuc_Category, levels = c("AT", "C", "G"))
  )

# 3. Create the Plot
pi_compart <- ggplot(plot_data, aes(x = Compartment, y = Pi_mean, fill = Nuc_Category)) +
  
  # A. Boxplots side-by-side (Dodged)
  # outlier.shape = NA hides the "duplicate" points since we add jitter below
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               alpha = 0.7) +
  
  # B. Jittered Points (The Chromosomes)
  # This shows the actual variance across the genome
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
             size = 1.2, alpha = 0.6, color = "black") +
  
  # C. Colors
  # Use Gray for AT (Background) and colors for C/G (Active)
  scale_fill_manual(values = c("AT" = "gray80", 
                               "C" = "#E41A1C", # Red (High Mutation)
                               "G" = "#377EB8"), # Blue
                    name = "Nucleotide") +
  
  # D. Aesthetics
  labs(title = "Nucleotide Diversity (Pi) by Genomic Compartment",
       subtitle = "Separation of C/G hypermutability from AT background",
       y = "Mean Nucleotide Diversity",
       x = NULL) +
  theme_custom() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "top"
  )

ggsave("./results/diversity_boxplot_improved.pdf",
       pi_compart, width = 10, height = 6)

# HUMP EFFECT TEST ----
# 1. Aggregate Data by "Selection Potential"
# We group C, G, and CG as "GC_Segregating" (Where selection acts)
# We keep AT as "AT_Only" (Where selection is absent/invisible)

hump_test_data <- pi_compartment |>
  dplyr::filter(Compartment %in% c("nonfirst_exon_4fold", "intron")) |>
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

# Create the summary dataframe directly from your results
plot_data_hump <- data.frame(
  Compartment = factor(c("Intron", "Intron", "Exon (4-fold)", "Exon (4-fold)"),
                       levels = c("Intron", "Exon (4-fold)")),
  Site_Type = c("AT_Only", "GC_Segregating", "AT_Only", "GC_Segregating"),
  Weighted_Pi = c(0.00504, 0.0244, 0.00651, 0.0365)
)

# Plot: The Interaction Effect
p_hump <- ggplot(plot_data_hump, aes(x = Compartment, y = Weighted_Pi, group = Site_Type, color = Site_Type)) +
  
  # Lines connecting the points highlight the differing slopes
  geom_line(size = 1.2) +
  geom_point(size = 5) +
  
  # Formatting
  scale_color_manual(values = c("AT_Only" = "gray60", "GC_Segregating" = "firebrick")) +
  labs(title = "Evidence for Weak Selection (The 'Hump' Effect)",
       subtitle = "Selection opposing Mutational Bias boosts diversity specifically at GC sites",
       y = "Nucleotide Diversity (Pi)",
       x = NULL,
       color = "Site Category") +
  theme_custom() +
  theme(axis.text = element_text(size = 12, face = "bold"),
        legend.position = "top")

ggsave("./results/Hump_Hypothesis_Confirmation.pdf", p_hump, width = 8, height = 6)

# Statistical test

# 1. Prepare Data with Aggregation step
paired_test_data <- pi_compartment |>
  dplyr::filter(Compartment %in% c("nonfirst_exon_4fold", "intron")) |>
  
  # Create the new categories
  dplyr::mutate(Site_Type = ifelse(Nuc_Category == "AT", "AT_Only", "GC_Segregating")) |>
  
  # CRITICAL FIX: Aggregate the multiple GC rows (C, G, CG) into one value per Chromosome
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

# ******************************************************************************
# 15) Testing the translational ramp hypothesis ----
# ______________________________________________________________________________

# 15.1) Reference based analysis (historic) ----
binary_preferred <- codons_to_preferred_state_bernoulli(trans, 
                                                        as.character(preferred_codons$Codon))

# 15.2) Polymorphism based (contemporaneous) ----