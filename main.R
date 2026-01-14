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
                        'FSA')

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
## 3) Load the data ----
## _____________________________________________________________________________

## 3.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnlyClean.fa", 
                                      format = 'fasta')

trans <- trans[check_canonical_start(trans)] |> check_cds()

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## 3.2) Clean the codon usage object (remove the STOP codon, alongside Trp and Met) ----

# codon_usage <- codon_usage |>
#   trim_uninformative(genetic_code = genetic_code_dna_long)
#   >> check_cds already does that

## 3.3) Load the expression data ----

# ORIGINAL (2 sources only)
# exp_data_bud <- read.table(file = "./data/bud_gene_expression_cpm_remapped.txt",
#                            header = T) |>
#   dplyr::rename(Exp_bud = Expression)
# 
# exp_data_leaf <- read.table(file = "./data/leaf_gene_expression_mean_cpm_renamed.txt",
#                             header = T) |>
#   dplyr::rename(Exp_leaf = Expression)

# Combining CUB metric with expression profiles for buds and leafs

# exp_complete <- dplyr::full_join(exp_data_leaf, exp_data_bud, by = "Gene")

# exp_complete <- exp_complete |> 
#   dplyr::rowwise() |>
#   dplyr::mutate(
#     High_exp = max(Exp_leaf, Exp_bud, na.rm = TRUE),
#     
#     Source_High_exp = case_when(
#       # If both are NA, source is NA
#       is.na(Exp_leaf) & is.na(Exp_bud) ~ NA_character_,
#       
#       # If bud is NA, leaf must be the max
#       is.na(Exp_bud) ~ "Leaf",
#       
#       # If leaf is NA, bud must be the max
#       is.na(Exp_leaf) ~ "Bud",
#       
#       # If leaf is greater or equal (handles ties)
#       Exp_leaf >= Exp_bud ~ "Leaf",
#       
#       # Otherwise, bud must be greater
#       Exp_bud > Exp_leaf ~ "Bud"
#     )
#   ) |>
#   # Fix the -Inf from max(NA, NA, na.rm=T)
#   dplyr::mutate(
#     High_exp = if_else(is.infinite(High_exp), NA_real_, High_exp)
#   ) |>
#   dplyr::ungroup() |>
#   dplyr::mutate(Source_High_exp = as.factor(Source_High_exp)) |>
#   na.exclude()

# ALTERNATIVE (multi-source)



## *****************************************************************************
## 4) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results",
                          aa_group = aa_chemistry_df)

# Creation of integrated data ----

integrated_data <- dplyr::left_join(exp_complete, cub_results$enc_results, 
                                    by = dplyr::join_by(Gene == Gene_name)) |>
  dplyr::rename(Gene_name = Gene) |>
  na.omit() |>
  distinct(Gene_name, .keep_all = TRUE)

# Add gene length (CDS length in codons and nucleotides)
codon_columns <- names(codon_usage)[names(codon_usage) != "Gene_name"]

gene_lengths <- codon_usage |>
  dplyr::mutate(
    Total_Codons = rowSums(across(all_of(codon_columns)), na.rm = TRUE),
    CDS_length_nt = Total_Codons * 3,  # nucleotides
    CDS_length_aa = Total_Codons        # amino acids (codons)
  ) |>
  dplyr::select(Gene_name, Total_Codons, CDS_length_nt, CDS_length_aa)

integrated_data <- integrated_data |>
  left_join(gene_lengths, by = "Gene_name")

# Adding GC content variables

integrated_data <- integrated_data |>
  left_join(cub_results$gc_results, by = "Gene_name")

## *****************************************************************************
## 5) CDC-based analysis ----
## _____________________________________________________________________________

# Full integration with your pipeline
cdc_results <- integrate_cdc_analysis(codon_usage, 
                                      genetic_code_dna_long, 
                                      integrated_data, 
                                      n_bootstrap = 10000,
                                      n_cores = 10)

# Re-plotting CDC-based neutrality plot highlighting the significant genes with CDC ----

# Extract just CDC columns we need
cdc_for_merge <- cdc_results |>
  dplyr::select(Gene_name, CDC, p_value, p_adj) |>
  dplyr::filter(!is.na(CDC))  # Remove genes without CDC

cat(sprintf("Valid CDC results: %d genes\n", nrow(cdc_for_merge)))

# Merge ENC, GC3s, and CDC results
enc_cdc_data <- cub_results$enc_results |>
  dplyr::left_join(cub_results$gc_results |> dplyr::select(Gene_name, GC3s), 
                   by = "Gene_name") |>
  dplyr::left_join(cdc_for_merge, by = "Gene_name") |>
  dplyr::filter(is.finite(ENC) & is.finite(GC3s) & ENC > 0 & ENC <= 61) |>
  dplyr::mutate(
    CDC_significant = !is.na(p_value) & p_value < 0.05,
    CDC_category = dplyr::case_when(
      is.na(p_value) ~ "No CDC data",
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "Not significant"
    )
  )

# Count significant genes
n_sig <- sum(enc_cdc_data$CDC_significant, na.rm = TRUE)
n_total <- sum(!is.na(enc_cdc_data$p_value))
pct_sig <- 100 * n_sig / n_total

cat(sprintf("Found %d / %d (%.1f%%) genes with significant CDC (p < 0.05)\n", 
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
  geom_point(data = enc_cdc_data |> filter(!CDC_significant | is.na(CDC_significant)),
             aes(x = GC3s, y = ENC), 
             color = "gray70", alpha = 0.3, size = 0.8) +
  # Foreground: CDC-significant genes
  geom_point(data = enc_cdc_data |> filter(CDC_significant),
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

cat("Enhanced ENC plot saved: ./results/ENC_plot_CDC_highlighted.pdf\n\n")

# Analyze CDC-significant genes: are they below the curve (under selection)?
cat("=== Position Analysis: CDC-Significant Genes Relative to Neutrality Curve ===\n")

# Calculate deviation from expected ENC
enc_cdc_data <- enc_cdc_data |>
  mutate(
    ENC_expected = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2),
    ENC_deviation = ENC - ENC_expected,
    Below_curve = ENC_deviation < 0
  )

# Compare CDC-significant vs non-significant genes
cdc_position_summary <- enc_cdc_data |>
  filter(!is.na(CDC_significant)) |>
  group_by(CDC_significant) |>
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
    data = enc_cdc_data |> filter(!is.na(CDC_significant))
  )
  cat(sprintf("ENC deviation (CDC-sig vs non-sig): W = %.0f, p = %.2e\n", 
              wilcox_enc$statistic, wilcox_enc$p.value))
  
  # Test if more CDC-significant genes are below the curve
  below_curve_table <- table(
    enc_cdc_data |> dplyr::filter(!is.na(CDC_significant)) |> dplyr::select(CDC_significant, Below_curve)
  )
  chi_test <- chisq.test(below_curve_table)
  cat(sprintf("Position relative to curve (chi-squared): X² = %.2f, p = %.2e\n", 
              chi_test$statistic, chi_test$p.value))
}

# Create density plot of ENC deviation
p_enc_deviation <- ggplot(enc_cdc_data |> filter(!is.na(CDC_significant)), 
                          aes(x = ENC_deviation, fill = CDC_significant)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("TRUE" = "#E41A1C", "FALSE" = "gray60"),
    labels = c("TRUE" = "CDC Significant", "FALSE" = "Not Significant"),
    name = ""
  ) +
  labs(
    title = "ENC Deviation from Expected: CDC-Significant vs Non-Significant Genes",
    subtitle = "Negative values = below neutrality curve (selection for codon bias)",
    x = "ENC Deviation (Observed - Expected)",
    y = "Density"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave("./results/ENC_deviation_by_CDC.pdf", p_enc_deviation, width = 9, height = 6)

## *****************************************************************************
## 6) Modeling relationship between CDC and Expression profiles ----
## _____________________________________________________________________________

# Box-Cox Transformation of expression data (lambda = 0.1) Log2 for simplicity
# Adding a small offset to the expression value (BoxCox works over positive numbers)

# integrated_data <- integrated_data |>
#   dplyr::mutate(High_exp = High_exp + 0.001)
# 
# box_cox_transformer <- preProcess(as.data.frame(integrated_data[, "High_exp"]), 
#                                   method = "BoxCox")
# 
# integrated_data <- integrated_data |>
#   dplyr::mutate(High_exp_BC = predict(box_cox_transformer, 
#                                       as.data.frame(integrated_data[, "High_exp"]))[[1]])

integrated_data <- integrated_data |>
  dplyr::mutate(High_exp_log2 = log2(High_exp + 1), # Adding 1 to avoid log2(0)
                High_exp_log10 = log10(High_exp + 1))  

# Generalized Linear Models ----

integrated_data <- integrated_data |> 
  left_join(enc_cdc_data |> dplyr::select(Gene_name, CDC), by = "Gene_name")
  
CDC_vs_exp <- glm(CDC ~ High_exp_log2 + CDS_length_nt, 
                 data = integrated_data, 
                 family = quasibinomial(link = "logit"))

summary(CDC_vs_exp)

# Check fitting results
integrated_data$predicted_CDC <- predict(CDC_vs_exp, type = "response")

# 2. Plot Observed vs Predicted
ggplot(integrated_data, aes(x = predicted_CDC, y = CDC)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_custom() +
  labs(title = "Quasibinomial GLM: Observed vs Predicted",
       x = "Predicted CDC",
       y = "Observed CDC")

ggsave("./results/CDC_observed_vs_predicted.pdf", 
       width = 8, height = 6)

# Density plots

# CDC ~ Exp
ggplot(data = integrated_data, 
       mapping = aes(x = High_exp_log2, y = CDC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/CDC_raw_vs_expression_density.pdf", 
       width = 10, height = 8)

# CDC ~ GC3
ggplot(data = integrated_data, 
       mapping = aes(x = GC3, y = CDC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/CDC_raw_vs_GC3_density.pdf", 
       width = 10, height = 8)

# Enc ~ Gene length
ggplot(data = integrated_data, 
       mapping = aes(x = CDS_length_nt, y = CDC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/CDC_raw_vs_gene_length_density.pdf", 
       width = 10, height = 8)

# GAM models ----

# Given the non-lineariry effect of the main confounder gene length
# we are going to fit a GAM model to account for this and assess effectively the
# effect of expression

cdc_model_quasibinom <- gam(CDC ~ High_exp_log2 + s(CDS_length_nt), 
                      data = integrated_data, family = quasibinomial(link = "logit"))

summary(cdc_model_quasibinom)

# Plotting detrended ENC against expression

confounder_model_gam <- gam(CDC ~ s(CDS_length_nt),
                            data = integrated_data,
                            family = quasibinomial(link = "logit"))

integrated_data$CDC_detrended <- residuals(confounder_model_gam)

summary(lm(CDC_detrended ~ High_exp_log2, data = integrated_data))

p_detrended <- ggplot(integrated_data, aes(x = High_exp_log2, y = CDC_detrended)) +
  # Use ggpointdensity for a clear view of the cluster
  geom_pointdensity(alpha = 0.5) + 
  
  # Add the linear regression line, which now shows the true effect
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  
  labs(
    title = "Detrended CDC vs. Gene Expression",
    subtitle = "Showing CDC after accounting for non-linear effects of gene length",
    y = "CDC Residuals (Detrended)",
    x = "log2(Expression + 1)"
  ) +
  theme_custom()

ggsave("./results/ENC_detrended_vs_expression.pdf", p_detrended, width = 8, height = 6)

# Define expression groups: Top 5% vs Bottom 5% (extreme comparison) ----

top_5_cutoff <- quantile(integrated_data$High_exp_log2, probs = 0.95)
bottom_5_cutoff <- quantile(integrated_data$High_exp_log2, probs = 0.05)

integrated_data$Expression_Group <- case_when(
  integrated_data$High_exp_log2 >= top_5_cutoff ~ "Top 5%",
  integrated_data$High_exp_log2 <= bottom_5_cutoff ~ "Bottom 5%",
  TRUE ~ "Middle 90%"
)

# Confounding out-based analysis (detendred CDC) ----

# Assesing significance of expression over the detrended residuals

cat("\n=== Kruskal-Wallis Test: Detrended ENC Residuals across Groups ===\n")

kw_detrended <- kruskal.test(CDC ~ Expression_Group, 
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

p_boxplot_detrended <- ggplot(integrated_data, aes(x = Expression_Group, y = CDC, fill = Expression_Group)) +
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
       y = "CDC",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/Detrended_ENC_by_expression_group.pdf", 
       p_boxplot_detrended, width = 8, height = 6)

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

cat("\n=== Step 10: Codon Adaptation Index (CAI) ===\n")
cat("CAI measures the degree of bias towards codons preferred in highly expressed genes\n")
cat("CAI ranges from 0 to 1, where higher values indicate stronger adaptation\n")
cat("Higher CAI = more similar to codon usage in highly expressed genes\n\n")

# Define reference set: Top 5% expressed genes
reference_genes <- read.table(file = 'data/CAI_Reference_Set_Mguttatus.txt')[, 1]

cat(sprintf("Using %d highly expressed genes as reference set with relevant functional annotations\n", 
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

ggsave("./results/optimal_codons_relative_adaptiveness.pdf", 
       width = 12, height = 10)

# Merge CAI with expression and ENC data
integrated_data <- integrated_data |>
  left_join(cai_values, by = "Gene_name")

cat("\n=== CAI vs Expression Level ===\n")
# Compare CAI across expression groups
cai_by_group <- integrated_data |>
  group_by(Expression_Group) |>
  summarise(
    n = n(),
    mean_CAI = mean(CAI, na.rm = TRUE),
    median_CAI = median(CAI, na.rm = TRUE),
    sd_CAI = sd(CAI, na.rm = TRUE),
    mean_ENC = mean(ENC, na.rm = TRUE)
  )

print(cai_by_group)

# Statistical tests for three groups
cat("\n=== Kruskal-Wallis Test: CAI across All Three Groups ===\n")
cat("H0: All three groups have the same median CAI\n")
kw_test <- kruskal.test(CAI ~ Expression_Group, data = integrated_data)
print(kw_test)

if (kw_test$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Install and load dunn.test if not available
  if (!require("dunn.test", quietly = TRUE)) {
    cat("Installing dunn.test package...\n")
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
    library(dunn.test)
  }
  
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
  
  cat("\nInterpretation of pairwise comparisons:\n")
  cat("  - Adjusted p-values account for multiple testing (FDR)\n")
  cat("  - p < 0.05 indicates significant difference between groups\n")
  
} else {
  cat("\nNo significant difference among groups (p >= 0.05)\n")
  cat("Post-hoc tests not necessary.\n")
}

# Additional pairwise effect sizes
cat("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

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
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle))
}

if (length(top_cai) > 0 && length(bottom_cai) > 0) {
  d_top_bottom <- cohens_d_calc(top_cai, bottom_cai)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom))
}

if (length(middle_cai) > 0 && length(bottom_cai) > 0) {
  d_middle_bottom <- cohens_d_calc(middle_cai, bottom_cai)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom))
}

cat("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

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
cat("\nBoxplot saved: ./results/CAI_by_expression_group.pdf\n")

# Correlation between CAI and other metrics
# Scatter plot: CAI vs ENC
p_cdc_enc <- ggplot(integrated_data, aes(x = CDC_detrended, y = CAI, color = Expression_Group)) +
  geom_point(alpha = 0.3, size = 1) +
  # Add lm per group
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                 "Bottom 5%" = "#377EB8",
                                 "Middle 90%" = "#999999")) +
  labs(title = "CAI vs CDC by Expression Level",
       subtitle = "Higher CDC and Higher CAI indicate stronger codon bias",
       x = "CDC",
       y = "CAI (Codon Adaptation Index)",
       color = "Expression Group") +
  theme_custom()

ggsave("./results/CAI_vs_ENC_scatter.pdf", p_cdc_enc, width = 10, height = 6)

# 7.2) Compare absolute codon frequencies: Top 5% vs Rest ----
# This shows that raw frequencies differ, but not all differences are due to selection
# Some codons are frequent simply because their amino acids are frequent
# This motivates the need for enrichment analysis to correct for amino acid composition

cat("\n=== 7.1: Absolute Codon Frequencies in Top 5% vs Rest ===\n")
cat("Comparing raw codon usage to motivate enrichment-based correction\n\n")

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
    y = "Absolute Frequency",
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

# Show codons with largest absolute differences
cat("\nCodons with largest frequency differences (Top 5% - Rest):\n")
top_diff <- freq_comparison |>
  arrange(desc(abs(Freq_Diff))) |>
  dplyr::select(Codon, AA, Freq_Top5, Freq_Rest, Freq_Diff) |>
  head(10)
print(top_diff)

cat("\nNote: These raw differences don't account for amino acid composition.\n")
cat("Enrichment analysis (via w-table) corrects for this by normalizing within amino acids.\n\n")

# 7.2.1) Enriched codons usage ----

# Get preferred codons (w = 1.0 from CAI)
enriched_codons_vec <- w_table |>
  dplyr::filter(relative_adaptiveness == 1.0) |>
  dplyr::pull(codon)

cat(sprintf("Using %d preferred codons (w = 1.0)\n\n", length(enriched_codons_vec)))

# Merge codon usage with expression groups
codon_usage_with_groups <- codon_usage |>
  dplyr::left_join(integrated_data |> dplyr::select(Gene_name, Expression_Group), 
                   by = "Gene_name")

# Filter to top 5% and rest (bottom 95%)
top5_genes <- codon_usage_with_groups |> dplyr::filter(Expression_Group == "Top 5%")
rest_genes <- codon_usage_with_groups |> dplyr::filter(Expression_Group != "Top 5%")

cat(sprintf("Top 5%% genes (selected): %d genes\n", nrow(top5_genes)))
cat(sprintf("Bottom 95%% genes (neutral/rest): %d genes\n\n", nrow(rest_genes)))

# Calculate for both groups
cat("Calculating preferred codon usage per amino acid...\n")

enriched_aa <- count_preferred_by_aa(top5_genes, enriched_codons_vec, genetic_code_dna_long)
enriched_aa$Group <- "Selected (Top 5%)"

rest_aa <- count_preferred_by_aa(rest_genes, enriched_codons_vec, genetic_code_dna_long)
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

cat("\n✓ Results saved: ./results/enriched_codon_usage.csv\n\n")

# Print table
cat("=== Enriched Codon Usage: Selected vs Rest ===\n\n")
cat(sprintf("%-4s %-4s %-15s %-12s %-12s %-12s %-8s\n",
            "AA", "Deg", "Preferred", "Top5%", "Rest95%", "Difference", "Fold"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:nrow(comparison_table)) {
  row <- comparison_table[i, ]
  cat(sprintf("%-4s %-4d %-15s %-12.4f %-12.4f %-12.4f %-8.2f\n",
              row$Amino_Acid,
              row$N_synonymous,
              substr(row$Preferred_codons, 1, 15),
              row$Selected_prop,
              row$Rest_prop,
              row$Difference,
              row$Fold_enrichment))
}
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# Statistical summary
cat("=== Summary Statistics ===\n\n")
cat(sprintf("Mean proportion enriched (Top 5%%): %.4f\n", 
            mean(comparison_table$Selected_prop, na.rm = TRUE)))
cat(sprintf("Mean proportion enriched (Rest 95%%): %.4f\n", 
            mean(comparison_table$Rest_prop, na.rm = TRUE)))
cat(sprintf("Mean difference: %.4f\n", 
            mean(comparison_table$Difference, na.rm = TRUE)))
cat(sprintf("Mean fold enrichment: %.2f\n\n", 
            mean(comparison_table$Fold_enrichment, na.rm = TRUE)))

# Wilcoxon test
wilcox_test <- wilcox.test(comparison_table$Selected_prop, 
                           comparison_table$Rest_prop,
                           paired = TRUE)

cat(sprintf("Wilcoxon signed-rank test (paired by amino acid):\n"))
cat(sprintf("  V = %.1f, p-value = %.2e\n", 
            wilcox_test$statistic, wilcox_test$p.value))

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
  labs(title = "Preferred Codon Usage: Top 5% vs Rest (Original w=1)",
       subtitle = "Proportion of preferred codons per amino acid (before enrichment correction)",
       x = "Amino Acid (ordered by difference)",
       y = "Proportion of Preferred Codons") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("./results/preferred_codon_usage_comparison.pdf", p_comparison,
       width = 10, height = 6)

cat("\n")

# 7.2.2) Statistical testing for codon-level differences ----

cat("=== Statistical Testing: Per-Codon Proportion Differences ===\n\n")

# Source the testing function
source("./src/test_codon_proportions.R")

# Test each codon individually
codon_test_results <- test_codon_proportions(
  selected_usage = top5_genes,
  neutral_usage = rest_genes,
  genetic_code = genetic_code_dna_long,
  method = "chisq",
  fdr_correction = TRUE
)

# Save results
write.csv(codon_test_results, 
          "./results/codon_proportion_test_results.csv",
          row.names = FALSE)

cat("\n✓ Results saved: ./results/codon_proportion_test_results.csv\n")

# Create summary plots
p_selection_summary <- plot_codon_selection_summary(
  codon_test_results,
  output_file = "./results/codon_selection_summary.pdf"
)

p_classification_heatmap <- plot_codon_classification_heatmap(
  codon_test_results,
  w_table,
  output_file = "./results/codon_classification_heatmap.pdf"
)

# Print key findings
cat("\n=== Key Findings ===\n\n")

sig_codons <- codon_test_results |> dplyr::filter(Significant)
sig_preferred <- sig_codons |> 
  dplyr::left_join(w_table |> dplyr::select(codon, relative_adaptiveness),
                   by = c("Codon" = "codon")) |>
  dplyr::filter(relative_adaptiveness == 1.0)

cat(sprintf("Codons with significant proportion differences (FDR < 0.05): %d / %d (%.1f%%)\n",
            nrow(sig_codons), 
            nrow(codon_test_results),
            100 * nrow(sig_codons) / nrow(codon_test_results)))

cat(sprintf("  - Enriched in Top 5%% (under selection): %d\n", 
            sum(sig_codons$Difference > 0)))
cat(sprintf("  - Depleted in Top 5%% (avoided): %d\n", 
            sum(sig_codons$Difference < 0)))

cat(sprintf("\nPreferred codons (w=1) that are significantly enriched: %d / %d (%.1f%%)\n",
            nrow(sig_preferred),
            length(enriched_codons_vec),
            100 * nrow(sig_preferred) / length(enriched_codons_vec)))

if (nrow(sig_preferred) > 0) {
  cat("\nThese 'under selection' preferred codons are:\n")
  sig_pref_sorted <- sig_preferred |> 
    dplyr::arrange(dplyr::desc(Difference)) |>
    dplyr::select(Codon, Amino_Acid, Selected_Prop, Neutral_Prop, 
                  Difference, p_adj)
  print(sig_pref_sorted)
}

cat("\n")

# 7.2.3) Diagnose CAI vs Proportion Discrepancies ----

cat("=== Diagnosing CAI w-values vs Proportion Differences ===\n\n")

# Source diagnostic function
source("./src/diagnose_cai_vs_proportion.R")

# Run diagnostic
diagnostic_results <- diagnose_cai_vs_proportion(
  w_table = w_table,
  test_results = codon_test_results,
  codon_usage = codon_usage,
  expression_groups = integrated_data,
  genetic_code = genetic_code_dna_long
)

# Create corrected classification
corrected_classification <- create_corrected_classification(
  w_table = w_table,
  test_results = codon_test_results
)

# Save corrected classification
write.csv(corrected_classification,
          "./results/codon_classification_corrected.csv",
          row.names = FALSE)

cat("\n✓ Corrected classification saved: ./results/codon_classification_corrected.csv\n")

# Update the codon_test_results with corrected classification for biplots
codon_test_results <- codon_test_results |>
  left_join(corrected_classification |> 
              dplyr::select(Codon, Combined_Classification),
            by = "Codon") |>
  dplyr::mutate(
    # Update Classification to be more accurate
    Classification_Original = Classification,
    Classification = Combined_Classification
  )

# 7.2.4) Enriched codons (corrected) ----
cat("\n=== 7.2.4: Defining Corrected Preferred Codons ===\n")
cat("Strategy: Preferentially choose w=1 codons, but use enrichment data when conflicts arise\n\n")

# Merge w-table with statistical test results
codon_combined <- w_table |>
  dplyr::left_join(
    codon_test_results |> 
      dplyr::select(Codon, Amino_Acid, Selected_Prop, Neutral_Prop, Difference, 
                    p_value, p_adj, Significant),
    by = c("codon" = "Codon")
  ) |>
  dplyr::rename(Codon = codon, AA = amino_acid)

# Decision rule for each amino acid:
# 1. If w=1 codon is significantly enriched (p_adj < 0.05 & Difference > 0): PREFERRED
# 2. If w=1 codon is NOT significantly enriched BUT another codon IS: Use the enriched one
# 3. If no codon is significantly enriched: Use w=1 as default (CAI definition)
# 4. If multiple codons enriched: Use the one with largest effect size

enriched_codons_corrected <- codon_combined |>
  dplyr::group_by(AA) |>
  dplyr::mutate(
    # Flag w=1 codons
    is_w1 = (relative_adaptiveness == 1.0),
    # Flag significantly enriched codons
    is_enriched = (Significant & Difference > 0),
    # Flag w=1 AND enriched (ideal case)
    is_w1_and_enriched = is_w1 & is_enriched
  ) |>
  dplyr::arrange(AA, desc(is_w1_and_enriched), desc(is_enriched), 
                 desc(relative_adaptiveness), desc(abs(Difference))) |>
  dplyr::slice(1) |>  # Take first (best) codon per amino acid
  dplyr::ungroup() |>
  dplyr::mutate(
    Selection_Rationale = dplyr::case_when(
      is_w1_and_enriched ~ "w=1 AND enriched (strong evidence)",
      is_w1 & !is_enriched & !any(is_enriched) ~ "w=1, no enrichment signal (CAI default)",
      is_w1 & !is_enriched & any(is_enriched) ~ "w=1 but NOT enriched (kept due to no better option)",
      !is_w1 & is_enriched ~ "NOT w=1 but significantly enriched (corrected)",
      TRUE ~ "Default (CAI w=1)"
    )
  )

# Summary of decisions
cat("Decision Summary:\n")
decision_summary <- enriched_codons_corrected |>
  dplyr::group_by(Selection_Rationale) |>
  dplyr::summarise(
    n_codons = n(),
    codons = paste(Codon, collapse = ", ")
  )
print(decision_summary)
cat("\n")

# Highlight cases where enrichment overrode w=1
conflicts <- codon_combined |>
  dplyr::group_by(AA) |>
  dplyr::filter(any(relative_adaptiveness == 1.0) & 
                any(Significant & Difference > 0)) |>
  dplyr::arrange(AA, desc(relative_adaptiveness)) |>
  dplyr::ungroup()

if (nrow(conflicts) > 0) {
  cat("=== Cases where w=1 and enrichment disagree ===\n")
  
  # For each AA with conflict, show w=1 vs enriched
  conflict_summary <- conflicts |>
    dplyr::group_by(AA) |>
    dplyr::summarise(
      w1_codon = Codon[relative_adaptiveness == 1.0][1],
      w1_enriched = Significant[relative_adaptiveness == 1.0][1] & 
                    Difference[relative_adaptiveness == 1.0][1] > 0,
      enriched_codons = paste(Codon[Significant & Difference > 0], collapse = ", "),
      n_enriched = sum(Significant & Difference > 0, na.rm = TRUE),
      final_choice = enriched_codons_corrected$Codon[enriched_codons_corrected$AA == AA[1]][1]
    )
  
  print(conflict_summary)
  cat("\n")
}

# Save corrected preferred codons
write.csv(enriched_codons_corrected, 
          "./results/enriched_codons_corrected.csv",
          row.names = FALSE)

cat("✓ Corrected preferred codons saved: ./results/enriched_codons_corrected.csv\n\n")

# Show final preferred codon set
cat("=== Final Enriched Codons (Corrected) ===\n")
final_table <- enriched_codons_corrected |>
  dplyr::select(AA, Codon, w = relative_adaptiveness, 
                Enriched = is_enriched, Diff = Difference, Selection_Rationale)
print(final_table, n = Inf)
cat("\n")

# Use corrected set for downstream analysis
enriched_codons_mg <- enriched_codons_corrected

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

# Access results:
# - dM_results$dM_introns          : dM data frame from introns
# - dM_results$dM_intergenic       : dM data frame from intergenic
# - dM_results$global_stats_introns: Nucleotide frequencies from introns
# - dM_results$global_stats_intergenic: Nucleotide frequencies from intergenic
# - dM_results$output_files        : Paths to generated CSV files
# - dM_results$intermediates       : Raw data (seq_data, nuc_composition, etc.)

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

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_naive_2/run_1",
  "./results/MCMC_results/results_naive_2/run_2",
  "./results/MCMC_results/results_naive_2/run_3",
  "./results/MCMC_results/results_naive_2/run_4",
  "./results/MCMC_results/results_naive_2/run_5",
  "./results/MCMC_results/results_naive_2/run_6"
)

Naive_conv <- GR_convergence(run_dirs)

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

# Convergence was achieved for dM_fixed

# 8.1.2.1) Checking the correlation between estimates of phi and the expression data ----

# From now on, we will work with chain 1, as an example

phi_hat_dM_fixed <- read.csv(file = "results/MCMC_results/results_dM_fixed/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi_dM_fixed <- exp_complete |>
  left_join(phi_hat_dM_fixed, by = join_by("Gene" == "GeneID")) |>
  dplyr::mutate(High_exp_log10 = log10(High_exp + 1)) |>
  na.exclude()

cor.test(phi_dM_fixed$Mean.log10.Phi, phi_dM_fixed$High_exp_log10)

# Visualization

ggplot(data = phi_dM_fixed, aes(x = Mean.log10.Phi,
                                y = High_exp_log10)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Estimated phi (log10)") +
  ylab("Empirical Max Expresion (log10)")

ggsave()

# There is no good correspondence with empirical data
# Next step is to pass expression data to the AnaCoDa

# 8.1.3) Preparing the expression data ----

# 1. Filter for complete cases (Intersection of Leaf and Bud)
# We strictly remove genes with 0 counts in either tissue
multi_tissue_phi <- exp_complete |>
  dplyr::select(Gene, Exp_leaf, Exp_bud) |>
  dplyr::filter(Exp_leaf > 0 & Exp_bud > 0) |>
  dplyr::rename(GeneID = Gene) |> # AnaCoDa expects "GeneID" as first col
  dplyr::filter(GeneID %in% names(trans)) # Ensures correspondence with transcriptome file
  
# 2. Calculate sphi (Global Prior)
# We estimate the "True Phi" shape by taking the mean of the log-expressions
# This gives the model the "width" of the overall distribution.
log_means <- rowMeans(log(multi_tissue_phi[, c("Exp_leaf", "Exp_bud")]))
sphi_init <- sd(log_means)

# 3. Calculate sepsilon (Noise per tissue)
# AnaCoDa needs a vector: c(noise_leaf, noise_bud)
# A good heuristic for initialization is the SD of the log-expression for that tissue.
# (The model will refine this during MCMC, but this puts it in the right ballpark)

sepsilon_leaf <- sd(log(multi_tissue_phi$Exp_leaf))
sepsilon_bud  <- sd(log(multi_tissue_phi$Exp_bud))

sepsilon_init <- c(sepsilon_leaf, sepsilon_bud)

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
  "./results/MCMC_results/results_dM_fixed_with_phi/run_1",
  "./results/MCMC_results/results_dM_fixed_with_phi/run_2",
  "./results/MCMC_results/results_dM_fixed_with_phi/run_3",
  "./results/MCMC_results/results_dM_fixed_with_phi/run_4",
  "./results/MCMC_results/results_dM_fixed_with_phi/run_5",
  "./results/MCMC_results/results_dM_fixed_with_phi/run_6"
)

dM_fixed_with_phi_conv <- GR_convergence(run_dirs, 
                                         parameter = 'selection') # Mutation is fixed

# 8.1.3.2) Codon frequency trajectories across expression levels ----

# This section visualizes whether the ROC multinomial model:
#   P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
# correctly predicts how codon frequencies change with expression.

# Load validation functions
source("./src/roc_model_validation.R")

cat("\n=== ROC Model: Codon Frequency Trajectories ===\n")

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
  dplyr::mutate(Exp_log10 = log10(High_exp + 1)) |>
  dplyr::select(Gene, Exp_log10)

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

cat("\n✓ Codon trajectory analysis complete!\n")
cat("  Plot saved to: ./results/ROC_codon_trajectories.pdf\n")

# 8.1.4) dM-fixed-intergenic ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed_intergenic/run_1",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_2",
  "./results/MCMC_results/results_dM_fixed_intergenic/run_3"
)

dM_fixed_intergenic <- GR_convergence(run_dirs, 
                                       parameter = 'selection') # Mutation is fixed

# Poor convergence (intergenic-based mutation bias is not adequate)

# 8.1.5) dM-fixed-with-phi-intergenic ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic/run_1",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic/run_2",
  "./results/MCMC_results/results_dM_fixed_with_phi_intergenic/run_3"
)

dM_fixed_with_phi_intergenic <- GR_convergence(run_dirs, 
                                       parameter = 'selection') # Mutation is fixed

# Poor convergence (intergenic-based mutation bias is not adequate)

# Exploring and overlapping traces to rule out bad mixing ----

parameters_objects <- list(run1 = loadParameterObject(file = paste(run_dirs[1], "R_objects/parameter.Rda", sep = '/')),
                           run2 = loadParameterObject(file = paste(run_dirs[2], "R_objects/parameter.Rda", sep = '/')),
                           run3 = loadParameterObject(file = paste(run_dirs[3], "R_objects/parameter.Rda", sep = '/')))

parameters_objects <- list(run1 = listRDA(paste(run_dirs[1], "R_objects/parameter.Rda", sep = '/')),
                           run2 = listRDA(paste(run_dirs[2], "R_objects/parameter.Rda", sep = '/')),
                           run3 = listRDA(paste(run_dirs[3], "R_objects/parameter.Rda", sep = '/')))

# 1. Select one of the high deviants codons
codon_index <- 36  

# 2. Extract Data function
# This pulls the specific codon trace from all 3 runs into one data frame
extract_traces <- function(objects, codon_idx) {
  df_list <- list()
  
  for (run_name in names(objects)) {
    # Access logic: Run -> selectionTrace -> Mixture 1 -> Codon Index
    trace_data <- objects[[run_name]]$selectionTrace[[1]][[codon_idx]]
    
    # Create temporary DF
    temp_df <- data.frame(
      Iteration = 1:length(trace_data),
      Value = trace_data,
      Run = run_name
    )
    df_list[[run_name]] <- temp_df
  }
  
  return(do.call(rbind, df_list))
}

# 3. Create the Master Data Frame
plot_data <- extract_traces(parameters_objects, codon_index)
plot_data <- subset(plot_data, Iteration > max(plot_data$Iteration) * 0.5)

# 4. PLOT 1: The Trace Overlay (The "Are we mixing?" Plot)
p1 <- ggplot(plot_data, aes(x = Iteration, y = Value, color = Run)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  theme_custom() +
  labs(title = paste("Trace Overlay: Codon Index", codon_index),
       subtitle = "Flat lines at different levels = Good Mixing + Multimodality",
       y = "Selection Cost (Delta Eta)",
       x = "MCMC Sample") +
  theme(legend.position = "top")

# 5. PLOT 2: Density Overlay (The "Distinct Realities" Plot)
p2 <- ggplot(plot_data, aes(x = Value, fill = Run)) +
  geom_density(alpha = 0.5) +
  theme_custom() +
  labs(title = paste("Posterior Density: Codon Index", codon_index),
       subtitle = "Non-overlapping peaks = Model Misspecification / Broken Energy Landscape",
       x = "Selection Cost (Delta Eta)",
       y = "Density") +
  theme(legend.position = "top")

# 6. PLOT 3: Autocorrelation (The "Stickiness" Check)
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

eta_data <- read.csv(file = "results/MCMC_results/results_dM_fixed_with_phi/run_1/Parameter_est/Cluster_1_Selection.csv")

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

cat(sprintf("✓ Preferred codons from ROC model: %d amino acids\n", nrow(preferred_codons_roc)))

# 8.3) Extracting selection estimates from the best model (dM-fixed-with_phi) ----

genome <- initializeGenomeObject(file = 'data/IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa',
                                 match.expression.by.id = TRUE,
                                 observed.expression.file = 'data/observed_expression_multitissue.csv') 

parameter_object <- loadParameterObject(file = "./results/MCMC_results/results_dM_fixed_with_phi/run_1/R_objects/parameter.Rda")

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
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
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

# Get total selection load per gene

counts_df <- as.data.frame(codon_usage)
rownames(counts_df) <- counts_df$Gene_name
counts_df$Gene_name <- NULL # Remove the ID column so it's all numeric

# Ensure 'selection_coeff' is a matrix
sel_mat <- as.matrix(selection_coeff)

# Find common genes (rows) and common codons (columns)
common_genes <- intersect(rownames(counts_df), rownames(sel_mat))
common_codons <- intersect(colnames(counts_df), colnames(sel_mat))

message(paste("Matching:", length(common_genes), "genes and", 
              length(common_codons), "codons."))

# Align genes and codons
counts_aligned <- counts_df[common_genes, common_codons]
sel_aligned <- sel_mat[common_genes, common_codons]
counts_aligned <- as.matrix(counts_aligned)

# Logic: Count * abs(Selection_Coefficient)
# We use abs() because our values are negative penalties (e.g., -0.06).
# A penalty of -0.06 is a "load" of 0.06.

# Element-wise multiplication (NOT matrix multiplication %*%)
gene_load_matrix <- counts_aligned * abs(sel_aligned)

# Sum across rows to get the total load per gene
total_selection_load <- rowSums(gene_load_matrix, na.rm = TRUE)

# Selection intensity ----
n_synonymous_codons <- rowSums(counts_aligned, na.rm = TRUE)

sel_intensity <- total_selection_load / n_synonymous_codons

# Unweighted selection intensity
sel_intensity_uw <- rowMeans(abs(sel_aligned), na.rm = TRUE)

selection_coeff_intensity <- data.frame(Gene_name = names(sel_intensity),
                                        S_coeff = as.vector(sel_intensity),
                                        S_coeff_uw = as.vector(sel_intensity_uw))

# 8.3.1) Analyzing the correlation between total selective pressure and CAI and CDC ----

integrated_data <- integrated_data |>
  dplyr::mutate(
    # Calculate Geometric Mean (add small epsilon if zeros exist)
    Geom_Exp = sqrt(Exp_leaf * Exp_bud),
    
    # Optional: Log transformation for plotting later
    Log_Geom_Exp = log10(Geom_Exp + 0.0001) 
  )

selection_coeff_intensity <- selection_coeff_intensity |>
  left_join(integrated_data |> dplyr::select(Gene_name, CAI, CDC, ENC, 
                                              Total_Codons, GC3s, Geom_Exp, Log_Geom_Exp,
                                              Pi_mean_4fold, TajimaD_4fold)) |>
  na.exclude()

# Relation between geom_expression and S

ggplot(selection_coeff_intensity, aes(x = Geom_Exp, y = S_coeff)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  theme_custom()

ggsave(filename = "./results/Selection_Coefficient_vs_Geom_Expression.pdf",
       width = 6, height = 4)

# Relation between geom_expression and S_uw

ggplot(selection_coeff_intensity, aes(x = Geom_Exp, y = S_coeff_uw)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  theme_custom()

ggsave(filename = "./results/Selection_Coefficient_Unweighted_vs_Geom_Expression.pdf",       width = 6, height = 4)

cor.test(selection_coeff_intensity$S_coeff_uw, selection_coeff_intensity$S_coeff)

# Correlation between selection metric and CUB metrics

cor_S_and_bias <- corrr::correlate(x = as.matrix(selection_coeff_intensity[, 2:5]),
                                   method = "spearman")

cor.test(selection_coeff_intensity$S_coeff, selection_coeff_intensity$CAI)

# Plot both selection measurements

ggplot(selection_coeff_intensity, aes(x = S_coeff_uw, y = S_coeff)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  theme_custom() +
  labs(x = "Selection Intensity (unweighted)", y = "Selection Intensity (weighted)")

ggsave(filename = "./results/Selection_Coefficient_Weighted_vs_Unweighted.pdf",
       width = 6, height = 4)

# 8.3.2) Final visualization ----
# Prepare Plot Data (Final Visualization)

plot_data <- selection_coeff_intensity |>
  mutate(
    # Log Transform Selection Load (Load = Total Cost per Gene)
    # Note: If you want Intensity (per codon), swap S_load for Selection_Intensity
    Log_S_coeff = log10(S_coeff_uw + 0.01), 
    
    # Log Transform Expression (using the new clear name)
    Log_Phi = log10(Geom_Exp + 0.0001), 
    
    # Log Transform Length
    Log_Length = log10(Total_Codons)
  ) |>
  filter(!is.na(ENC), !is.na(Total_Codons), !is.na(CAI))

# 4. Visualization Setup

# Define common color scale limits
phi_range <- range(plot_data$Log_Phi, na.rm = TRUE)

# Fit linear model for Panel C annotation (Load vs Length)
tail_data <- plot_data |> 
  dplyr::filter(Log_S_coeff > -0.5)

# Fit LM specifically on the tail
lm_tail <- lm(GC3s ~ Log_S_coeff, data = tail_data)
lm_tail_eq <- sprintf("Tail: y = %.2f + %.2fx\nR² = %.3f, p = %.4f",
                      coef(lm_tail)[1], 
                      coef(lm_tail)[2], 
                      summary(lm_tail)$r.squared,
                      summary(lm_tail)$coefficients["Log_S_coeff","Pr(>|t|)"])

# Panel A: Selection Load Distribution
drift_thresh <- log10(1 + 0.01)   
strong_thresh <- log10(5 + 0.01)
y_max_anno <- 10000

p1 <- ggplot(plot_data, aes(x = Log_S_coeff)) +
  
  # --- Background Shading ---
  annotate("rect", xmin = -Inf, xmax = drift_thresh, 
           ymin = 0, ymax = Inf, fill = "gray95", alpha = 0.8) +
  annotate("rect", xmin = strong_thresh, xmax = Inf, 
           ymin = 0, ymax = Inf, fill = "#ffe5e5", alpha = 0.5) +
  
  # --- Histogram ---
  geom_histogram(bins = 100, fill = "#69b3a2", color = "white", linewidth = 0.05) +
  
  # --- Vertical Threshold Lines ---
  geom_vline(xintercept = drift_thresh, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = strong_thresh, linetype = "dotted", color = "red") +
  
  # --- Vertical Text Annotations (No Ne*s text) ---
  # Drift Label (Left side, Vertical)
  annotate("text", x = log10(0.05), y = y_max_anno/6, 
           label = "Drift Dominated", 
           color = "gray50", fontface = "bold", size = 4, 
           angle = 0) +
  
  # Strong Selection Label (Right side, Vertical)
  annotate("text", x = log10(12), y = y_max_anno/15, 
           label = "Strong Selection", 
           color = "red", fontface = "bold", size = 4, 
           angle = 90) +
  
  # --- Custom "Separated" Rug ---
  # We draw segments explicitly at y = -0.5 (below axis) to create the gap
  geom_segment(aes(x = Log_S_coeff, xend = Log_S_coeff, 
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
  
  labs(x = expression(Log[10]("Selection Intensity" ~ (S[avg]))), 
       y = "Gene Count (Log1p Scale)") +
  
  theme_custom() +
  # Add margin at bottom to ensure the new rug doesn't get cut off
  theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))

# Panel B: CAI vs Selection Load
p2 <- ggplot(plot_data, aes(x = Log_S_coeff, y = CAI)) +
  geom_point(aes(color = Log_Phi), alpha = 0.6, size = 1) +
  scale_color_viridis_c(option = "plasma", name = expression(Log[10](Phi[geom])), 
                        limits = phi_range, direction = 1) +
  geom_smooth(color = "black") +
  labs(x = expression(Log[10](S[avg])), y = "CAI") +
  theme_custom()

# Panel C: Gene Length Effect
p3_main <- ggplot(plot_data, aes(x = Log_S_coeff, y = GC3s)) +
  geom_point(aes(color = Log_Phi), alpha = 0.3, size = 0.8) +
  
  # Use GAM to show the true shape (Flat -> Rising)
  geom_smooth(method = "gam", color = "black", se = TRUE, size = 0.8) +
  
  scale_color_viridis_c(option = "plasma", name = expression(Log[10](Phi[geom])), 
                        limits = phi_range, direction = 1) +
  
  # Add a box to show where the inset comes from
  annotate("rect", xmin = -0.5, xmax = max(plot_data$Log_S_coeff), 
           ymin = min(tail_data$GC3s), ymax = max(tail_data$GC3s),
           fill = NA, color = "red", linetype = "dashed", alpha = 0.5) +
  
  labs(x = expression(Log[10]("Selection Intensity" ~ (S[avg]))), 
       y = expression(GC3s)) +
  theme_custom() +
  theme(legend.position = "none") # Remove legend (shared in combined plot)

# --- 3. Create the Inset Plot (The Linear Tail) ---
p3_zoom_in <- ggplot(tail_data, aes(x = Log_S_coeff, y = GC3s)) +
  geom_point(aes(color = Log_Phi), alpha = 0.6, size = 1) +
  # Linear model for the tail
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  scale_color_viridis_c(option = "plasma", limits = phi_range, direction = 1) +
  
  # Add the equation inside the inset
  annotate("text", x = min(tail_data$Log_S_coeff), y = max(tail_data$GC3s), 
           label = lm_tail_eq, hjust = 0, vjust = 1, size = 3, color = "red") +
  
  theme_custom() +
  labs(x = expression(Log[10]("Selection Intensity" ~ (S[avg]))), 
       y = expression(GC3s)) +
  theme(legend.position = "none") 

# Combine
combined_plot <- (p1 | p2) / p3_main / p3_zoom_in + plot_annotation(tag_levels = 'A')
ggsave("results/Selection_Landscape_Final.pdf", combined_plot, width = 11, height = 14)

# 8.4) GO-enrichment analysis of genes with a massive selection load ----

thr_sel <- 1

subset_strongly_shaped_by_s <- selection_coeff_intensity |>
  dplyr::filter(S_coeff > thr_sel) |>
  dplyr::pull(Gene_name)

custom_bag <- selection_coeff_intensity |> dplyr::pull(Gene_name)

GO_results <- gost(query = subset_strongly_shaped_by_s,
                   organism = 'gp__q7VP_EAck_dZk',
                   multi_query = F,
                   significant = T,
                   correction_method = 'fdr',
                   domain_scope = "custom",
                   custom_bg = custom_bag,
                   user_threshold = 0.05)
  
GO1_plot <- gostplot(GO_results, capped = TRUE, interactive = FALSE)

ggsave(filename = "./results/Manhattan_like_GO.pdf", plot = GO1_plot, 
       width = 10, height = 8)

# Export results

write.csv(x = GO_results$result |> dplyr::select(-parents), 
          file = "./results/Go_enrichment.csv", quote = T, 
          row.names = F)

# 8.5) Getting top 10 genes in terms of S_load ----

subset_strongly_shaped_by_s <- selection_coeff_intensity |>
  dplyr::filter(S_coeff > thr_sel) |>
  dplyr::arrange(desc(S_coeff)) |>
  dplyr::slice(1:10) |>
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
  dplyr::filter(locusName %in% subset_strongly_shaped_by_s)

# Export information about top 10 genes

write.csv(x = detailed_annotation, 
          file = "./results/Top10_genes_strong_selection_load.csv", 
          quote = T, row.names = F)

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

merge_2_and_4_to_6_fold <- function(preference_df, AA_family_col)
{
  #' This function will condense together the 4 and 2 fold families from a 6-fold
  #' aminoacid family back to the one preferred codon per amino acid. It will take
  #' the aminoacid with the greater adaptiveness.
  #' 
  #' FOR COMPATIBILITY WITH OTHER PLANT STUDIES
  #' 
  #' Args:
  #' preference_df: Data frame with columns for Amino Acid, Codon_RNA, relative_adaptiveness
  #' AA_family_col: Column name indicating the root amino acid family (e.g., "Leu" for "Leu_2" and "Leu_4")
  #' 
  #' ___________________________________________________________________________
  
  condensed_preferences <- preference_df |>
    dplyr::group_by(!!sym(AA_family_col)) |>
    dplyr::arrange(!!sym(AA_family_col), 
                   eta) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(Amino_Acid = !!sym(AA_family_col), 
                  Codon_RNA, eta)
  
  return(condensed_preferences)
}

preferred_codons_mg <- merge_2_and_4_to_6_fold(
  preferred_codons_comparative,
  "AA_root"
)

cat(sprintf("Found %d preferred codons for M. guttatus\n\n", 
            nrow(preferred_codons_mg)))

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

cat("Extended comparison table saved: ./results/plant_preferred_codons_comparison.csv\n\n")

# Print summary
cat("=== M. guttatus Preferred Codons ===\n")
print(preferred_codons_mg |> dplyr::select(Amino_Acid, Codon = Codon_RNA, Weight = eta))

# Calculate codon preference similarity between species
cat("\n\n=== Cross-Species Codon Preference Analysis ===\n\n")

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

cat("Jaccard Similarity Matrix (codon preference overlap):\n")
print(round(similarity_matrix, 3))
cat("\n")

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

cat("Dendrogram saved: ./results/plant_codon_preference_dendrogram.pdf\n\n")

# Create unrooted phylogram using ape package
tree <- as.phylo(hc)

# Save unrooted tree
pdf("./results/plant_codon_preference_unrooted.pdf", width = 10, height = 10)
par(mar = c(1, 1, 3, 1))
plot(tree, type = "unrooted", main = "Unrooted Tree: Codon Preference Similarity",
     cex = 1.2, lab4ut = "axial", edge.width = 2)
dev.off()

cat("Unrooted tree saved: ./results/plant_codon_preference_unrooted.pdf\n\n")

# Create heatmap of codon preferences
cat("Creating codon preference heatmap...\n")

# Prepare data for heatmap
codon_df <- as.data.frame(t(codon_matrix))
codon_df$Codon <- rownames(codon_df)
codon_df$AA <- genetic_code_dna_long[gsub("U", "T", codon_df$Codon)]

# Reshape for ggplot
codon_long <- codon_df |>
  pivot_longer(cols = all_of(species), names_to = "Species", values_to = "Preferred") |>
  dplyr::mutate(Species = gsub("_", " ", Species),
                Species = factor(Species, levels = gsub("_", " ", species)))

# Create heatmap
p_heatmap <- ggplot(codon_long, aes(x = Species, y = Codon, fill = factor(Preferred))) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("0" = "gray90", "1" = "#E41A1C"),
                    labels = c("Not Preferred", "Preferred")) +
  facet_grid(AA ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Preferred Codon Usage Across Plant Species",
       subtitle = "Based on highest expression genes",
       x = "", y = "Codon", fill = "") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "bottom")

ggsave("./results/plant_codon_preference_heatmap.pdf", p_heatmap, 
       width = 10, height = 18)

cat("Heatmap saved: ./results/plant_codon_preference_heatmap.pdf\n\n")

# Create color-coded comparison plot showing M. guttatus sharing patterns

cat("Creating color-coded codon preference comparison plot...\n")

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

cat("✓ Color-coded comparison plot saved: ./results/plant_codon_preference_comparison_colored.pdf\n\n")

# Print summary of M. guttatus sharing patterns
cat("=== M. guttatus Codon Preference Sharing Patterns ===\n\n")

mg_summary <- plot_data |> 
  dplyr::filter(Species == "M. guttatus") |>
  dplyr::count(Color_Category) |>
  dplyr::arrange(dplyr::desc(n))

total_aa <- nrow(mg_summary |> dplyr::summarise(total = sum(n)))

for (i in 1:nrow(mg_summary)) {
  cat_name <- mg_summary$Color_Category[i]
  count <- mg_summary$n[i]
  pct <- 100 * count / sum(mg_summary$n)
  
  cat_label <- switch(cat_name,
                      "All_three" = "Shares with all three species",
                      "Two_species" = "Shares with two species",
                      "Only_Arabidopsis" = "Shares only with A. thaliana",
                      "Only_Populus" = "Shares only with P. trichocarpa",
                      "Only_Physcomitrella" = "Shares only with P. patens",
                      "Unique" = "Unique to M. guttatus",
                      cat_name)
  
  cat(sprintf("  %-40s: %2d amino acids (%.1f%%)\n", cat_label, count, pct))
}

cat("\n")

## *****************************************************************************
## 10) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________
## Simplified version: One biplot per analysis showing preferred vs non-preferred codons

source("./src/enhanced_biplot.R")

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

# Filter to extreme expression groups only
codon_usage_CA_coord_extremes <- codon_usage_CA_coord |>
  dplyr::filter(Expression_Group %in% c("Top 5%", "Bottom 5%"))

# Prepare gene data for biplot
gene_data_ca <- codon_usage_CA_coord_extremes |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\n--- CA Analysis: Preferred vs Non-preferred Codons ---\n")
p_ca <- create_preference_biplot(
  ordination_result = codon_usage_CA,
  gene_data = gene_data_ca,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 2),
  arrow_scale = 1.0,
  title = "CA Biplot: Codon Preference (ROC Model)",
  subtitle = "Top 5% vs Bottom 5% expressed genes",
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

# Filter to extreme expression groups
rscu_PCA_coord_extremes <- rscu_PCA_coord |>
  dplyr::filter(Expression_Group %in% c("Top 5%", "Bottom 5%"))

# Prepare gene data for biplot
gene_data_pca <- rscu_PCA_coord_extremes |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create single enhanced biplot: preferred vs non-preferred codons
cat("\n--- PCA Analysis: Preferred vs Non-preferred Codons ---\n")
p_pca <- create_preference_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 2),
  arrow_scale = 1.5,
  title = "PCA Biplot: Codon Preference (ROC Model)",
  subtitle = "Top 5% vs Bottom 5% expressed genes",
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

# 10.3) Selection Load (S_load) Based Analysis ----

# Test whether genes under strong vs weak selection show distinct codon patterns

if (exists("selection_coeff_intensity") && "S_load" %in% names(selection_coeff_intensity)) {
  
  cat("\n--- CA/PCA Analysis by Selection Load ---\n")
  
  # Create S_load-based groups
  s_quantiles <- quantile(selection_coeff_intensity$S_load, 
                          probs = c(0.05, 0.95), na.rm = TRUE)
  
  selection_groups <- selection_coeff_intensity |>
    dplyr::mutate(
      S_Group = dplyr::case_when(
        S_load >= s_quantiles[2] ~ "High Selection (Top 5%)",
        S_load <= s_quantiles[1] ~ "Low Selection (Bottom 5%)",
        TRUE ~ "Intermediate"
      )
    ) |>
    dplyr::select(Gene_name, S_load, S_Group)
  
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

# 10.5) 3D visuals for PCA results ----

cat("\n=== 10.5: Creating 3D PCA Visualizations ===\n")

source("./src/create_3d_pca_video.R")

# 10.3.1) Generate dynamics 3D videos for presentation ----

cat("\n10.3.1: Creating interactive 3D PCA plot...\n")

# Ensure preferred_codons_roc has the required columns for 3D functions
preferred_codons_roc <- preferred_codons_roc |>
  dplyr::mutate(
    Codon = Preferred_Codons,
    relative_adaptiveness = 1
  )

# Create interactive 3D plot (HTML)
pca_3d_interactive <- create_3d_pca_plot(
  pca_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 2, 3),
  color_by = "expression",
  show_loadings = TRUE,
  loading_scale = 5.0,
  title = "3D PCA: RSCU Analysis with Codon Loadings (ROC)"
)

# Save interactive plot
htmlwidgets::saveWidget(
  widget = pca_3d_interactive,
  file = "./results/PCA_3D_interactive.html",
  selfcontained = TRUE
)

cat("✓ Interactive 3D plot saved: ./results/PCA_3D_interactive.html\n")
cat("  Open in browser to explore (rotate, zoom, hover for details)\n\n")

# Create rotating animation (HTML with auto-rotation)
cat("10.3.2: Creating rotating 3D animation...\n")

pca_3d_animation <- create_3d_pca_animation(
  pca_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  preferred_codons = preferred_codons_roc,
  dims = c(1, 2, 3),
  color_by = "expression",
  show_loadings = TRUE,
  loading_scale = 5.0,
  title = "3D PCA Animation - RSCU Analysis (ROC)",
  output_file = "./results/PCA_3D_animation.html",
  n_frames = 360,
  frame_duration = 50
)

# Create simple GIF animation (uses ggplot2, no heavy dependencies)
cat("8.3.3: Creating simple GIF video...\n")
if (requireNamespace("magick", quietly = TRUE)) {
  
  source("./src/create_simple_3d_gif.R")
  
  create_simple_3d_gif(
    pca_result = rscu_PCA,
    gene_data = gene_data_pca,
    codon_test_results = codon_test_results,
    preferred_codons = preferred_codons_roc,
    dims = c(1, 2, 3),
    color_by = "expression",
    show_loadings = TRUE,
    loading_scale = 5.0,
    title = "3D PCA - RSCU Analysis",
    output_file = "./results/PCA_3D_rotation.gif",
    n_frames = 60,
    width = 1000,
    height = 800,
    point_size = 1.5,
    resolution = 120
  )
  
} else {
  cat("  Skipping GIF creation (requires 'magick' package)\n")
  cat("  Install with: install.packages('magick')\n\n")
}

cat("✓ 3D visualizations complete\n\n")

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

cat("\n✓ tRNA copy number correlation analysis complete!\n")

# Analysis 2: By tRNA expression levels (all genes)
cat("\n=== Analysis 2: tRNA Expression Levels (All Genes) ===\n")
cat("Analyzing correlation using tRNA gene expression from RNA-seq\n\n")

# Prepare expression data
expression_df <- integrated_data |>
  dplyr::select(Gene_name, Expression = High_exp)

tRNA_expression_all_results <- tRNA_codon_correlation(
  codon_counts = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_expression_all",
  test_method = "spearman",
  mode = "by.expression",
  ann = "./data/Mguttatusvar_IM767_887_v2.1.gene.gff3",
  expression_data = expression_df
)

cat("\n✓ tRNA expression correlation analysis (all genes) complete!\n")

# Analysis 3: By tRNA expression levels (top 5% expressed genes only)
cat("\n=== Analysis 3: tRNA Expression Levels (Top 5% Genes) ===\n")
cat("Analyzing correlation for highly expressed genes only\n\n")

# Filter to top 5% genes
top5_threshold <- quantile(integrated_data$High_exp, probs = 0.95, na.rm = TRUE)
top5_genes <- integrated_data |>
  filter(High_exp >= top5_threshold) |>
  pull(Gene_name)

codon_usage_top5 <- codon_usage |>
  filter(Gene_name %in% top5_genes)

cat(sprintf("Analyzing %d genes in top 5%% (expression >= %.2f)\n", 
            nrow(codon_usage_top5), top5_threshold))

tRNA_expression_top5_results <- tRNA_codon_correlation(
  codon_counts = codon_usage_top5,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis_expression_top5",
  test_method = "spearman",
  mode = "by.expression",
  ann = "./data/Mguttatusvar_IM767_887_v2.1.gene.gff3",
  expression_data = expression_df
)

cat("\n✓ tRNA expression correlation analysis (top 5%) complete!\n")

# Sanity check: Does amino acid frequency match tRNA supply?
cat("\n=== Sanity Check: Amino Acid Frequency vs tRNA Supply ===\n")
cat("Testing if amino acids with more tRNA genes are used more frequently\n")
cat("This validates the tRNA adaptation hypothesis at the amino acid level\n\n")

source("./src/check_aa_trna_supply.R")

aa_trna_check <- check_aa_frequency_vs_tRNA_supply(
  codon_usage = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/aa_trna_sanity_check"
)

cat("✓ Amino acid vs tRNA supply sanity check complete!\n\n")

# ---- Final Analysis: Translational Accuracy Hypothesis ----
cat("\n=== TRANSLATIONAL ACCURACY HYPOTHESIS TEST ===\n")
cat("Testing if selection favors Watson-Crick pairing (high fidelity) over wobble pairing\n")
cat("Hypothesis: Preferred codons should use accurate Watson-Crick pairs at wobble position\n\n")

source("./src/classify_codon_anticodon_pairing.R")

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

cat("\n✓ Translational accuracy hypothesis test complete!\n")
cat("  Results saved to: ./results/tRNA_analysis_pairing/\n\n")

## *****************************************************************************
## 12) Polymorphism data integration ----
## _____________________________________________________________________________

pi_data <- fread(input = "data/all_chromosomes.bygene.pi.txt")

# Homogenizing gene names to match the previous convention

pi_data <- pi_data |>
  dplyr::select(Chr, Gene, contains("Tajima"), contains("mean")) |>
  dplyr::mutate(Gene = paste0("MgIM767.", pi_data[['Gene']])) |>
  dplyr::rename(Gene_name = Gene)

# Join polymorphism data to integrated_data
integrated_data <- integrated_data |>
  dplyr::left_join(pi_data, by = "Gene_name")

# Differences in pi at synonymous sites as a function of expression

pi_per_expression <- anova(lm(Pi_mean_4fold ~ Expression_Group, 
                           data = integrated_data))

# Post-hoc test

pi_posthoc <- TukeyHSD(aov(lm(Pi_mean_4fold ~ Expression_Group, 
                              data = integrated_data)))

# Getting required information

plot_data_pi_1 <- integrated_data |>
  dplyr::select(Gene_name, Pi_mean_4fold, Expression_Group) |>
  na.exclude()

alpha <- 0.05

plot_data_pi_1_ready <- plot_data_pi_1 |>
  dplyr::group_by(Expression_Group) |>
  dplyr::summarize(Mean_pi_4_fold = mean(Pi_mean_4fold),
                   LL = mean(Pi_mean_4fold) - qt(1 - alpha/2, (n() - 1)) * sd(Pi_mean_4fold) / sqrt(n()),
                   UL = mean(Pi_mean_4fold) + qt(1 - alpha/2, (n() - 1)) * sd(Pi_mean_4fold) / sqrt(n()))

# Visual (Mean + CI)

ggplot(data = plot_data_pi_1_ready, 
       mapping = aes(x = Expression_Group,
                     y = Mean_pi_4_fold)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 0.2) +
  theme_custom() +
  labs(title = "Nucleotide Diversity (Pi) at 4-fold Sites by Expression Group",
       x = "Expression Group",
       y = "Mean Pi at 4-fold Sites") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./results/diversity_modeling/Pi_by_expression_group_Mean_CI.pdf",
       width = 6, height = 4)



# 12.1) Tracking frequency of preferred allele as a function of expression ----

preferred_data <- read.delim("./data/all_chromosomes.codon_frequencies_preferred.txt", 
                             stringsAsFactors = FALSE) |>
  dplyr::mutate(Gene = paste0("MgIM767.", Gene))

preferred_data <- preferred_data |>
  dplyr::select(Gene, Preferred_Freq) |>
  dplyr::rename(Gene_name = Gene) |>
  dplyr::group_by(Gene_name) |>
  summarize(
    Mean_preferred_freq = mean(Preferred_Freq),
    n_codons = n()
  ) |>
  ungroup()

integrated_data <- integrated_data |>
  left_join(preferred_data) |>
  na.exclude()

summary(lm(Mean_preferred_freq ~ High_exp_log10, data = integrated_data))

# Assesing significance of expression over the detrended residuals

cat("\n=== Kruskal-Wallis Test: Frequency of preferred codons across Groups ===\n")

kw_preferred_freq <- kruskal.test(Mean_preferred_freq ~ Expression_Group, 
                                  data = integrated_data)

# Plotting and assessing significance using Dunn

print(kw_preferred_freq)
if (kw_detrended$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Perform Dunn's test with FDR correction
  dunn_result_detrended <- dunn.test::dunn.test(
    x = integrated_data$Mean_preferred_freq,
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
  dplyr::pull(Mean_preferred_freq)

middle_preferred <- integrated_data |>
  dplyr::filter(Expression_Group == "Middle 90%") |>
  dplyr::pull(Mean_preferred_freq)

bottom5_preferred <- integrated_data |>
  dplyr::filter(Expression_Group == "Bottom 5%") |>
  dplyr::pull(Mean_preferred_freq)

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
                                                   y = Mean_preferred_freq, 
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

# 12.2) Relationship between CDC-detrended and freq_preferred ----

summary(lm(CDC ~ Mean_preferred_freq, 
           data = integrated_data))

CDC_f_Mean_preferred_freq <- ggplot(integrated_data, aes(x = Mean_preferred_freq, 
                                                         y = CDC)) +
  # Use ggpointdensity for a clear view of the cluster
  geom_pointdensity(alpha = 0.5) + 
  
  # Add the linear regression line, which now shows the true effect
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  
  labs(
    title = "CDC vs. Mean_preferred_freq",
    y = "CDC",
    x = "Mean_preferred_freq"
  ) +
  theme_custom()

ggsave("./results/CDC_vs_Mean_preferred_freq.pdf", 
       CDC_f_Mean_preferred_freq, width = 8, height = 6)

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

cutoffs_exp <- quantile(integrated_data$Geom_Exp, seq(0, 1, by = 0.05),
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
    dplyr::filter(Geom_Exp > low_thr & Geom_Exp <= high_thr) |>
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

cutoffs_exp <- quantile(selection_coeff_intensity$Geom_Exp, seq(0, 1, by = 0.05),
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
  
  sub_data_genes <- selection_coeff_intensity |>
    dplyr::filter(Geom_Exp > low_thr & Geom_Exp <= high_thr)
  
  list(
    mean_S = mean(sub_data_genes$S_coeff, na.rm = TRUE),
    se_S = sd(sub_data_genes$S_coeff, na.rm = TRUE) / sqrt(sum(!is.na(sub_data_genes$S_coeff))),
    n_genes = nrow(sub_data_genes),
    mean_exp = mean(sub_data_genes$Geom_Exp, na.rm = TRUE)
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
compute_mi_analysis <- function(x, y, n_bins = 10) {
  #' Compute mutual information between two continuous variables
  #' Uses discretization with equal-frequency binning
  #' 
  #' @param x First variable (e.g., S_ROC)
  #' @param y Second variable (e.g., Gamma)
  #' @param n_bins Number of bins for discretization
  #' @return List with MI, normalized MI, and diagnostics
  
  # Remove NA values
  valid_idx <- complete.cases(x, y)
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  n <- length(x)
  
  # Equal-frequency binning (quantile-based)
  x_bins <- cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = n_bins + 1)), 
                include.lowest = TRUE, labels = FALSE)
  y_bins <- cut(y, breaks = quantile(y, probs = seq(0, 1, length.out = n_bins + 1)), 
                include.lowest = TRUE, labels = FALSE)
  
  # Compute joint and marginal probabilities
  joint_table <- table(x_bins, y_bins)
  p_xy <- joint_table / sum(joint_table)
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  
  # Compute mutual information: I(X;Y) = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
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
  nmi <- mi / sqrt(H_x * H_y)
  
  # Uncertainty coefficient (asymmetric): fraction of Y explained by X
  u_y_given_x <- mi / H_y
  u_x_given_y <- mi / H_x
  
  return(list(
    MI = mi,
    NMI = nmi,
    H_X = H_x,
    H_Y = H_y,
    U_Y_given_X = u_y_given_x,
    U_X_given_Y = u_x_given_y,
    n_samples = n,
    n_bins = n_bins
  ))
}

# Compute MI for different numbers of bins (sensitivity analysis)
mi_results <- lapply(c(5, 10, 15, 20), function(nb) {
  res <- compute_mi_analysis(scaling_df$S_ROC, scaling_df$Gamma_avg, n_bins = nb)
  res$n_bins <- nb
  res
})

cat("=== Mutual Information Results (varying bin sizes) ===\n")
for (res in mi_results) {
  cat(sprintf("  Bins=%2d: MI=%.4f bits, NMI=%.4f, U(Gamma|S)=%.3f\n", 
              res$n_bins, res$MI, res$NMI, res$U_Y_given_X))
}

# Use optimal binning (Sturges' rule approximation)
optimal_bins <- max(5, min(15, ceiling(1 + log2(nrow(scaling_df)))))
mi_optimal <- compute_mi_analysis(scaling_df$S_ROC, scaling_df$Gamma_avg, n_bins = optimal_bins)

cat(sprintf("\nOptimal Binning (%d bins):\n", optimal_bins))
cat(sprintf("  Mutual Information: %.4f bits\n", mi_optimal$MI))
cat(sprintf("  Normalized MI:      %.4f (0=independent, 1=perfect dependence)\n", mi_optimal$NMI))
cat(sprintf("  U(Gamma|S_ROC):     %.3f (fraction of Gamma entropy explained by S_ROC)\n", mi_optimal$U_Y_given_X))
cat(sprintf("  U(S_ROC|Gamma):     %.3f (fraction of S_ROC entropy explained by Gamma)\n\n", mi_optimal$U_X_given_Y))

# Permutation test for significance
cat("=== Permutation Test for MI Significance ===\n")
n_permutations <- 1000
observed_mi <- mi_optimal$MI

permuted_mi <- replicate(n_permutations, {
  y_shuffled <- sample(scaling_df$Gamma_avg)
  compute_mi_analysis(scaling_df$S_ROC, y_shuffled, n_bins = optimal_bins)$MI
})

mi_p_value <- mean(permuted_mi >= observed_mi)
mi_zscore <- (observed_mi - mean(permuted_mi)) / sd(permuted_mi)

cat(sprintf("Observed MI: %.4f bits\n", observed_mi))
cat(sprintf("Permuted MI: %.4f ± %.4f bits (mean ± SD)\n", mean(permuted_mi), sd(permuted_mi)))
cat(sprintf("Z-score:     %.2f\n", mi_zscore))
cat(sprintf("P-value:     %.4f (from %d permutations)\n", mi_p_value, n_permutations))

if (mi_p_value < 0.05) {
  cat("RESULT: Significant non-random association between S_ROC and Gamma\n\n")
} else {
  cat("RESULT: Association not significantly different from random\n\n")
}

# Compare linear vs non-linear information
cat("=== Linear vs Non-Linear Information ===\n")
spearman_cor <- cor(scaling_df$S_ROC, scaling_df$Gamma_avg, method = "spearman")
pearson_cor <- cor(scaling_df$S_ROC, scaling_df$Gamma_avg, method = "pearson")

# Spearman captures monotonic (possibly non-linear) relationships
# If NMI >> R², there's substantial non-linear information
cat(sprintf("Pearson R:   %.4f (linear correlation)\n", pearson_cor))
cat(sprintf("Pearson R²:  %.4f (linear variance explained)\n", pearson_cor^2))
cat(sprintf("Spearman ρ:  %.4f (monotonic correlation)\n", spearman_cor))
cat(sprintf("Spearman ρ²: %.4f (monotonic variance explained)\n", spearman_cor^2))
cat(sprintf("NMI:         %.4f (total dependence including non-linear)\n", mi_optimal$NMI))
cat(sprintf("\nNon-linearity indicator (NMI - R²): %.4f\n", mi_optimal$NMI - pearson_cor^2))

if (mi_optimal$NMI - pearson_cor^2 > 0.1) {
  cat("  → Substantial non-linear component in the relationship\n")
} else if (mi_optimal$NMI - pearson_cor^2 > 0.05) {
  cat("  → Moderate non-linear component\n")
} else {
  cat("  → Relationship is predominantly linear\n")
}

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
cat("\n✓ MI heatmap saved: ./results/MI_heatmap_Sroc_vs_Gamma.pdf\n")

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
cat("✓ Permutation test plot saved: ./results/MI_permutation_test.pdf\n")

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
cat("✓ MI summary saved: ./results/MI_analysis_summary.csv\n\n")

## *****************************************************************************
## 14) Diversity Hump Validation (EMPIRICAL Pi vs. Selection Intensity) ----
## _____________________________________________________________________________
## This section tests the hump pattern using OBSERVED Pi_mean_4fold values
## binned by S_coeff from AnaCoDa. This is the primary empirical validation.

# Kruskal-Wallis with Post-Hoc Dunn's Test (Mean + CI Visualization)

# Prepare Data
# Prepare Data
clean_data_4cat <- selection_coeff_intensity |>
  dplyr::filter(
    !is.na(S_coeff), !is.na(Geom_Exp), !is.na(Pi_mean_4fold),
    is.finite(S_coeff), is.finite(Geom_Exp), is.finite(Pi_mean_4fold),
    Pi_mean_4fold < 0.05
  ) |>
  dplyr::mutate(
    # Use safe internal codes (A-D) to prevent regex errors in cldList
    S_Bin_Code = cut(S_coeff, 
                     breaks = c(-Inf, 0.5, 1, 2, Inf), 
                     labels = c("A", "B", "C", "D"), 
                     include.lowest = TRUE),
    
    # Expression remains as Deciles
    Exp_Bin = cut(Geom_Exp, 
                  breaks = quantile(Geom_Exp, probs = seq(0, 1, 0.1)),
                  include.lowest = TRUE, labels = 1:10)
  )

# Define pretty labels for mapping back later
s_bin_labels <- c("A" = "0-0.5", "B" = "0.5-1", "C" = "1-2", "D" = ">2")

# Analysis Function
analyze_and_summarize <- function(data, group_col, response_col = "Pi_mean_4fold") {
  
  f <- as.formula(paste(response_col, "~", group_col))
  dunn_res <- FSA::dunnTest(f, data = data, method = "bh")
  
  # Generate Letters
  cld <- rcompanion::cldList(P.adj ~ Comparison, data = dunn_res$res, threshold = 0.05)
  
  # Clean CLD
  cld_clean <- cld |>
    dplyr::rename(!!group_col := Group, Letter = Letter) |>
    dplyr::mutate(!!group_col := trimws(as.character(!!sym(group_col)))) |>
    dplyr::select(!!group_col, Letter)
  
  # Summary Stats
  summary_stats <- data |>
    dplyr::group_by(!!sym(group_col)) |>
    dplyr::summarise(
      Mean = mean(!!sym(response_col), na.rm = TRUE),
      SD = sd(!!sym(response_col), na.rm = TRUE),
      N = n(),
      SE = SD / sqrt(N),
      CI_95 = 1.96 * SE,
      Upper = Mean + CI_95,
      Lower = Mean - CI_95,
      .groups = "drop"
    ) |>
    dplyr::mutate(!!group_col := as.character(!!sym(group_col))) |>
    dplyr::left_join(cld_clean, by = group_col) |>
    dplyr::mutate(!!group_col := factor(!!sym(group_col), 
                                        levels = levels(data[[group_col]])))
  
  return(summary_stats)
}

# Run Analysis
stats_S <- analyze_and_summarize(clean_data_4cat, "S_Bin_Code")
stats_Exp <- analyze_and_summarize(clean_data_4cat, "Exp_Bin")

# Map codes back to pretty labels for S_ROC
stats_S$Pretty_Label <- s_bin_labels[as.character(stats_S$S_Bin_Code)]
stats_S$Pretty_Label <- factor(stats_S$Pretty_Label, levels = c("0-0.5", "0.5-1", "1-2", ">2"))

# For Expression, the labels are already correct
stats_Exp$Pretty_Label <- stats_Exp$Exp_Bin

# Plotting Function
plot_mean_ci <- function(summary_data, title, subtitle, x_lab) {
  
  y_max <- max(summary_data$Upper, na.rm = TRUE)
  y_min <- min(summary_data$Lower, na.rm = TRUE)
  text_offset <- (y_max - y_min) * 0.15
  
  ggplot(summary_data, aes(x = Pretty_Label, y = Mean)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 0.8) +
    geom_point(aes(fill = Pretty_Label), shape = 21, size = 5, color = "black", stroke = 1) +
    geom_text(aes(y = Upper + text_offset, label = Letter), size = 6, fontface = "bold") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = "Nucleotide Diversity (Pi) ± 95% CI"
    ) +
    theme_custom() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
}

# Generate and Save
p_select <- plot_mean_ci(stats_S, 
                         "Diversity vs Selection (4 Categories)", 
                         "Mean ± 95% CI | Letters: Post-hoc Dunn's Test (p < 0.05)",
                         "Selection Intensity (S_ROC)")

p_exp <- plot_mean_ci(stats_Exp, 
                      "Diversity vs Expression (Deciles)", 
                      "Mean ± 95% CI | Letters: Post-hoc Dunn's Test (p < 0.05)",
                      "Expression Decile")

ggsave("results/Results_Pi_vs_S_4Cat.pdf", p_select, width = 7, height = 6)
ggsave("results/Results_Pi_vs_Exp_Deciles.pdf", p_exp, width = 8, height = 6)

cat("\nPlots saved successfully using the 4-category S_ROC binning.\n")

## *****************************************************************************
## 15) Diversity across different genomic compartment ----
## _____________________________________________________________________________