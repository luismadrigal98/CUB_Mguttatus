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
                        'admisc')

set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = T, n_cores = 10)

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

exp_data_bud <- read.table(file = "./data/bud_gene_expression_cpm_remapped.txt",
                           header = T) |>
  dplyr::rename(Exp_bud = Expression)

exp_data_leaf <- read.table(file = "./data/leaf_gene_expression_mean_cpm_renamed.txt",
                            header = T) |>
  dplyr::rename(Exp_leaf = Expression)

# Combining CUB metric with expression profiles for buds and leafs

exp_complete <- dplyr::full_join(exp_data_leaf, exp_data_bud, by = "Gene")

exp_complete <- exp_complete |> 
  dplyr::rowwise() |>
  dplyr::mutate(
    High_exp = max(Exp_leaf, Exp_bud, na.rm = TRUE),
    
    Source_High_exp = case_when(
      # If both are NA, source is NA
      is.na(Exp_leaf) & is.na(Exp_bud) ~ NA_character_,
      
      # If bud is NA, leaf must be the max
      is.na(Exp_bud) ~ "Leaf",
      
      # If leaf is NA, bud must be the max
      is.na(Exp_leaf) ~ "Bud",
      
      # If leaf is greater or equal (handles ties)
      Exp_leaf >= Exp_bud ~ "Leaf",
      
      # Otherwise, bud must be greater
      Exp_bud > Exp_leaf ~ "Bud"
    )
  ) |>
  # Fix the -Inf from max(NA, NA, na.rm=T)
  dplyr::mutate(
    High_exp = if_else(is.infinite(High_exp), NA_real_, High_exp)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(Source_High_exp = as.factor(Source_High_exp)) |>
  na.exclude()

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
    Gene_name_clean = sub("\\.1$", "", Gene_name),  # Remove .1 suffix
    Total_Codons = rowSums(across(all_of(codon_columns)), na.rm = TRUE),
    CDS_length_nt = Total_Codons * 3,  # nucleotides
    CDS_length_aa = Total_Codons        # amino acids (codons)
  ) |>
  dplyr::select(Gene_name_clean, Total_Codons, CDS_length_nt, CDS_length_aa) |>
  dplyr::rename(Gene_name = Gene_name_clean)

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

# Re-plotting ENC-based neutrality plot highlighting the significant genes with CDC ----

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

# Linear models ----

integrated_data <- integrated_data |> 
  left_join(enc_cdc_data |> dplyr::select(Gene_name, CDC), by = "Gene_name")
  
CDC_vs_exp <- lm(CDC ~ High_exp_log2 + CDS_length_nt, 
                 data = integrated_data)
summary(CDC_vs_exp)

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

cdc_model_beta <- gam(CDC ~ High_exp_log2 + s(CDS_length_nt), 
                      data = integrated_data, family = betar(link = "logit"))

summary(cdc_model_beta)

# Plotting detrended ENC against expression

confounder_model_gam <- gam(CDC ~ s(CDS_length_nt),
                            data = integrated_data)

integrated_data$CDC_detrended <- residuals(confounder_model_gam)

p_detrended <- ggplot(integrated_data, aes(x = High_exp_log2, y = CDC_detrended)) +
  # Use ggpointdensity for a clear view of the cluster
  geom_pointdensity(alpha = 0.5) + 
  
  # Add the linear regression line, which now shows the true effect
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  
  labs(
    title = "Detrended CDC vs. Gene Expression",
    subtitle = "Showing CDC after accounting for non-linear effects of gene length",
    y = "ENC Residuals (Detrended)",
    x = "log2(Expression + 1)"
  ) +
  theme_custom()

print(p_detrended)
ggsave("./results/ENC_detrended_vs_expression.pdf", p_detrended, width = 8, height = 6)

# Define expression groups: Top 5% vs Bottom 5% (extreme comparison) ----

top_5_cutoff <- quantile(integrated_data$High_exp_log2, probs = 0.95)
bottom_5_cutoff <- quantile(integrated_data$High_exp_log2, probs = 0.05)

integrated_data$Expression_Group <- case_when(
  integrated_data$High_exp_log2 >= top_5_cutoff ~ "Top 5%",
  integrated_data$High_exp_log2 <= bottom_5_cutoff ~ "Bottom 5%",
  TRUE ~ "Middle 90%"
)

# Boxplot comparison

p_boxplot <- ggplot(integrated_data, aes(x = Expression_Group, y = CDC, fill = Expression_Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "CDC by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR",
       y = "CDC",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/CDC_by_expression_group.pdf", p_boxplot, width = 8, height = 6)

# Statistical tests for three groups
cat("\n=== Kruskal-Wallis Test: ENC across All Three Groups ===\n")
cat("H0: All three groups have the same median ENC\n")
kw_test_cdc <- kruskal.test(CDC ~ Expression_Group, data = integrated_data)
print(kw_test_cdc)

if (kw_test_cdc$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Install and load dunn.test if not available
  if (!require("dunn.test", quietly = TRUE)) {
    cat("Installing dunn.test package...\n")
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
    library(dunn.test)
  }
  
  # Perform Dunn's test with FDR correction
  dunn_result_cdc <- dunn.test::dunn.test(
    x = integrated_data$CDC,
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

cat("\n=== Summary Statistics ===\n")
summary_stats <- integrated_data |>
  group_by(Expression_Group) |>
  summarise(
    n = n(),
    mean_CDC = mean(CDC, na.rm = TRUE),
    median_CDC = median(CDC, na.rm = TRUE),
    sd_CDC = sd(CDC, na.rm = TRUE),
    mean_Expression = mean(High_exp, na.rm = TRUE)
  )
print(summary_stats)

# Effect sizes for pairwise comparisons
cat("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

# Get ENC values for each group
top5_cdc <- integrated_data |> filter(Expression_Group == "Top 5%") |> pull(CDC)
middle_cdc <- integrated_data |> filter(Expression_Group == "Middle 90%") |> pull(CDC)
bottom5_cdc <- integrated_data |> filter(Expression_Group == "Bottom 5%") |> pull(CDC)

# Calculate effect sizes
if (length(top5_cdc) > 0 && length(middle_cdc) > 0) {
  d_top_middle <- cohens_d_calc(top5_cdc, middle_cdc)
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle))
}

if (length(top5_cdc) > 0 && length(bottom5_cdc) > 0) {
  d_top_bottom <- cohens_d_calc(top5_cdc, bottom5_cdc)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom))
}

if (length(middle_cdc) > 0 && length(bottom5_cdc) > 0) {
  d_middle_bottom <- cohens_d_calc(middle_cdc, bottom5_cdc)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom))
}

cat("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

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

p_boxplot_detrended <- ggplot(integrated_data, aes(x = Expression_Group, y = CDC_detrended, fill = Expression_Group)) +
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
       y = "CDC Residuals (detrended)",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/Detrended_ENC_by_expression_group.pdf", 
       p_boxplot_detrended, width = 8, height = 6)

# Get ENC values for each group
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
reference_genes <- integrated_data |>
  filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

cat(sprintf("Using %d highly expressed genes as reference set\n", length(reference_genes)))

# Remove .1 suffix from codon_usage gene names to match gene-level IDs
# (codon_usage has transcript IDs like MgIM767.10G127000.1,
#  expression data has gene IDs like MgIM767.10G127000)
codon_usage <- codon_usage |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", Gene_name))

cat(sprintf("Converted codon usage transcript IDs to gene IDs\n"))

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

ggplot(data = w_table, mapping = aes(x = reorder(codon, relative_adaptiveness), 
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
  labs(y = "Relative Adaptiveness (w)", x = "Codon",
       title = "Codon Preference in Highly Expressed Genes")

ggsave("./results/optimal_codons_relative_adaptiveness.pdf", 
       width = 12, height = 10)

# Merge CAI with expression and ENC data
integrated_data <- integrated_data |>
  left_join(cai_values, by = "Gene_name")

# Save results
write.csv(integrated_data, "./results/expression_enc_cai.csv", row.names = FALSE)
write.csv(w_table, "./results/optimal_codons_weights.csv", row.names = FALSE)

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
                               "Freq_Rest" = "Rest"))

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

# 7.2.1) Preferred codon usage: Selected vs Neutral genes ---- ----

# Get preferred codons (w = 1.0 from CAI)
preferred_codons_vec <- w_table |>
  dplyr::filter(relative_adaptiveness == 1.0) |>
  dplyr::pull(codon)

cat(sprintf("Using %d preferred codons (w = 1.0)\n\n", length(preferred_codons_vec)))

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

selected_aa <- count_preferred_by_aa(top5_genes, preferred_codons_vec, genetic_code_dna_long)
selected_aa$Group <- "Selected (Top 5%)"

rest_aa <- count_preferred_by_aa(rest_genes, preferred_codons_vec, genetic_code_dna_long)
rest_aa$Group <- "Rest (Bottom 95%)"

# Combine for comparison
comparison_table <- selected_aa |>
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
write.csv(comparison_table, "./results/preferred_codon_usage_selected_vs_neutral.csv",
          row.names = FALSE)

cat("\n✓ Results saved: ./results/preferred_codon_usage_selected_vs_neutral.csv\n\n")

# Print table
cat("=== Preferred Codon Usage: Selected vs Rest ===\n\n")
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
cat(sprintf("Mean proportion preferred (Top 5%%): %.4f\n", 
            mean(comparison_table$Selected_prop, na.rm = TRUE)))
cat(sprintf("Mean proportion preferred (Rest 95%%): %.4f\n", 
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

if (wilcox_test$p.value < 0.001) {
  cat("  *** Highly significant (p < 0.001)\n")
  cat("  → Top 5%% genes use MORE preferred codons than rest\n")
} else if (wilcox_test$p.value < 0.05) {
  cat("  * Significant (p < 0.05)\n")
} else {
  cat("  Not significant (p >= 0.05)\n")
}

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

# Now create the same plot with CORRECTED preferred codons (after section 7.2.4)
# This will be added in section 7.2.5 for comparison

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
            length(preferred_codons_vec),
            100 * nrow(sig_preferred) / length(preferred_codons_vec)))

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

# 7.2.4) Preferred codons (corrected) ----
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

preferred_codons_corrected <- codon_combined |>
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
decision_summary <- preferred_codons_corrected |>
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
      final_choice = preferred_codons_corrected$Codon[preferred_codons_corrected$AA == AA[1]][1]
    )
  
  print(conflict_summary)
  cat("\n")
}

# Save corrected preferred codons
write.csv(preferred_codons_corrected, 
          "./results/preferred_codons_corrected.csv",
          row.names = FALSE)

cat("✓ Corrected preferred codons saved: ./results/preferred_codons_corrected.csv\n\n")

# Show final preferred codon set
cat("=== Final Preferred Codons (Corrected) ===\n")
final_table <- preferred_codons_corrected |>
  dplyr::select(AA, Codon, w = relative_adaptiveness, 
                Enriched = is_enriched, Diff = Difference, Selection_Rationale)
print(final_table, n = Inf)
cat("\n")

# Use corrected set for downstream analysis
preferred_codons_mg <- preferred_codons_corrected

# 7.2.5) Recalculate CAI with corrected preferred codons ----
cat("\n=== 7.2.5: Recalculating CAI with Corrected Preferred Codons ===\n\n")

# Create corrected w-table from preferred_codons_corrected
w_table_corrected <- preferred_codons_corrected |>
  dplyr::select(amino_acid = AA, codon = Codon, relative_adaptiveness)

# For codons not marked as preferred, calculate their relative adaptiveness
# based on their frequency relative to the preferred codon
all_codons_by_aa <- codon_combined |>
  dplyr::select(AA, Codon, Selected_Prop) |>
  dplyr::group_by(AA) |>
  dplyr::mutate(
    max_prop = max(Selected_Prop, na.rm = TRUE),
    w_corrected = Selected_Prop / max_prop
  ) |>
  dplyr::ungroup() |>
  dplyr::select(amino_acid = AA, codon = Codon, relative_adaptiveness = w_corrected)

# Get Top 5% gene names (matching codon_usage IDs - without .1 suffix)
top5_genes_for_cai <- integrated_data |>
  dplyr::filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

cat(sprintf("Using %d Top 5%% genes as reference for corrected CAI\n", 
            length(top5_genes_for_cai)))

# Recalculate CAI using the corrected w-table
cat("Recalculating CAI with corrected w-values...\n")
cai_results_corrected <- calculate_cai(
  codon_counts = codon_usage,
  reference_genes = top5_genes_for_cai,
  genetic_code = genetic_code_dna_long
)

# Extract corrected CAI values
cai_values_corrected <- cai_results_corrected$cai_values |>
  dplyr::rename(CAI_corrected = CAI)

# Merge with integrated data
integrated_data <- integrated_data |>
  dplyr::left_join(cai_values_corrected, by = "Gene_name")

# Compare original vs corrected CAI
cat("\n=== CAI Comparison: Original vs Corrected ===\n")
cai_comparison <- integrated_data |>
  dplyr::select(Gene_name, Expression_Group, CAI, CAI_corrected) |>
  dplyr::filter(!is.na(CAI) & !is.na(CAI_corrected))

# Correlation
cor_cai <- cor(cai_comparison$CAI, cai_comparison$CAI_corrected, 
               use = "complete.obs")
cat(sprintf("Correlation between original and corrected CAI: %.4f\n\n", cor_cai))

# Summary by expression group
cai_comparison_summary <- cai_comparison |>
  dplyr::group_by(Expression_Group) |>
  dplyr::summarise(
    n = n(),
    mean_CAI_original = mean(CAI, na.rm = TRUE),
    mean_CAI_corrected = mean(CAI_corrected, na.rm = TRUE),
    median_CAI_original = median(CAI, na.rm = TRUE),
    median_CAI_corrected = median(CAI_corrected, na.rm = TRUE),
    diff_mean = mean_CAI_corrected - mean_CAI_original
  )

print(cai_comparison_summary)
cat("\n")

# Plot comparison
p_cai_comparison <- ggplot(cai_comparison, 
                           aes(x = CAI, y = CAI_corrected, color = Expression_Group)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                 "Bottom 5%" = "#377EB8",
                                 "Middle 90%" = "#999999")) +
  labs(
    title = "CAI: Original vs Corrected",
    subtitle = sprintf("Correlation: %.3f | Dashed line = perfect agreement", cor_cai),
    x = "Original CAI (w=1 from raw frequencies)",
    y = "Corrected CAI (w=1 from enrichment-corrected preferences)",
    color = "Expression Group"
  ) +
  theme_custom() +
  theme(legend.position = "right")

ggsave("./results/CAI_original_vs_corrected.pdf", p_cai_comparison, 
       width = 10, height = 7)

cat("✓ Plot saved: ./results/CAI_original_vs_corrected.pdf\n\n")

# Analyze preferred codon usage: Top 5% vs Bottom 5%
cat("=== Preferred Codon Usage: Top 5% vs Bottom 5% ===\n\n")

# Get gene lists
top5_genes_list <- integrated_data |>
  dplyr::filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

bottom5_genes_list <- integrated_data |>
  dplyr::filter(Expression_Group == "Bottom 5%") |>
  pull(Gene_name)

# Count preferred codon usage
cat("Counting preferred codon usage in Top 5% genes...\n")
preferred_usage_top5 <- count_preferred_codons(
  top5_genes_list, codon_usage, preferred_codons_corrected
)

cat("Counting preferred codon usage in Bottom 5% genes...\n")
preferred_usage_bottom5 <- count_preferred_codons(
  bottom5_genes_list, codon_usage, preferred_codons_corrected
)

# Summary statistics
cat("\n=== Summary: Preferred Codon Usage ===\n")
cat(sprintf("\nTop 5%% genes (n=%d):\n", nrow(preferred_usage_top5)))
cat(sprintf("  Mean proportion of preferred codons: %.4f (SD = %.4f)\n",
            mean(preferred_usage_top5$Preferred_Proportion, na.rm = TRUE),
            sd(preferred_usage_top5$Preferred_Proportion, na.rm = TRUE)))
cat(sprintf("  Median proportion: %.4f\n",
            median(preferred_usage_top5$Preferred_Proportion, na.rm = TRUE)))

cat(sprintf("\nBottom 5%% genes (n=%d):\n", nrow(preferred_usage_bottom5)))
cat(sprintf("  Mean proportion of preferred codons: %.4f (SD = %.4f)\n",
            mean(preferred_usage_bottom5$Preferred_Proportion, na.rm = TRUE),
            sd(preferred_usage_bottom5$Preferred_Proportion, na.rm = TRUE)))
cat(sprintf("  Median proportion: %.4f\n",
            median(preferred_usage_bottom5$Preferred_Proportion, na.rm = TRUE)))

# Statistical test
cat("\n=== Statistical Test ===\n")
wilcox_test_pref <- wilcox.test(
  preferred_usage_top5$Preferred_Proportion,
  preferred_usage_bottom5$Preferred_Proportion,
  alternative = "greater"
)

cat(sprintf("Wilcoxon rank-sum test (Top 5%% > Bottom 5%%):\n"))
cat(sprintf("  W = %.2f, p-value = %.2e\n", 
            wilcox_test_pref$statistic, wilcox_test_pref$p.value))

# Effect size
d_preferred <- cohens_d_calc(
  preferred_usage_top5$Preferred_Proportion,
  preferred_usage_bottom5$Preferred_Proportion
)
cat(sprintf("  Cohen's d = %.3f\n", d_preferred))
cat(sprintf("  Interpretation: %s\n",
            ifelse(abs(d_preferred) < 0.2, "negligible",
                   ifelse(abs(d_preferred) < 0.5, "small",
                          ifelse(abs(d_preferred) < 0.8, "medium", "large")))))

# Visualization
preferred_usage_combined <- rbind(
  preferred_usage_top5 |> dplyr::mutate(Group = "Top 5%"),
  preferred_usage_bottom5 |> dplyr::mutate(Group = "Bottom 5%")
)

p_preferred_usage <- ggplot(preferred_usage_combined, 
                            aes(x = Group, y = Preferred_Proportion, fill = Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.3, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", "Bottom 5%" = "#377EB8")) +
  labs(
    title = "Preferred Codon Usage: Top 5% vs Bottom 5%",
    subtitle = sprintf("Wilcoxon p = %.2e, Cohen's d = %.3f", 
                       wilcox_test_pref$p.value, d_preferred),
    x = "Expression Group",
    y = "Proportion of Preferred Codons",
    caption = sprintf("%d corrected preferred codons (enrichment-based)", 
                      nrow(preferred_codons_corrected))
  ) +
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/preferred_codon_usage_top5_vs_bottom5.pdf", 
       p_preferred_usage, width = 8, height = 6)

cat("\n✓ Plot saved: ./results/preferred_codon_usage_top5_vs_bottom5.pdf\n\n")

# Save preferred codon usage data
preferred_usage_summary <- preferred_usage_combined |>
  dplyr::select(Gene_name, Group, Total_Codons, Preferred_Codons, Preferred_Proportion)

write.csv(preferred_usage_summary, 
          "./results/preferred_codon_usage_by_expression.csv",
          row.names = FALSE)

cat("✓ Preferred codon usage data saved: ./results/preferred_codon_usage_by_expression.csv\n\n")

# Create amino acid-level comparison plot with corrected codons
cat("Creating amino acid-level comparison plot with corrected preferred codons...\n")

# Calculate preferred codon usage per amino acid for corrected set
# Reuse the count_preferred_by_aa function but with corrected codons
preferred_codons_corrected_vec <- preferred_codons_corrected$Codon

selected_aa_corrected <- count_preferred_by_aa(top5_genes, 
                                              preferred_codons_corrected_vec, 
                                              genetic_code_dna_long)
selected_aa_corrected$Group <- "Selected (Top 5%)"

rest_aa_corrected <- count_preferred_by_aa(rest_genes, 
                                          preferred_codons_corrected_vec, 
                                          genetic_code_dna_long)
rest_aa_corrected$Group <- "Rest (Bottom 95%)"

# Combine for comparison
comparison_table_corrected <- selected_aa_corrected |>
  dplyr::select(Amino_Acid, N_synonymous, Preferred_codons, 
                Selected_count = Preferred_count, 
                Selected_prop = Prop_preferred) |>
  dplyr::left_join(
    rest_aa_corrected |> dplyr::select(Amino_Acid, 
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
write.csv(comparison_table_corrected, 
          "./results/preferred_codon_usage_selected_vs_neutral_CORRECTED.csv",
          row.names = FALSE)

cat("✓ Results saved: ./results/preferred_codon_usage_selected_vs_neutral_CORRECTED.csv\n\n")

# Create corrected visualization
p_comparison_corrected <- ggplot(comparison_table_corrected, 
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
  labs(title = "Preferred Codon Usage: Top 5% vs Rest (CORRECTED)",
       subtitle = "Proportion of corrected preferred codons per amino acid (after enrichment correction)",
       x = "Amino Acid (ordered by difference)",
       y = "Proportion of Preferred Codons") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("./results/preferred_codon_usage_comparison_CORRECTED.pdf", 
       p_comparison_corrected,
       width = 10, height = 6)

cat("✓ Plot saved: ./results/preferred_codon_usage_comparison_CORRECTED.pdf\n\n")

# Statistical test for corrected version
wilcox_test_corrected <- wilcox.test(comparison_table_corrected$Selected_prop, 
                                     comparison_table_corrected$Rest_prop,
                                     paired = TRUE)

cat("=== Corrected Preferred Codons: Statistical Test ===\n")
cat(sprintf("Wilcoxon signed-rank test (paired by amino acid):\n"))
cat(sprintf("  V = %.1f, p-value = %.2e\n", 
            wilcox_test_corrected$statistic, 
            wilcox_test_corrected$p.value))

if (wilcox_test_corrected$p.value < 0.001) {
  cat("  *** Highly significant (p < 0.001)\n")
  cat("  → Top 5%% genes use MORE corrected preferred codons than rest\n\n")
} else if (wilcox_test_corrected$p.value < 0.05) {
  cat("  * Significant (p < 0.05)\n\n")
} else {
  cat("  Not significant (p >= 0.05)\n\n")
}

## 7.2.6) Comparing preferred codon of Mimulus guttatus to other plants ----

# Use w_table from CAI analysis (already calculated preferred codons)
cat("Using optimal codons from corrected reference set...\n")

# Get preferred codons (those with relative_adaptiveness == 1.0)
preferred_codons_comparative <- preferred_codons_mg |>
  dplyr::mutate(Codon_RNA = gsub("T", "U", Codon)) |>
  dplyr::select(Amino_Acid = AA, Codon_RNA, relative_adaptiveness)

# Collapse amino acids with six codons back into six, based on relative adaptiveness

preferred_codons_comparative <- preferred_codons_comparative |>
  dplyr::mutate(AA_root = sapply(preferred_codons_comparative$Amino_Acid, 
                                 function(x) 
                                 {
                                   unlist(strsplit(x, "_"))[1]
                                 }))

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
                   dplyr::desc(relative_adaptiveness)) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(Amino_Acid = !!sym(AA_family_col), 
                  Codon_RNA, relative_adaptiveness)
  
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
print(preferred_codons_mg |> dplyr::select(Amino_Acid, Codon = Codon_RNA, Weight = relative_adaptiveness))

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
  theme_minimal(base_size = 11) +
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
  theme_minimal(base_size = 12) +
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
## 8) Binomial/Multinomial modeling ----
## _____________________________________________________________________________

integrated_data <- as.data.table(integrated_data)

# 1. Get all unique, synonymous amino acid families from your genetic code
all_families <- unique(genetic_code_dna_long[!(genetic_code_dna_long %in% c("STOP", "Met", "Trp"))])
# [1] "Phe" "Leu_2" "Ser_4" "Tyr" "Cys" "Leu_4" "Pro" "His" ...

# 2. Run the model for every family
#    This will take some time to run
all_results_list <- lapply(all_families, function(fam) {
  cat("Fitting model for family:", fam, "\n")
  tryCatch(
    fit_bi_multinom_family_model(
      family_name = fam,
      genetic_code = genetic_code_dna_long,
      usage_dt = codon_usage,
      meta_dt = integrated_data,
      preferred_codons_df = preferred_codons_mg
    ),
    error = function(e) {
      cat("ERROR fitting", fam, ":", e$message, "\n")
      return(NULL)
    }
  )
})

names(all_results_list) <- all_families

# Example
Ala_data <- integrated_data |>
  dplyr::select(Gene_name, GC3s, CDS_length_nt, High_exp_log2)

Ala_data <- Ala_data |>
  left_join(codon_usage[, c("Gene_name", "GCT", "GCC", "GCA", "GCG")])

Ala_data_clean <- Ala_data |>
  filter(
    !is.na(GCT) & !is.na(GCC) & !is.na(GCA) & !is.na(GCG) &
      !is.na(High_exp_log2) & !is.na(GC3s) & !is.na(CDS_length_nt)
  ) |>
  mutate(
    total_Ala = GCT + GCC + GCA + GCG
  ) |>
  filter(total_Ala > 0)  # Remove genes with no Ala codons

cat("Original rows:", nrow(Ala_data), "\n")
cat("After cleaning:", nrow(Ala_data_clean), "\n")
cat("Removed:", nrow(Ala_data) - nrow(Ala_data_clean), "rows\n\n")

# Calculate proportions for each codon
Ala_data_clean <- Ala_data_clean |>
  mutate(
    prop_GCT = GCT / total_Ala,
    prop_GCC = GCC / total_Ala,
    prop_GCA = GCA / total_Ala,
    prop_GCG = GCG / total_Ala
  )

# Use cleaned data for analysis
Ala_data <- Ala_data_clean

# Visualize: How do codon proportions change with expression?
plot_data <- Ala_data |>
  dplyr::select(High_exp_log2, prop_GCT, prop_GCC, prop_GCA, prop_GCG) |>
  pivot_longer(cols = starts_with("prop_"), 
               names_to = "Codon", 
               values_to = "Proportion") |>
  dplyr::mutate(Codon = gsub("prop_", "", Codon))

ggplot(plot_data, aes(x = High_exp_log2, y = Proportion, color = Codon)) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Alanine Codon Usage vs Expression",
    subtitle = "Do certain codons increase in frequency at high expression?",
    x = "Gene Expression (log2)",
    y = "Proportion of Alanine Codons"
  ) +
  theme_custom()

cat("Fitting multinomial GAM...\n")

# Try with explicit na.action
model_ala <- vgam(
  cbind(GCT, GCC, GCA, GCG) ~ High_exp_log2 + s(GC3s, df = 4) + s(CDS_length_nt, df = 4),
  family = multinomial(refLevel = 1),  # GCT is reference
  data = Ala_data,
  weights = total_Ala,
  na.action = na.exclude
)

plotvgam(model_ala)


cat("Model fitted successfully!\n\n")

# View summary
summary(model_ala)

# Method 2: Get standard errors from the vcov matrix
coef_est <- coef(model_ala)
se <- sqrt(diag(vcov(model_ala)))

# Create a coefficient table manually
coef_table <- data.frame(
  Estimate = coef_est,
  Std_Error = se,
  z_value = coef_est / se,
  p_value = 2 * pnorm(-abs(coef_est / se))
)

# Show just the High_exp_log2 terms (selection effects)
coef_table[grep("High_exp_log2", rownames(coef_table)), ]

## *****************************************************************************
## 8.1) Pairwise Binomial Regression Analysis ----
## _____________________________________________________________________________

cat("\n\n===============================================================\n")
cat("SECTION 8.1: PAIRWISE BINOMIAL REGRESSION ANALYSIS\n")
cat("===============================================================\n\n")

cat("This section fits pairwise binomial regressions for each amino acid family\n")
cat("to model codon choice as a function of expression level (selection),\n")
cat("while controlling for confounding effects of gene length and GC content.\n\n")

cat("Two parallel approaches are implemented:\n")
cat("  1. GAM approach: Uses smoothers s() for non-linear confounders\n")
cat("  2. GLM approach: Uses Box-Cox transformed confounders for interpretability\n\n")

# Source the necessary functions
source("./src/fit_pairwise_binomial_models.R")
source("./src/plot_scurves.R")
source("./src/check_concordance.R")

# Get all synonymous amino acid families (exclude Met, Trp, STOP)
all_families <- unique(genetic_code_dna_long[
  !(genetic_code_dna_long %in% c("STOP", "Met", "Trp"))
])

cat(sprintf("Analyzing %d amino acid families with synonymous codons\n\n", 
            length(all_families)))

cat("Using existing preferred codons as baseline for interpretability.\n")
cat("This ensures all non-preferred codons get negative slopes.\n\n")

# 8.1.1) GAM-based approach ----

cat("=== 8.1.1: GAM Approach (Smoothers for Confounders) ===\n\n")

all_gam_results <- lapply(all_families, function(fam) {
  cat(sprintf("Fitting GAM models for family: %s\n", fam))
  tryCatch(
    fit_pairwise_gams(
      family_name = fam,
      genetic_code = genetic_code_dna_long,
      usage_dt = codon_usage,
      meta_dt = integrated_data,
      preferred_codons_df = preferred_codons_mg  # Use existing preferred codons
    ),
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NULL)
    }
  )
})

names(all_gam_results) <- all_families

# Remove NULL results
all_gam_results <- all_gam_results[!sapply(all_gam_results, is.null)]

cat(sprintf("\n✓ Successfully fitted models for %d families\n\n", 
            length(all_gam_results)))

# Aggregate coefficient tables
master_gam_table <- data.table::rbindlist(
  lapply(all_gam_results, function(x) x$coefficients)
)

# Save results
write.csv(master_gam_table, 
          "./results/section_8.1_GAM_selection_coefficients.csv",
          row.names = FALSE)

cat("✓ GAM coefficients saved: ./results/section_8.1_GAM_selection_coefficients.csv\n\n")

# Find data-driven preferred codons (highest selection slope per family)
gam_preferred_codons <- master_gam_table |>
  dplyr::group_by(Family) |>
  dplyr::filter(Selection_Slope == max(Selection_Slope)) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(Family, Codon, Selection_Slope, p_value)

cat("=== Data-Driven Preferred Codons (GAM approach) ===\n")
print(gam_preferred_codons, n = Inf)
cat("\n")

write.csv(gam_preferred_codons,
          "./results/section_8.1_GAM_preferred_codons.csv",
          row.names = FALSE)

cat("Note: GAM results saved for comparison, but GLM is the gold standard\n\n")

# 8.1.2) GLM approach with Box-Cox transformation ----

cat("\n=== 8.1.2: GLM Approach (Box-Cox Transformed Confounders) ===\n\n")
cat("This approach transforms confounders to achieve linearity,\n")
cat("avoiding the need for GAM smoothers and improving interpretability.\n\n")

all_glm_results <- lapply(all_families, function(fam) {
  cat(sprintf("Fitting GLM models with Box-Cox for family: %s\n", fam))
  tryCatch(
    fit_pairwise_glms(
      family_name = fam,
      genetic_code = genetic_code_dna_long,
      usage_dt = codon_usage,
      meta_dt = integrated_data,
      boxcox_confounders = TRUE,
      preferred_codons_df = preferred_codons_mg  # Use existing preferred codons
    ),
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NULL)
    }
  )
})

names(all_glm_results) <- all_families

# Remove NULL results
all_glm_results <- all_glm_results[!sapply(all_glm_results, is.null)]

cat(sprintf("\n✓ Successfully fitted GLM models for %d families\n\n", 
            length(all_glm_results)))

# Aggregate coefficient tables
master_glm_table <- data.table::rbindlist(
  lapply(all_glm_results, function(x) x$coefficients)
)

# Save results
write.csv(master_glm_table, 
          "./results/section_8.1_GLM_BoxCox_selection_coefficients.csv",
          row.names = FALSE)

cat("✓ GLM coefficients saved: ./results/section_8.1_GLM_BoxCox_selection_coefficients.csv\n\n")

# Find data-driven preferred codons (GLM approach)
glm_preferred_codons <- master_glm_table |>
  dplyr::group_by(Family) |>
  dplyr::filter(Selection_Slope == max(Selection_Slope)) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(Family, Codon, Selection_Slope, p_value)

cat("=== Data-Driven Preferred Codons (GLM Box-Cox approach) ===\n")
print(glm_preferred_codons, n = Inf)
cat("\n")

write.csv(glm_preferred_codons,
          "./results/section_8.1_GLM_preferred_codons.csv",
          row.names = FALSE)

cat("Note: GLM detected 3 additional significant cases vs GAM - using as gold standard\n\n")

# 8.1.2.1) Update preferred codons based on GLM results (Gold Standard) ----

cat("\n=== 8.1.2.1: Updating Preferred Codons Based on GLM Significance Patterns ===\n\n")
cat("GLM models are used as gold standard (detected 3 additional significant cases vs GAM)\n")
cat("This analysis identifies:\n")
cat("  - Families with single clear preference\n")
cat("  - Families with multiple co-optimal codons\n")
cat("  - Families with no clear preference (all non-significant)\n\n")

source("./src/update_preferred_codons.R")

# Update preferences based on GLM significance patterns
preferences_updated_glm <- update_preferred_codons_from_models(
  model_results = all_glm_results,
  existing_preferred_codons = preferred_codons_mg
)

# Save updated preferences
write.csv(preferences_updated_glm,
          "./results/section_8.1_GLM_updated_preferences_patterns.csv",
          row.names = FALSE)

cat("✓ Updated preference patterns saved\n\n")

# Create formatted table for downstream use
preferred_codons_updated <- create_preferred_codons_table(preferences_updated_glm)

write.csv(preferred_codons_updated,
          "./results/section_8.1_GLM_preferred_codons_updated.csv",
          row.names = FALSE)

cat("✓ Updated preferred codons table saved: ./results/section_8.1_GLM_preferred_codons_updated.csv\n\n")

# Show summary
cat("Summary of preference patterns:\n")
summary_table <- preferences_updated_glm |>
  dplyr::group_by(Preference_Pattern) |>
  dplyr::summarise(Count = n(), .groups = "drop")
print(summary_table)
cat("\n")

# 8.1.3) Compare GAM vs GLM approaches ----

cat("\n=== 8.1.3: Comparing GAM vs GLM Approaches ===\n\n")

comparison_df <- compare_gam_vs_glm(
  gam_results = all_gam_results,
  glm_results = all_glm_results,
  output_file = "./results/section_8.1_GAM_vs_GLM_comparison.pdf"
)

# Save comparison table
write.csv(comparison_df,
          "./results/section_8.1_GAM_vs_GLM_comparison_table.csv",
          row.names = FALSE)

cat("✓ Comparison table saved: ./results/section_8.1_GAM_vs_GLM_comparison_table.csv\n\n")

# 8.1.4) Visualize S-curves using GLM models (Gold Standard) ----

cat("\n=== 8.1.4: Creating S-Curve Visualizations (GLM Models) ===\n\n")
cat("Using GLM models as gold standard (more sensitive, detected 3 additional cases)\n")
cat("Visualizations reflect updated preference patterns:\n")
cat("  - ★ marks preferred codon(s)\n")
cat("  - Multiple ★ for co-optimal codons\n")
cat("  - No ★ for families without clear preference\n")
cat("  - Solid lines = significant, Dashed lines = non-significant\n\n")

# Create individual plots for interesting families
selected_families <- c("Ala", "Leu_4", "Ser_4", "Arg_4", "Gly", "Val", "Asp")

cat("Creating individual S-curve plots for selected families...\n")
for (fam in selected_families) {
  if (fam %in% names(all_glm_results)) {
    plot_family_scurves(
      model_result = all_glm_results[[fam]],
      meta_dt = integrated_data,
      output_file = sprintf("./results/section_8.1_scurve_GLM_%s.pdf", fam),
      alpha_significance = 0.05,
      preferred_codons_updated = preferences_updated_glm,
      n_points = 10000
    )
  }
}

cat("\n")

# Create multi-panel plot with all families
cat("Creating multi-panel plot with all families...\n")
plot_all_families_panel(
  all_model_results = all_glm_results,
  meta_dt = integrated_data,
  output_file = "./results/section_8.1_all_families_scurves_GLM.pdf",
  ncol = 4,
  preferred_codons_updated = preferences_updated_glm
)

cat("\n")

# 8.1.5) Statistical summary ----

cat("\n=== 8.1.5: Statistical Summary ===\n\n")

# How many codons show significant selection (using FDR-corrected p-values)?
gam_significant <- master_gam_table |>
  dplyr::filter(Codon != Baseline, 
                !is.na(p_adj),
                p_adj < 0.05)

glm_significant <- master_glm_table |>
  dplyr::filter(Codon != Baseline, 
                !is.na(p_adj),
                p_adj < 0.05)

cat(sprintf("GAM approach: %d / %d codons show significant selection (p < 0.05)\n",
            nrow(gam_significant), 
            nrow(master_gam_table |> dplyr::filter(Codon != Baseline))))

cat(sprintf("GLM approach: %d / %d codons show significant selection (p < 0.05)\n\n",
            nrow(glm_significant),
            nrow(master_glm_table |> dplyr::filter(Codon != Baseline))))

# Direction of selection (positive = increases with expression)
cat("Direction of selection:\n")
cat(sprintf("  GAM: %d positive, %d negative\n",
            sum(gam_significant$Selection_Slope > 0),
            sum(gam_significant$Selection_Slope < 0)))
cat(sprintf("  GLM: %d positive, %d negative\n\n",
            sum(glm_significant$Selection_Slope > 0),
            sum(glm_significant$Selection_Slope < 0)))

# Families with strongest selection signal
cat("Families with strongest selection signal (GAM approach):\n")
family_summary_gam <- master_gam_table |>
  dplyr::filter(Codon != Baseline) |>
  dplyr::group_by(Family) |>
  dplyr::summarise(
    Max_Slope = max(Selection_Slope),
    Min_p = min(p_value),
    N_Significant = sum(p_value < 0.05)
  ) |>
  dplyr::arrange(dplyr::desc(N_Significant), Min_p)

print(family_summary_gam, n = 10)
cat("\n")

# Compare with previous CAI-based preferred codons
cat("=== Comparing with CAI-Based Preferred Codons ===\n\n")

# Merge GAM preferred with CAI preferred
cai_vs_gam <- preferred_codons_corrected |>
  dplyr::select(Family = AA, CAI_Preferred = Codon) |>
  dplyr::left_join(
    gam_preferred_codons |> dplyr::select(Family, GAM_Preferred = Codon),
    by = "Family"
  ) |>
  dplyr::mutate(Agreement = (CAI_Preferred == GAM_Preferred))

cat(sprintf("Agreement between CAI and GAM preferred codons: %d / %d (%.1f%%)\n\n",
            sum(cai_vs_gam$Agreement, na.rm = TRUE),
            nrow(cai_vs_gam),
            100 * mean(cai_vs_gam$Agreement, na.rm = TRUE)))

# Show disagreements
disagreements <- cai_vs_gam |> dplyr::filter(!Agreement)
if (nrow(disagreements) > 0) {
  cat("Families where CAI and GAM disagree:\n")
  print(disagreements)
  cat("\n")
}

write.csv(cai_vs_gam,
          "./results/section_8.1_CAI_vs_GAM_comparison.csv",
          row.names = FALSE)

# Get preferred codons for export
preferred_codons <- preferred_codons_updated |>
  dplyr::filter(Preference_Pattern != "Neutral_Family") |>
  dplyr::select(Codon)

write.table(x = preferred_codons, file = 'results/preferred_codons.txt', 
            col.names = F, row.names = F, quote = F)

## *****************************************************************************
## 9) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________

# 9.1) CA analysis ---- 

codon_usage_m <- as.matrix(codon_usage[, -1])
rownames(codon_usage_m) <- codon_usage[[1]]
colnames(codon_usage_m) <- names(codon_usage)[-1]

codon_usage_CA <- CA(X = codon_usage_m, graph = F)
codon_usage_CA_coord <- as.data.frame(codon_usage_CA$row$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", row.names(codon_usage_CA$row$coord)))

codon_usage_CA_coord <- integrated_data |>
  left_join(y = codon_usage_CA_coord, by = "Gene_name")

# Rename dimensions to match plotting function expectations
names(codon_usage_CA_coord)[names(codon_usage_CA_coord) %in% c("Dim 1", "Dim 2", "Dim 3", "Dim 4", "Dim 5")] <- 
  c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5")

# Ensure Expression_Group is character (not factor) for color matching
codon_usage_CA_coord$Expression_Group <- as.character(codon_usage_CA_coord$Expression_Group)

# Create version with only extreme groups
codon_usage_CA_coord_extremes <- codon_usage_CA_coord |>
  filter(Expression_Group %in% c("Top 5%", "Bottom 5%")) |>  # Updated group names
  mutate(Expression_Group = as.character(Expression_Group))  # Ensure character type

# Define colors for all groups
# Note: Update these if you re-run Step 7 with new group names
colors_all <- c("Top 5%" = "#E41A1C",  # Current name in data
                "Bottom 5%" = "#377EB8",         # Current name in data
                "Middle 90%" = "#CCCCCC")         # Gray for middle (if exists)

colors_extremes <- c("Top 5%" = "#E41A1C", 
                     "Bottom 5%" = "#377EB8")

# Plot variance explained
plot_variance_explained(codon_usage_CA, 
                        analysis_type = "CA",
                        n_dims = 10,
                        output_file = "./results/CA_variance_explained.pdf")

cat("\n--- CA Analysis: Comparing All Three Groups ---\n")

# Bivariate ellipse plot - All groups (CA Dim 1 vs Dim 2)
plot_multivariate_analysis(codon_usage_CA_coord,
                           dims = c("Dim.1", "Dim.2"),
                           group_var = "Expression_Group",
                           analysis_type = "CA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_all,
                           output_file = "./results/CA_all_groups_D1_D2.pdf")

# Biplot with all groups
create_biplot(codon_usage_CA,
              coord_data = codon_usage_CA_coord,
              group_var = "Expression_Group",
              analysis_type = "CA",
              dims = c(1, 2),
              n_loadings = 20,
              colors = colors_all,
              output_file = "./results/CA_biplot_all_groups.pdf")

cat("\n--- CA Analysis: Top 5% vs Bottom 5% Only (Clearer Contrast) ---\n")

# Bivariate ellipse plot - Extremes only (CA Dim 1 vs Dim 2)
plot_multivariate_analysis(codon_usage_CA_coord_extremes,
                           dims = c("Dim.1", "Dim.2"),
                           group_var = "Expression_Group",
                           analysis_type = "CA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_extremes,
                           output_file = "./results/CA_extremes_only_D1_D2.pdf")

# Bivariate ellipse plot - Extremes (CA Dim 1 vs Dim 3)
plot_multivariate_analysis(codon_usage_CA_coord_extremes,
                           dims = c("Dim.1", "Dim.3"),
                           group_var = "Expression_Group",
                           analysis_type = "CA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_extremes,
                           output_file = "./results/CA_extremes_only_D1_D3.pdf")

# 3D static plot - Extremes only
plot_multivariate_analysis(codon_usage_CA_coord_extremes,
                           dims = c("Dim.1", "Dim.2", "Dim.3"),
                           group_var = "Expression_Group",
                           analysis_type = "CA",
                           plot_type = "3D_static",
                           colors = colors_extremes,
                           output_file = "./results/CA_extremes_3D.pdf")

# Biplot - Extremes only
create_biplot(codon_usage_CA,
              coord_data = codon_usage_CA_coord_extremes,
              group_var = "Expression_Group",
              analysis_type = "CA",
              dims = c(1, 3),
              n_loadings = 20,
              colors = colors_extremes,
              output_file = "./results/CA_biplot_extremes_only.pdf")

cat("\n--- Enhanced Biplots with Codon Classification ---\n")

# Source enhanced biplot functions
source("./src/enhanced_biplot.R")

# Prepare data for enhanced biplots
# Need to convert CA result to format expected by enhanced_biplot
ca_for_biplot <- list(
  li = codon_usage_CA$row$coord,  # Gene scores
  co = codon_usage_CA$col$coord,  # Codon loadings
  eig = codon_usage_CA$eig         # Eigenvalues
)
class(ca_for_biplot) <- "coa"

gene_data_ca <- codon_usage_CA_coord_extremes |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create enhanced biplots with different color schemes
cat("\n1. Creating CA biplot colored by selection status...\n")
p_ca_selection <- create_enhanced_biplot(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "selection",
  show_only_significant = FALSE,
  arrow_scale = 1.0,
  title = "CA Biplot: Codon Selection Status",
  output_file = "./results/CA_enhanced_biplot_selection.pdf"
)

cat("\n2. Creating CA biplot colored by preference (w=1)...\n")
p_ca_preference <- create_enhanced_biplot(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "preference",
  show_only_significant = FALSE,
  arrow_scale = 1.0,
  title = "CA Biplot: CAI Preferred Codons (w=1)",
  output_file = "./results/CA_enhanced_biplot_preference.pdf"
)

cat("\n3. Creating CA biplot colored by AT vs GC ending...\n")
p_ca_ending <- create_enhanced_biplot(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "ending",
  show_only_significant = FALSE,
  arrow_scale = 1.0,
  title = "CA Biplot: AT vs GC Codon Ending",
  output_file = "./results/CA_enhanced_biplot_ending.pdf"
)

cat("\n4. Creating CA biplot with combined classification...\n")
p_ca_combined <- create_enhanced_biplot(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "combined",
  show_only_significant = FALSE,
  arrow_scale = 1.0,
  title = "CA Biplot: Combined Codon Classification",
  output_file = "./results/CA_enhanced_biplot_combined.pdf"
)

cat("\n5. Creating CA biplot with significant codons only...\n")
p_ca_significant <- create_enhanced_biplot(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  preferred_codons = preferences_updated_glm,
  dims = c(1, 3),
  color_by = "combined",
  show_only_significant = TRUE,
  arrow_scale = 1.0,
  title = "CA Biplot: Significant Codons Only",
  output_file = "./results/CA_enhanced_biplot_significant_only.pdf"
)

# Analyze codon loading patterns
cat("\n6. Analyzing codon loading patterns...\n")
ca_loading_analysis <- analyze_codon_loading_patterns(
  ordination_result = ca_for_biplot,
  codon_test_results = codon_test_results,
  dims = c(1, 2),
  preferred_codons = preferred_codons_corrected
)

write.csv(ca_loading_analysis, 
          "./results/CA_codon_loading_analysis.csv",
          row.names = FALSE)

cat("\n✓ CA loading analysis saved: ./results/CA_codon_loading_analysis.csv\n")

# Statistical test for CA dimension separation
ca_manova <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ Expression_Group, 
                    data = codon_usage_CA_coord_extremes)

# Univariate tests for each dimension
cat("\n=== Univariate Tests for Each CA Dimension ===\n")
for (dim in c("Dim.1", "Dim.2", "Dim.3")) {
  wtest <- wilcox.test(as.formula(paste(dim, "~ Expression_Group")), 
                       data = codon_usage_CA_coord_extremes)
  cat(sprintf("%s: W = %.2f, p-value = %.4f %s\n", 
              dim, wtest$statistic, wtest$p.value,
              ifelse(wtest$p.value < 0.05, "***", "")))
}

# 8.2) PCA analysis ----

rscu_m <- as.matrix(cub_results$rscu_results[, -1])
rownames(rscu_m) <- cub_results$rscu_results[[1]]
colnames(rscu_m) <- names(cub_results$rscu_results)[-1]

rscu_PCA <- PCA(rscu_m, graph = F)

rscu_PCA_coord <- as.data.frame(rscu_PCA$ind$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", row.names(rscu_PCA$ind$coord)))

rscu_PCA_coord <- integrated_data |>
  left_join(y = rscu_PCA_coord, by = "Gene_name")

# Ensure Expression_Group is character (not factor) for color matching
rscu_PCA_coord$Expression_Group <- as.character(rscu_PCA_coord$Expression_Group)

# Create version with only extreme groups
rscu_PCA_coord_extremes <- rscu_PCA_coord |>
  filter(Expression_Group %in% c("Top 5%", "Bottom 5%")) |>  # Updated group names
  mutate(Expression_Group = as.character(Expression_Group))  # Ensure character type

# Plot variance explained
plot_variance_explained(rscu_PCA, 
                        analysis_type = "PCA",
                        n_dims = 10,
                        output_file = "./results/PCA_variance_explained.pdf")

cat("\n--- PCA Analysis: Comparing All Three Groups ---\n")

# Bivariate ellipse plot - All groups (PC1 vs PC2)
plot_multivariate_analysis(rscu_PCA_coord,
                           dims = c("Dim.1", "Dim.2"),
                           group_var = "Expression_Group",
                           analysis_type = "PCA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_all,
                           output_file = "./results/PCA_all_groups_PC1_PC2.pdf")

# Biplot with all groups
create_biplot(rscu_PCA,
              coord_data = rscu_PCA_coord,
              group_var = "Expression_Group",
              analysis_type = "PCA",
              dims = c(1, 2),
              n_loadings = 20,
              colors = colors_all,
              output_file = "./results/PCA_biplot_all_groups.pdf")

cat("\n--- PCA Analysis: Top 5% vs Bottom 5% Only (Clearer Contrast) ---\n")

# Bivariate ellipse plot - Extremes only (PC1 vs PC2)
plot_multivariate_analysis(rscu_PCA_coord_extremes,
                           dims = c("Dim.1", "Dim.2"),
                           group_var = "Expression_Group",
                           analysis_type = "PCA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_extremes,
                           output_file = "./results/PCA_extremes_only_PC1_PC2.pdf")

# Bivariate ellipse plot - Extremes (PC1 vs PC3)
plot_multivariate_analysis(rscu_PCA_coord_extremes,
                           dims = c("Dim.1", "Dim.3"),
                           group_var = "Expression_Group",
                           analysis_type = "PCA",
                           plot_type = "bivariate_ellipse",
                           confidence_level = 0.95,
                           colors = colors_extremes,
                           output_file = "./results/PCA_extremes_only_PC1_PC3.pdf")

# 3D static plot - Extremes only
plot_multivariate_analysis(rscu_PCA_coord_extremes,
                           dims = c("Dim.1", "Dim.2", "Dim.3"),
                           group_var = "Expression_Group",
                           analysis_type = "PCA",
                           plot_type = "3D_static",
                           colors = colors_extremes,
                           output_file = "./results/PCA_extremes_3D.pdf")

# Biplot - Extremes only
create_biplot(rscu_PCA,
              coord_data = rscu_PCA_coord_extremes,
              group_var = "Expression_Group",
              analysis_type = "PCA",
              dims = c(1, 2),
              n_loadings = 20,
              colors = colors_extremes,
              output_file = "./results/PCA_biplot_extremes_only.pdf")

cat("\n--- Enhanced PCA Biplots with Codon Classification ---\n")

# Prepare data for enhanced PCA biplots
gene_data_pca <- rscu_PCA_coord_extremes |>
  dplyr::select(Gene_name, expression_group = Expression_Group)

# Create enhanced PCA biplots with different color schemes
cat("\n1. Creating PCA biplot colored by selection status...\n")
p_pca_selection <- create_enhanced_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "selection",
  show_only_significant = FALSE,
  arrow_scale = 1.5,
  title = "PCA Biplot: Codon Selection Status",
  output_file = "./results/PCA_enhanced_biplot_selection.pdf"
)

cat("\n2. Creating PCA biplot colored by preference (w=1)...\n")
p_pca_preference <- create_enhanced_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "preference",
  show_only_significant = FALSE,
  arrow_scale = 1.5,
  title = "PCA Biplot: CAI Preferred Codons (w=1)",
  output_file = "./results/PCA_enhanced_biplot_preference.pdf"
)

cat("\n3. Creating PCA biplot colored by AT vs GC ending...\n")
p_pca_ending <- create_enhanced_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "ending",
  show_only_significant = FALSE,
  arrow_scale = 1.5,
  title = "PCA Biplot: AT vs GC Codon Ending",
  output_file = "./results/PCA_enhanced_biplot_ending.pdf"
)

cat("\n4. Creating PCA biplot with combined classification...\n")
p_pca_combined <- create_enhanced_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "combined",
  show_only_significant = FALSE,
  arrow_scale = 1.5,
  title = "PCA Biplot: Combined Codon Classification",
  output_file = "./results/PCA_enhanced_biplot_combined.pdf"
)

cat("\n5. Creating PCA biplot with significant codons only...\n")
p_pca_significant <- create_enhanced_biplot(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  color_by = "combined",
  show_only_significant = TRUE,
  arrow_scale = 1.5,
  title = "PCA Biplot: Significant Codons Only",
  output_file = "./results/PCA_enhanced_biplot_significant_only.pdf"
)

# Analyze PCA codon loading patterns
cat("\n6. Analyzing PCA codon loading patterns...\n")
pca_loading_analysis <- analyze_codon_loading_patterns(
  ordination_result = rscu_PCA,
  codon_test_results = codon_test_results,
  dims = c(1, 2),
  preferred_codons = preferred_codons_corrected
)

write.csv(pca_loading_analysis, 
          "./results/PCA_codon_loading_analysis.csv",
          row.names = FALSE)

cat("\n✓ PCA loading analysis saved: ./results/PCA_codon_loading_analysis.csv\n")

# Create combined panel for publication
cat("\n7. Creating combined biplot panels...\n")

# CA panel
create_biplot_panel(
  ordination_result = ca_for_biplot,
  gene_data = gene_data_ca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  output_file = "./results/CA_biplot_panel_4plots.pdf"
)

# PCA panel
create_biplot_panel(
  ordination_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  w_table = w_table,
  dims = c(1, 2),
  output_file = "./results/PCA_biplot_panel_4plots.pdf"
)

cat("\n✓ All enhanced biplots created successfully!\n\n")

# Statistical test for PCA dimension separation
cat("\n=== MANOVA Test: PCA Dimensions by Expression Group (Extremes) ===\n")
pca_manova <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ Expression_Group, 
                     data = rscu_PCA_coord_extremes)
print(summary(pca_manova))

# Univariate tests for each dimension
cat("\n=== Univariate Tests for Each PC Dimension ===\n")
for (dim in c("Dim.1", "Dim.2", "Dim.3")) {
  wtest <- wilcox.test(as.formula(paste(dim, "~ Expression_Group")), 
                       data = rscu_PCA_coord_extremes)
  cat(sprintf("%s: W = %.2f, p-value = %.4f %s\n", 
              dim, wtest$statistic, wtest$p.value,
              ifelse(wtest$p.value < 0.05, "***", "")))
}

# 8.3) 3D visuals for PCA results ----

cat("\n=== 8.3: Creating 3D PCA Visualizations ===\n")

source("./src/create_3d_pca_video.R")

# 8.3.1) Generate dynamics 3D videos for presentation ----

cat("\n8.3.1: Creating interactive 3D PCA plot...\n")

preferences_updated_glm <- preferences_updated_glm |>
  dplyr::mutate(Codon = Preferred_Codons) |>
  dplyr::mutate(relative_adaptiveness = 1)

# Create interactive 3D plot (HTML)
pca_3d_interactive <- create_3d_pca_plot(
  pca_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  preferred_codons = preferences_updated_glm,
  dims = c(1, 2, 3),
  color_by = "expression",
  show_loadings = TRUE,
  loading_scale = 5.0,
  title = "3D PCA: RSCU Analysis with Codon Loadings"
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
cat("8.3.2: Creating rotating 3D animation...\n")

pca_3d_animation <- create_3d_pca_animation(
  pca_result = rscu_PCA,
  gene_data = gene_data_pca,
  codon_test_results = codon_test_results,
  preferred_codons = preferred_codons_corrected,
  dims = c(1, 2, 3),
  color_by = "expression",
  show_loadings = TRUE,
  loading_scale = 5.0,
  title = "3D PCA Animation - RSCU Analysis",
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
    preferred_codons = preferred_codons_corrected,
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
## 10) Analyze codon loading patterns (AT vs GC bias) ----
## _____________________________________________________________________________

cat("\n=== Analyzing Codon Loading Patterns ===\n")
cat("Understanding mutational bias: AT-ending vs GC-ending codons\n\n")

source("./src/analyze_codon_loadings.R")

# Analyze CA loadings
ca_loadings <- analyze_codon_loadings(
  analysis_result = codon_usage_CA,
  analysis_type = "CA",
  dims = c(1, 2, 3),
  genetic_code = genetic_code_dna_long,
  output_file = "./results/CA_codon_loadings_AT_vs_GC.pdf"
)

# Analyze PCA loadings
pca_loadings <- analyze_codon_loadings(
  analysis_result = rscu_PCA,
  analysis_type = "PCA",
  dims = c(1, 2, 3),
  genetic_code = genetic_code_dna_long,
  output_file = "./results/PCA_codon_loadings_AT_vs_GC.pdf"
)

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

# Use CA loadings to define "preferred" codons
# Codons with high positive loadings on Dim 1 (which separates high vs low expression)
# Get codon loadings from CA
ca_codon_loadings <- as.data.frame(codon_usage_CA$col$coord)
ca_codon_loadings$Codon <- rownames(ca_codon_loadings)

# Get top codons on Dim 1 (one per amino acid family for non-circular test)
ca_codon_loadings$Amino_Acid <- genetic_code_dna_long[ca_codon_loadings$Codon]
ca_codon_loadings <- ca_codon_loadings[ca_codon_loadings$Amino_Acid != "STOP", ]

# Create full codon status table (all codons, not just top ones)
# Rank all codons by Dim 1 within each amino acid family
ca_codon_status <- ca_codon_loadings |>
  group_by(Amino_Acid) |>
  mutate(
    Rank = rank(-`Dim 1`, ties.method = "first"),
    Status = ifelse(Rank == 1, "Preferred", "Non-Preferred")
  ) |>
  ungroup() |>
  select(Codon, Amino_Acid, Status, `Dim 1`, Rank) |>
  as.data.frame()

# Convert codons to RNA format (T -> U) for compatibility with pairing function
ca_codon_status$Codon <- gsub("T", "U", ca_codon_status$Codon)

n_preferred <- sum(ca_codon_status$Status == "Preferred")
n_total <- nrow(ca_codon_status)

cat("Using CA Dimension 1 to classify all codons:\n")
cat("  Preferred:", n_preferred, "codons (highest Dim 1 in each AA family)\n")
cat("  Non-Preferred:", n_total - n_preferred, "codons\n")
cat("  Total:", n_total, "codons\n\n")

# Run the translational accuracy test
pairing_analysis <- classify_codon_anticodon_pairing(
  tRNA_data = tRNA_data,
  codon_supply = codon_supply,
  preferred_codons = ca_codon_status[, c("Codon", "Amino_Acid", "Status")],
  output_dir = "./results/tRNA_analysis_pairing",
  save_results = TRUE
)

cat("\n✓ Translational accuracy hypothesis test complete!\n")
cat("  Results saved to: ./results/tRNA_analysis_pairing/\n\n")

# ---- Parallel Analysis: Using Expression-Based Preferred Codons ----
cat("\n=== PARALLEL ANALYSIS: Expression-Based Preferred Codons ===\n")
cat("Testing if results change when using CAI-derived preferred codons (w = 1.0)\n")
cat("This provides an independent test using a different criterion\n\n")

# Run parallel analysis with expression-based preferred codons
preferences_updated_glm <- preferences_updated_glm |>
  dplyr::rename(Amino_Acid = Family)

pairing_analysis_expression <- classify_pairing_with_expression_preferred(
  tRNA_data = tRNA_data,
  codon_supply = codon_supply,
  preferred_codons_corrected = preferences_updated_glm,
  output_dir = "./results/tRNA_analysis_pairing_expression",
  save_results = TRUE
)

## 12) Polymorphism data integration ----

pi_data <- fread(input = "data/all_chromosomes.bygene.pi.txt")

# Homogenizing gene names to match the previous convention

pi_data <- pi_data |>
  dplyr::select(Chr, Gene, contains("Tajima"), contains("mean")) |>
  dplyr::mutate(Gene = paste0("MgIM767.", pi_data[['Gene']])) |>
  dplyr::rename(Gene_name = Gene)

# Analyzing the diversity of synonymous sites 

integrated_data <- integrated_data |>
  left_join(y = pi_data, by = "Gene_name")

ggplot(integrated_data, aes(x = Expression_Group,
                            y = Pi_mean_4fold)) +
  geom_boxplot() +
  theme_custom()

kruskal.test(Pi_mean_4fold ~ Expression_Group, data = integrated_data)

# Dunn test

dunn_result_cdc <- dunn.test::dunn.test(
  x = integrated_data$Pi_mean_4fold,
  g = integrated_data$Expression_Group,
  method = "bh",
  kw = TRUE,
  label = TRUE,
  wrap = FALSE,
  table = TRUE,
  list = FALSE,
  altp = TRUE
)

library(ggplot2)

# Ensure 'Expression_Group' is a factor with the correct order
integrated_data$Expression_Group <- factor(
  integrated_data$Expression_Group, 
  levels = c("Bottom 5%", "Middle 90%", "Top 5%")
)

detrending_pi_4fold <- residuals(gam(Pi_mean_4fold ~ s(GC3) + s(CDS_length_nt),
                               data = integrated_data,
                               na.action = na.exclude))

integrated_data$Pi_4fold_detrended <- detrending_pi_4fold

# Create the plot with Mean and 95% Confidence Intervals
p_pi_mean_ci <- ggplot(integrated_data, aes(x = Expression_Group, 
                                            y = Pi_4fold_detrended,
                                            fill = Expression_Group)) +
  
  # 2. Add the error bars for the 95% Confidence Interval
  stat_summary(
    fun.data = "mean_cl_normal",  # Calculates mean and 95% CI
    geom = "errorbar",
    width = 0.25,                # Width of the error bar caps
    color = "black",
    linewidth = 1
  ) +
  
  # 3. Add the point for the Mean
  stat_summary(
    fun = "mean",                # Calculates the mean
    geom = "point",
    size = 4,
    color = "black",
    shape = 23,                  # Diamond shape
    fill = "white"               # White fill for the diamond
  ) +
  
  # Use your group colors for the jitter points
  # scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
  #                              "Bottom 5%" = "#377EB8",
  #                              "Middle 90%" = "#999999")) +
  
  labs(
    title = "Mean Nucleotide Diversity (π) at 4-fold Sites",
    subtitle = "Showing mean and 95% confidence interval by expression group",
    y = "Pi (4-fold Synonymous Sites)",
    x = "Expression Group"
  ) +
  theme_custom() + # Use your custom theme
  theme(legend.position = "none") # Hide the fill legend

print(p_pi_mean_ci)

# Save the plot
ggsave("./results/pi_4fold_by_expression_mean_ci.pdf", p_pi_mean_ci, width = 8, height = 6)

# Exploring the relationship between diversity (4-fold) and ENC

plot(x = integrated_data$Pi_mean_4fold, 
           y = integrated_data$High_exp_log2)
lm(High_exp_log2 ~ Pi_mean_4fold, data = integrated_data)
plot(x = integrated_data$Pi_mean_all, 
     y = integrated_data$High_exp_log2)
lm(High_exp_log2 ~ Pi_mean_all, data = integrated_data)

# Does 4-fold differs from background?

t.test(integrated_data$Pi_mean_4fold, integrated_data$Pi_mean_all)

# Boxplot for expression groups and 4-fold pi

p_pi_4fold <- ggplot(integrated_data, aes(x = Expression_Group, y = Pi_mean_4fold)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightblue") +
  labs(title = "Nucleotide Diversity at 4-fold Synonymous Sites by Expression Group",
       x = "Expression Group",
       y = "Pi (4-fold Synonymous Sites)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

overall_pi <- ggplot(integrated_data, aes(x = Expression_Group, y = Pi_mean_all)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightblue") +
  labs(title = "Nucleotide Diversity by Expression Group",
       x = "Expression Group",
       y = "Pi (overall)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

int_variables <- integrated_data |>
  dplyr::select("CDC", "CDC_detrended", "Pi_mean_0fold", "Pi_mean_2fold", "Pi_mean_3fold",
                "Pi_mean_4fold", "Pi_mean_all", "p_adj", "ENC") |>
  as.matrix()

int_cor <- corrr::correlate(x = int_variables)

# 12.1) Searching for purifying selection signatures ----

# Run the GAM (to control for confounders)
# This asks: "Does expression predict Tajima's D at 4-fold sites?"
model_tajimaD <- gam(TajimaD_4fold ~ High_exp_log2 + s(CDS_length_nt) + s(GC12),
                     data = integrated_data)

summary(model_tajimaD)

# Plot it
ggplot(integrated_data, aes(x = High_exp_log2, y = TajimaD_all)) +
  geom_pointdensity() +
  geom_smooth(method = "gam", show.legend = T) +
  theme_custom() +
  labs(title = "Selection Signature vs. Expression",
       y = "Tajima's D (4-fold Synonymous Sites)",
       x = "log2(Expression)") +
  ylim(c(0.25, -0.25))

# 12.2) TajimaD vs CDC ----

model_cdc_tajima <- gam(TajimaD_Overall ~ CDC + s(CDS_length_nt) + s(GC3s),
                        data = integrated_data)

summary(model_cdc_tajima)

# 12.2) TajimaD vs expression

model_expression_tajima <- gam(TajimaD_4fold ~ High_exp_log2 + s(CDS_length_nt) + s(GC3s),
                        data = integrated_data)

summary(model_expression_tajima)

plot_tajimaD <- data.frame(detrended_TajimaD = residuals(gam(TajimaD_4fold ~ s(CDS_length_nt) + s(GC3s),
                                                              data = integrated_data,
                                                             na.action = 'na.exclude')))
plot_tajimaD$High_exp_log2 <- integrated_data$High_exp_log2

ggplot(plot_tajimaD, aes(x = High_exp_log2, y = detrended_TajimaD)) +
  geom_pointdensity() +
  geom_smooth(method = "lm") +
  theme_custom() +
  ylim(c(0.05, -0.05))

ggplot(integrated_data, aes(x = High_exp_log2, y = TajimaD_4fold)) +
  
  # Add points, maybe with a little less alpha
  geom_point(alpha = 0.2, size = 1) +
  
  # *** THIS IS THE KEY ***
  # Add a separate, linear trend line FOR EACH FACET
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  
  # *** THIS IS THE OTHER KEY ***
  # Split the plot into 5 panels, one for each GC3s bin
  facet_wrap(~ gc3s_label) +
  
  # Make it look clean
  theme_bw(base_size = 14) +
  labs(
    title = "Controlling for GC3s Reveals Negative Selection",
    subtitle = "Relationship between Expression and Tajima's D (π), binned by GC3s",
    x = "Gene Expression (log2)",
    y = "Tajima's D at 4-fold Synonymous Sites"
  )

# 12.3) pi vs expression ----

model_pi <- gam(Pi_mean_4fold ~ High_exp_log2 + s(CDS_length_nt) + s(GC3s),
                data = integrated_data)

summary(model_pi)

ggplot(integrated_data, aes(x = High_exp_log2, y = Pi_mean_4fold)) +
  geom_pointdensity() +
  geom_smooth(method = "lm") +
  theme_custom()

ggplot(integrated_data, aes(x = High_exp_log2, y = Pi_mean_4fold)) +
  
  # Add points, mapping GC3s to color
  geom_point(aes(color = GC3s), alpha = 0.5, size = 1.5) +
  
  # This adds the "wrong" marginal trend line
  geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") +
  
  # Use a nice color scale
  scale_color_viridis_c(name = "GC3s") +
  
  # Make it look clean
  theme_minimal(base_size = 14) +
  labs(
    title = "Nucleotide Diversity (π) vs. Gene Expression",
    subtitle = "Apparent positive trend (red) is a confound of GC3s (color)",
    x = "Gene Expression (log2)",
    y = "π at 4-fold Synonymous Sites",
    caption = "Each point is one gene."
  )

# We'll create 5 quantile groups (quintiles) for GC3s
integrated_data <- integrated_data |>
  mutate(gc3s_bin = ntile(GC3s, 5)) |>
  
  # Optional: Make the labels nicer for the plot
  mutate(gc3s_label = factor(gc3s_bin, 
                             labels = c("Lowest 20%", "20-40%", "40-60%", "60-80%", "Highest 80%")))

ggplot(integrated_data, aes(x = High_exp_log2, y = Pi_mean_4fold)) +
  
  # Add points, maybe with a little less alpha
  geom_point(alpha = 0.2, size = 1) +
  
  # *** THIS IS THE KEY ***
  # Add a separate, linear trend line FOR EACH FACET
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  
  # *** THIS IS THE OTHER KEY ***
  # Split the plot into 5 panels, one for each GC3s bin
  facet_wrap(~ gc3s_label) +
  
  # Make it look clean
  theme_bw(base_size = 14) +
  labs(
    title = "Controlling for GC3s Reveals Negative Selection",
    subtitle = "Relationship between Expression and Diversity (π), binned by GC3s",
    x = "Gene Expression (log2)",
    y = "π at 4-fold Synonymous Sites"
  )

# 12.4) Loading the site specific data ----

sfp_data <- fread("data/all_chromosomes.site_freq_by_preference.txt")

# 12.5) Loading the codon specific data ----

csf_data <- fread("data/all_chromosomes.codon_frequencies_preferred.txt")
csf_data <- csf_data |>
  dplyr::mutate(Gene_name = paste0("MgIM767.", Gene))

# 12.5.1) F_preferred vs Expression ----

all_codon_data <- integrated_data |>
  left_join(csf_data, by = "Gene_name") |>
  na.omit()

# Standardizing the gene position from 0 to 1 per gene

all_codon_data <- all_codon_data |>
  dplyr::group_by(Gene_name) |>
  dplyr::mutate(
    Codon_Pos_Rel = Codon_Pos / max(Codon_Pos)
    ) |>
  dplyr::ungroup()

n <- nrow(all_codon_data)
all_codon_data$Preferred_Freq_beta <- (all_codon_data$Preferred_Freq * (n - 1) + 0.5) / n

model_selection <- gam(
  Preferred_Freq_beta ~ High_exp_log2 + s(CDS_length_nt) + s(GC3s),
  family = betar(link = "logit"),
  data = all_codon_data
)

summary(model_selection)

# 12.5.2) 5' towards 3' idea ----

model_position <- gam(
  Preferred_Freq_beta ~ s(Codon_Pos, k = 20),
  family = betar(link = "logit"),
  data = all_codon_data
)

summary(model_position)

p_pos_trend <- ggplot(all_codon_data, aes(x = Codon_Pos_Rel, y = Preferred_Freq)) +
  # geom_smooth() will draw the average trend line
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  geom_pointdensity() +
  labs(
    title = "Frequence of preferred codons vs position",
    y = "Preferred freq",
    x = "Codon position (relative)"
  ) +
  theme_bw()

ggsave(filename = "results/FPvsPosition.pdf", plot = p_pos_trend)

## *****************************************************************************
## 13) Selection Coefficient Analysis (Mutation-Selection-Drift Balance) ----
## _____________________________________________________________________________

cat("\n=== SELECTION COEFFICIENT ANALYSIS ===\n")
cat("Estimating population-scaled selection (S = 4Nes) using Hershberg & Petrov model\n\n")

# Load selection coefficient functions
source("./src/selection_coefficient_analysis.R")

# Get preferred codons from CAI analysis (w = 1.0)
preferred_codons <- preferences_updated_glm |>
  pull(Codon)

cat(sprintf("Using %d optimal codons as 'preferred' codons:\n", length(preferred_codons)))
cat(paste(preferred_codons, collapse = ", "), "\n")

# Prepare expression data (using High_exp from bud tissue as primary metric)
expression_df <- integrated_data |>
  dplyr::select(Gene_name, Expression = High_exp_log2)

# Calculate selection coefficients for all genes
selection_results <- run_selection_analysis(
  codon_usage = codon_usage,
  expression_data = expression_df,
  preferred_codons = preferred_codons,
  genetic_code = genetic_code_dna_long,
  low_expr_quantile = 0.10,  # Use bottom 10% to estimate mutation bias
  run_diagnostics = TRUE
)

selection_results_long <- selection_results |>
  pivot_longer(cols = c("P_preferred", "S"), names_to = "Metric",
               values_to = "Values")

ggplot(data = selection_results_long, mapping = aes(x = Expression,
                                                    y = Values,
                                                    color = Metric,
                                                    fill = Metric)) +
  geom_smooth() +
  theme_custom()

# Sensitivity analysis: How does s vary with different Ne values?
# Literature estimates for M. guttatus Ne: ~200,000 - 500,000
Ne_sensitivity <- analyze_Ne_sensitivity(
  selection_results,
  Ne_values = c(1e5, 2e5, 3e5, 5e5, 1e6)
)

# Create visualizations
cat("\nGenerating plots...\n")
plot_S_vs_expression(selection_results, "./results/S_vs_expression.pdf")
plot_S_distribution(selection_results, "./results/S_distribution.pdf")

# Compare S between expression groups
selection_with_groups <- selection_results |>
  left_join(integrated_data |> 
              dplyr::select(Gene_name, Expression_Group), by = "Gene_name")

cat("\n=== S by Expression Group ===\n")
S_by_group <- selection_with_groups |>
  group_by(Expression_Group) |>
  summarize(
    n = n(),
    mean_S = mean(S, na.rm = TRUE),
    median_S = median(S, na.rm = TRUE),
    sd_S = sd(S, na.rm = TRUE)
  )
print(S_by_group)

# Statistical test
kw_S <- kruskal.test(S ~ Expression_Group, data = selection_with_groups)
cat("\nKruskal-Wallis test for S across expression groups:\n")
print(kw_S)

# Key biological interpretation
cat("\n=== BIOLOGICAL INTERPRETATION ===\n")
M <- attr(selection_results, "mutation_bias")
cat(sprintf("Mutation bias (M = μ_p/μ_u): %.4f\n", M))
if (M > 1) {
  cat("  → Mutation pressure OPPOSES preferred codons (selection maintains them)\n")
} else {
  cat("  → Mutation pressure FAVORS preferred codons (selection reinforced by mutation)\n")
}

median_S <- median(selection_results$S, na.rm = TRUE)
cat(sprintf("\nMedian S = 4Nes: %.4f\n", median_S))

# For Ne = 300,000 (midpoint estimate)
Ne_midpoint <- 3e5
s_midpoint <- calculate_s_from_S(median_S, Ne_midpoint)
cat(sprintf("\nAssuming Ne = %s:\n", format(Ne_midpoint, scientific = FALSE)))
cat(sprintf("  Median s = %.2e\n", s_midpoint))
cat(sprintf("  s·Ne = %.2f\n", s_midpoint * Ne_midpoint))

if (abs(s_midpoint * Ne_midpoint) < 2) {
  cat("\n✓ This confirms WEAK SELECTION regime (s·Ne ~ 1)\n")
  cat("  → Selection and drift are of comparable magnitude\n")
  cat("  → Consistent with intermediate codon bias patterns\n")
} else {
  cat("\n  Selection is relatively strong (s·Ne >> 1)\n")
}

# Save results
write.table(selection_results, 
            "./results/selection_coefficients.csv",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(Ne_sensitivity,
            "./results/Ne_sensitivity_analysis.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✓ Selection coefficient analysis complete!\n")
cat("  Results saved to ./results/selection_coefficients.csv\n")

## *****************************************************************************
## Estimate mutation rates ----
## _____________________________________________________________________________

# Extract intronic sequences
introns_list <- get_intron_sequences(fasta_file = "./data/Mguttatusvar_IM767_887_v2.0.hardmasked.fa",
                                     ann_file = "./data/Mguttatusvar_IM767_887_v2.1.gene.gff3",
                                     organism = "Mimulus guttatus")

# Calculate nucleotide composicion per window
nuc_composition <- get_base_composition_per_windows(genome_seqinfo = introns_list$genome_seqinfo, 
                                                    trimmed_introns = introns_list$trimmed_introns,
                                                    intron_seqs = introns_list$intron_seqs,
                                                    window_size = 100000)

windows_thinned <- refine_windows_for_genes(nuc_composition, 1000)

nuc_composition_filtered <- nuc_composition |>
  filter(total_bp >= 1000) |>
  dplyr::mutate(mid_point = (start + end ) / 2)

# Calculate Q matrix

Q_matrices <- apply_q_matrix_to_windows(nuc_composition_filtered)

# Plot frequency sprectrum across windows

freq_pi_nuc_long <- nuc_composition_filtered |>
  pivot_longer(cols = contains("pi"),
               names_to = "Pi_nuc",
               values_to = "Freq")

ggplot(data = freq_pi_nuc_long, 
       mapping = aes(x = mid_point, 
                     y = Freq, 
                     color = Pi_nuc)) +
  
  # 1. Light lines for raw data (shows the noise/variance)
  geom_line(alpha = 0.2, linewidth = 0.3) + 
  
  # 2. Smooth trend lines (shows the mutational pressure signal)
  geom_smooth(se = FALSE, span = 0.2, linewidth = 1) +
  
  # 3. Facet by Chromosome to separate genomic contexts
  # scales = "free_x" ensures chromosomes with different lengths fit well
  # space = "free_x" keeps the physical scale consistent across panels
  facet_grid(. ~ seqnames, scales = "free_x", space = "free_x") +
  
  # 4. Formatting
  scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6), 
                     breaks = pretty_breaks(n = 3)) +
  labs(
    x = "Genomic Position (Mb)",
    y = "Nucleotide Frequency (Pi)",
    title = "Mutational Spectrum Across the Genome"
  ) +
  theme_custom() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate labels if crowded
    panel.spacing = unit(0.1, "lines") # Tighten the gap between chromosomes
  )

ggsave(filename = "results/pi_spectrum_across_windows.pdf", width = 16,
       height = 8)

# Extract the list of matrices from the data frame
q_list <- Q_matrices$Q_matrix
names(q_list) <- Q_matrices$window_idx

# Use abind to stack the matrices along the third dimension (windows)
# The dimensions will be [from base, to base, window index]
Q_array <- abind::abind(q_list, along = 3)

# Assign dimnames for clarity (optional, but good practice)
dimnames(Q_array) <- list(
  From = c("A", "C", "G", "T"), 
  To = c("A", "C", "G", "T"), 
  Window = Q_matrices$window_idx
)

message(paste("Created Q_array with dimensions:", 
              paste(dim(Q_array), collapse = " x ")))

# Searching for evidence of variation in M

df_analysis <- Q_matrices |>
  # Add the extracted rate as a new column
  dplyr::mutate(Q_AG_rate = Q_array["A", "G", ]) |>
  dplyr::mutate(Q_AT_rate = Q_array["A", "T", ]) |>
  # Ensure chromosome names are a factor for ANOVA
  dplyr::mutate(seqnames = as.factor(seqnames))

# Test if the mean A->G transition rate differs significantly by chromosome
rate_anova <- aov(Q_AT_rate ~ seqnames, data = df_analysis)
summary(rate_anova)

df_analysis <- df_analysis |>
  dplyr::mutate(GC_content = pi_G + pi_C)

GC_content <- aov(GC_content ~ seqnames, data = df_analysis)
summary(GC_content)

df_plot <- df_analysis |>
  dplyr::mutate(GC_content = pi_G + pi_C) |>
  
  # B. Define the X-axis position (midpoint of the window)
  dplyr::mutate(midpoint = (start + end) / 2) |>
  
  # C. Select only the columns needed for the plot
  #    (Chromosome, Position, and the two metrics)
  dplyr::select(seqnames, midpoint, GC_content, Q_AG_rate) |>
  
  # D. Reshape to "Long" format for ggplot
  #    This stacks 'GC_content' and 'Q_AG_rate' into a single column
  pivot_longer(
    cols = c(GC_content, Q_AG_rate),
    names_to = "Variable",
    values_to = "Rate_Value"
  )

plot_genomic_rate_variation(df_plot)

ggsave(filename = "results/genomic_rate_variation.pdf", width = 16,
       height = 8)

# These analysis suggest that the mutational pressure is not shared across all
# genes. Let's cluster windows as a function of these matrices and exmploy Gaussian
# mixed models to get putative clusters.

# Getting rates out for each nucleotide and normalized directional mutation spectrum

window_data <- data.frame(window_idx = Q_matrices$window_idx)
out_rates <- base::do.call("rbind", base::lapply(X = q_list, FUN = function(x)
                    {
                      diag(x)
                    }))
window_data <- window_data |> cbind(as.data.frame(out_rates))

# Getting six representative instantaneous transition rates

trans_rates <- base::do.call("rbind", base::lapply(X = q_list, FUN = function(x)
{
  c(
    "A>C" = x["A", "C"], # A --> C
    "A>G" = x["A", "G"], # A --> G
    "A>T" = x["A", "T"], # A --> T
    "C>G" = x["C", "G"], # C --> G
    "C>T" = x["C", "T"], # C --> T
    "G>T" = x["G", "T"] # G --> T
  )
}))

window_data <- window_data |> cbind(trans_rates)

# Getting summary variables

widow_data <- window_data |>
  dplyr::select(window_idx) |>
  cbind(prcomp(x = as.matrix(window_data[, -1]), center = T,
               scale = T)$x[, paste0("PC", 1:4)])

# Getting the clusters

clusters_localM <- make_clusters(data = window_data[, -1], G = 1:10)

# NOTE: GMM does not find evidence for multiple clusters, and the model with greater
# BIC was EEE with one component.

# Calcultate dM

# 1. Calculate Global Average Nucleotide Frequencies (weighted by total_bp)
global_stats <- nuc_composition_filtered |>
  summarize(
    total_genome_bp = sum(total_bp),
    avg_pi_A = sum(pi_A * total_bp) / total_genome_bp,
    avg_pi_C = sum(pi_C * total_bp) / total_genome_bp,
    avg_pi_G = sum(pi_G * total_bp) / total_genome_bp,
    avg_pi_T = sum(pi_T * total_bp) / total_genome_bp
  )

# 2. Generate the file
dM_data <- generate_anacoda_dM(
  pi_A = global_stats$avg_pi_A,
  pi_C = global_stats$avg_pi_C,
  pi_G = global_stats$avg_pi_G,
  pi_T = global_stats$avg_pi_T,
  output_file = "./data/Mguttatus_intron_derived_dM.csv"
)

## *****************************************************************************
## 14) AnaCoDa-based analysis ----
## _____________________________________________________________________________

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
# --max_num_runs 3 \
# --fix_dM \
# --dM ./data/Mguttatus_intron_derived_dM.csv
# 
# echo "Job finished on $(date)"

# 14.1) Retrieving AnaCoDa results to analyze congruence between runs ----

# 14.1.1) Naive model ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./MCMC_results/results_naive_2/run_1",
  "./MCMC_results/results_naive_2/run_2",
  "./MCMC_results/results_naive_2/run_3"
)

Naive_conv <- GR_convergence(run_dirs)

# 14.1.2) dM-fixed model ----

# Setup paths for the 3 runs
run_dirs <- c(
  "./MCMC_results/results_dM_fixed/run_1",
  "./MCMC_results/results_dM_fixed/run_2",
  "./MCMC_results/results_dM_fixed/run_3"
)

dM_fixed_conv <- GR_convergence(run_dirs, parameter = 'selection') # Mutation is fixed

# Convergence was achieved for dM_fixed

# 14.2) Checking the correlation between estimates of phi and the expression data ----

# From now on, we will work with chain 1, given that all chains are statistically
# equal.

phi_hat <- read.csv(file = "results/MCMC_results/results_dM_fixed/run_1/Parameter_est/gene_expression.txt") |>
  dplyr::select(GeneID, Mean, Mean.log10) |>
  dplyr::rename(MeanPhi = Mean, Mean.log10.Phi = Mean.log10)

phi <- exp_complete |>
  left_join(phi_hat, by = join_by("Gene" == "GeneID")) |>
  dplyr::mutate(High_exp_log10 = log10(High_exp + 1)) |>
  na.exclude()

cor.test(phi$Mean.log10.Phi, phi$High_exp_log10)

# There is no good correspondence with empirical data
# Next step is to pass expression data to the AnaCoDa

# 14.3) Preparing the expression data ----

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