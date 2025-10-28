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
                        'factoextra', 'dplyr')

set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = T, n_cores = 10)

# 1.1) Definition of globals ----
# Look-up table

genetic_code_dna_long <- c(
  "TTT"="Phe", "TTC"="Phe", "TTA"="Leu", "TTG"="Leu",
  "TCT"="Ser", "TCC"="Ser", "TCA"="Ser", "TCG"="Ser",
  "TAT"="Tyr", "TAC"="Tyr", "TAA"="STOP", "TAG"="STOP",
  "TGT"="Cys", "TGC"="Cys", "TGA"="STOP", "TGG"="Trp",
  "CTT"="Leu", "CTC"="Leu", "CTA"="Leu", "CTG"="Leu",
  "CCT"="Pro", "CCC"="Pro", "CCA"="Pro", "CCG"="Pro",
  "CAT"="His", "CAC"="His", "CAA"="Gln", "CAG"="Gln",
  "CGT"="Arg", "CGC"="Arg", "CGA"="Arg", "CGG"="Arg",
  "ATT"="Ile", "ATC"="Ile", "ATA"="Ile", "ATG"="Met",
  "ACT"="Thr", "ACC"="Thr", "ACA"="Thr", "ACG"="Thr",
  "AAT"="Asn", "AAC"="Asn", "AAA"="Lys", "AAG"="Lys",
  "AGT"="Ser", "AGC"="Ser", "AGA"="Arg", "AGG"="Arg",
  "GTT"="Val", "GTC"="Val", "GTA"="Val", "GTG"="Val",
  "GCT"="Ala", "GCC"="Ala", "GCA"="Ala", "GCG"="Ala",
  "GAT"="Asp", "GAC"="Asp", "GAA"="Glu", "GAG"="Glu",
  "GGT"="Gly", "GGC"="Gly", "GGA"="Gly", "GGG"="Gly"
)

## *****************************************************************************
## 3) Load the data ----
## _____________________________________________________________________________

## 3.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
                                      format = 'fasta')

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## 3.2) Clean the codon usage object (remove the STOP codon)

codon_usage <- codon_usage |>
  trim_uninformative(genetic_code = genetic_code_dna_long)

## *****************************************************************************
## 4) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results")

# Create amino acid-specific codon logos
create_aa_specific_logos(codon_usage, genetic_code_dna_long,
                         output_dir = "./results/codon_logos")

## *****************************************************************************
## 5) tRNA abundance correlation analysis ----
## _____________________________________________________________________________

message("Analyzing correlation between codon usage and tRNA abundance...")

# Perform tRNA-codon correlation analysis using filtered tRNA data and genes
# with a ENC < 0.35

tRNA_correlation_results <- tRNA_codon_correlation(
  codon_counts = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis",
  test_method = "spearman"  # Can also use "pearson" or "kendall"
)

message("tRNA correlation analysis complete!")

## *****************************************************************************
## 6) Additional analyses (optional) ----
## _____________________________________________________________________________

## Individual metric calculations (if needed separately):

# Calculate RSCU
rscu_values <- calculate_rscu(codon_usage, genetic_code_dna_long)

# Calculate ENC
enc_values <- calculate_enc(codon_usage, genetic_code_dna_long)

# Calculate RF
rf_values <- calculate_rf(codon_usage, genetic_code_dna_long)

# Get he PSPM
pspm_overall <- calculate_overall_PSPM(rf_values, genetic_code_dna_long)

# Create logos
create_aa_logo(pspm_overall)

# Calculate GC content
gc_content <- calculate_gc_content(codon_usage)

# Create specific visualizations
visualize_codon_usage(codon_usage, genetic_code_dna_long,
                     "custom_heatmap.pdf", type = "heatmap")
neutrality_plot(gc_content, "custom_neutrality.pdf")
enc_plot(enc_values, gc_content, "custom_enc.pdf")
pr2_bias_plot(codon_usage, "custom_pr2.pdf")

message("\nAnalysis complete! Check the './results' directory for all outputs.")

## *****************************************************************************
## 7) Modeling relationship between ENC and Expression profiles ----
## _____________________________________________________________________________

# Trimming suffix from ENC table in gene names
enc_values[, Gene_name := sub("\\.1$", "", Gene_name)]

exp_data_bud <- read.table(file = "./data/bud_gene_expression_cpm_remapped.txt",
                       header = T)

# Combining CUB metric with expression profiles for buds

exp_bud_enc <- dplyr::left_join(exp_data_bud, enc_values, 
                                by = dplyr::join_by(Remapped_Gene == Gene_name)) |>
  dplyr::select(Remapped_Gene, Expression, ENC) |>
  dplyr::rename(Gene_name = Remapped_Gene) |>
  na.omit()

# Add gene length (CDS length in codons and nucleotides)
cat("\n=== Adding Gene Length Information ===\n")

# First, clean gene names in codon_usage to match exp_bud_enc
# Calculate total codons from numeric columns only
codon_columns <- names(codon_usage)[names(codon_usage) != "Gene_name"]

gene_lengths <- codon_usage |>
  mutate(
    Gene_name_clean = sub("\\.1$", "", Gene_name),  # Remove .1 suffix
    Total_Codons = rowSums(across(all_of(codon_columns)), na.rm = TRUE),
    CDS_length_nt = Total_Codons * 3,  # nucleotides
    CDS_length_aa = Total_Codons        # amino acids (codons)
  ) |>
  select(Gene_name_clean, Total_Codons, CDS_length_nt, CDS_length_aa) |>
  rename(Gene_name = Gene_name_clean)

exp_bud_enc <- exp_bud_enc |>
  left_join(gene_lengths, by = "Gene_name")

cat(sprintf("Added length for %d genes\n", sum(!is.na(exp_bud_enc$CDS_length_nt))))
cat(sprintf("Mean CDS length: %.0f nt (%.0f codons)\n", 
            mean(exp_bud_enc$CDS_length_nt, na.rm = TRUE),
            mean(exp_bud_enc$Total_Codons, na.rm = TRUE)))

cor(exp_bud_enc$ENC, exp_bud_enc$Expression)
plot(exp_bud_enc$ENC, exp_bud_enc$Expression)
lm(Expression ~ ENC, data = exp_bud_enc)

# Define expression groups: Top 5% vs Bottom 5% (extreme comparison)

top_5_cutoff <- quantile(exp_bud_enc$Expression, probs = 0.95)
bottom_5_cutoff <- quantile(exp_bud_enc$Expression, probs = 0.05)

exp_bud_enc$Expression_Group <- case_when(
  exp_bud_enc$Expression >= top_5_cutoff ~ "Top 5%",
  exp_bud_enc$Expression <= bottom_5_cutoff ~ "Bottom 5%",
  TRUE ~ "Middle 90%"
)

# Also create a version with only extremes for clearer visualization
exp_bud_enc_extremes <- exp_bud_enc |>
  filter(Expression_Group %in% c("Top 5%", "Bottom 5%"))

# Print the actual group names for verification
cat("\n=== WARNING: Check Expression Group Names ===\n")
cat("If you see different group names in Step 8 plots, re-run this Step 7!\n")
cat("Expected names: Top 5%, Bottom 5%, Middle 90%\n")
cat("Actual names in data:\n")
print(table(exp_bud_enc$Expression_Group))
cat("\n")

# Statistical comparison
cat("\n=== Expression Group Statistics ===\n")
cat(sprintf("Top 5%% threshold: %.2f CPM\n", top_5_cutoff))
cat(sprintf("Bottom 5%% threshold: %.2f CPM\n", bottom_5_cutoff))
cat(sprintf("Top 5%% genes: %d\n", sum(exp_bud_enc$Expression_Group == "Top 5%")))
cat(sprintf("Bottom 5%% genes: %d\n", sum(exp_bud_enc$Expression_Group == "Bottom 5%")))
cat(sprintf("Middle 90%% genes: %d\n", sum(exp_bud_enc$Expression_Group == "Middle 90%")))

# Boxplot comparison
library(ggplot2)
p_boxplot <- ggplot(exp_bud_enc, aes(x = Expression_Group, y = ENC, fill = Expression_Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "ENC by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR",
       y = "Effective Number of Codons (ENC)",
       x = "Expression Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("./results/ENC_by_expression_group.pdf", p_boxplot, width = 8, height = 6)

# Statistical tests for three groups
cat("\n=== Kruskal-Wallis Test: ENC across All Three Groups ===\n")
cat("H0: All three groups have the same median ENC\n")
kw_test_enc <- kruskal.test(ENC ~ Expression_Group, data = exp_bud_enc)
print(kw_test_enc)

if (kw_test_enc$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with Bonferroni Correction ===\n")
  
  # Install and load dunn.test if not available
  if (!require("dunn.test", quietly = TRUE)) {
    cat("Installing dunn.test package...\n")
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
    library(dunn.test)
  }
  
  # Perform Dunn's test with Bonferroni correction
  dunn_result_enc <- dunn.test::dunn.test(
    x = exp_bud_enc$ENC,
    g = exp_bud_enc$Expression_Group,
    method = "bonferroni",
    kw = TRUE,
    label = TRUE,
    wrap = FALSE,
    table = TRUE,
    list = FALSE,
    altp = TRUE
  )
  
  cat("\nInterpretation of pairwise comparisons:\n")
  cat("  - Adjusted p-values account for multiple testing (Bonferroni)\n")
  cat("  - p < 0.05 indicates significant difference between groups\n")
  
} else {
  cat("\nNo significant difference among groups (p >= 0.05)\n")
  cat("Post-hoc tests not necessary.\n")
}

cat("\n=== Summary Statistics ===\n")
summary_stats <- exp_bud_enc %>%
  group_by(Expression_Group) %>%
  summarise(
    n = n(),
    mean_ENC = mean(ENC, na.rm = TRUE),
    median_ENC = median(ENC, na.rm = TRUE),
    sd_ENC = sd(ENC, na.rm = TRUE),
    mean_Expression = mean(Expression, na.rm = TRUE)
  )
print(summary_stats)

# Effect sizes for pairwise comparisons
cat("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

# Helper function to calculate Cohen's d
cohens_d_calc <- function(x1, x2) {
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- var(x1, na.rm = TRUE)
  s2 <- var(x2, na.rm = TRUE)
  pooled_sd <- sqrt((s1 + s2) / 2)
  d <- (m1 - m2) / pooled_sd
  return(d)
}

# Get ENC values for each group
top5_enc <- exp_bud_enc |> filter(Expression_Group == "Top 5%") |> pull(ENC)
middle_enc <- exp_bud_enc |> filter(Expression_Group == "Middle 90%") |> pull(ENC)
bottom5_enc <- exp_bud_enc |> filter(Expression_Group == "Bottom 5%") |> pull(ENC)

# Calculate effect sizes
if (length(top5_enc) > 0 && length(middle_enc) > 0) {
  d_top_middle <- cohens_d_calc(top5_enc, middle_enc)
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle))
}

if (length(top5_enc) > 0 && length(bottom5_enc) > 0) {
  d_top_bottom <- cohens_d_calc(top5_enc, bottom5_enc)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom))
}

if (length(middle_enc) > 0 && length(bottom5_enc) > 0) {
  d_middle_bottom <- cohens_d_calc(middle_enc, bottom5_enc)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom))
}

cat("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

# Gene length analysis
cat("\n=== Gene Length by Expression Group ===\n")
cat("Checking if gene length explains ENC patterns\n\n")
length_stats <- exp_bud_enc %>%
  group_by(Expression_Group) %>%
  summarise(
    n = n(),
    mean_length_aa = mean(CDS_length_aa, na.rm = TRUE),
    median_length_aa = median(CDS_length_aa, na.rm = TRUE),
    sd_length_aa = sd(CDS_length_aa, na.rm = TRUE),
    mean_length_nt = mean(CDS_length_nt, na.rm = TRUE)
  )
print(length_stats)

# Test for length differences
cat("\n=== Kruskal-Wallis Test: Gene Length across Groups ===\n")
kw_length <- kruskal.test(CDS_length_aa ~ Expression_Group, data = exp_bud_enc)
print(kw_length)

# Correlation between length and ENC
cat("\n=== Correlation: Gene Length vs ENC ===\n")
cor_length_enc <- cor(exp_bud_enc$CDS_length_aa, exp_bud_enc$ENC, use = "complete.obs")
cat(sprintf("Pearson r = %.4f\n", cor_length_enc))
cat("Note: Negative correlation = shorter genes have lower ENC (artifact)\n")

## *****************************************************************************
## 8) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________

# Source plotting functions
source("./src/plot_multivariate_analysis.R")

cat("\n=== Running Multivariate Analyses ===\n")
cat("This will create plots for:\n")
cat("  1. All three groups (Top 5%, Middle 90%, Bottom 5%)\n")
cat("  2. Only extremes (Top 5% vs Bottom 5% - clearer contrast)\n\n")

# 8.1) CA analysis ---- 

codon_usage_m <- as.matrix(codon_usage[, -1])
rownames(codon_usage_m) <- codon_usage[[1]]
colnames(codon_usage_m) <- names(codon_usage)[-1]

codon_usage_CA <- CA(X = codon_usage_m, graph = F)
codon_usage_CA_coord <- as.data.frame(codon_usage_CA$row$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", row.names(codon_usage_CA$row$coord)))

codon_usage_CA_coord <- exp_bud_enc |>
  left_join(y = codon_usage_CA_coord, by = "Gene_name")

# Rename dimensions to match plotting function expectations
names(codon_usage_CA_coord)[names(codon_usage_CA_coord) %in% c("Dim 1", "Dim 2", "Dim 3", "Dim 4", "Dim 5")] <- 
  c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5")

# Ensure Expression_Group is character (not factor) for color matching
codon_usage_CA_coord$Expression_Group <- as.character(codon_usage_CA_coord$Expression_Group)

# Debug: check that Expression_Group exists and has correct values
cat("Checking Expression_Group column:\n")
cat("  Column exists:", "Expression_Group" %in% names(codon_usage_CA_coord), "\n")
cat("  Unique values:", paste(unique(codon_usage_CA_coord$Expression_Group), collapse = ", "), "\n")
cat("  Counts:", paste(table(codon_usage_CA_coord$Expression_Group), collapse = ", "), "\n")

# Create version with only extreme groups
codon_usage_CA_coord_extremes <- codon_usage_CA_coord |>
  filter(Expression_Group %in% c("Top 5% Expressed", "Bottom 95%")) |>  # Updated group names
  mutate(Expression_Group = as.character(Expression_Group))  # Ensure character type

# Define colors for all groups
# Note: Update these if you re-run Step 7 with new group names
colors_all <- c("Top 5% Expressed" = "#E41A1C",  # Current name in data
                "Bottom 95%" = "#377EB8",         # Current name in data
                "Middle 90%" = "#CCCCCC")         # Gray for middle (if exists)

colors_extremes <- c("Top 5% Expressed" = "#E41A1C", 
                     "Bottom 95%" = "#377EB8")

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
             dims = c(1, 2),
             n_loadings = 20,
             colors = colors_extremes,
             output_file = "./results/CA_biplot_extremes_only.pdf")

# Statistical test for CA dimension separation
cat("\n=== MANOVA Test: CA Dimensions by Expression Group (Extremes) ===\n")
ca_manova <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ Expression_Group, 
                    data = codon_usage_CA_coord_extremes)
print(summary(ca_manova))

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

rscu_m <- as.matrix(rscu_values[, -1])
rownames(rscu_m) <- rscu_values[[1]]
colnames(rscu_m) <- names(rscu_values)[-1]

rscu_PCA <- PCA(rscu_m, graph = F)

rscu_PCA_coord <- as.data.frame(rscu_PCA$ind$coord) |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", row.names(rscu_PCA$ind$coord)))

rscu_PCA_coord <- exp_bud_enc |>
  left_join(y = rscu_PCA_coord, by = "Gene_name")

# Ensure Expression_Group is character (not factor) for color matching
rscu_PCA_coord$Expression_Group <- as.character(rscu_PCA_coord$Expression_Group)

# Create version with only extreme groups
rscu_PCA_coord_extremes <- rscu_PCA_coord |>
  filter(Expression_Group %in% c("Top 5% Expressed", "Bottom 95%")) |>  # Updated group names
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

message("\n=== Multivariate Analysis Complete ===")
message("\nPlots saved to ./results/:")
message("  - Variance explained plots (CA and PCA)")
message("  - All groups comparison (Top 5%, Middle 90%, Bottom 5%)")
message("  - Extremes only comparison (Top 5% vs Bottom 5%)")
message("  - Biplots showing codon loadings (AT vs GC pattern)")
message("\nCheck the statistical tests above for significance of group separation.")

## *****************************************************************************
## 9) Analyze codon loading patterns (AT vs GC bias) ----
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

cat("\n=== Interpretation ===\n")
cat("If AT-ending and GC-ending codons load in opposite directions:\n")
cat("  → This indicates MUTATIONAL BIAS is the dominant force\n")
cat("  → Genome-wide GC content variation drives the main axis\n")
cat("\nIf selection for expression is important, look for:\n")
cat("  → Small but significant separation in extremes (Top 5% vs Bottom 5%)\n")
cat("  → Lower ENC in highly expressed genes (already observed!)\n")
cat("  → Specific codons preferred in highly expressed genes\n")

## *****************************************************************************
## 10) Calculate Codon Adaptation Index (CAI) ----
## _____________________________________________________________________________

cat("\n=== Step 10: Codon Adaptation Index (CAI) ===\n")
cat("CAI measures the degree of bias towards codons preferred in highly expressed genes\n")
cat("CAI ranges from 0 to 1, where higher values indicate stronger adaptation\n")
cat("Higher CAI = more similar to codon usage in highly expressed genes\n\n")

source("./src/calculate_cai.R")

# Define reference set: Top 5% expressed genes
reference_genes <- exp_bud_enc |>
  filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

cat(sprintf("Using %d highly expressed genes as reference set\n", length(reference_genes)))

# Calculate CAI for all genes
cai_results <- calculate_cai(
  codon_counts = codon_usage,
  reference_genes = reference_genes,
  genetic_code = genetic_code_dna_long
)

# Extract CAI values and merge with expression data
cai_values <- cai_results$cai_values
w_table <- cai_results$w_table

# Merge CAI with expression and ENC data
exp_bud_enc_cai <- exp_bud_enc |>
  left_join(cai_values, by = "Gene_name")

# Save results
write.csv(exp_bud_enc_cai, "./results/expression_enc_cai.csv", row.names = FALSE)
write.csv(w_table, "./results/optimal_codons_weights.csv", row.names = FALSE)

cat("\n=== CAI vs Expression Level ===\n")
# Compare CAI across expression groups
cai_by_group <- exp_bud_enc_cai |>
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
kw_test <- kruskal.test(CAI ~ Expression_Group, data = exp_bud_enc_cai)
print(kw_test)

if (kw_test$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with Bonferroni Correction ===\n")
  
  # Install and load dunn.test if not available
  if (!require("dunn.test", quietly = TRUE)) {
    cat("Installing dunn.test package...\n")
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
    library(dunn.test)
  }
  
  # Perform Dunn's test with Bonferroni correction
  dunn_result <- dunn.test::dunn.test(
    x = exp_bud_enc_cai$CAI,
    g = exp_bud_enc_cai$Expression_Group,
    method = "bonferroni",
    kw = TRUE,
    label = TRUE,
    wrap = FALSE,
    table = TRUE,
    list = FALSE,
    altp = TRUE
  )
  
  cat("\nInterpretation of pairwise comparisons:\n")
  cat("  - Adjusted p-values account for multiple testing (Bonferroni)\n")
  cat("  - p < 0.05 indicates significant difference between groups\n")
  
} else {
  cat("\nNo significant difference among groups (p >= 0.05)\n")
  cat("Post-hoc tests not necessary.\n")
}

# Additional pairwise effect sizes
cat("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

# Helper function to calculate Cohen's d
cohens_d_calc <- function(x1, x2) {
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- var(x1, na.rm = TRUE)
  s2 <- var(x2, na.rm = TRUE)
  pooled_sd <- sqrt((s1 + s2) / 2)
  d <- (m1 - m2) / pooled_sd
  return(d)
}

# Get CAI values for each group
top_cai <- exp_bud_enc_cai |> filter(Expression_Group == "Top 5%") |> pull(CAI)
middle_cai <- exp_bud_enc_cai |> filter(Expression_Group == "Middle 90%") |> pull(CAI)
bottom_cai <- exp_bud_enc_cai |> filter(Expression_Group == "Bottom 5%") |> pull(CAI)

# If groups don't exist with new names, try old names
if (length(top_cai) == 0) {
  top_cai <- exp_bud_enc_cai |> filter(Expression_Group == "Top 5% Expressed") |> pull(CAI)
}
if (length(bottom_cai) == 0) {
  bottom_cai <- exp_bud_enc_cai |> filter(Expression_Group == "Bottom 95%") |> pull(CAI)
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
library(ggplot2)
p_cai_boxplot <- ggplot(exp_bud_enc_cai, aes(x = Expression_Group, y = CAI, fill = Expression_Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5% Expressed" = "#E41A1C", 
                                "Bottom 95%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "Codon Adaptation Index by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR. Higher CAI = more adapted to highly expressed genes",
       y = "CAI (Codon Adaptation Index)",
       x = "Expression Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("./results/CAI_by_expression_group.pdf", p_cai_boxplot, width = 8, height = 6)
cat("\nBoxplot saved: ./results/CAI_by_expression_group.pdf\n")

# Correlation between CAI and other metrics
cat("\n=== Correlations ===\n")
cat(sprintf("CAI vs Expression: r = %.4f\n", cor(exp_bud_enc_cai$CAI, exp_bud_enc_cai$Expression, use = "complete.obs")))
cat(sprintf("CAI vs ENC: r = %.4f\n", cor(exp_bud_enc_cai$CAI, exp_bud_enc_cai$ENC, use = "complete.obs")))
cat(sprintf("ENC vs Expression: r = %.4f\n", cor(exp_bud_enc_cai$ENC, exp_bud_enc_cai$Expression, use = "complete.obs")))

# Scatter plot: CAI vs ENC
p_cai_enc <- ggplot(exp_bud_enc_cai, aes(x = ENC, y = CAI, color = Expression_Group)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_color_manual(values = c("Top 5% Expressed" = "#E41A1C", 
                                 "Bottom 95%" = "#377EB8",
                                 "Middle 90%" = "#999999")) +
  labs(title = "CAI vs ENC by Expression Level",
       subtitle = "Lower ENC and Higher CAI indicate stronger codon bias",
       x = "ENC (Effective Number of Codons)",
       y = "CAI (Codon Adaptation Index)",
       color = "Expression Group") +
  theme_minimal(base_size = 12)

ggsave("./results/CAI_vs_ENC_scatter.pdf", p_cai_enc, width = 10, height = 6)