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
                        'ape', 'tidyr')

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

# Preferred codons in three additional model plants

model_plants_PC <- read.table(file = "data/plant_preferred_codons.txt", 
                              header = T, sep = ',')

## *****************************************************************************
## 3) Load the data ----
## _____________________________________________________________________________

## 3.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
                                      format = 'fasta')

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## 3.2) Clean the codon usage object (remove the STOP codon, alongside Trp and Met) ----

codon_usage <- codon_usage |>
  trim_uninformative(genetic_code = genetic_code_dna_long)

## *****************************************************************************
## 4) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results")

## *****************************************************************************
## 5) Additional analyses (optional) ----
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
                     "results/custom_heatmap.pdf", type = "heatmap")
neutrality_plot(gc_content, "results/custom_neutrality.pdf")
enc_plot(enc_values, gc_content, "results/custom_enc.pdf")
pr2_bias_plot(codon_usage, "results/custom_pr2.pdf")

message("\nAnalysis complete! Check the './results' directory for all outputs.")

## *****************************************************************************
## 6) Modeling relationship between ENC and Expression profiles ----
## _____________________________________________________________________________

# Trimming suffix from ENC table in gene names
enc_values[, Gene_name := sub("\\.1$", "", Gene_name)]

exp_data_bud <- read.table(file = "./data/bud_gene_expression_cpm_remapped.txt",
                       header = T) |>
  dplyr::select(Remapped_Gene, Expression) |>
  dplyr::rename(Gene = Remapped_Gene,
                Exp_bud = Expression)

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
  dplyr::mutate(Source_High_exp = as.factor(Source_High_exp))

# Add the ENC values per gene

exp_enc_data <- dplyr::left_join(exp_complete, enc_values, 
                                by = dplyr::join_by(Gene == Gene_name)) |>
  dplyr::rename(Gene_name = Gene) |>
  na.omit() |>
  # Remove duplicates - keep first occurrence of each gene
  # (duplicates come from bud expression having multiple entries per gene)
  distinct(Gene_name, .keep_all = TRUE)

cat(sprintf("Removed %d duplicate gene entries\n", 
            nrow(dplyr::left_join(exp_complete, enc_values, by = dplyr::join_by(Gene == Gene_name))) - nrow(exp_enc_data)))

# Add gene length (CDS length in codons and nucleotides)
cat("\n=== Adding Gene Length Information ===\n")

# First, clean gene names in codon_usage to match exp_enc_data
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

exp_enc_data <- exp_enc_data |>
  left_join(gene_lengths, by = "Gene_name")

cat(sprintf("Added length for %d genes\n", sum(!is.na(exp_enc_data$CDS_length_nt))))
cat(sprintf("Mean CDS length: %.0f nt (%.0f codons)\n", 
            mean(exp_enc_data$CDS_length_nt, na.rm = TRUE),
            mean(exp_enc_data$Total_Codons, na.rm = TRUE)))

log2exp <- log2(exp_enc_data$High_exp + 1)
cor(x = exp_enc_data$ENC, y = log2exp)
plot(exp_enc_data$ENC, log2exp)
lm(log2(High_exp+1) ~ ENC, data = exp_enc_data)

# Define expression groups: Top 5% vs Bottom 5% (extreme comparison)

top_5_cutoff <- quantile(exp_enc_data$High_exp, probs = 0.95)
bottom_5_cutoff <- quantile(exp_enc_data$High_exp, probs = 0.05)

exp_enc_data$Expression_Group <- case_when(
  exp_enc_data$High_exp >= top_5_cutoff ~ "Top 5%",
  exp_enc_data$High_exp <= bottom_5_cutoff ~ "Bottom 5%",
  TRUE ~ "Middle 90%"
)

# Boxplot comparison

p_boxplot <- ggplot(exp_enc_data, aes(x = Expression_Group, y = ENC, fill = Expression_Group)) +
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
  theme_custom() +
  theme(legend.position = "none")

ggsave("./results/ENC_by_expression_group.pdf", p_boxplot, width = 8, height = 6)

# Statistical tests for three groups
cat("\n=== Kruskal-Wallis Test: ENC across All Three Groups ===\n")
cat("H0: All three groups have the same median ENC\n")
kw_test_enc <- kruskal.test(ENC ~ Expression_Group, data = exp_enc_data)
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
    x = exp_enc_data$ENC,
    g = exp_enc_data$Expression_Group,
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
summary_stats <- exp_enc_data |>
  group_by(Expression_Group) |>
  summarise(
    n = n(),
    mean_ENC = mean(ENC, na.rm = TRUE),
    median_ENC = median(ENC, na.rm = TRUE),
    sd_ENC = sd(ENC, na.rm = TRUE),
    mean_Expression = mean(High_exp, na.rm = TRUE)
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
top5_enc <- exp_enc_data |> filter(Expression_Group == "Top 5%") |> pull(ENC)
middle_enc <- exp_enc_data |> filter(Expression_Group == "Middle 90%") |> pull(ENC)
bottom5_enc <- exp_enc_data |> filter(Expression_Group == "Bottom 5%") |> pull(ENC)

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
length_stats <- exp_enc_data |>
  group_by(Expression_Group) |>
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
kw_length <- kruskal.test(CDS_length_aa ~ Expression_Group, data = exp_enc_data)
print(kw_length)

# Correlation between length and ENC
cat("\n=== Correlation: Gene Length vs ENC ===\n")
cor_length_enc <- cor(exp_enc_data$CDS_length_aa, exp_enc_data$ENC, use = "complete.obs")
cat(sprintf("Pearson r = %.4f\n", cor_length_enc))

# Detrending the ENC and the length of genes

exp_enc_data <- exp_enc_data |>
  mutate(
    ENC_residuals = resid(lm(ENC ~ CDS_length_aa, data = exp_enc_data))
  )

# Assesing significance of expression over the detrended residuals

cat("\n=== Kruskal-Wallis Test: Detrended ENC Residuals across Groups ===\n")

kw_detrended <- kruskal.test(ENC_residuals ~ Expression_Group, 
                             data = exp_enc_data)

# Plotting and assessing significance using Dunn

print(kw_detrended)
if (kw_detrended$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with Bonferroni Correction ===\n")
  
  # Perform Dunn's test with Bonferroni correction
  dunn_result_detrended <- dunn.test::dunn.test(
    x = exp_enc_data$ENC_residuals,
    g = exp_enc_data$Expression_Group,
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

p_boxplot_detrended <- ggplot(exp_enc_data, aes(x = Expression_Group, y = ENC_residuals, fill = Expression_Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
                                "Middle 90%" = "#999999")) +
  labs(title = "Detrended ENC Residuals by Expression Level",
       subtitle = "Diamond = mean, box = median ± IQR",
       y = "ENC Residuals (detrended)",
       x = "Expression Group") +
  theme_custom() +
  theme(legend.position = "none")

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

codon_usage_CA_coord <- exp_enc_data |>
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

rscu_PCA_coord <- exp_enc_data |>
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

# Define reference set: Top 5% expressed genes
reference_genes <- exp_enc_data |>
  filter(Expression_Group == "Top 5%") |>
  pull(Gene_name)

cat(sprintf("Using %d highly expressed genes as reference set\n", length(reference_genes)))

# Remove .1 suffix from codon_usage gene names to match gene-level IDs
# (codon_usage has transcript IDs like MgIM767.10G127000.1,
#  expression data has gene IDs like MgIM767.10G127000)
codon_usage <- codon_usage |>
  mutate(Gene_name = sub("\\.1$", "", Gene_name))

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

# 10.1) Making bar plot per aminoacid to show differences in use

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

# Merge CAI with expression and ENC data
exp_enc_data_cai <- exp_enc_data |>
  left_join(cai_values, by = "Gene_name")

# Save results
write.csv(exp_enc_data_cai, "./results/expression_enc_cai.csv", row.names = FALSE)
write.csv(w_table, "./results/optimal_codons_weights.csv", row.names = FALSE)

cat("\n=== CAI vs Expression Level ===\n")
# Compare CAI across expression groups
cai_by_group <- exp_enc_data_cai |>
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
kw_test <- kruskal.test(CAI ~ Expression_Group, data = exp_enc_data_cai)
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
    x = exp_enc_data_cai$CAI,
    g = exp_enc_data_cai$Expression_Group,
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
top_cai <- exp_enc_data_cai |> filter(Expression_Group == "Top 5%") |> pull(CAI)
middle_cai <- exp_enc_data_cai |> filter(Expression_Group == "Middle 90%") |> pull(CAI)
bottom_cai <- exp_enc_data_cai |> filter(Expression_Group == "Bottom 5%") |> pull(CAI)

# If groups don't exist with new names, try old names
if (length(top_cai) == 0) {
  top_cai <- exp_enc_data_cai |> filter(Expression_Group == "Top 5% Expressed") |> pull(CAI)
}
if (length(bottom_cai) == 0) {
  bottom_cai <- exp_enc_data_cai |> filter(Expression_Group == "Bottom 95%") |> pull(CAI)
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
p_cai_boxplot <- ggplot(exp_enc_data_cai, aes(x = Expression_Group, y = CAI, fill = Expression_Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("Top 5%" = "#E41A1C", 
                                "Bottom 5%" = "#377EB8",
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
cat(sprintf("CAI vs Expression: r = %.4f\n", cor(exp_enc_data_cai$CAI, exp_enc_data_cai$High_exp, use = "complete.obs")))
cat(sprintf("CAI vs ENC: r = %.4f\n", cor(exp_enc_data_cai$CAI, exp_enc_data_cai$ENC, use = "complete.obs")))
cat(sprintf("ENC vs Expression: r = %.4f\n", cor(exp_enc_data_cai$ENC, exp_enc_data_cai$High_exp, use = "complete.obs")))

# Scatter plot: CAI vs ENC
p_cai_enc <- ggplot(exp_enc_data_cai, aes(x = ENC, y = CAI, color = Expression_Group)) +
  geom_point(alpha = 0.3, size = 1) +
  # Add lm per group
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                 "Bottom 5%" = "#377EB8",
                                 "Middle 90%" = "#999999")) +
  labs(title = "CAI vs ENC by Expression Level",
       subtitle = "Lower ENC and Higher CAI indicate stronger codon bias",
       x = "ENC (Effective Number of Codons)",
       y = "CAI (Codon Adaptation Index)",
       color = "Expression Group") +
  theme_custom()

ggsave("./results/CAI_vs_ENC_scatter.pdf", p_cai_enc, width = 10, height = 6)

# 10.2) Comparing preferred codon of Mimulus guttatus to other plants ----

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  10.2 Preferred Codons: Cross-Species Comparison                ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Use w_table from CAI analysis (already calculated preferred codons)
cat("Using optimal codons from CAI reference set...\n")

# Get preferred codons (those with relative_adaptiveness == 1.0)
preferred_codons_mg <- w_table |>
  dplyr::filter(relative_adaptiveness == 1.0) |>
  dplyr::mutate(Codon_RNA = gsub("T", "U", codon)) |>
  dplyr::select(Amino_Acid = amino_acid, Codon_RNA, relative_adaptiveness)

cat(sprintf("Found %d preferred codons for M. guttatus\n\n", nrow(preferred_codons_mg)))

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

# Calculate Jaccard similarity
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
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
  pivot_longer(cols = all_of(species), names_col = "Species", values_to = "Preferred") |>
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

# Summary statistics
cat("=== Summary Statistics ===\n\n")

# How many codons does M. guttatus share with each species?
for (sp in species) {
  if (sp != "Mimulus_guttatus") {
    shared <- sum(codon_matrix["Mimulus_guttatus", ] & codon_matrix[sp, ])
    total <- length(all_codons_rna)
    pct <- 100 * shared / total
    cat(sprintf("M. guttatus shares %d/%d (%.1f%%) preferred codons with %s\n",
                shared, total, pct, gsub("_", " ", sp)))
  }
}

cat("\n")

## *****************************************************************************
## xx) tRNA abundance correlation analysis ----
## _____________________________________________________________________________

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  tRNA Abundance and Codon Usage Correlation Analysis            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Load tRNA correlation functions
source("./src/tRNA_codon_correlation.R")
source("./src/get_codon_supply_map.R")

# Analysis 1: By tRNA gene copy number (traditional approach)
cat("=== Analysis 1: tRNA Gene Copy Number ===\n")
cat("Analyzing correlation using tRNA gene counts (traditional approach)\n\n")

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
expression_df <- exp_enc_data |>
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
top5_threshold <- quantile(exp_enc_data$High_exp, probs = 0.95, na.rm = TRUE)
top5_genes <- exp_enc_data |>
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

# Summary of all three analyses
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  Summary of tRNA-Codon Correlation Results                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

cat("1. tRNA Gene Copy Number (all genes):\n")
cat("   RSCU vs tRNA Supply:\n")
if (!is.null(tRNA_copynumber_results$correlation_results$overall)) {
  cor_val <- tRNA_copynumber_results$correlation_results$overall$estimate
  p_val <- tRNA_copynumber_results$correlation_results$overall$p.value
  cat(sprintf("     Spearman ρ = %.4f (p = %.2e)\n", cor_val, p_val))
  if (!is.null(tRNA_copynumber_results$significant_amino_acids) && 
      nrow(tRNA_copynumber_results$significant_amino_acids) > 0) {
    cat(sprintf("     %d amino acids significant (p < 0.05)\n", 
                nrow(tRNA_copynumber_results$significant_amino_acids)))
  }
} else {
  cat("     No correlation computed\n")
}

cat("\n2. tRNA Expression (all genes):\n")
cat("   RSCU vs tRNA Expression:\n")
if (!is.null(tRNA_expression_all_results$correlation_results$overall)) {
  cor_val <- tRNA_expression_all_results$correlation_results$overall$estimate
  p_val <- tRNA_expression_all_results$correlation_results$overall$p.value
  cat(sprintf("     Spearman ρ = %.4f (p = %.2e)\n", cor_val, p_val))
  if (!is.null(tRNA_expression_all_results$significant_amino_acids) && 
      nrow(tRNA_expression_all_results$significant_amino_acids) > 0) {
    cat(sprintf("     %d amino acids significant (p < 0.05)\n", 
                nrow(tRNA_expression_all_results$significant_amino_acids)))
  }
}
if (!is.null(tRNA_expression_all_results$tAI_analysis)) {
  cat("   tAI vs Gene Expression:\n")
  cor_val <- tRNA_expression_all_results$tAI_analysis$spearman$estimate
  p_val <- tRNA_expression_all_results$tAI_analysis$spearman$p.value
  cat(sprintf("     Spearman ρ = %.4f (p = %.2e)\n", cor_val, p_val))
}

cat("\n3. tRNA Expression (top 5% genes):\n")
cat("   RSCU vs tRNA Expression:\n")
if (!is.null(tRNA_expression_top5_results$correlation_results$overall)) {
  cor_val <- tRNA_expression_top5_results$correlation_results$overall$estimate
  p_val <- tRNA_expression_top5_results$correlation_results$overall$p.value
  cat(sprintf("     Spearman ρ = %.4f (p = %.2e)\n", cor_val, p_val))
  if (!is.null(tRNA_expression_top5_results$significant_amino_acids) && 
      nrow(tRNA_expression_top5_results$significant_amino_acids) > 0) {
    cat(sprintf("     %d amino acids significant (p < 0.05)\n", 
                nrow(tRNA_expression_top5_results$significant_amino_acids)))
  }
}
if (!is.null(tRNA_expression_top5_results$tAI_analysis)) {
  cat("   tAI vs Gene Expression:\n")
  cor_val <- tRNA_expression_top5_results$tAI_analysis$spearman$estimate
  p_val <- tRNA_expression_top5_results$tAI_analysis$spearman$p.value
  cat(sprintf("     Spearman ρ = %.4f (p = %.2e)\n", cor_val, p_val))
}

cat("\nResults saved to:\n")
cat("  - ./results/tRNA_analysis_copynumber/\n")
cat("  - ./results/tRNA_analysis_expression_all/\n")
cat("    • tRNA_codon_correlations.csv (per-AA correlations)\n")
cat("    • tAI_vs_expression.pdf (gene-level adaptation)\n")
cat("  - ./results/tRNA_analysis_expression_top5/\n\n")

## *****************************************************************************
## xx) CDC-based analysis ----
## _____________________________________________________________________________

# Full integration with your pipeline
cdc_results <- integrate_cdc_analysis(codon_usage, genetic_code_dna_long, 
                                      exp_enc_data, n_bootstrap = 10000,
                                      n_cores = 10)

# Re-plotting ENC-based neutrality plot highlighting the significant genes with CDC ----

cat("\n=== Creating Enhanced ENC Plot with CDC-Significant Genes ===\n")
cat("Highlighting genes deviating from neutral codon usage (significant CDC)\n\n")

# Check what columns cdc_results has
cat("CDC results columns:", paste(names(cdc_results), collapse = ", "), "\n")
cat(sprintf("CDC results has %d rows\n", nrow(cdc_results)))

# Clean gene names: remove .1 suffix from all data frames
enc_values_clean <- enc_values |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", Gene_name))

gc_content_clean <- gc_content |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", Gene_name))

# Extract just CDC columns we need
cdc_for_merge <- cdc_results |>
  dplyr::select(Gene_name, CDC, p_value, p_adj) |>
  dplyr::filter(!is.na(CDC))  # Remove genes without CDC

cat(sprintf("Valid CDC results: %d genes\n", nrow(cdc_for_merge)))

# Merge ENC, GC3s, and CDC results
enc_cdc_data <- enc_values_clean |>
  dplyr::left_join(gc_content_clean |> dplyr::select(Gene_name, GC3s), by = "Gene_name") |>
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
  geom_point(data = enc_cdc_data %>% filter(!CDC_significant | is.na(CDC_significant)),
             aes(x = GC3s, y = ENC), 
             color = "gray70", alpha = 0.3, size = 0.8) +
  # Foreground: CDC-significant genes
  geom_point(data = enc_cdc_data %>% filter(CDC_significant),
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
enc_cdc_data <- enc_cdc_data %>%
  mutate(
    ENC_expected = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2),
    ENC_deviation = ENC - ENC_expected,
    Below_curve = ENC_deviation < 0
  )

# Compare CDC-significant vs non-significant genes
cdc_position_summary <- enc_cdc_data %>%
  filter(!is.na(CDC_significant)) %>%
  group_by(CDC_significant) %>%
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
    data = enc_cdc_data %>% filter(!is.na(CDC_significant))
  )
  cat(sprintf("ENC deviation (CDC-sig vs non-sig): W = %.0f, p = %.2e\n", 
              wilcox_enc$statistic, wilcox_enc$p.value))
  
  # Test if more CDC-significant genes are below the curve
  below_curve_table <- table(
    enc_cdc_data %>% filter(!is.na(CDC_significant)) %>% select(CDC_significant, Below_curve)
  )
  chi_test <- chisq.test(below_curve_table)
  cat(sprintf("Position relative to curve (chi-squared): X² = %.2f, p = %.2e\n", 
              chi_test$statistic, chi_test$p.value))
}

# Create density plot of ENC deviation
p_enc_deviation <- ggplot(enc_cdc_data %>% filter(!is.na(CDC_significant)), 
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

cat("\nENC deviation density plot saved: ./results/ENC_deviation_by_CDC.pdf\n")

cat("\n✓ Enhanced ENC plot with CDC analysis complete!\n\n")



## *****************************************************************************
## xx) Selection Coefficient Analysis (Mutation-Selection-Drift Balance) ----
## _____________________________________________________________________________

cat("\n=== SELECTION COEFFICIENT ANALYSIS ===\n")
cat("Estimating population-scaled selection (S = 4Nes) using Hershberg & Petrov model\n\n")

# Load selection coefficient functions
source("./src/selection_coefficient_analysis.R")

# Get preferred codons from CAI analysis (w = 1.0)
preferred_codons <- cai_results$w_table |>
  filter(relative_adaptiveness == 1.0, amino_acid != "STOP") |>
  pull(codon)

cat(sprintf("Using %d optimal codons as 'preferred' codons:\n", length(preferred_codons)))
cat(paste(preferred_codons, collapse = ", "), "\n")

# Prepare expression data (using High_exp from bud tissue as primary metric)
expression_df <- exp_enc_data |>
  select(Gene_name, Expression = High_exp)

# Calculate selection coefficients for all genes
selection_results <- calculate_selection_coefficients(
  codon_usage = codon_usage,
  expression_data = expression_df,
  preferred_codons = preferred_codons,
  genetic_code = genetic_code_dna_long,
  low_expr_quantile = 0.10  # Use bottom 10% to estimate mutation bias
)

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
  left_join(exp_enc_data |> select(Gene_name, Expression_Group), by = "Gene_name")

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

save.image('Env')
