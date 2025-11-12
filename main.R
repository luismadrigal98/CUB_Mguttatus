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
                        'DescTools', 'mgcv')

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
                          output_dir = "./results",
                          aa_group = aa_chemistry_df)

## *****************************************************************************
## 5) Modeling relationship between ENC and Expression profiles ----
## _____________________________________________________________________________

# Trimming suffix from ENC table in gene names
cub_results$enc_results[, Gene_name := sub("\\.1$", "", Gene_name)]

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
  dplyr::mutate(Source_High_exp = as.factor(Source_High_exp))

# Add the ENC values per gene
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
  dplyr::mutate(High_exp_log2 = log2(High_exp + 1))  # Adding 1 to avoid log2(0)

# Linear models ----

ENC_vs_exp <- lm(ENC ~ High_exp_log2 + GC3 + CDS_length_nt, 
                 data = integrated_data)
summary(ENC_vs_exp)

# Density plots

# Enc ~ Exp
ggplot(data = integrated_data, 
       mapping = aes(x = High_exp_log2, y = ENC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/ENC_raw_vs_expression_density.pdf", 
       width = 10, height = 8)

# Enc ~ GC3
ggplot(data = integrated_data, 
       mapping = aes(x = GC3, y = ENC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/ENC_raw_vs_GC3_density.pdf", 
       width = 10, height = 8)

# Enc ~ Gene length
ggplot(data = integrated_data, 
       mapping = aes(x = CDS_length_nt, y = ENC)) +
  geom_pointdensity() +
  geom_smooth(method = lm, color = 'red') +
  theme_custom()

ggsave("./results/ENC_raw_vs_gene_length_density.pdf", 
       width = 10, height = 8)

# GAM models ----

# Given the non-lineariry effect of main confounders gene length and GC3 content
# we are going to fit a GAM model to account for those and assess effectively the
# effect of expression

ENC_exp_and_conf_gam <- gam(ENC ~ High_exp_log2 + s(CDS_length_nt) + s(GC3),
                            data = integrated_data)

summary(ENC_exp_and_conf_gam)

# Plotting detrended ENC against expression

confounder_model_gam <- gam(ENC ~ s(CDS_length_nt) + s(GC3s),
                            data = integrated_data)

integrated_data$ENC_detrended <- residuals(confounder_model_gam)

p_detrended <- ggplot(integrated_data, aes(x = High_exp_log2, y = ENC_detrended)) +
  # Use ggpointdensity for a clear view of the cluster
  geom_pointdensity(alpha = 0.5) + 
  
  # Add the linear regression line, which now shows the true effect
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  
  labs(
    title = "Detrended ENC vs. Gene Expression",
    subtitle = "Showing ENC after accounting for non-linear effects of GC3 and gene length",
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

p_boxplot <- ggplot(integrated_data, aes(x = Expression_Group, y = ENC, fill = Expression_Group)) +
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
kw_test_enc <- kruskal.test(ENC ~ Expression_Group, data = integrated_data)
print(kw_test_enc)

if (kw_test_enc$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Install and load dunn.test if not available
  if (!require("dunn.test", quietly = TRUE)) {
    cat("Installing dunn.test package...\n")
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
    library(dunn.test)
  }
  
  # Perform Dunn's test with FDR correction
  dunn_result_enc <- dunn.test::dunn.test(
    x = integrated_data$ENC,
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
    mean_ENC = mean(ENC, na.rm = TRUE),
    median_ENC = median(ENC, na.rm = TRUE),
    sd_ENC = sd(ENC, na.rm = TRUE),
    mean_Expression = mean(High_exp, na.rm = TRUE)
  )
print(summary_stats)

# Effect sizes for pairwise comparisons
cat("\n=== Effect Sizes (Cohen's d) for Pairwise Comparisons ===\n")

# Get ENC values for each group
top5_enc <- integrated_data |> filter(Expression_Group == "Top 5%") |> pull(ENC)
middle_enc <- integrated_data |> filter(Expression_Group == "Middle 90%") |> pull(ENC)
bottom5_enc <- integrated_data |> filter(Expression_Group == "Bottom 5%") |> pull(ENC)

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

# Confounding out-based analysis (detendred ENC) ----

# Assesing significance of expression over the detrended residuals

cat("\n=== Kruskal-Wallis Test: Detrended ENC Residuals across Groups ===\n")

kw_detrended <- kruskal.test(ENC_detrended ~ Expression_Group, 
                             data = integrated_data)

# Plotting and assessing significance using Dunn

print(kw_detrended)
if (kw_detrended$p.value < 0.05) {
  cat("\nSignificant difference detected! Performing post-hoc pairwise comparisons...\n")
  cat("\n=== Dunn's Test: Pairwise Comparisons with FDR Correction ===\n")
  
  # Perform Dunn's test with FDR correction
  dunn_result_detrended <- dunn.test::dunn.test(
    x = integrated_data$ENC_detrended,
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

p_boxplot_detrended <- ggplot(integrated_data, aes(x = Expression_Group, y = ENC_detrended, fill = Expression_Group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(outlier.alpha = 0.3) +
  # geom_boxplot(outlier.alpha = 0.3) +
  # geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
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

ggsave("./results/Detrended_ENC_by_expression_group.pdf", 
       p_boxplot_detrended, width = 8, height = 6)

# Get ENC values for each group
top5_enc_de <- integrated_data |> filter(Expression_Group == "Top 5%") |> pull(ENC_detrended)
middle_enc_de <- integrated_data |> filter(Expression_Group == "Middle 90%") |> pull(ENC_detrended)
bottom5_enc_de <- integrated_data |> filter(Expression_Group == "Bottom 5%") |> pull(ENC_detrended)

# Calculate effect sizes
if (length(top5_enc_de) > 0 && length(middle_enc_de) > 0) {
  d_top_middle_de <- cohens_d_calc(top5_enc_de, middle_enc_de)
  cat(sprintf("Top 5%% vs Middle 90%%: d = %.3f\n", d_top_middle_de))
}

if (length(top5_enc_de) > 0 && length(bottom5_enc_de) > 0) {
  d_top_bottom_de <- cohens_d_calc(top5_enc_de, bottom5_enc_de)
  cat(sprintf("Top 5%% vs Bottom 5%%: d = %.3f\n", d_top_bottom_de))
}

if (length(middle_enc_de) > 0 && length(bottom5_enc_de) > 0) {
  d_middle_bottom_de <- cohens_d_calc(middle_enc_de, bottom5_enc_de)
  cat(sprintf("Middle 90%% vs Bottom 5%%: d = %.3f\n", d_middle_bottom_de))
}

cat("\nInterpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large\n")

## *****************************************************************************
## 6) Calculate Codon Adaptation Index (CAI) ----
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
p_cai_enc <- ggplot(integrated_data, aes(x = ENC, y = CAI, color = Expression_Group)) +
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

# 6.2) Comparing preferred codon of Mimulus guttatus to other plants ----

# Use w_table from CAI analysis (already calculated preferred codons)
cat("Using optimal codons from CAI reference set...\n")

# Get preferred codons for export
preferred_codons_mg <- w_table |>
  dplyr::filter(relative_adaptiveness == 1.0) |>
  dplyr::select(codon)

write.table(x = preferred_codons_mg, file = 'results/preferred_codons.txt', 
            col.names = F, row.names = F, quote = F)

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

## 6.3) Preferred codon usage: Selected vs Neutral genes ----

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
               color = "gray70", size = 1) +
  geom_point(aes(y = Selected_prop, color = "Top 5%"), 
             size = 3, shape = 16) +
  geom_point(aes(y = Rest_prop, color = "Rest 95%"), 
             size = 3, shape = 16) +
  scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                "Rest 95%" = "#377EB8"),
                    name = "") +
  labs(title = "Preferred Codon Usage: Top 5% vs Rest",
       subtitle = "Proportion of preferred codons per amino acid",
       x = "Amino Acid (ordered by difference)",
       y = "Proportion of Preferred Codons") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("./results/preferred_codon_usage_comparison.pdf", p_comparison,
       width = 10, height = 6)

cat("\n")

# 6.3.1) Statistical testing for codon-level differences ----

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

# 6.3.2) Diagnose CAI vs Proportion Discrepancies ----

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
codon_test_results <- codon_test_results %>%
  left_join(corrected_classification %>% 
              dplyr::select(Codon, Combined_Classification),
            by = "Codon") %>%
  dplyr::mutate(
    # Update Classification to be more accurate
    Classification_Original = Classification,
    Classification = Combined_Classification
  )

# 6.4) Split 6-codon amino acids by degeneracy ----

# Define codon families for 6-codon amino acids
# Arg: CGN (4-fold) + AGA/AGG (2-fold)
# Leu: CTN (4-fold) + TTA/TTG (2-fold)  
# Ser: TCN (4-fold) + AGT/AGC (2-fold)

codon_families <- list(
  Arg_4fold = c("CGT", "CGC", "CGA", "CGG"),
  Arg_2fold = c("AGA", "AGG"),
  Leu_4fold = c("CTT", "CTC", "CTA", "CTG"),
  Leu_2fold = c("TTA", "TTG"),
  Ser_4fold = c("TCT", "TCC", "TCA", "TCG"),
  Ser_2fold = c("AGT", "AGC")
)

cat("Codon families:\n")
for (family in names(codon_families)) {
  cat(sprintf("  %s: %s\n", family, paste(codon_families[[family]], collapse = ", ")))
}
cat("\n")

# Calculate for both groups
cat("Calculating preferred codon usage by degeneracy family...\n")

selected_families <- count_preferred_by_family(top5_genes, preferred_codons_vec, codon_families)
selected_families$Group <- "Top 5%"

rest_families <- count_preferred_by_family(rest_genes, preferred_codons_vec, codon_families)
rest_families$Group <- "Rest 95%"

# Combine for comparison
family_comparison <- selected_families |>
  dplyr::select(Amino_Acid, Degeneracy, Family, N_codons, Preferred_codons,
                Selected_prop = Prop_preferred) |>
  dplyr::left_join(
    rest_families |> dplyr::select(Family, Rest_prop = Prop_preferred),
    by = "Family"
  ) |>
  dplyr::mutate(
    Difference = Selected_prop - Rest_prop,
    Fold_enrichment = Selected_prop / Rest_prop
  ) |>
  dplyr::arrange(Amino_Acid, Degeneracy)

# Save table
write.csv(family_comparison, "./results/preferred_codon_usage_by_degeneracy.csv",
          row.names = FALSE)

cat("\n✓ Results saved: ./results/preferred_codon_usage_by_degeneracy.csv\n\n")

# Print table
cat("=== Preferred Codon Usage by Degeneracy Level ===\n\n")
cat(sprintf("%-4s %-8s %-4s %-20s %-12s %-12s %-12s %-8s\n",
            "AA", "Degen", "N", "Preferred", "Top5%", "Rest95%", "Difference", "Fold"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (i in 1:nrow(family_comparison)) {
  row <- family_comparison[i, ]
  cat(sprintf("%-4s %-8s %-4d %-20s %-12.4f %-12.4f %-12.4f %-8.2f\n",
              row$Amino_Acid,
              row$Degeneracy,
              row$N_codons,
              substr(row$Preferred_codons, 1, 20),
              row$Selected_prop,
              row$Rest_prop,
              row$Difference,
              row$Fold_enrichment))
}
cat(paste(rep("-", 90), collapse = ""), "\n\n")

# Statistical comparison: 2-fold vs 4-fold
cat("=== Comparison: 2-fold vs 4-fold Degeneracy ===\n\n")

fold2 <- family_comparison |> dplyr::filter(Degeneracy == "2fold")
fold4 <- family_comparison |> dplyr::filter(Degeneracy == "4fold")

cat(sprintf("2-fold degenerate families (n=%d):\n", nrow(fold2)))
cat(sprintf("  Mean difference (Top5%% - Rest95%%): %.4f\n", mean(fold2$Difference, na.rm = TRUE)))
cat(sprintf("  Mean fold enrichment: %.2f\n\n", mean(fold2$Fold_enrichment, na.rm = TRUE)))

cat(sprintf("4-fold degenerate families (n=%d):\n", nrow(fold4)))
cat(sprintf("  Mean difference (Top5%% - Rest95%%): %.4f\n", mean(fold4$Difference, na.rm = TRUE)))
cat(sprintf("  Mean fold enrichment: %.2f\n\n", mean(fold4$Fold_enrichment, na.rm = TRUE)))

# Wilcoxon test
if (nrow(fold2) > 0 && nrow(fold4) > 0) {
  wilcox_degen <- wilcox.test(fold2$Difference, fold4$Difference)
  
  cat(sprintf("Wilcoxon test (2-fold vs 4-fold difference):\n"))
  cat(sprintf("  W = %.1f, p-value = %.4f\n", 
              wilcox_degen$statistic, wilcox_degen$p.value))
  
  if (wilcox_degen$p.value < 0.05) {
    cat("  * Significant difference between degeneracy levels\n")
  } else {
    cat("  Not significant (p >= 0.05)\n")
  }
}

# Create visualization
p_degeneracy <- ggplot(family_comparison, 
                       aes(x = Family, fill = Degeneracy)) +
  geom_segment(aes(xend = Family, y = Rest_prop, yend = Selected_prop),
               color = "gray70", size = 1) +
  geom_point(aes(y = Selected_prop, shape = "Top 5%"), 
             size = 3, color = "#E41A1C") +
  geom_point(aes(y = Rest_prop, shape = "Rest 95%"), 
             size = 3, color = "#377EB8") +
  scale_fill_manual(values = c("2fold" = "#FDB462", "4fold" = "#8DD3C7"),
                    name = "Degeneracy") +
  scale_shape_manual(values = c("Top 5%" = 16, "Rest 95%" = 17),
                     name = "Group") +
  facet_wrap(~ Amino_Acid, scales = "free_x", ncol = 3) +
  labs(title = "Preferred Codon Usage: 2-fold vs 4-fold Degenerate Families",
       subtitle = "6-codon amino acids (Arg, Leu, Ser) split by degeneracy level",
       x = "", y = "Proportion of Preferred Codons") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 11))

ggsave("./results/preferred_codon_usage_degeneracy.pdf", p_degeneracy,
       width = 12, height = 6)

cat("\n✓ Plot saved: ./results/preferred_codon_usage_degeneracy.pdf\n\n")

## *****************************************************************************
## 7) Correspondence analysis over counts and PCA over RSCU ----
## _____________________________________________________________________________

# 7.1) CA analysis ---- 

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
  w_table = w_table,
  dims = c(1, 2),
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
  dims = c(1, 2)
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

# 7.2) PCA analysis ----

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
  dims = c(1, 2)
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

## *****************************************************************************
## xx) tRNA abundance correlation analysis ----
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

## *****************************************************************************
## xx) CDC-based analysis ----
## _____________________________________________________________________________

# Full integration with your pipeline
cdc_results <- integrate_cdc_analysis(codon_usage, genetic_code_dna_long, 
                                      integrated_data, n_bootstrap = 10000,
                                      n_cores = 10)

# Re-plotting ENC-based neutrality plot highlighting the significant genes with CDC ----

cat("\n=== Creating Enhanced ENC Plot with CDC-Significant Genes ===\n")
cat("Highlighting genes deviating from neutral codon usage (significant CDC)\n\n")

# Check what columns cdc_results has
cat("CDC results columns:", paste(names(cdc_results), collapse = ", "), "\n")
cat(sprintf("CDC results has %d rows\n", nrow(cdc_results)))

# Clean gene names: remove .1 suffix from all data frames
enc_values_clean <- cub_results$enc_results |>
  dplyr::mutate(Gene_name = sub("\\.1$", "", Gene_name))

gc_content_clean <- cub_results$gc_results |>
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

# Relationship between CDC and expression levels ----

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
expression_df <- integrated_data |>
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
  left_join(integrated_data |> select(Gene_name, Expression_Group), by = "Gene_name")

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

# Add polymorphism data integration here if available (not implemented)

## xx) Polymorphism data integration ----

pi_data <- fread(input = "data/all_chromosomes.bygene.pi.txt")

# Homogenizing gene names to match the previous convention

pi_data <- pi_data |>
  dplyr::select(Chr, Gene, contains("Tajima"), contains("mean")) |>
  dplyr::mutate(Gene = paste0("MgIM767.", pi_data[['Gene']])) |>
  dplyr::rename(Gene_name = Gene)

# Analyzing the diversity of synonymous sites 

integrated_data <- integrated_data |>
  left_join(y = pi_data, by = "Gene_name")

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

# Bringing in the CDC results

integrated_data <- integrated_data |>
  left_join(cdc_results[, c("Gene_name", "CDC", "p_adj")], by = 'Gene_name')

int_variables <- integrated_data |>
  dplyr::select("CDC", "Pi_mean_0fold", "Pi_mean_2fold", "Pi_mean_3fold",
                "Pi_mean_4fold", "Pi_mean_all", "p_adj", "ENC") |>
  as.matrix()

int_cor <- corrr::correlate(x = int_variables)

save.image('Env')