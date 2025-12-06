check_aa_frequency_vs_tRNA_supply <- function(codon_usage,
                                              tRNA_file,
                                              genetic_code,
                                              output_dir = "./results/aa_trna_sanity_check") {
  #' Check if amino acid usage frequency matches tRNA gene supply
  #' 
  #' This is a sanity check to verify that amino acids with more tRNA genes
  #' are used more frequently in the genome, as expected under the tRNA
  #' adaptation hypothesis.
  #' 
  #' @param codon_usage Codon usage data frame (genes x codons)
  #' @param tRNA_file Path to tRNA gene data file
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_dir Directory to save results
  #' 
  #' @return List with correlation results and plot
  
  suppressPackageStartupMessages({
    require(dplyr)
    require(ggplot2)
    require(data.table)
  })
  
  cat("\n=== Amino Acid Frequency vs tRNA Supply Sanity Check ===\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load tRNA data
  tRNA_data <- fread(tRNA_file)
  cat(sprintf("Loaded %d tRNA genes\n", nrow(tRNA_data)))
  
  # Count tRNA genes per amino acid
  # Column name is tRNA_type (lowercase) based on the data structure
  tRNA_counts <- tRNA_data %>%
    dplyr::group_by(tRNA_type) %>%
    dplyr::summarise(
      tRNA_gene_count = n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(Amino_Acid = tRNA_type) %>%
    dplyr::filter(Amino_Acid != "iMet")  # Remove initiator Met (keep regular Met)
  
  cat(sprintf("Found %d amino acids with tRNA genes\n", nrow(tRNA_counts)))
  
  # Calculate amino acid usage frequency across all genes
  codon_cols <- names(codon_usage)[names(codon_usage) %in% names(genetic_code)]
  
  # Sum codon counts across all genes
  total_codon_counts <- codon_usage %>%
    dplyr::select(all_of(codon_cols)) %>%
    summarise(across(everything(), ~sum(., na.rm = TRUE))) %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "Codon",
      values_to = "Total_Count"
    )
  
  # Map codons to amino acids
  total_codon_counts$Amino_Acid <- genetic_code[total_codon_counts$Codon]
  
  # Sum by amino acid
  aa_usage <- total_codon_counts %>%
    dplyr::filter(!is.na(Amino_Acid) & Amino_Acid != "STOP") %>%
    dplyr::group_by(Amino_Acid) %>%
    dplyr::summarise(
      Total_AA_Count = sum(Total_Count),
      N_Codons = n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      AA_Frequency = Total_AA_Count / sum(Total_AA_Count)
    )
  
  cat(sprintf("Calculated usage for %d amino acids\n", nrow(aa_usage)))
  
  # Merge tRNA counts with AA usage
  aa_trna_data <- aa_usage %>%
    left_join(tRNA_counts, by = "Amino_Acid") %>%
    dplyr::filter(!is.na(tRNA_gene_count))
  
  cat(sprintf("Matched %d amino acids with both usage and tRNA data\n\n", 
              nrow(aa_trna_data)))
  
  # Calculate correlations
  cor_pearson <- cor.test(aa_trna_data$AA_Frequency, 
                          aa_trna_data$tRNA_gene_count,
                          method = "pearson")
  
  cor_spearman <- cor.test(aa_trna_data$AA_Frequency, 
                           aa_trna_data$tRNA_gene_count,
                           method = "spearman",
                           exact = F)
  
  # Print results
  cat("=== Correlation Results ===\n")
  cat(sprintf("Pearson r = %.3f, p-value = %.2e\n", 
              cor_pearson$estimate, cor_pearson$p.value))
  cat(sprintf("Spearman rho = %.3f, p-value = %.2e\n\n", 
              cor_spearman$estimate, cor_spearman$p.value))
  
  if (cor_pearson$p.value < 0.05) {
    cat("✓ PASS: Significant positive correlation between AA frequency and tRNA supply\n")
    cat("  This supports the tRNA adaptation hypothesis\n\n")
  } else {
    cat("⚠ WARNING: No significant correlation detected\n")
    cat("  This may indicate other factors influencing AA usage\n\n")
  }
  
  # Create visualization
  p <- ggplot(aa_trna_data, aes(x = tRNA_gene_count, y = AA_Frequency)) +
    geom_point(aes(size = Total_AA_Count), alpha = 0.6, color = "#1F77B4") +
    geom_smooth(method = "lm", se = TRUE, color = "#D62728", linewidth = 1) +
    geom_text(aes(label = Amino_Acid), vjust = -0.5, hjust = 0.5, size = 3.5) +
    labs(
      title = "Amino Acid Usage vs tRNA Gene Supply",
      subtitle = sprintf("Pearson r = %.3f (p = %.2e), Spearman rho = %.3f (p = %.2e)",
                         cor_pearson$estimate, cor_pearson$p.value,
                         cor_spearman$estimate, cor_spearman$p.value),
      x = "Number of tRNA Genes",
      y = "Amino Acid Frequency in Genome",
      size = "Total Count",
      caption = "Point size represents total amino acid usage across all genes"
    ) +
    theme_custom() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, "aa_frequency_vs_tRNA_supply.pdf"),
    plot = p,
    width = 10,
    height = 8
  )
  
  cat(sprintf("✓ Plot saved: %s/aa_frequency_vs_tRNA_supply.pdf\n", output_dir))
  
  # Save data
  write.csv(
    aa_trna_data,
    file = file.path(output_dir, "aa_frequency_vs_tRNA_supply.csv"),
    row.names = FALSE
  )
  
  cat(sprintf("✓ Data saved: %s/aa_frequency_vs_tRNA_supply.csv\n\n", output_dir))
  
  # Return results
  return(list(
    data = aa_trna_data,
    pearson = cor_pearson,
    spearman = cor_spearman,
    plot = p
  ))
}
