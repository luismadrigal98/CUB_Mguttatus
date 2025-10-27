tRNA_codon_correlation <- function(codon_counts, tRNA_file, genetic_code,
                                   output_dir = "./results", 
                                   test_method = "spearman")
{
  #' Analyze correlation between codon usage and tRNA abundance
  #' 
  #' @description Performs statistical tests to examine if codon usage bias
  #' can be explained by tRNA gene abundance. Calculates correlations between
  #' codon usage frequencies and tRNA gene copy numbers.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param tRNA_file Path to filtered tRNA annotation file (tab-separated)
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_dir Directory for output files
  #' @param test_method Correlation method: "spearman", "pearson", or "kendall"
  #' 
  #' @return List with correlation results, plots, and statistics
  #' ___________________________________________________________________________
  
  require(data.table)
  require(ggplot2)
  
  # Read tRNA data and get GCN for each anticodon
  trna_data <- fread(tRNA_file)
  trna_counts <- trna_data[, .(tRNA_count = .N), by = Anticodon]
  
  # Calculate codon supply based on wobble rules
  codon_supply <- get_codon_supply_map(trna_counts)
  
  # Calculate genome-wide codon usage (from the input subset)
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  total_codon_usage <- colSums(codon_counts[, codon_cols, with = FALSE])
  
  # Create analysis dataset
  analysis_data <- data.frame(
    Codon = names(total_codon_usage),
    Codon_usage = as.numeric(total_codon_usage),
    AA = genetic_code[names(total_codon_usage)],
    stringsAsFactors = FALSE
  )
  
  # Remove STOP, Met, and Trp
  analysis_data <- analysis_data[analysis_data$AA != "STOP" &
                                   analysis_data$AA != "Trp" &
                                   analysis_data$AA != "Met", ]
  
  # Merge with tRNA supply data
  analysis_data <- merge(analysis_data, codon_supply[, .(Codon, tRNA_supply)], 
                         by = "Codon", all.x = TRUE)
  analysis_data$tRNA_supply[is.na(analysis_data$tRNA_supply)] <- 0
  
  # Calculate relative frequencies (RSCU)
  analysis_dt <- as.data.table(analysis_data)
  analysis_dt[, Codon_freq := Codon_usage / sum(Codon_usage), by = AA]
  analysis_dt[, RSCU := Codon_freq / (1 / .N), by = AA]
  
  # Perform correlation tests
  correlation_results <- list()
  
  # Overall correlation
  if (sum(analysis_dt$tRNA_supply > 0) > 2) {
    overall_cor <- cor.test(analysis_dt$RSCU, analysis_dt$tRNA_supply, 
                            method = test_method, exact = FALSE)
    correlation_results$overall <- overall_cor
  }
  
  # Per amino acid correlations
  aa_correlations <- list()
  for (aa in unique(analysis_dt$AA)) {
    aa_data <- analysis_dt[AA == aa]
    # Check for variance in both variables to avoid errors
    if (nrow(aa_data) > 2 && sum(aa_data$tRNA_supply > 0) >= 2 && 
        var(aa_data$RSCU) > 0 && var(aa_data$tRNA_supply) > 0) {
      
      aa_cor <- cor.test(aa_data$RSCU, aa_data$tRNA_supply, 
                         method = test_method, exact = FALSE)
      aa_correlations[[aa]] <- aa_cor
    }
  }
  correlation_results$per_amino_acid <- aa_correlations
  
  # --- Create Visualizations ---
  
  # 1. Overall scatter plot
  p1 <- ggplot(analysis_dt, aes(x = tRNA_supply, y = RSCU, color = AA)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Codon Usage (RSCU) vs tRNA Supply (GCN)",
         subtitle = paste("Overall correlation:", 
                          ifelse(exists("overall_cor", correlation_results), 
                                 round(correlation_results$overall$estimate, 3), "N/A")),
         x = "tRNA Supply (Gene Copy Number)", 
         y = "Relative Synonymous Codon Usage (RSCU)",
         color = "Amino Acid") +
    theme(legend.position = "none")
  
  # 2. Per amino acid correlations bar plot
  aa_cor_data <- data.frame()
  if (length(aa_correlations) > 0) {
    for (aa in names(aa_correlations)) {
      aa_cor_data <- rbind(aa_cor_data, data.frame(
        AA = aa,
        Correlation = aa_correlations[[aa]]$estimate,
        P_value = aa_correlations[[aa]]$p.value,
        Significant = aa_correlations[[aa]]$p.value < 0.05
      ))
    }
  }
  
  p2 <- NULL # Initialize as NULL
  if (nrow(aa_cor_data) > 0) {
    p2 <- ggplot(aa_cor_data, aes(x = reorder(AA, Correlation), y = Correlation, 
                                  fill = Significant)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "red")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Codon Usage - tRNA Abundance Correlations by Amino Acid",
           subtitle = paste("Method:", test_method),
           x = "Amino Acid", 
           y = paste(tools::toTitleCase(test_method), "Correlation"),
           fill = "Significant (p < 0.05)")
  }
  
  # 3. Faceted scatter plots for significant correlations
  p3 <- NULL # Initialize as NULL
  if (nrow(aa_cor_data) > 0) {
    significant_aas <- aa_cor_data$AA[aa_cor_data$Significant]
    
    if (length(significant_aas) > 0) {
      sig_data <- analysis_dt[AA %in% significant_aas]
      
      p3 <- ggplot(sig_data, aes(x = tRNA_supply, y = RSCU)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", se = TRUE) +
        facet_wrap(~ AA, scales = "free") +
        theme_minimal() +
        labs(title = "Significant Correlations: Codon Usage vs tRNA Abundance",
             x = "tRNA Supply (Gene Copy Number)",
             y = "RSCU")
    }
  }
  
  # --- Save plots and results ---
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(output_dir, "tRNA_codon_correlation_overall.pdf"), 
         p1, width = 10, height = 8)
  
  if (!is.null(p2)) {
    ggsave(file.path(output_dir, "tRNA_codon_correlation_by_AA.pdf"), 
           p2, width = 12, height = 6)
  }
  
  if (!is.null(p3)) {
    ggsave(file.path(output_dir, "tRNA_codon_correlation_significant.pdf"), 
           p3, width = 12, height = 8)
  }
  
  # Save correlation results as CSV
  if (exists("overall_cor", correlation_results)) {
    overall_results <- data.frame(
      Test = "Overall",
      Correlation = correlation_results$overall$estimate,
      P_value = correlation_results$overall$p.value,
      Method = test_method
    )
    
    if (nrow(aa_cor_data) > 0) {
      aa_results <- data.frame(
        Test = paste("AA:", aa_cor_data$AA),
        Correlation = aa_cor_data$Correlation,
        P_value = aa_cor_data$P_value,
        Method = test_method
      )
      
      all_results <- rbind(overall_results, aa_results)
    } else {
      all_results <- overall_results
    }
    
    fwrite(all_results, file.path(output_dir, "tRNA_codon_correlations.csv"))
  }
  
  # Save analysis data
  fwrite(analysis_dt, file.path(output_dir, "tRNA_codon_analysis_data.csv"))
  
  message("tRNA-codon correlation analysis complete!")
  message(paste("Results saved to:", output_dir))
  
  # Return results
  result_list <- list(
    correlation_results = correlation_results,
    analysis_data = analysis_dt,
    plots = list(overall = p1)
  )
  
  if (!is.null(p2)) result_list$plots$by_aa <- p2
  if (!is.null(p3)) result_list$plots$significant <- p3
  
  return(result_list)
}