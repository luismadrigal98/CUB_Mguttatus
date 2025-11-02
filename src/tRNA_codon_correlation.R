# Helper Functions for tRNA Expression Analysis

map_tRNA_to_geneIDs <- function(trna_data, gff3_file) {
  #' Map tRNA coordinates to Gene IDs from GFF3 annotation
  #' 
  #' @param trna_data Data table with tRNA coordinates
  #' @param gff3_file Path to GFF3 file
  #' @return Data table with Gene_ID column added
  
  require(data.table)
  
  cat(sprintf("  Loading %d tRNA predictions\n", nrow(trna_data)))
  
  # Read GFF3 - standard 9 columns, skip comment lines
  gff3 <- fread(cmd = paste("grep -v '^#'", gff3_file), 
                sep = "\t", header = FALSE)
  
  # Set standard GFF3 column names
  setnames(gff3, c("Chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
  
  # Filter for gene features only
  gff3_genes <- gff3[type == "gene"]
  
  # Extract gene ID from attributes column
  gff3_genes[, Gene_ID := sub("ID=([^;]+).*", "\\1", attributes)]
  
  # Convert coordinates to numeric
  gff3_genes[, Start := as.numeric(start)]
  gff3_genes[, End := as.numeric(end)]
  
  cat(sprintf("  Loaded %d gene annotations\n", nrow(gff3_genes)))
  
  # Add Gene_ID column
  trna_data[, Gene_ID := NA_character_]
  
  # Find overlaps
  for (i in seq_len(nrow(trna_data))) {
    chr <- trna_data$Sequence_name[i]
    start <- min(trna_data$tRNA_start[i], trna_data$tRNA_end[i])
    end <- max(trna_data$tRNA_start[i], trna_data$tRNA_end[i])
    
    overlaps <- gff3_genes[Chr == chr & Start <= end & End >= start]
    
    if (nrow(overlaps) > 0) {
      trna_data$Gene_ID[i] <- overlaps$Gene_ID[1]
    }
  }
  
  matched <- sum(!is.na(trna_data$Gene_ID))
  cat(sprintf("  Matched %d / %d tRNAs to gene IDs (%.1f%%)\n", 
              matched, nrow(trna_data), 100 * matched / nrow(trna_data)))
  
  return(trna_data)
}

get_tRNA_expression <- function(trna_with_ids, expression_data) {
  #' Get tRNA expression levels from expression data
  #'
  #' @param trna_with_ids Data table with Gene_ID column
  #' @param expression_data Data frame with Gene_name and Expression
  #' @return Data table with Expression column added
  
  # Clean Gene_ID formatting
  trna_with_ids[, Gene_ID_clean := sub("\\.v2\\.1$", "", Gene_ID)]
  expression_dt <- as.data.table(expression_data)
  expression_dt[, Gene_name_clean := sub("\\.v2\\.1$", "", Gene_name)]
  expression_dt[, Gene_name_clean := sub("\\.1$", "", Gene_name_clean)]
  
  # Merge
  trna_expr <- merge(trna_with_ids, expression_dt, 
                     by.x = "Gene_ID_clean", by.y = "Gene_name_clean",
                     all.x = TRUE)
  
  with_expr <- sum(!is.na(trna_expr$Expression))
  cat(sprintf("  Found expression for %d / %d tRNAs (%.1f%%)\n",
              with_expr, nrow(trna_expr), 100 * with_expr / nrow(trna_expr)))
  
  if (with_expr > 0) {
    cat(sprintf("  Expression range: %.2f - %.2f (mean: %.2f)\n",
                min(trna_expr$Expression, na.rm = TRUE),
                max(trna_expr$Expression, na.rm = TRUE),
                mean(trna_expr$Expression, na.rm = TRUE)))
  }
  
  return(trna_expr)
}

calculate_codon_supply_from_expression <- function(trna_expr) {
  #' Calculate codon supply from tRNA expression levels
  #'
  #' @param trna_expr Data table with Anticodon and Expression columns
  #' @return Data table with Codon and tRNA_abundance columns
  
  require(data.table)
  
  # Filter tRNAs with expression data
  trna_with_expr <- trna_expr[!is.na(Expression)]
  
  # Standard wobble rules
  wobble_rules <- list(
    "G" = c("T", "C"),
    "C" = c("G"),
    "A" = c("T"),
    "T" = c("A", "G")
  )
  
  complement <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  
  # Map anticodons to codons with expression weights
  codon_expr_list <- list()
  
  for (i in seq_len(nrow(trna_with_expr))) {
    anticodon <- trna_with_expr$Anticodon[i]
    expr <- trna_with_expr$Expression[i]
    
    ac_1 <- substr(anticodon, 1, 1)
    ac_2 <- substr(anticodon, 2, 2)
    ac_3 <- substr(anticodon, 3, 3)
    
    codon_1 <- complement[ac_3]
    codon_2 <- complement[ac_2]
    codon_3_list <- wobble_rules[[ac_1]]
    
    if (!is.null(codon_3_list)) {
      for (c3 in codon_3_list) {
        codon <- paste0(codon_1, codon_2, c3)
        codon_expr_list[[length(codon_expr_list) + 1]] <- data.table(
          Codon = codon,
          tRNA_expr = expr
        )
      }
    }
  }
  
  codon_expr_dt <- rbindlist(codon_expr_list)
  
  # Sum expression for codons recognized by multiple tRNAs
  codon_supply <- codon_expr_dt[, .(tRNA_abundance = sum(tRNA_expr)), by = Codon]
  
  cat(sprintf("  Calculated expression for %d codons\n", nrow(codon_supply)))
  
  return(codon_supply)
}

# Main Function

tRNA_codon_correlation <- function(codon_counts, tRNA_file, genetic_code,
                                   output_dir = "./results", 
                                   test_method = "spearman",
                                   mode = "by.expression",
                                   ann = NULL,
                                   expression_data = NULL)
{
  #' Analyze correlation between codon usage and tRNA abundance
  #' 
  #' @description Performs statistical tests to examine if codon usage bias
  #' can be explained by tRNA gene abundance. Calculates correlations between
  #' codon usage frequencies and tRNA gene copy numbers (mode="by.copy.number")
  #' or tRNA gene expression levels (mode="by.expression").
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param tRNA_file Path to filtered tRNA annotation file (tab-separated)
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_dir Directory for output files
  #' @param test_method Correlation method: "spearman", "pearson", or "kendall"
  #' @param mode "by.copy.number" (uses tRNA gene counts) or "by.expression" 
  #' (uses tRNA expression levels from RNA-seq)
  #' @param ann Path to GFF3 file (required if mode is "by.expression") for 
  #' mapping tRNA coordinates to gene IDs
  #' @param expression_data Data frame with Gene_name and Expression columns
  #' (required if mode is "by.expression")
  #' 
  #' @return List with correlation results, plots, and statistics
  #' ___________________________________________________________________________
  
  require(data.table)
  require(ggplot2)
  
  # Read tRNA data
  trna_data <- fread(tRNA_file)
  
  # Determine tRNA abundance metric based on mode
  if (mode == "by.expression") {
    cat("\n=== Using tRNA Expression Levels ===\n")
    
    # Check required inputs
    if (is.null(ann) || is.null(expression_data)) {
      stop("Mode 'by.expression' requires both 'ann' (GFF3 path) and 'expression_data'")
    }
    
    # Map tRNA coordinates to Gene IDs from GFF3
    cat("Mapping tRNA coordinates to gene IDs...\n")
    trna_with_ids <- map_tRNA_to_geneIDs(trna_data, ann)
    
    # Get tRNA expression levels
    cat("Extracting tRNA expression levels...\n")
    trna_expr <- get_tRNA_expression(trna_with_ids, expression_data)
    
    # Calculate codon supply based on expression
    codon_supply <- calculate_codon_supply_from_expression(trna_expr)
    
  } else {
    cat("\n=== Using tRNA Gene Copy Numbers ===\n")
    
    # Original behavior: count tRNA genes per anticodon
    trna_counts <- trna_data[, .(tRNA_count = .N), by = Anticodon]
    
    # Calculate codon supply based on wobble rules
    codon_supply <- get_codon_supply_map(trna_counts)
    setnames(codon_supply, "tRNA_supply", "tRNA_abundance")
  }
  
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
  
  # Merge with tRNA abundance data (either from expression or copy number)
  abundance_col <- if (mode == "by.expression") "tRNA_abundance" else "tRNA_abundance"
  merge_cols <- c("Codon", abundance_col)
  analysis_data <- merge(analysis_data, codon_supply[, ..merge_cols], 
                         by = "Codon", all.x = TRUE)
  analysis_data[[abundance_col]][is.na(analysis_data[[abundance_col]])] <- 0
  
  # Rename for consistency
  names(analysis_data)[names(analysis_data) == abundance_col] <- "tRNA_supply"
  
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