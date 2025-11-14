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
  
  # Extended wobble rules including inosine (I) modifications
  # Inosine at wobble position can pair with U, C, or A
  # A at wobble is often modified to I in eukaryotes
  wobble_rules <- list(
    "G" = c("T", "C"),           # G pairs with U/C
    "C" = c("G"),                # C pairs with G only
    "A" = c("T", "C", "A"),      # A (often modified to I) pairs with U/C/A
    "T" = c("A", "G"),           # U pairs with A/G
    "I" = c("T", "C", "A")       # Inosine pairs with U/C/A
  )
  
  complement <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G", "I" = "I")
  
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

#' Calculate tRNA Adaptation Index (tAI) for genes
#'
#' @param codon_counts Data table with gene codon counts
#' @param codon_supply Data table with Codon and tRNA_abundance columns
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Data table with Gene_name and tAI
calculate_tAI <- function(codon_counts, codon_supply, genetic_code) {
  #' tAI = geometric mean of tRNA abundances weighted by codon usage
  #' Based on dos Reis et al. 2004, Nucleic Acids Research
  
  require(data.table)
  
  cat("\n=== Calculating tRNA Adaptation Index (tAI) ===\n")
  
  # Create tRNA weights (relative abundance)
  # Add pseudocount to avoid zeros
  codon_weights <- codon_supply
  codon_weights$tRNA_weight <- codon_weights$tRNA_abundance + 1
  
  # Normalize weights (0-1 scale)
  max_weight <- max(codon_weights$tRNA_weight, na.rm = TRUE)
  codon_weights$tRNA_weight <- codon_weights$tRNA_weight / max_weight
  
  # Create lookup table
  weight_lookup <- setNames(codon_weights$tRNA_weight, codon_weights$Codon)
  
  # Calculate tAI for each gene
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  tAI_results <- data.table(
    Gene_name = codon_counts$Gene_name,
    tAI = numeric(nrow(codon_counts))
  )
  
  cat(sprintf("Calculating tAI for %d genes...\n", nrow(codon_counts)))
  
  for (i in seq_len(nrow(codon_counts))) {
    gene_codons <- as.numeric(codon_counts[i, codon_cols, with = FALSE])
    names(gene_codons) <- codon_cols
    
    # Get sense codons (exclude STOP, Met, Trp)
    sense_codons <- names(genetic_code)[genetic_code != "STOP" & 
                                         genetic_code != "Met" & 
                                         genetic_code != "Trp"]
    gene_codons <- gene_codons[names(gene_codons) %in% sense_codons]
    
    # Get used codons
    used_codons <- names(gene_codons)[gene_codons > 0]
    
    if (length(used_codons) == 0) {
      tAI_results$tAI[i] <- NA
      next
    }
    
    # Get weights for used codons
    weights <- weight_lookup[used_codons]
    counts <- gene_codons[used_codons]
    
    # Handle missing weights
    weights[is.na(weights)] <- 0.01  # Small pseudocount
    
    # tAI = geometric mean weighted by codon counts
    # log(tAI) = sum(n_i * log(w_i)) / sum(n_i)
    log_tAI <- sum(counts * log(weights)) / sum(counts)
    tAI_results$tAI[i] <- exp(log_tAI)
  }
  
  cat(sprintf("tAI range: %.4f - %.4f (mean: %.4f)\n",
              min(tAI_results$tAI, na.rm = TRUE),
              max(tAI_results$tAI, na.rm = TRUE),
              mean(tAI_results$tAI, na.rm = TRUE)))
  
  return(tAI_results)
}


#' Analyze tAI vs gene expression correlation
#'
#' @param tAI_results Data table from calculate_tAI()
#' @param expression_data Data frame with Gene_name and Expression
#' @param output_dir Directory for output files
#' @return List with correlation results and plot
analyze_tAI_expression <- function(tAI_results, expression_data, output_dir = "./results") {
  
  require(ggplot2)
  require(data.table)
  
  cat("\n=== Analyzing tAI vs Gene Expression ===\n")
  
  # Merge tAI with expression
  analysis_data <- merge(tAI_results, expression_data, by = "Gene_name")
  analysis_data <- analysis_data[!is.na(tAI) & !is.na(Expression)]
  
  cat(sprintf("Analyzing %d genes with both tAI and expression\n", nrow(analysis_data)))
  
  # Log-transform expression for better correlation
  analysis_data$log_Expression <- log10(analysis_data$Expression + 1)
  
  # Correlations
  cor_pearson <- cor.test(analysis_data$tAI, analysis_data$log_Expression, method = "pearson")
  cor_spearman <- cor.test(analysis_data$tAI, analysis_data$log_Expression, method = "spearman")
  
  cat(sprintf("Pearson r = %.4f (p = %.2e)\n", cor_pearson$estimate, cor_pearson$p.value))
  cat(sprintf("Spearman ρ = %.4f (p = %.2e)\n", cor_spearman$estimate, cor_spearman$p.value))
  
  # Create plot
  p <- ggplot(analysis_data, aes(x = tAI, y = log_Expression)) +
    geom_hex(bins = 50) +
    geom_smooth(method = "lm", color = "red", linewidth = 1.5, se = TRUE) +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = "tRNA Adaptation Index vs Gene Expression",
      subtitle = sprintf("Pearson r = %.3f (p = %.2e), Spearman ρ = %.3f (p = %.2e)\n%d genes",
                        cor_pearson$estimate, cor_pearson$p.value,
                        cor_spearman$estimate, cor_spearman$p.value,
                        nrow(analysis_data)),
      x = "tRNA Adaptation Index (tAI)",
      y = "log10(Expression + 1)",
      fill = "Gene\nDensity"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  # Save plot
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, "tAI_vs_expression.pdf"), p, width = 10, height = 8)
  
  # Save data
  fwrite(analysis_data, file.path(output_dir, "tAI_expression_data.csv"))
  
  cat(sprintf("Results saved to %s\n", output_dir))
  
  return(list(
    data = analysis_data,
    pearson = cor_pearson,
    spearman = cor_spearman,
    plot = p
  ))
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
    Codon_count = as.numeric(total_codon_usage),
    AA = genetic_code[names(total_codon_usage)],
    stringsAsFactors = FALSE
  )
  
  # Remove STOP, Met, and Trp
  analysis_data <- analysis_data[analysis_data$AA != "STOP" &
                                   analysis_data$AA != "Trp" &
                                   analysis_data$AA != "Met", ]
  
  # Calculate simple genome-wide frequency
  total_codons <- sum(analysis_data$Codon_count)
  analysis_data$Codon_frequency <- analysis_data$Codon_count / total_codons
  
  # Merge with tRNA abundance data (either from expression or copy number)
  abundance_col <- if (mode == "by.expression") "tRNA_abundance" else "tRNA_abundance"
  merge_cols <- c("Codon", abundance_col)
  analysis_data <- merge(analysis_data, codon_supply[, ..merge_cols], 
                         by = "Codon", all.x = TRUE)
  analysis_data[[abundance_col]][is.na(analysis_data[[abundance_col]])] <- 0
  
  # Rename for consistency
  names(analysis_data)[names(analysis_data) == abundance_col] <- "tRNA_supply"
  
  # Convert to data.table
  analysis_dt <- as.data.table(analysis_data)
  
  # Add some useful metrics
  analysis_dt[, N_synonyms := .N, by = AA]
  
  # For within-family comparisons, also calculate proportions within each AA
  analysis_dt[, Total_AA_count := sum(Codon_count), by = AA]
  analysis_dt[, Codon_proportion_in_AA := Codon_count / Total_AA_count]
  analysis_dt[, Total_tRNA_per_AA := sum(tRNA_supply), by = AA]
  analysis_dt[, tRNA_proportion_in_AA := ifelse(Total_tRNA_per_AA > 0, 
                                                 tRNA_supply / Total_tRNA_per_AA, 
                                                 0)]
  
  # Perform correlation tests - SIMPLE APPROACH
  correlation_results <- list()
  
  # Overall correlation: genome-wide codon frequency vs tRNA supply
  if (sum(analysis_dt$tRNA_supply > 0) > 2) {
    # Main test: simple codon frequency vs tRNA supply
    overall_cor <- cor.test(analysis_dt$Codon_frequency, 
                           analysis_dt$tRNA_supply, 
                           method = test_method, exact = FALSE)
    correlation_results$overall <- overall_cor
  }
  
  # Per amino acid correlations
  # Within each AA family: codon proportion vs tRNA proportion
  aa_correlations <- list()
  
  for (aa in unique(analysis_dt$AA)) {
    aa_data <- analysis_dt[AA == aa]
    
    # Skip if not enough codons
    if (nrow(aa_data) < 3) next  # Need at least 3 for meaningful correlation
    
    # Check for actual variation (handle edge cases)
    codon_var <- var(aa_data$Codon_proportion_in_AA, na.rm = TRUE)
    trna_var <- var(aa_data$tRNA_proportion_in_AA, na.rm = TRUE)
    
    # Skip if no variation or NA variance
    if (is.na(codon_var) || is.na(trna_var) || 
        codon_var < 1e-10 || trna_var < 1e-10) next
    
    # Test within-family: does codon proportion match tRNA proportion?
    tryCatch({
      aa_cor <- cor.test(aa_data$Codon_proportion_in_AA, 
                        aa_data$tRNA_proportion_in_AA,
                        method = test_method, exact = FALSE)
      # Only store if we got a valid result
      if (!is.na(aa_cor$estimate) && !is.na(aa_cor$p.value)) {
        aa_correlations[[aa]] <- aa_cor
      }
    }, error = function(e) {
      # Skip this AA if correlation fails
    })
  }
  
  correlation_results$per_amino_acid <- aa_correlations
  
  # --- Create Visualizations ---
  
  # 1. Overall scatter plot: genome-wide codon frequency vs tRNA supply
  p1 <- ggplot(analysis_dt, aes(x = tRNA_supply, y = Codon_frequency, color = AA)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(title = "Genome-wide Codon Frequency vs tRNA Supply",
         subtitle = sprintf("%s r = %.3f (p = %.2e) | %d codons", 
                          tools::toTitleCase(test_method),
                          ifelse(!is.null(correlation_results$overall), 
                                 correlation_results$overall$estimate, NA),
                          ifelse(!is.null(correlation_results$overall), 
                                 correlation_results$overall$p.value, NA),
                          nrow(analysis_dt)),
         x = "tRNA Supply (gene count or expression)", 
         y = "Codon Frequency (proportion of all codons)",
         color = "Amino Acid") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 2. Per amino acid correlations bar plot (using proportion-based)
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
    
    # Add info about number of codons per AA
    aa_info <- analysis_dt[, .(N_codons = .N), by = AA]
    aa_cor_data <- merge(aa_cor_data, aa_info, by = "AA")
  }
  
  p2 <- NULL # Initialize as NULL
  if (nrow(aa_cor_data) > 0) {
    p2 <- ggplot(aa_cor_data, aes(x = reorder(AA, Correlation), y = Correlation, 
                                  fill = Significant)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("FALSE" = "#CCCCCC", "TRUE" = "#E74C3C"),
                       name = "p < 0.05") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      labs(title = "Codon-tRNA Correlation by Amino Acid",
           subtitle = sprintf("Using proportion-based measures | Method: %s", test_method),
           x = "Amino Acid (sorted by correlation strength)", 
           y = paste(tools::toTitleCase(test_method), "Correlation (r)"),
           caption = sprintf("%d / %d families show significant correlation",
                           sum(aa_cor_data$Significant), nrow(aa_cor_data)))
  }
  
  # 3. Faceted scatter plots for significant correlations
  p3 <- NULL # Initialize as NULL
  if (nrow(aa_cor_data) > 0) {
    significant_aas <- aa_cor_data$AA[aa_cor_data$Significant]
    
    if (length(significant_aas) > 0 && !any(is.na(significant_aas))) {
      sig_data <- analysis_dt[AA %in% significant_aas]
      
      if (nrow(sig_data) > 0) {
        # Add correlation info to facet labels
        facet_labels <- sapply(significant_aas, function(aa) {
          cor_val <- aa_cor_data$Correlation[aa_cor_data$AA == aa]
          p_val <- aa_cor_data$P_value[aa_cor_data$AA == aa]
          sprintf("%s (r=%.2f, p=%.3f)", aa, cor_val, p_val)
        })
        names(facet_labels) <- significant_aas
        
        p3 <- ggplot(sig_data, aes(x = tRNA_proportion_in_AA, y = Codon_proportion_in_AA)) +
          geom_point(size = 3, alpha = 0.7, color = "#3498DB") +
          geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#E74C3C") +
          geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
          facet_wrap(~ AA, scales = "free", labeller = labeller(AA = facet_labels)) +
          theme_minimal(base_size = 11) +
          theme(strip.text = element_text(face = "bold", size = 9),
                plot.title = element_text(face = "bold", hjust = 0.5)) +
          labs(title = "Within-Family Correlations: Codon vs tRNA Proportions",
               subtitle = sprintf("%d amino acid families with p < 0.05", length(significant_aas)),
               x = "tRNA Proportion (within family)",
               y = "Codon Proportion (within family)",
               caption = "Dotted line = perfect agreement")
      }
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
  all_results <- list()
  
  if (!is.null(correlation_results$overall)) {
    overall <- data.frame(
      Test = "Overall",
      Correlation = correlation_results$overall$estimate,
      P_value = correlation_results$overall$p.value,
      Method = test_method,
      Measure = "Genome-wide frequency vs tRNA supply",
      N_codons = sum(analysis_dt$Codon_count > 0)
    )
    all_results[[1]] <- overall
  }
  
  # Check if per_amino_acid exists and is a list with results
  if (!is.null(correlation_results$per_amino_acid) && 
      length(correlation_results$per_amino_acid) > 0) {
    # Convert list to data frame
    aa_results <- data.frame(
      Test = paste0("AA_", names(correlation_results$per_amino_acid)),
      Correlation = sapply(correlation_results$per_amino_acid, function(x) x$estimate),
      P_value = sapply(correlation_results$per_amino_acid, function(x) x$p.value),
      Method = test_method,
      Measure = "Within-family proportion",
      N_codons = sapply(names(correlation_results$per_amino_acid), function(aa) {
        sum(analysis_dt$AA == aa)
      })
    )
    all_results[[2]] <- aa_results
  }
  
  if (length(all_results) > 0) {
    all_results_df <- do.call(rbind, all_results)
    fwrite(all_results_df, file.path(output_dir, "tRNA_codon_correlations.csv"))
  }
  
  # Save analysis data
  fwrite(analysis_dt, file.path(output_dir, "tRNA_codon_analysis_data.csv"))
  
  # Calculate tAI if expression data provided
  tAI_results <- NULL
  if (!is.null(expression_data) && mode == "by.expression") {
    cat("\n")
    tAI_results <- calculate_tAI(codon_counts, codon_supply, genetic_code)
    
    # Analyze tAI vs expression
    tAI_analysis <- analyze_tAI_expression(tAI_results, expression_data, output_dir)
  }
  
  # Print summary statistics
  cat("\n=== Correlation Summary ===\n")
  if (!is.null(correlation_results$overall)) {
    cat(sprintf("Overall genome-wide frequency: r = %.3f, p = %.2e\n",
                correlation_results$overall$estimate,
                correlation_results$overall$p.value))
  }
  
  # Extract significant per-AA correlations for summary
  significant_aa <- NULL
  if (!is.null(correlation_results$per_amino_acid) && length(correlation_results$per_amino_acid) > 0) {
    # per_amino_acid is a list, convert to data frame
    aa_summary <- data.frame(
      amino_acid = names(correlation_results$per_amino_acid),
      correlation = sapply(correlation_results$per_amino_acid, function(x) x$estimate),
      p_value = sapply(correlation_results$per_amino_acid, function(x) x$p.value),
      n_codons = sapply(names(correlation_results$per_amino_acid), function(aa) {
        sum(analysis_dt$AA == aa)
      })
    )
    
    significant_aa <- aa_summary[aa_summary$p_value < 0.05, ]
    if (nrow(significant_aa) > 0) {
      cat(sprintf("\n%d / %d amino acid families show significant within-family correlation (p < 0.05):\n", 
                  nrow(significant_aa), nrow(aa_summary)))
      significant_aa <- significant_aa[order(-abs(significant_aa$correlation)), ]
      for (i in seq_len(min(10, nrow(significant_aa)))) {
        cat(sprintf("  %s (%d codons): r = %.3f (p = %.2e)\n", 
                    significant_aa$amino_acid[i],
                    significant_aa$n_codons[i],
                    significant_aa$correlation[i],
                    significant_aa$p_value[i]))
      }
      if (nrow(significant_aa) > 10) {
        cat(sprintf("  ... and %d more\n", nrow(significant_aa) - 10))
      }
    } else {
      cat("\nNo amino acid families show significant within-family correlation\n")
    }
  }
  
  message("\ntRNA-codon correlation analysis complete!")
  message(paste("Results saved to:", output_dir))
  
  # Return results
  result_list <- list(
    correlation_results = correlation_results,
    analysis_data = analysis_dt,
    plots = list(overall = p1),
    significant_amino_acids = significant_aa,
    codon_supply = codon_supply
  )
  
  if (!is.null(p2)) result_list$plots$by_aa <- p2
  if (!is.null(p3)) result_list$plots$significant <- p3
  if (!is.null(tAI_results)) {
    result_list$tAI_results = tAI_results
    result_list$tAI_analysis = tAI_analysis
  }
  
  return(result_list)
}