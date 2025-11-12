##' Test significance of codon usage differences between groups
##' 
##' @description Tests whether codon usage proportions differ significantly
##' between selected (highly expressed) and neutral (rest) genes using
##' chi-squared or Fisher's exact test for each codon.
##' 
##' @param selected_usage Codon usage data for selected genes (data.table)
##' @param neutral_usage Codon usage data for neutral genes (data.table)
##' @param genetic_code Named vector mapping codons to amino acids
##' @param method Test method: "chisq" or "fisher" (default: "chisq")
##' @param fdr_correction Apply FDR correction? (default: TRUE)
##' 
##' @return Data table with test results per codon
##' 
##' @author Luis J. Madrigal-Roca
##' @date November 12, 2025

test_codon_proportions <- function(selected_usage, neutral_usage, genetic_code,
                                   method = "chisq", fdr_correction = TRUE) {
  
  cat("\n=== Testing Codon Usage Differences ===\n")
  cat(sprintf("Selected genes: %d\n", nrow(selected_usage)))
  cat(sprintf("Neutral genes: %d\n", nrow(neutral_usage)))
  cat(sprintf("Test method: %s\n", method))
  
  # Get codon columns (exclude Gene_name)
  codon_cols <- setdiff(names(selected_usage), "Gene_name")
  codon_cols <- codon_cols[codon_cols %in% names(genetic_code)]
  
  cat(sprintf("Testing %d codons...\n\n", length(codon_cols)))
  
  # Initialize results
  results <- data.table(
    Codon = character(),
    Amino_Acid = character(),
    Selected_Count = numeric(),
    Neutral_Count = numeric(),
    Selected_Prop = numeric(),
    Neutral_Prop = numeric(),
    Difference = numeric(),
    Test_Statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Test each codon
  for (codon in codon_cols) {
    
    # Get counts
    selected_count <- sum(selected_usage[[codon]], na.rm = TRUE)
    neutral_count <- sum(neutral_usage[[codon]], na.rm = TRUE)
    
    # Get total counts for this amino acid
    aa <- genetic_code[codon]
    aa_codons <- names(genetic_code)[genetic_code == aa]
    aa_codons <- intersect(aa_codons, codon_cols)
    
    selected_aa_total <- sum(sapply(aa_codons, function(c) sum(selected_usage[[c]], na.rm = TRUE)))
    neutral_aa_total <- sum(sapply(aa_codons, function(c) sum(neutral_usage[[c]], na.rm = TRUE)))
    
    # Calculate proportions
    selected_prop <- selected_count / selected_aa_total
    neutral_prop <- neutral_count / neutral_aa_total
    difference <- selected_prop - neutral_prop
    
    # Create contingency table
    # Rows: codon vs other codons for this AA
    # Cols: selected vs neutral
    cont_table <- matrix(
      c(selected_count, selected_aa_total - selected_count,
        neutral_count, neutral_aa_total - neutral_count),
      nrow = 2, byrow = FALSE
    )
    
    # Perform test
    if (method == "fisher") {
      test_result <- fisher.test(cont_table)
      test_stat <- NA  # Fisher's test doesn't have a simple test statistic
      p_value <- test_result$p.value
    } else {
      # Chi-squared test
      if (any(cont_table < 5)) {
        # Use Fisher's exact for small counts
        test_result <- fisher.test(cont_table)
        test_stat <- NA
        p_value <- test_result$p.value
      } else {
        test_result <- chisq.test(cont_table, correct = TRUE)
        test_stat <- test_result$statistic
        p_value <- test_result$p.value
      }
    }
    
    # Add to results
    results <- rbind(results, data.table(
      Codon = codon,
      Amino_Acid = aa,
      Selected_Count = selected_count,
      Neutral_Count = neutral_count,
      Selected_Prop = selected_prop,
      Neutral_Prop = neutral_prop,
      Difference = difference,
      Test_Statistic = test_stat,
      p_value = p_value
    ))
  }
  
  # Apply FDR correction
  if (fdr_correction) {
    results$p_adj <- p.adjust(results$p_value, method = "BH")
    results$Significant <- results$p_adj < 0.05
  } else {
    results$p_adj <- results$p_value
    results$Significant <- results$p_value < 0.05
  }
  
  # Classify codons
  results$Classification <- with(results, case_when(
    !Significant ~ "Neutral",
    Difference > 0 ~ "Under Selection (Preferred in High Expr)",
    Difference < 0 ~ "Under Selection (Avoided in High Expr)",
    TRUE ~ "Neutral"
  ))
  
  # Add AT/GC ending classification
  results$Ending <- ifelse(substr(results$Codon, 3, 3) %in% c("A", "T"), "AT", "GC")
  
  # Summary
  cat("\n=== Test Results Summary ===\n")
  cat(sprintf("Total codons tested: %d\n", nrow(results)))
  cat(sprintf("Significant (FDR < 0.05): %d (%.1f%%)\n", 
              sum(results$Significant), 
              100 * sum(results$Significant) / nrow(results)))
  
  sig_results <- results[results$Significant, ]
  if (nrow(sig_results) > 0) {
    cat(sprintf("  - Preferred in high expression: %d\n", 
                sum(sig_results$Difference > 0)))
    cat(sprintf("  - Avoided in high expression: %d\n", 
                sum(sig_results$Difference < 0)))
  }
  
  cat("\n=== Top 10 Most Significant Codons ===\n")
  top10 <- results[order(results$p_value), ][1:min(10, nrow(results)), ]
  print(top10[, .(Codon, Amino_Acid, Selected_Prop, Neutral_Prop, 
                  Difference, p_value, p_adj, Classification)])
  
  return(results)
}


##' Create summary plot of codon selection classification
##' 
##' @param test_results Output from test_codon_proportions()
##' @param output_file Path to save plot
##' 
##' @return ggplot object

plot_codon_selection_summary <- function(test_results, output_file = NULL) {
  
  require(ggplot2)
  
  # Order by difference
  test_results$Codon <- reorder(test_results$Codon, test_results$Difference)
  
  # Create color scheme
  colors <- c(
    "Under Selection (Preferred in High Expr)" = "#E41A1C",
    "Under Selection (Avoided in High Expr)" = "#377EB8", 
    "Neutral" = "gray70"
  )
  
  p <- ggplot(test_results, aes(x = Codon, y = Difference, 
                                 fill = Classification,
                                 alpha = Significant)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = colors, name = "") +
    scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.3),
                      name = "FDR < 0.05") +
    facet_wrap(~ Amino_Acid, scales = "free_x", ncol = 4) +
    labs(
      title = "Codon Usage Difference: Top 5% vs Rest",
      subtitle = "Positive = preferred in highly expressed genes",
      x = "Codon",
      y = "Proportion Difference (Selected - Neutral)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 14, height = 10)
    cat(sprintf("\n✓ Plot saved: %s\n", output_file))
  }
  
  return(p)
}


##' Create heatmap showing codon classification
##'
##' @param test_results Output from test_codon_proportions()
##' @param w_table CAI weight table with preferred codons
##' @param output_file Path to save plot
##'
##' @return ggplot object

plot_codon_classification_heatmap <- function(test_results, w_table, 
                                              output_file = NULL) {
  
  require(ggplot2)
  require(dplyr)
  
  # Merge with w values
  plot_data <- test_results %>%
    left_join(w_table %>% dplyr::select(codon, relative_adaptiveness), 
              by = c("Codon" = "codon")) %>%
    dplyr::mutate(
      Preferred = relative_adaptiveness == 1.0,
      Selection_Status = case_when(
        !Significant ~ "Neutral",
        Difference > 0 & Preferred ~ "Selection + Preferred",
        Difference > 0 & !Preferred ~ "Selection (non-pref)",
        Difference < 0 ~ "Avoided",
        TRUE ~ "Neutral"
      )
    )
  
  # Create tile plot
  p <- ggplot(plot_data, aes(x = Ending, y = Codon, 
                              fill = Selection_Status,
                              alpha = abs(Difference))) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Selected_Prop)), 
              size = 2.5, color = "black") +
    scale_fill_manual(
      values = c(
        "Selection + Preferred" = "#d73027",
        "Selection (non-pref)" = "#fc8d59",
        "Avoided" = "#4575b4",
        "Neutral" = "gray85"
      ),
      name = "Classification"
    ) +
    scale_alpha_continuous(range = c(0.3, 1.0), name = "|Difference|") +
    facet_grid(Amino_Acid ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "Codon Classification: Selection, Preference, and Ending",
      subtitle = "Numbers = proportion in highly expressed genes",
      x = "Codon Ending (AT vs GC)",
      y = "Codon"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 16)
    cat(sprintf("\n✓ Heatmap saved: %s\n", output_file))
  }
  
  return(p)
}