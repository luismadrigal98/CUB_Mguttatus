##' Diagnose discrepancies between CAI w-values and proportion tests
##' 
##' @description Helps understand cases where a codon has w=1.0 (CAI preferred)
##' but negative proportion difference (avoided in highly expressed genes).
##' This happens when ALL codons for an amino acid shift in the same direction.
##' 
##' @param w_table CAI weight table with relative_adaptiveness values
##' @param test_results Output from test_codon_proportions()
##' @param codon_usage Codon usage data
##' @param expression_groups Data frame with Gene_name and Expression_Group
##' @param genetic_code Named vector mapping codons to amino acids
##' 
##' @return List with diagnostic information
##' 
##' @author Luis J. Madrigal-Roca
##' @date November 12, 2025

diagnose_cai_vs_proportion <- function(w_table, test_results, 
                                       codon_usage, expression_groups,
                                       genetic_code) {
  
  cat("\n=== Diagnosing CAI vs Proportion Discrepancies ===\n\n")
  
  # Merge data
  comparison <- test_results %>%
    left_join(w_table %>% select(codon, relative_adaptiveness), 
              by = c("Codon" = "codon")) %>%
    mutate(
      Is_Preferred = relative_adaptiveness == 1.0,
      Is_Avoided = Significant & Difference < 0,
      Discrepancy = Is_Preferred & Is_Avoided
    )
  
  # Find discrepancies
  discrepancies <- comparison %>% filter(Discrepancy)
  
  if (nrow(discrepancies) == 0) {
    cat("✓ No discrepancies found - all preferred codons are enriched in high expression\n")
    return(list(discrepancies = NULL, amino_acid_shifts = NULL))
  }
  
  cat(sprintf("⚠ Found %d discrepant codons (w=1.0 but avoided in high expression):\n\n", 
              nrow(discrepancies)))
  
  # Print discrepancies
  for (i in 1:nrow(discrepancies)) {
    row <- discrepancies[i, ]
    cat(sprintf("  %s (%s): w=%.3f, Diff=%.4f, p_adj=%.4e\n",
                row$Codon, row$Amino_Acid, 
                row$relative_adaptiveness, row$Difference, row$p_adj))
  }
  
  cat("\n--- Detailed Analysis ---\n\n")
  
  # Analyze each discrepant amino acid
  discrepant_aas <- unique(discrepancies$Amino_Acid)
  aa_analysis <- list()
  
  for (aa in discrepant_aas) {
    cat(sprintf("\n=== %s (Arginine) ===\n", aa))
    
    # Get all codons for this AA
    aa_codons <- names(genetic_code)[genetic_code == aa]
    aa_codons <- aa_codons[aa_codons != "STOP"]
    
    # Get their data
    aa_data <- comparison %>% filter(Amino_Acid == aa)
    
    # Calculate absolute frequencies (not just proportions within AA)
    codon_usage_merged <- codon_usage %>%
      left_join(expression_groups %>% select(Gene_name, Expression_Group),
                by = "Gene_name")
    
    top5 <- codon_usage_merged %>% filter(Expression_Group == "Top 5%")
    rest <- codon_usage_merged %>% filter(Expression_Group != "Top 5%")
    
    cat("\nAbsolute codon counts:\n")
    cat(sprintf("%-6s %10s %10s %12s %12s %10s %8s\n",
                "Codon", "Top5_cnt", "Rest_cnt", "Top5_prop", "Rest_prop", "w_value", "Status"))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    # Get total codon counts (all amino acids)
    top5_total <- sum(sapply(aa_codons, function(c) {
      if (c %in% names(top5)) sum(top5[[c]], na.rm = TRUE) else 0
    }))
    rest_total <- sum(sapply(aa_codons, function(c) {
      if (c %in% names(rest)) sum(rest[[c]], na.rm = TRUE) else 0
    }))
    
    for (codon in aa_codons) {
      if (!codon %in% aa_data$Codon) next
      
      codon_row <- aa_data %>% filter(Codon == codon)
      
      top5_cnt <- if (codon %in% names(top5)) sum(top5[[codon]], na.rm = TRUE) else 0
      rest_cnt <- if (codon %in% names(rest)) sum(rest[[codon]], na.rm = TRUE) else 0
      
      top5_prop <- codon_row$Selected_Prop
      rest_prop <- codon_row$Neutral_Prop
      w_val <- codon_row$relative_adaptiveness
      
      status <- ""
      if (w_val == 1.0) status <- "PREFERRED"
      if (codon_row$Significant & codon_row$Difference < 0) status <- paste(status, "AVOIDED")
      
      cat(sprintf("%-6s %10d %10d %12.4f %12.4f %10.3f %8s\n",
                  codon, top5_cnt, rest_cnt, top5_prop, rest_prop, w_val, status))
    }
    
    cat("\n** Interpretation **\n")
    cat("This amino acid shows a general shift in codon usage between groups.\n")
    cat("CAI w=1.0 means 'most common within this AA in Top 5%'\n")
    cat("But if ALL codons for this AA decline in Top 5%, the 'preferred' one\n")
    cat("can still have LOWER absolute proportion than in rest.\n\n")
    
    cat("Biological meaning:\n")
    cat("- This amino acid is used LESS in highly expressed genes overall\n")
    cat("- Within this AA, certain codons are relatively preferred\n")
    cat("- But the 'preference' is only relative, not absolute enrichment\n\n")
    
    aa_analysis[[aa]] <- aa_data
  }
  
  # Summary recommendations
  cat("\n=== Recommendations ===\n\n")
  cat("1. CAI identifies codons that are RELATIVELY preferred WITHIN each amino acid\n")
  cat("2. Statistical tests identify codons ABSOLUTELY enriched in high expression\n")
  cat("3. For selection analysis, focus on codons that are BOTH:\n")
  cat("   - Preferred (w=1.0) AND\n")
  cat("   - Significantly enriched (positive difference, p_adj < 0.05)\n\n")
  
  cat("4. Codons with w=1.0 but negative difference indicate:\n")
  cat("   - The amino acid overall is avoided in high expression\n")
  cat("   - When this AA is used, this codon is preferred\n")
  cat("   - But overall usage is still lower than in low expression genes\n\n")
  
  # Calculate better metric: absolute enrichment of preferred codons
  preferred_codons <- comparison %>% filter(Is_Preferred)
  truly_selected <- preferred_codons %>% 
    filter(Significant, Difference > 0)
  
  cat(sprintf("Summary:\n"))
  cat(sprintf("- Total preferred codons (w=1.0): %d\n", nrow(preferred_codons)))
  cat(sprintf("- Preferred + significantly enriched: %d (%.1f%%)\n", 
              nrow(truly_selected),
              100 * nrow(truly_selected) / nrow(preferred_codons)))
  cat(sprintf("- Preferred but avoided: %d (%.1f%%)\n",
              nrow(discrepancies),
              100 * nrow(discrepancies) / nrow(preferred_codons)))
  cat(sprintf("- Preferred but neutral: %d (%.1f%%)\n",
              sum(preferred_codons$Is_Preferred & !preferred_codons$Significant),
              100 * sum(preferred_codons$Is_Preferred & !preferred_codons$Significant) / 
                nrow(preferred_codons)))
  
  return(list(
    discrepancies = discrepancies,
    amino_acid_analysis = aa_analysis,
    summary = data.frame(
      Total_preferred = nrow(preferred_codons),
      Truly_selected = nrow(truly_selected),
      Avoided = nrow(discrepancies),
      Neutral = sum(preferred_codons$Is_Preferred & !preferred_codons$Significant)
    )
  ))
}


##' Create corrected classification combining CAI and proportion tests
##' 
##' @description Creates a corrected classification that properly handles
##' the difference between relative (CAI) and absolute (proportion) measures
##' 
##' @param w_table CAI weight table
##' @param test_results Proportion test results
##' 
##' @return Data frame with corrected classifications

create_corrected_classification <- function(w_table, test_results) {
  
  cat("\n=== Creating Corrected Codon Classification ===\n\n")
  
  classification <- test_results %>%
    left_join(w_table %>% select(codon, relative_adaptiveness, amino_acid),
              by = c("Codon" = "codon")) %>%
    mutate(
      CAI_Status = case_when(
        relative_adaptiveness == 1.0 ~ "CAI Preferred (w=1)",
        relative_adaptiveness >= 0.5 ~ "CAI Intermediate",
        TRUE ~ "CAI Rare"
      ),
      Proportion_Status = case_when(
        !Significant ~ "Neutral (no sig diff)",
        Difference > 0 ~ "Enriched in High Expr",
        Difference < 0 ~ "Avoided in High Expr",
        TRUE ~ "Neutral"
      ),
      Combined_Classification = case_when(
        # True selection: preferred AND enriched
        relative_adaptiveness == 1.0 & Significant & Difference > 0 ~ 
          "Under Selection (pref + enriched)",
        
        # Relative preference but not absolutely enriched
        relative_adaptiveness == 1.0 & Significant & Difference < 0 ~ 
          "Relatively Preferred (but AA avoided)",
        
        relative_adaptiveness == 1.0 & !Significant ~ 
          "Relatively Preferred (neutral)",
        
        # Enriched but not CAI preferred
        relative_adaptiveness < 1.0 & Significant & Difference > 0 ~ 
          "Enriched (not CAI pref)",
        
        # Avoided
        relative_adaptiveness < 1.0 & Significant & Difference < 0 ~ 
          "Avoided in High Expr",
        
        # Neutral
        TRUE ~ "Neutral"
      )
    )
  
  # Summary table
  cat("Classification summary:\n")
  summary_table <- classification %>%
    group_by(Combined_Classification) %>%
    summarise(
      Count = n(),
      Mean_w = mean(relative_adaptiveness, na.rm = TRUE),
      Mean_Diff = mean(Difference, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Count))
  
  print(summary_table)
  
  cat("\n** Key Insights **\n")
  cat("- 'Under Selection': codons truly under selection (w=1 AND enriched)\n")
  cat("- 'Relatively Preferred (but AA avoided)': CAI artifact (w=1 but AA declines)\n")
  cat("- 'Relatively Preferred (neutral)': stable within-AA preference\n")
  cat("- 'Enriched (not CAI pref)': compensatory changes within AA\n\n")
  
  return(classification)
}
