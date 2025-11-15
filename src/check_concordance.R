check_preferred_codon_concordance <- function(model_results, 
                                            existing_preferred_codons,
                                            alpha_significance = 0.05) {
  #' Compares model-driven preferred codons with existing preferred codons
  #' and provides recommendations for updates
  #'
  #' @param model_results List of model results from fit_pairwise_gams/glms
  #' @param existing_preferred_codons Data frame with Codon and AA/Amino_Acid columns
  #' @param alpha_significance Significance threshold for selection slopes
  #' @return Data frame with concordance analysis
  
  suppressPackageStartupMessages({
    require(dplyr)
    require(data.table)
  })
  
  # Extract all coefficients
  all_coefs <- data.table::rbindlist(
    lapply(model_results, function(x) {
      if (is.null(x)) return(NULL)
      x$coefficients
    })
  )
  
  # For each family, determine the preferred codon(s) based on significance patterns
  # This follows the same logic as update_preferred_codons_from_models()
  model_preferred <- all_coefs %>%
    dplyr::group_by(Family) %>%
    dplyr::mutate(
      # Count significant patterns
      n_sig_positive = sum(Significant & Selection_Slope > 0, na.rm = TRUE),
      n_sig_negative = sum(Significant & Selection_Slope < 0, na.rm = TRUE),
      
      # Determine preference pattern (simplified logic)
      Preference_Pattern = dplyr::case_when(
        # Case 1: At least one positive slope = non-baseline codon beats baseline
        n_sig_positive >= 1 ~ "NonBaseline_Preferred",
        
        # Case 2: No positive slopes AND at least one significantly negative = baseline is best
        n_sig_positive == 0 & n_sig_negative >= 1 ~ "Baseline_Preferred",
        
        # Case 3: No significant differences = all codons equivalent
        n_sig_positive == 0 & n_sig_negative == 0 ~ "No_Preference",
        
        # Default
        TRUE ~ "Review_Needed"
      )
    ) %>%
    dplyr::ungroup()
  
  # For each family, select the preferred codon based on pattern
  model_preferred_summary <- model_preferred %>%
    dplyr::group_by(Family, Baseline, Preference_Pattern) %>%
    dplyr::summarise(
      Model_Preferred = dplyr::case_when(
        # If non-baseline preferred, pick the one with highest positive slope
        unique(Preference_Pattern) == "NonBaseline_Preferred" ~ 
          Codon[which.max(Selection_Slope)],
        
        # If baseline preferred, use baseline
        unique(Preference_Pattern) == "Baseline_Preferred" ~ 
          unique(Baseline),
        
        # For no preference, indicate none
        unique(Preference_Pattern) == "No_Preference" ~ 
          "None",
        
        # Otherwise needs review
        TRUE ~ "REVIEW_NEEDED"
      ),
      Model_Slope = ifelse(
        unique(Preference_Pattern) == "Baseline_Preferred",
        0,  # Baseline has slope of 0 by definition
        Selection_Slope[which.max(abs(Selection_Slope))]
      ),
      Model_p_adj = min(p_adj, na.rm = TRUE),
      Model_Significant = any(Significant, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(
      Family, 
      Model_Preferred,
      Model_Slope,
      Model_p_adj,
      Model_Significant,
      Baseline,
      Preference_Pattern
    )
  
  # Prepare existing preferred codons
  if ("AA" %in% names(existing_preferred_codons)) {
    existing_prep <- existing_preferred_codons %>%
      dplyr::select(Family = AA, Existing_Preferred = Codon)
  } else if ("Amino_Acid" %in% names(existing_preferred_codons)) {
    existing_prep <- existing_preferred_codons %>%
      dplyr::select(Family = Amino_Acid, Existing_Preferred = Codon)
  } else {
    stop("existing_preferred_codons must have either 'AA' or 'Amino_Acid' column")
  }
  
  # Merge
  concordance <- model_preferred_summary %>%
    dplyr::left_join(existing_prep, by = "Family")
  
  # Check concordance
  concordance <- concordance %>%
    dplyr::mutate(
      # Check if existing preferred is in the model preferred (handles multiple)
      Concordant = dplyr::case_when(
        grepl("/", Model_Preferred) ~ 
          grepl(Existing_Preferred, Model_Preferred, fixed = TRUE),
        TRUE ~ (Model_Preferred == Existing_Preferred)
      ),
      Concordant_with_Baseline = (Existing_Preferred == Baseline),
      
      Recommendation = dplyr::case_when(
        # No clear preference
        Preference_Pattern == "No_Preference" ~ 
          "Keep - No significant selection detected",
        
        # Pattern needs manual review
        Preference_Pattern == "Review_Needed" | Model_Preferred == "REVIEW_NEEDED" ~ 
          "Review manually - Unclear pattern",
        
        # Perfect agreement: existing preferred matches model
        Concordant & Preference_Pattern == "Baseline_Preferred" ~ 
          "Keep - Concordant (Baseline preferred)",
        
        Concordant & Preference_Pattern == "NonBaseline_Preferred" ~ 
          "Keep - Concordant",
        
        # Model clearly prefers baseline but existing doesn't match
        !Concordant & Preference_Pattern == "Baseline_Preferred" ~ 
          sprintf("Update to %s (Significantly better than alternatives)", Baseline),
        
        # Model clearly prefers non-baseline and it's different from existing
        !Concordant & Preference_Pattern == "NonBaseline_Preferred" & Model_Significant ~ 
          sprintf("Update to %s (Significantly better than baseline)", Model_Preferred),
        
        # Default
        TRUE ~ "Review manually"
      )
    )
  
  # Add summary statistics
  n_total <- nrow(concordance)
  n_concordant <- sum(concordance$Concordant, na.rm = TRUE)
  n_update_recommended <- sum(grepl("Update|Consider", concordance$Recommendation), 
                             na.rm = TRUE)
  
  cat("\n=== Preferred Codon Concordance Analysis ===\n\n")
  cat(sprintf("Total families: %d\n", n_total))
  cat(sprintf("Concordant: %d (%.1f%%)\n", 
              n_concordant, 100 * n_concordant / n_total))
  cat(sprintf("Updates recommended: %d (%.1f%%)\n", 
              n_update_recommended, 100 * n_update_recommended / n_total))
  cat(sprintf("Keep existing: %d (%.1f%%)\n\n", 
              n_total - n_update_recommended, 
              100 * (n_total - n_update_recommended) / n_total))
  
  # Show families needing review
  needs_review <- concordance %>%
    dplyr::filter(!Concordant | grepl("Review|Update", Recommendation))
  
  if (nrow(needs_review) > 0) {
    cat("Families needing attention:\n")
    print(needs_review %>% 
            dplyr::select(Family, Existing_Preferred, Model_Preferred, 
                         Model_Slope, Model_p_adj, Preference_Pattern, Recommendation))
    cat("\n")
  }
  
  return(concordance)
}
