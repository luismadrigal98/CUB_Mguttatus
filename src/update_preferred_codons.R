update_preferred_codons_from_models <- function(model_results,
                                                existing_preferred_codons = NULL) 
{
  #' Updates preferred codons based on model significance patterns
  #' This ROBUST version correctly identifies "co-optimal" (multiple winners)
  #' and "neutral" families.
  #'
  #' @param model_results List of model results from fit_pairwise_gams/glms
  #' @param existing_preferred_codons Optional data frame with existing preferences
  #' @return Data frame with updated preferred codons and preference patterns
  
  suppressPackageStartupMessages({
    require(dplyr)
    require(data.table)
  })
  
  # Extract all coefficients
  all_coefs <- data.table::rbindlist(
    lapply(model_results, function(x) {
      if (is.null(x)) return(NULL)
      x$coefficients %>%
        dplyr::mutate(Family = x$family, Baseline = x$baseline)
    })
  )
  
  # Analyze each family
  family_preferences <- lapply(unique(all_coefs$Family), function(fam) {
    
    fam_coefs <- all_coefs %>%
      dplyr::filter(Family == fam) %>%
      dplyr::arrange(dplyr::desc(Selection_Slope))
    
    baseline <- unique(fam_coefs$Baseline)[1]
    n_codons <- nrow(fam_coefs)
    
    # --- NEW ROBUST LOGIC ---
    
    # 1. Find the "Best" codon (the one with the maximum slope)
    best_codon_row <- fam_coefs[1, ] # Already sorted
    best_codon <- best_codon_row$Codon
    best_slope <- best_codon_row$Selection_Slope
    
    # 2. Find all "Loser" codons (significantly worse than baseline)

    sig_losers <- fam_coefs %>%
      dplyr::filter(Significant == TRUE & Selection_Slope < 0)
    
    # 3. Find all "Winner" codons (significantly better than baseline)
    sig_winners <- fam_coefs %>%
      dplyr::filter(Significant == TRUE & Selection_Slope > 0)
    
    n_winners <- nrow(sig_winners)
    n_losers <- nrow(sig_losers)
    n_neutral <- n_codons - n_winners - n_losers - 1
    
    # --- New Decision Logic ---
    preference_pattern <- "Review_Needed" # Default
    preferred_codons_str <- "Review"
    
    if (n_winners == 0 & n_losers == 0) {
      # --- Case 1: Neutral Family ---
      # No codon is significantly different from the baseline.

      preference_pattern <- "Neutral_Family"
      preferred_codons_str <- "None"
      
    } else if (n_winners > 0) {
      # --- Case 2: At least one codon is significantly *better* than the baseline ---
      # These are the "winners". There might be one, or multiple.

      preference_pattern <- "Preferred_NonBaseline"
      
      # Pick the codon with the higher slope
      preferred_codons_str <- best_codon_row |>
        pull(Codon)
      
    } else if (n_winners == 0 & n_losers > 0) {
      # --- Case 3: No codon is better, but at least one is worse ---
      # This means the baseline is the sole, unchallenged winner.
      # e.g., Ala, Gly, Lys, etc.
      preference_pattern <- "Preferred_is_Baseline"
      preferred_codons_str <- baseline
      
    } else {
      # This case (n_winners == 0, n_losers == 0, n_neutral > 0) is covered by Case 1
      # But we'll add a fallback
      preference_pattern <- "Neutral_Family"
      preferred_codons_str <- "None"
    }
    
    # Return summary
    data.frame(
      Family = fam,
      Baseline = baseline,
      N_Codons = n_codons,
      N_Sig_Negative = n_losers,
      N_Sig_Positive = n_winners,
      N_Non_Significant = n_neutral,
      Preference_Pattern = preference_pattern,
      Preferred_Codons = preferred_codons_str,
      stringsAsFactors = FALSE
    )
  })
  
  preferences_df <- data.table::rbindlist(family_preferences)
  
  # Print summary
  cat("\n=== Updated Preferred Codons Analysis ===\n\n")
  cat(sprintf("Total families: %d\n", nrow(preferences_df)))
  cat(sprintf("Preferred is Baseline: %d\n", 
              sum(preferences_df$Preference_Pattern == "Preferred_is_Baseline")))
  cat(sprintf("Preferred is Non-Baseline: %d\n", 
              sum(preferences_df$Preference_Pattern == "Preferred_NonBaseline")))
  cat(sprintf("Neutral Family (no preference): %d\n", 
              sum(preferences_df$Preference_Pattern == "Neutral_Family")))
  cat(sprintf("Needs review: %d\n\n", 
              sum(preferences_df$Preference_Pattern == "Review_Needed")))
  
  # Show patterns
  cat("Preference patterns by family:\n")
  print(preferences_df %>% 
          dplyr::select(Family, Preference_Pattern, Preferred_Codons, Baseline))
  cat("\n")
  
  return(preferences_df)
}

create_preferred_codons_table <- function(preferences_df) {
  #' Converts preferences_df to format compatible with existing pipeline
  #' Handles multiple preferred codons by creating multiple rows
  #'
  #' @param preferences_df Output from update_preferred_codons_from_models
  #' @return Data frame in preferred_codons_mg format
  
  suppressPackageStartupMessages({
    require(dplyr)
  })
  
  # Expand families with multiple preferred codons
  expanded <- lapply(1:nrow(preferences_df), function(i) {
    row <- preferences_df[i, ]
    
    if (row$Preferred_Codons == "None") {
      # No preference - use baseline but mark as non-significant
      data.frame(
        AA = row$Family,
        Codon = row$Baseline,
        relative_adaptiveness = 0.5,  # Neutral
        Significant = FALSE,
        Selection_Rationale = "No clear preference",
        Preference_Pattern = row$Preference_Pattern,
        stringsAsFactors = FALSE
      )
      
    } else if (row$Preferred_Codons == "Review") {
      # Needs review
      data.frame(
        AA = row$Family,
        Codon = row$Baseline,
        relative_adaptiveness = NA,
        Significant = FALSE,
        Selection_Rationale = "Manual review needed",
        Preference_Pattern = row$Preference_Pattern,
        stringsAsFactors = FALSE
      )
      
    } else if (grepl(";", row$Preferred_Codons)) {
      # Multiple preferred codons
      codons <- strsplit(row$Preferred_Codons, ";")[[1]]
      data.frame(
        AA = rep(row$Family, length(codons)),
        Codon = codons,
        relative_adaptiveness = rep(1.0, length(codons)),
        Significant = rep(TRUE, length(codons)),
        Selection_Rationale = rep("Co-optimal codons", length(codons)),
        Preference_Pattern = rep(row$Preference_Pattern, length(codons)),
        stringsAsFactors = FALSE
      )
      
    } else {
      # Single preferred codon
      data.frame(
        AA = row$Family,
        Codon = row$Preferred_Codons,
        relative_adaptiveness = 1.0,
        Significant = TRUE,
        Selection_Rationale = ifelse(
          row$Preference_Pattern == "Single_Preferred_Baseline",
          "Selection (baseline)",
          "Selection (model-driven)"
        ),
        Preference_Pattern = row$Preference_Pattern,
        stringsAsFactors = FALSE
      )
    }
  })
  
  preferred_table <- data.table::rbindlist(expanded)
  
  return(preferred_table)
}