#!/usr/bin/env Rscript
# classify_codon_anticodon_pairing.R
# 
# This function classifies codon-anticodon pairing types (Watson-Crick vs Wobble)
# and tests the translational accuracy hypothesis: that selection favors codons
# with high-fidelity Watson-Crick pairing over error-prone wobble pairing.
#
# Key biological concepts:
# - Watson-Crick pairs (G-C, C-G, A-U, U-A): High fidelity, accurate translation
# - Wobble pairs (G-U, U-G, and Inosine pairs): Faster but less accurate
# - Hypothesis: "Preferred" codons use Watson-Crick pairs at codon position 3

classify_codon_anticodon_pairing <- function(
  tRNA_data,
  codon_supply,
  preferred_codons,
  output_dir = "results/tRNA_analysis",
  save_results = TRUE
) {
  #' Classify Codon-Anticodon Pairing and Test Translational Accuracy Hypothesis
  #'
  #' @param tRNA_data data.frame with columns: tRNA_type (AA), Anticodon, and abundance measure
  #' @param codon_supply data.frame with columns: Codon, tRNA_gene_count (or other abundance)
  #' @param preferred_codons data.frame with columns: Codon, Amino_Acid, Status (Preferred/Non-Preferred)
  #' @param output_dir character, directory to save results
  #' @param save_results logical, whether to save outputs
  #'
  #' @return list with:
  #'   - pairing_classification: table of all codons with pairing types
  #'   - contingency_table: Status vs Pairing_Type cross-tabulation
  #'   - chi_squared_test: statistical test result
  #'   - plot: visualization of pairing type by codon status
  #'
  #' @details
  #' Wobble base pairing rules (Crick, 1966):
  #' Position 1 of anticodon (pairs with position 3 of codon):
  #' - Standard pairs: G-C, C-G, A-U, U-A (Watson-Crick)
  #' - Wobble pairs: G-U, U-G (wobble)
  #' - Inosine (I) pairs: I-U, I-C, I-A (wobble, very permissive)
  #'
  #' The third codon position is where wobble occurs. This function classifies
  #' each codon based on the pairing at this critical position.
  
  require(dplyr)
  require(data.table)
  require(ggplot2)
  
  cat("\n=== Translational Accuracy Hypothesis Test ===\n")
  cat("Testing if selection favors Watson-Crick over wobble pairing\n\n")
  
  # ---- Step 1: Create codon-anticodon pairing map ----
  cat("Step 1: Mapping codons to anticodons and classifying pairing...\n")
  
  # Get reverse complement function
  reverse_complement <- function(seq) {
    # Handle both DNA and RNA
    bases <- c("A" = "U", "U" = "A", "G" = "C", "C" = "G", "T" = "A", 
               "a" = "u", "u" = "a", "g" = "c", "c" = "g", "t" = "a")
    seq_bases <- strsplit(seq, "")[[1]]
    rev_seq <- rev(seq_bases)
    comp <- bases[rev_seq]
    # Handle any NA (unknown bases)
    comp[is.na(comp)] <- "N"
    paste(comp, collapse = "")
  }
  
  # Classify pairing type at wobble position (codon pos 3, anticodon pos 1)
  classify_wobble_pair <- function(codon_base, anticodon_base) {
    # Watson-Crick pairs
    if ((codon_base == "A" && anticodon_base == "U") ||
        (codon_base == "U" && anticodon_base == "A") ||
        (codon_base == "G" && anticodon_base == "C") ||
        (codon_base == "C" && anticodon_base == "G")) {
      return("Watson-Crick")
    }
    # Wobble G-U pairs
    else if ((codon_base == "G" && anticodon_base == "U") ||
             (codon_base == "U" && anticodon_base == "G")) {
      return("Wobble_GU")
    }
    # Inosine wobble (very permissive)
    else if (anticodon_base == "I") {
      if (codon_base %in% c("A", "C", "U")) {
        return("Wobble_Inosine")
      } else {
        return("Unknown")
      }
    }
    else {
      return("Unknown")
    }
  }
  
  # Prepare tRNA data
  tRNA_dt <- data.table::as.data.table(tRNA_data)
  
  # Get anticodon column name (flexible)
  anticodon_col <- names(tRNA_dt)[grep("anticodon", names(tRNA_dt), ignore.case = TRUE)]
  if (length(anticodon_col) == 0) {
    anticodon_col <- "Anticodon"
  } else {
    anticodon_col <- anticodon_col[1]
  }
  
  # Get amino acid column
  aa_col <- names(tRNA_dt)[grep("tRNA_type|amino.acid|aa", names(tRNA_dt), ignore.case = TRUE)]
  if (length(aa_col) == 0) {
    stop("Cannot find amino acid column in tRNA_data")
  }
  aa_col <- aa_col[1]
  
  # Standardize column names
  tRNA_dt <- tRNA_dt[, .(
    Amino_Acid = get(aa_col),
    Anticodon = get(anticodon_col)
  )]
  
  # Convert T to U in anticodons (keep DNA format for now, convert later)
  tRNA_dt[, Anticodon := toupper(Anticodon)]
  
  # Get unique anticodons per amino acid
  anticodon_map <- tRNA_dt[, .(Anticodons = list(unique(Anticodon))), by = Amino_Acid]
  
  cat("  Found", nrow(tRNA_dt), "tRNA genes with", length(unique(tRNA_dt$Anticodon)), "unique anticodons\n")
  
  # ---- Step 2: Classify each codon using wobble base pairing rules ----
  cat("Step 2: Applying wobble rules to predict codon-anticodon pairing...\n")
  
  # Get genetic code
  genetic_code <- c(
    "UUU" = "Phe", "UUC" = "Phe", "UUA" = "Leu", "UUG" = "Leu",
    "UCU" = "Ser", "UCC" = "Ser", "UCA" = "Ser", "UCG" = "Ser",
    "UAU" = "Tyr", "UAC" = "Tyr", "UAA" = "Stop", "UAG" = "Stop",
    "UGU" = "Cys", "UGC" = "Cys", "UGA" = "Stop", "UGG" = "Trp",
    "CUU" = "Leu", "CUC" = "Leu", "CUA" = "Leu", "CUG" = "Leu",
    "CCU" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
    "CAU" = "His", "CAC" = "His", "CAA" = "Gln", "CAG" = "Gln",
    "CGU" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg",
    "AUU" = "Ile", "AUC" = "Ile", "AUA" = "Ile", "AUG" = "Met",
    "ACU" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
    "AAU" = "Asn", "AAC" = "Asn", "AAA" = "Lys", "AAG" = "Lys",
    "AGU" = "Ser", "AGC" = "Ser", "AGA" = "Arg", "AGG" = "Arg",
    "GUU" = "Val", "GUC" = "Val", "GUA" = "Val", "GUG" = "Val",
    "GCU" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
    "GAU" = "Asp", "GAC" = "Asp", "GAA" = "Glu", "GAG" = "Glu",
    "GGU" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly"
  )
  
  # Create codon classification table
  codon_classification <- data.table::data.table(
    Codon = names(genetic_code),
    Amino_Acid = as.character(genetic_code)
  )
  
  # Remove stop codons
  codon_classification <- codon_classification[Amino_Acid != "Stop"]
  
  # Add pairing classification
  codon_classification[, Pairing_Type := NA_character_]
  codon_classification[, Matching_Anticodon := NA_character_]
  codon_classification[, Wobble_Position_Pair := NA_character_]
  
  # Create a map of which anticodons can read which codons (using wobble rules)
  # For each codon, find ALL anticodons that could potentially read it
  
  for (i in seq_len(nrow(codon_classification))) {
    codon <- codon_classification$Codon[i]
    aa <- codon_classification$Amino_Acid[i]
    
    # Get anticodons for this amino acid
    aa_match <- anticodon_map[Amino_Acid == aa]
    if (nrow(aa_match) == 0) {
      # No tRNA for this AA - infer pairing from standard wobble rules
      # This shouldn't happen for standard amino acids, but handle it
      next
    }
    anticodons <- aa_match$Anticodons[[1]]
    
    if (length(anticodons) == 0) {
      next
    }
    
    # Wobble base pairing rules (Crick 1966):
    # Codon (mRNA):      5'-A  -B  -C  -3'
    # Anticodon (tRNA):  3'-A' -B' -C' -5' (antiparallel)
    #
    # In your tRNA file, anticodons are written 5'->3' as X-Y-Z
    # During translation they flip to 3'-X-Y-Z-5'
    # So pairing is:
    #   Codon position 1 (A) pairs with anticodon position 3 (Z)
    #   Codon position 2 (B) pairs with anticodon position 2 (Y)  
    #   Codon position 3 (C) pairs with anticodon position 1 (X) <- WOBBLE POSITION
    
    pairing_types <- c()
    matching_anticodons <- c()
    wobble_pairs <- c()
    
    for (anticodon in anticodons) {
      # Extract bases
      # Anticodon as written in file: 5'-X-Y-Z-3'
      # During pairing (antiparallel): 3'-X-Y-Z-5'
      ac_x <- substr(anticodon, 1, 1)  # Pairs with codon pos 3 (WOBBLE)
      ac_y <- substr(anticodon, 2, 2)  # Pairs with codon pos 2
      ac_z <- substr(anticodon, 3, 3)  # Pairs with codon pos 1
      
      codon_a <- substr(codon, 1, 1)
      codon_b <- substr(codon, 2, 2)
      codon_c <- substr(codon, 3, 3)
      
      # Convert DNA to RNA
      ac_x <- gsub("T", "U", ac_x)
      ac_y <- gsub("T", "U", ac_y)
      ac_z <- gsub("T", "U", ac_z)
      codon_a <- gsub("T", "U", codon_a)
      codon_b <- gsub("T", "U", codon_b)
      codon_c <- gsub("T", "U", codon_c)
      
      # Check if positions 1 and 2 match (standard Watson-Crick)
      complement <- c("A" = "U", "U" = "A", "G" = "C", "C" = "G")
      
      if (complement[codon_a] == ac_z && complement[codon_b] == ac_y) {
        # First two positions match, now classify wobble position
        pair_type <- classify_wobble_pair(codon_c, ac_x)
        
        if (pair_type != "Unknown") {
          pairing_types <- c(pairing_types, pair_type)
          matching_anticodons <- c(matching_anticodons, anticodon)
          wobble_pairs <- c(wobble_pairs, paste0(codon_c, "-", ac_x))
        }
      }
    }
    
    if (length(pairing_types) > 0) {
      # Prioritize Watson-Crick if available (most accurate)
      if ("Watson-Crick" %in% pairing_types) {
        best_idx <- which(pairing_types == "Watson-Crick")[1]
        codon_classification$Pairing_Type[i] <- "Watson-Crick"
      } else {
        # Otherwise use the first wobble type
        best_idx <- 1
        codon_classification$Pairing_Type[i] <- pairing_types[best_idx]
      }
      codon_classification$Matching_Anticodon[i] <- matching_anticodons[best_idx]
      codon_classification$Wobble_Position_Pair[i] <- wobble_pairs[best_idx]
    }
  }
  
  # Simplify pairing types for analysis
  codon_classification[, Pairing_Category := ifelse(
    Pairing_Type == "Watson-Crick", 
    "Watson-Crick", 
    "Wobble"
  )]
  
  cat("  Found", sum(codon_classification$Pairing_Category == "Watson-Crick", na.rm = TRUE), 
      "Watson-Crick codons\n")
  cat("  Found", sum(codon_classification$Pairing_Category == "Wobble", na.rm = TRUE), 
      "Wobble codons\n")
  cat("  Found", sum(is.na(codon_classification$Pairing_Category)), 
      "codons with unknown pairing\n\n")
  
  # ---- Step 3: Merge with preferred codon status ----
  cat("Step 3: Merging with preferred codon classification...\n")
  
  # Ensure preferred_codons has Status column
  if (!"Status" %in% names(preferred_codons)) {
    # Assume all codons in this table are preferred
    preferred_codons$Status <- "Preferred"
    
    # Mark all other codons as non-preferred
    all_codons <- unique(codon_classification$Codon)
    non_preferred <- all_codons[!all_codons %in% preferred_codons$Codon]
    
    non_preferred_df <- data.frame(
      Codon = non_preferred,
      Status = "Non-Preferred"
    )
    
    preferred_codons <- rbind(
      preferred_codons[, c("Codon", "Status")],
      non_preferred_df
    )
  }
  
  # Merge
  analysis_dt <- merge(
    codon_classification,
    preferred_codons[, c("Codon", "Status")],
    by = "Codon",
    all.x = TRUE
  )
  
  # Remove codons without classification
  analysis_dt <- analysis_dt[!is.na(Pairing_Category) & !is.na(Status)]
  
  cat("  Total codons in analysis:", nrow(analysis_dt), "\n\n")
  
  # ---- Step 4: Statistical test ----
  cat("Step 4: Running Fisher's exact test...\n")
  
  contingency <- table(analysis_dt$Status, analysis_dt$Pairing_Category)
  cat("\nContingency table:\n")
  print(contingency)
  cat("\n")
  
  # Use Fisher's exact test (better for small expected frequencies)
  fisher_test <- fisher.test(contingency)
  cat("Fisher's exact test:\n")
  cat("  Odds ratio =", round(fisher_test$estimate, 3), "\n")
  cat("  p-value =", format.pval(fisher_test$p.value, digits = 3), "\n")
  cat("  95% CI: [", round(fisher_test$conf.int[1], 3), ",", 
      round(fisher_test$conf.int[2], 3), "]\n\n")
  
  # Also calculate Cramer's V for effect size
  # Using chi-squared statistic for effect size calculation
  chi_stat <- sum((contingency - rowSums(contingency) %*% t(colSums(contingency)) / sum(contingency))^2 / 
                    (rowSums(contingency) %*% t(colSums(contingency)) / sum(contingency)))
  cramers_v <- sqrt(chi_stat / (sum(contingency) * (min(dim(contingency)) - 1)))
  cat("  Cramer's V =", round(cramers_v, 3), "(effect size)\n\n")
  
  # Calculate proportions
  prop_table <- prop.table(contingency, margin = 1)
  cat("Proportions (by row):\n")
  print(round(prop_table, 3))
  cat("\n")
  
  # ---- Step 5: Visualization ----
  cat("Step 5: Creating visualization...\n")
  
  # Prepare data for plotting
  plot_data <- analysis_dt %>%
    group_by(Status, Pairing_Category) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Status) %>%
    mutate(
      Total = sum(Count),
      Proportion = Count / Total,
      Percentage = Proportion * 100
    )
  
  # Create bar plot
  p1 <- ggplot(plot_data, aes(x = Status, y = Percentage, fill = Pairing_Category)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold") +
    scale_fill_manual(
      values = c("Watson-Crick" = "#2E7D32", "Wobble" = "#F57C00"),
      labels = c("Watson-Crick\n(High Fidelity)", "Wobble\n(Low Fidelity)")
    ) +
    labs(
      title = "Translational Accuracy Hypothesis Test",
      subtitle = sprintf("Fisher's exact: OR = %.2f, p = %s, Cramer's V = %.3f",
                        fisher_test$estimate, 
                        format.pval(fisher_test$p.value, digits = 2),
                        cramers_v),
      x = "Codon Status",
      y = "Percentage of Codons",
      fill = "Pairing Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    )
  
  # Create grouped bar plot for easier comparison
  p2 <- ggplot(plot_data, aes(x = Pairing_Category, y = Percentage, fill = Status)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 4) +
    scale_fill_manual(
      values = c("Preferred" = "#1976D2", "Non-Preferred" = "#757575")
    ) +
    labs(
      title = "Codon Preference by Pairing Type",
      subtitle = "Do preferred codons use more accurate Watson-Crick pairing?",
      x = "Pairing Type at Wobble Position",
      y = "Percentage of Codons",
      fill = "Codon Status"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    ) +
    ylim(0, max(plot_data$Percentage) * 1.1)
  
  # Create detailed table plot showing each amino acid family
  family_summary <- analysis_dt %>%
    group_by(Amino_Acid, Status, Pairing_Category) %>%
    summarise(N_Codons = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = c(Status, Pairing_Category),
      values_from = N_Codons,
      values_fill = 0
    )
  
  # ---- Step 6: Save results ----
  if (save_results) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    cat("Step 6: Saving results...\n")
    
    # Save classification table
    write.csv(
      analysis_dt,
      file.path(output_dir, "codon_anticodon_pairing_classification.csv"),
      row.names = FALSE
    )
    
    # Save contingency table
    write.csv(
      as.data.frame.matrix(contingency),
      file.path(output_dir, "pairing_contingency_table.csv")
    )
    
    # Save test results
    test_summary <- data.frame(
      Test = "Fisher_exact",
      Odds_ratio = fisher_test$estimate,
      p_value = fisher_test$p.value,
      CI_lower = fisher_test$conf.int[1],
      CI_upper = fisher_test$conf.int[2],
      Cramers_V = cramers_v
    )
    write.csv(
      test_summary,
      file.path(output_dir, "pairing_hypothesis_test.csv"),
      row.names = FALSE
    )
    
    # Save plots
    ggsave(
      file.path(output_dir, "pairing_type_by_status.png"),
      plot = p1,
      width = 10, height = 7, dpi = 300
    )
    
    ggsave(
      file.path(output_dir, "status_by_pairing_type.png"),
      plot = p2,
      width = 10, height = 7, dpi = 300
    )
    
    cat("  Results saved to:", output_dir, "\n")
  }
  
  # ---- Step 7: Interpretation ----
  cat("\n=== INTERPRETATION ===\n")
  
  if (fisher_test$p.value < 0.05) {
    cat("✓ SIGNIFICANT association between codon status and pairing type (p < 0.05)\n\n")
    
    # Check direction of effect
    pref_wc <- prop_table["Preferred", "Watson-Crick"]
    nonpref_wc <- prop_table["Non-Preferred", "Watson-Crick"]
    
    if (pref_wc > nonpref_wc) {
      cat("✓ SUPPORTS Translational Accuracy Hypothesis:\n")
      cat("  - Preferred codons have MORE Watson-Crick pairing (", 
          round(pref_wc * 100, 1), "% vs ", round(nonpref_wc * 100, 1), "%)\n", sep = "")
      cat("  - Selection appears to favor HIGH-FIDELITY translation\n")
      cat("  - This explains why simple tRNA abundance didn't correlate:\n")
      cat("    Selection is for ACCURACY, not just SPEED\n")
    } else {
      cat("✗ CONTRADICTS Translational Accuracy Hypothesis:\n")
      cat("  - Preferred codons have LESS Watson-Crick pairing\n")
      cat("  - Selection may favor wobble pairing for some reason\n")
    }
  } else {
    cat("✗ NO SIGNIFICANT association (p >= 0.05)\n")
    cat("  - Pairing type does not explain codon preference\n")
    cat("  - Selection may be driven by other factors\n")
  }
  
  cat("\n=== Analysis Complete ===\n\n")
  
  # Return results
  return(list(
    pairing_classification = analysis_dt,
    contingency_table = contingency,
    fisher_test = fisher_test,
    cramers_v = cramers_v,
    proportions = prop_table,
    family_summary = family_summary,
    plot_stacked = p1,
    plot_grouped = p2
  ))
}


classify_pairing_with_expression_preferred <- function(
  tRNA_data,
  codon_supply,
  preferred_codons_corrected,
  output_dir = "results/tRNA_analysis_pairing_expression",
  save_results = TRUE
) {
  #' Classify Codon-Anticodon Pairing Using Expression-Based Preferred Codons
  #'
  #' @param tRNA_data data.frame with tRNA information
  #' @param codon_supply data.frame with codon supply information
  #' @param preferred_codons_corrected data.frame with corrected preferred codons (from expression analysis)
  #' @param output_dir character, directory to save results
  #' @param save_results logical, whether to save outputs
  #'
  #' @return list with analysis results
  #'
  #' @details
  #' This function parallels the CA-based analysis but uses expression-based
  #' preferred codons (those identified as preferred from expression analysis).
  #' This tests whether the translational accuracy hypothesis holds when using 
  #' a different criterion for defining "preferred" codons.
  
  require(dplyr)
  
  cat("\n=== Translational Accuracy Test (Expression-Based Preferred Codons) ===\n")
  cat("Using expression-based preferred codons (w=1 AND/OR enriched)\n\n")
  
  # Extract preferred codons from the corrected table
  # These are codons identified as preferred through expression analysis
  preferred_codons_expr <- preferred_codons_corrected %>%
    mutate(
      Codon = toupper(gsub("T", "U", Codon)),
      Amino_Acid = Amino_Acid,
      Status = "Preferred"
    ) %>%
    select(Codon, Amino_Acid, Status)
  
  cat("Found", nrow(preferred_codons_expr), "preferred codons (w = 1.0)\n")
  
  # Get all other codons as non-preferred
  genetic_code <- c(
    "UUU" = "Phe", "UUC" = "Phe", "UUA" = "Leu", "UUG" = "Leu",
    "UCU" = "Ser", "UCC" = "Ser", "UCA" = "Ser", "UCG" = "Ser",
    "UAU" = "Tyr", "UAC" = "Tyr", "UAA" = "Stop", "UAG" = "Stop",
    "UGU" = "Cys", "UGC" = "Cys", "UGA" = "Stop", "UGG" = "Trp",
    "CUU" = "Leu", "CUC" = "Leu", "CUA" = "Leu", "CUG" = "Leu",
    "CCU" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
    "CAU" = "His", "CAC" = "His", "CAA" = "Gln", "CAG" = "Gln",
    "CGU" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg",
    "AUU" = "Ile", "AUC" = "Ile", "AUA" = "Ile", "AUG" = "Met",
    "ACU" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
    "AAU" = "Asn", "AAC" = "Asn", "AAA" = "Lys", "AAG" = "Lys",
    "AGU" = "Ser", "AGC" = "Ser", "AGA" = "Arg", "AGG" = "Arg",
    "GUU" = "Val", "GUC" = "Val", "GUA" = "Val", "GUG" = "Val",
    "GCU" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
    "GAU" = "Asp", "GAC" = "Asp", "GAA" = "Glu", "GAG" = "Glu",
    "GGU" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly"
  )
  
  all_codons <- data.frame(
    Codon = names(genetic_code),
    Amino_Acid = as.character(genetic_code)
  ) %>%
    filter(Amino_Acid != "Stop")
  
  # Add status to all codons
  all_codons <- all_codons %>%
    left_join(preferred_codons_expr, by = c("Codon", "Amino_Acid")) %>%
    mutate(Status = ifelse(is.na(Status), "Non-Preferred", Status)) %>%
    select(Codon, Amino_Acid, Status)
  
  cat("Total:", sum(all_codons$Status == "Preferred"), "preferred,",
      sum(all_codons$Status == "Non-Preferred"), "non-preferred\n\n")
  
  # Run the main classification function
  results <- classify_codon_anticodon_pairing(
    tRNA_data = tRNA_data,
    codon_supply = codon_supply,
    preferred_codons = all_codons,
    output_dir = output_dir,
    save_results = save_results
  )
  
  cat("\n=== Expression-Based Analysis Complete ===\n")
  cat("Comparing with CA-based results will show if criterion matters\n\n")
  
  return(results)
}
