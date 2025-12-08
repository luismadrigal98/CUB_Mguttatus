#!/usr/bin/env Rscript
# cfs_analysis.R
# 
# Codon Frequency Spectrum (CFS) Analysis for Preferred Codon Frequencies
# 
# NOTE: We use "Codon Frequency Spectrum" rather than "Site Frequency Spectrum"
# because we lack ancestral state information and outgroup data. This is a
# folded spectrum based on preferred vs non-preferred codon classification.
#
# This script analyzes the distribution of preferred codon frequencies
# to detect signatures of selection on codon usage. Neutral expectations
# depend on family degeneracy:
#   - 2-fold: expected freq = 0.50 (1/2 codons preferred)
#   - 3-fold: expected freq = 0.33 (1/3 codons preferred)  
#   - 4-fold: expected freq = 0.25 (1/4 codons preferred)
#   - 6-fold: expected freq = 0.17 (1/6 codons preferred)
#
# Deviations from these expectations indicate selection:
# - Excess above expectation: selection FOR preferred codons
# - Deficit below expectation: selection AGAINST preferred codons
#
# Author: Luis J. Madrigal-Roca
# Date: December 2025


# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

#' Load site frequency spectrum data
#' 
#' @param raw_freq_file Path to raw frequencies file
#' @param full_dist_file Path to full distribution file
#' @param data_dir Directory containing data files
#' 
#' @return List with raw_frequencies and binned_distribution data frames
load_cfs_data <- function(
    raw_freq_file = "all_chromosomes.raw_frequencies.txt",
    full_dist_file = "all_chromosomes.full_distribution.txt",
    data_dir = "./data"
) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("LOADING SITE FREQUENCY SPECTRUM DATA\n")
  cat(strrep("=", 70), "\n\n")
  
  # Load raw frequencies (polymorphic sites only)
  raw_path <- file.path(data_dir, raw_freq_file)
  if (!file.exists(raw_path)) {
    # Try results directory
    raw_path <- file.path("./results", raw_freq_file)
  }
  
  raw_freq <- fread(raw_path, skip = 3)  # Skip comment lines
  colnames(raw_freq) <- c("Amino_Acid", "Family", "Preferred_Freq")
  
  cat("Raw frequencies loaded:", nrow(raw_freq), "polymorphic sites\n")
  cat("  Amino acids:", paste(unique(raw_freq$Amino_Acid), collapse = ", "), "\n")
  cat("  Frequency range:", round(min(raw_freq$Preferred_Freq), 4), "-", 
      round(max(raw_freq$Preferred_Freq), 4), "\n\n")
  
  # Load full distribution (including fixed sites)
  dist_path <- file.path(data_dir, full_dist_file)
  if (!file.exists(dist_path)) {
    dist_path <- file.path("./results", full_dist_file)
  }
  
  full_dist <- fread(dist_path, skip = 8)  # Skip comment lines
  
  cat("Full distribution loaded:", nrow(full_dist), "amino acids\n")
  cat("  Total sites per AA (example A):", 
      sum(as.numeric(full_dist[1, -c(1,2)]), na.rm = TRUE), "\n\n")
  
  return(list(
    raw_frequencies = raw_freq,
    binned_distribution = full_dist
  ))
}


# =============================================================================
# SITE FREQUENCY SPECTRUM VISUALIZATION
# =============================================================================

#' Create suspension bridge plot (classic CFS visualization)
#' 
#' This plot shows ALL SITES including both polymorphic AND fixed sites.
#' The characteristic U-shape ("suspension bridge") reflects:
#' - Many sites fixed for non-preferred codons (mutation bias)
#' - Many sites fixed for preferred codons (selection)
#' - Fewer polymorphic sites in between
#' 
#' @param full_dist Full distribution data from load_cfs_data()
#' @param by_family Facet by degeneracy family?
#' @param output_dir Directory to save plots
#' 
#' @return List with ggplot objects and data
plot_suspension_bridge <- function(
    full_dist,
    by_family = TRUE,
    output_dir = "./results/sfs_analysis"
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Creating suspension bridge plot...\n")
  
  # Reshape to long format
  dist_long <- full_dist %>%
    pivot_longer(
      cols = -c(Amino_Acid, Family),
      names_to = "Bin",
      values_to = "Count"
    ) %>%
    mutate(
      # Extract bin number
      Bin_Num = case_when(
        grepl("Bin_0", Bin) ~ 0,
        grepl("Bin_11", Bin) ~ 11,
        TRUE ~ as.numeric(gsub("Bin_(\\d+)_.*", "\\1", Bin))
      ),
      # Create frequency midpoint for x-axis
      Freq_Midpoint = case_when(
        Bin_Num == 0 ~ 0,
        Bin_Num == 11 ~ 1,
        TRUE ~ (Bin_Num - 1) / 10 + 0.05  # Midpoint of each 0.1 bin
      ),
      # Bin labels
      Bin_Label = case_when(
        Bin_Num == 0 ~ "Fixed\nNon-Pref",
        Bin_Num == 11 ~ "Fixed\nPref",
        TRUE ~ sprintf("%.1f-%.1f", (Bin_Num - 1) / 10, Bin_Num / 10)
      )
    )
  
  # Aggregate by family for faceting
  dist_by_family <- dist_long %>%
    group_by(Family, Bin_Num, Freq_Midpoint, Bin_Label) %>%
    summarise(
      Total_Count = sum(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Family) %>%
    mutate(
      Proportion = Total_Count / sum(Total_Count),
      Log_Count = log10(Total_Count + 1)
    )
  
  # Total across all amino acids
  dist_total <- dist_long %>%
    group_by(Bin_Num, Freq_Midpoint, Bin_Label) %>%
    summarise(
      Total_Count = sum(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Proportion = Total_Count / sum(Total_Count),
      Log_Count = log10(Total_Count + 1)
    )
  
  # Plot 1: Overall suspension bridge (log scale)
  p_bridge <- ggplot(dist_total, aes(x = Freq_Midpoint, y = Total_Count)) +
    geom_col(fill = "steelblue", color = "black", width = 0.08, position = "identity") +
    geom_line(color = "darkred", linewidth = 1, group = 1) +
    geom_point(color = "darkred", size = 2) +
    scale_y_log10(labels = scales::comma) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0\n(Fixed Non-Pref)", "0.25", "0.5", "0.75", "1\n(Fixed Pref)")
    ) +
    labs(
      title = "Codon Frequency Spectrum: ALL SITES (Fixed + Polymorphic)",
      subtitle = "Classic 'Suspension Bridge' pattern - U-shape indicates selection",
      x = "Preferred Codon Frequency",
      y = "Number of Sites (log scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = 0.1, y = max(dist_total$Total_Count) * 0.5,
             label = "Non-preferred\nfixed/common", hjust = 0, size = 3.5) +
    annotate("text", x = 0.9, y = max(dist_total$Total_Count) * 0.5,
             label = "Preferred\nfixed/common", hjust = 1, size = 3.5)
  
  ggsave(file.path(output_dir, "cfs_suspension_bridge.pdf"),
         p_bridge, width = 10, height = 6)
  
  # Plot 2: By degeneracy family
  if (by_family) {
    p_by_family <- ggplot(dist_by_family, aes(x = Freq_Midpoint, y = Total_Count)) +
      geom_col(aes(fill = Family), color = "black", width = 0.08, position = "identity") +
      geom_line(color = "darkred", linewidth = 0.8, group = 1) +
      scale_y_log10(labels = scales::comma) +
      scale_fill_viridis_d(option = "plasma") +
      facet_wrap(~ Family, scales = "free_y") +
      labs(
        title = "CFS by Degeneracy Family: ALL SITES (Fixed + Polymorphic)",
        x = "Preferred Codon Frequency",
        y = "Number of Sites (log scale)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )
    
    ggsave(file.path(output_dir, "cfs_by_family.pdf"),
           p_by_family, width = 12, height = 8)
  }
  
  # Plot 3: Proportion view (normalized)
  p_proportion <- ggplot(dist_total, aes(x = Freq_Midpoint, y = Proportion)) +
    geom_col(fill = "coral", color = "black", width = 0.08, position = "identity") +
    geom_line(color = "darkred", linewidth = 1, group = 1) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Normalized CFS: ALL SITES (Fixed + Polymorphic)",
      subtitle = "Proportion of sites in each frequency bin",
      x = "Preferred Codon Frequency",
      y = "Proportion of Sites"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "cfs_proportions.pdf"),
         p_proportion, width = 10, height = 6)
  
  cat("Suspension bridge plots saved to:", output_dir, "\n\n")
  
  return(list(
    overall = p_bridge,
    by_family = if(by_family) p_by_family else NULL,
    proportions = p_proportion,
    data_total = dist_total,
    data_by_family = dist_by_family
  ))
}


#' Plot frequency distributions for POLYMORPHIC SITES ONLY
#' 
#' All plots in this function use only polymorphic sites (0 < freq < 1).
#' Fixed sites are excluded. Neutral expectation lines are added where appropriate.
#' 
#' @param raw_freq Raw frequency data from load_cfs_data() - polymorphic sites only
#' @param output_dir Directory to save plots
#' 
#' @return List of ggplot objects
plot_polymorphic_frequency_distributions <- function(
    raw_freq,
    output_dir = "./results/cfs_analysis"
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Creating polymorphic site frequency distribution plots...\n")
  cat("  NOTE: All plots show POLYMORPHIC SITES ONLY (0 < freq < 1)\n\n")
  
 # Define neutral expectations by family
  neutral_expectations <- data.frame(
    Family = c("2-fold", "3-fold", "4-fold", "6-fold"),
    Neutral_Expectation = c(1/2, 1/3, 1/4, 1/6),
    N_Codons = c(2, 3, 4, 6)
  )
  
  # Add neutral expectation to raw_freq
  raw_freq_with_null <- raw_freq %>%
    left_join(neutral_expectations, by = "Family")
  
  # -------------------------------------------------------------------------
  # Plot 1: Overall density with annotation
  # -------------------------------------------------------------------------
  overall_mean <- mean(raw_freq$Preferred_Freq)
  
  p_density_overall <- ggplot(raw_freq, aes(x = Preferred_Freq)) +
    geom_density(fill = "steelblue", alpha = 0.6, color = "darkblue") +
    geom_vline(xintercept = overall_mean, color = "red", linewidth = 1, linetype = "solid") +
    geom_vline(xintercept = 0.5, color = "gray40", linewidth = 1, linetype = "dashed") +
    annotate("text", x = overall_mean + 0.02, y = Inf, vjust = 2, hjust = 0,
             label = sprintf("Observed mean = %.3f", overall_mean), color = "red", size = 4) +
    annotate("text", x = 0.48, y = Inf, vjust = 4, hjust = 1,
             label = "Neutral (if all 2-fold)", color = "gray40", size = 3.5) +
    labs(
      title = "Preferred Codon Frequency Distribution (Polymorphic Sites Only)",
      subtitle = sprintf("n = %s polymorphic sites | Red = observed mean, Dashed = neutral expectation",
                        format(nrow(raw_freq), big.mark = ",")),
      x = "Preferred Codon Frequency",
      y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "polymorphic_density_overall.pdf"),
         p_density_overall, width = 10, height = 6)
  
  # -------------------------------------------------------------------------
  # Plot 2: Density by FAMILY with family-specific neutral expectations
  # -------------------------------------------------------------------------
  p_density_family <- ggplot(raw_freq_with_null, aes(x = Preferred_Freq)) +
    geom_density(aes(fill = Family), alpha = 0.6, color = "black") +
    # Add family-specific neutral expectation as dashed vertical line
    geom_vline(aes(xintercept = Neutral_Expectation), 
               linetype = "dashed", linewidth = 1, color = "red") +
    facet_wrap(~ Family, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Preferred Codon Frequency by Degeneracy Family (Polymorphic Sites)",
      subtitle = "Red dashed line = neutral expectation (1/n codons)",
      x = "Preferred Codon Frequency",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(file.path(output_dir, "polymorphic_density_by_family.pdf"),
         p_density_family, width = 10, height = 8)
  
  # -------------------------------------------------------------------------
  # Plot 3: Density by AMINO ACID with family-specific neutral expectations
  # -------------------------------------------------------------------------
  p_density_aa <- ggplot(raw_freq_with_null, aes(x = Preferred_Freq)) +
    geom_density(aes(fill = Family), alpha = 0.6, color = "black") +
    # Add neutral expectation line for each AA (based on its family)
    geom_vline(aes(xintercept = Neutral_Expectation), 
               linetype = "dashed", linewidth = 0.8, color = "red") +
    facet_wrap(~ Amino_Acid, scales = "free_y", ncol = 6) +
    scale_fill_brewer(palette = "Set2", name = "Degeneracy") +
    labs(
      title = "Preferred Codon Frequency by Amino Acid (Polymorphic Sites)",
      subtitle = "Red dashed line = neutral expectation (1/n codons for each family)",
      x = "Preferred Codon Frequency",
      y = "Density"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, "polymorphic_density_by_aminoacid.pdf"),
         p_density_aa, width = 14, height = 10)
  
  # -------------------------------------------------------------------------
  # Plot 4: Overlaid densities by family (single panel for comparison)
  # -------------------------------------------------------------------------
  # Calculate family means for annotation
  family_means <- raw_freq_with_null %>%
    group_by(Family, Neutral_Expectation) %>%
    summarise(Mean_Freq = mean(Preferred_Freq), .groups = "drop")
  
  p_overlay <- ggplot(raw_freq_with_null, aes(x = Preferred_Freq, fill = Family, color = Family)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    # Add neutral expectation lines
    geom_vline(data = neutral_expectations, 
               aes(xintercept = Neutral_Expectation, color = Family),
               linetype = "dashed", linewidth = 1, show.legend = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    labs(
      title = "Comparing CFS Across Degeneracy Families (Polymorphic Sites)",
      subtitle = "Dashed lines = neutral expectations (1/n codons)",
      x = "Preferred Codon Frequency",
      y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "polymorphic_density_overlay.pdf"),
         p_overlay, width = 10, height = 6)
  
  # -------------------------------------------------------------------------
  # Plot 5: Violin plot by amino acid with neutral expectation
  # -------------------------------------------------------------------------
  p_violin <- ggplot(raw_freq_with_null, 
                     aes(x = reorder(Amino_Acid, Preferred_Freq, median), 
                         y = Preferred_Freq, fill = Family)) +
    geom_violin(scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3, outlier.alpha = 0.3) +
    # Add horizontal lines for neutral expectations
    geom_hline(yintercept = 1/2, linetype = "dashed", color = "#66C2A5", alpha = 0.8) +
    geom_hline(yintercept = 1/3, linetype = "dashed", color = "#FC8D62", alpha = 0.8) +
    geom_hline(yintercept = 1/4, linetype = "dashed", color = "#8DA0CB", alpha = 0.8) +
    geom_hline(yintercept = 1/6, linetype = "dashed", color = "#E78AC3", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", name = "Degeneracy") +
    labs(
      title = "Preferred Codon Frequency by Amino Acid (Polymorphic Sites)",
      subtitle = "Dashed lines = neutral expectations (0.50, 0.33, 0.25, 0.17 for 2/3/4/6-fold)",
      x = "Amino Acid (ordered by median frequency)",
      y = "Preferred Codon Frequency"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold", size = 10),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, "polymorphic_violin_by_aminoacid.pdf"),
         p_violin, width = 14, height = 7)
  
  # -------------------------------------------------------------------------
  # Plot 6: Excess over neutral by family (showing selection strength)
  # -------------------------------------------------------------------------
  excess_data <- raw_freq_with_null %>%
    mutate(Excess = Preferred_Freq - Neutral_Expectation)
  
  p_excess <- ggplot(excess_data, aes(x = Excess, fill = Family)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    facet_wrap(~ Family, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Excess Preferred Codon Frequency Over Neutral Expectation",
      subtitle = "Values > 0 indicate selection FOR preferred codons",
      x = "Excess (Observed - Expected)",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(file.path(output_dir, "polymorphic_excess_over_neutral.pdf"),
         p_excess, width = 10, height = 8)
  
  cat("Polymorphic site plots saved to:", output_dir, "\n\n")
  
  return(list(
    density_overall = p_density_overall,
    density_family = p_density_family,
    density_aa = p_density_aa,
    density_overlay = p_overlay,
    violin = p_violin,
    excess = p_excess
  ))
}


# =============================================================================
# STATISTICAL TESTS FOR SELECTION
# =============================================================================

#' Test for departure from neutrality
#' 
#' Compare observed SFS to neutral expectation. Under neutrality,
#' polymorphic sites should follow a specific distribution based on
#' mutation-drift balance.
#' 
#' @param raw_freq Raw frequency data
#' @param full_dist Full distribution data
#' 
#' @return List with test results
test_selection_signature <- function(raw_freq, full_dist) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("TESTING FOR SELECTION SIGNATURES\n")
  cat(strrep("=", 70), "\n\n")
  
  results <- list()
  
  # 1. Summary statistics for polymorphic sites
  cat("=== Polymorphic Site Statistics ===\n\n")
  
  poly_summary <- raw_freq %>%
    group_by(Amino_Acid, Family) %>%
    summarise(
      n_sites = n(),
      mean_freq = mean(Preferred_Freq),
      median_freq = median(Preferred_Freq),
      sd_freq = sd(Preferred_Freq),
      skewness = (mean(Preferred_Freq) - median(Preferred_Freq)) / sd(Preferred_Freq),
      prop_high = mean(Preferred_Freq > 0.5),  # Proportion with preferred > 50%
      prop_very_high = mean(Preferred_Freq > 0.9),  # Proportion with preferred > 90%
      .groups = "drop"
    )
  
  print(poly_summary)
  results$polymorphic_summary <- poly_summary
  
  # 2. Fixation index: ratio of fixed preferred to fixed non-preferred
  cat("\n\n=== Fixation Index (Fixed Preferred / Fixed Non-Preferred) ===\n")
  cat("Values > 1 indicate net selection for preferred codons\n\n")
  
  fixation_index <- full_dist %>%
    mutate(
      Fixed_Pref = as.numeric(Bin_11_Fixed_Pref),
      Fixed_NonPref = as.numeric(Bin_0_Fixed_NonPref),
      Fixation_Index = Fixed_Pref / Fixed_NonPref,
      Log_Fixation_Index = log2(Fixation_Index)
    ) %>%
    dplyr::select(Amino_Acid, Family, Fixed_Pref, Fixed_NonPref, 
           Fixation_Index, Log_Fixation_Index)
  
  print(fixation_index)
  results$fixation_index <- fixation_index
  
 # 3. Asymmetry test: are frequencies skewed toward preferred?
  # IMPORTANT: Neutral expectation depends on family size!
  # 2-fold: 1/2 = 0.50, 3-fold: 1/3 = 0.333, 4-fold: 1/4 = 0.25, 6-fold: 1/6 = 0.167
  cat("\n\n=== Frequency Asymmetry Test (Family-Specific Null) ===\n")
  cat("H0: Mean preferred frequency = 1/n_codons (neutral expectation)\n")
  cat("H1: Mean preferred frequency > 1/n_codons (selection for preferred)\n\n")
  
  # Define neutral expectations by family
  neutral_expectations <- data.frame(
    Family = c("2-fold", "3-fold", "4-fold", "6-fold"),
    Neutral_Expectation = c(1/2, 1/3, 1/4, 1/6),
    N_Codons = c(2, 3, 4, 6)
  )
  
  cat("Neutral expectations by degeneracy:\n")
  print(neutral_expectations)
  cat("\n")
  
  # By family with correct null hypothesis
  asymmetry_by_family <- raw_freq %>%
    left_join(neutral_expectations, by = "Family") %>%
    group_by(Family, Neutral_Expectation, N_Codons) %>%
    summarise(
      n = n(),
      mean_freq = mean(Preferred_Freq),
      median_freq = median(Preferred_Freq),
      # Use family-specific neutral expectation
      t_stat = t.test(Preferred_Freq, mu = first(Neutral_Expectation), 
                      alternative = "greater")$statistic,
      p_value = t.test(Preferred_Freq, mu = first(Neutral_Expectation),
                       alternative = "greater")$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      Excess_Over_Neutral = mean_freq - Neutral_Expectation,
      Fold_Enrichment = mean_freq / Neutral_Expectation,
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05,
      interpretation = case_when(
        significant & Fold_Enrichment > 1 ~ "Selection FOR preferred",
        significant & Fold_Enrichment < 1 ~ "Selection AGAINST preferred",
        TRUE ~ "Not significant"
      )
    )
  
  cat("Results by degeneracy family:\n")
  print(asymmetry_by_family)
  results$asymmetry_by_family <- asymmetry_by_family
  
  # Overall weighted test (accounting for different expectations)
  cat("\n\n=== Overall Test (Excess Over Neutral) ===\n")
  
  # Calculate excess for each site relative to its family expectation
  excess_data <- raw_freq %>%
    left_join(neutral_expectations, by = "Family") %>%
    mutate(Excess = Preferred_Freq - Neutral_Expectation)
  
  overall_excess_test <- t.test(excess_data$Excess, mu = 0, alternative = "greater")
  
  cat("Testing if observed frequencies exceed neutral expectation:\n")
  cat("  Mean excess over neutral:", round(mean(excess_data$Excess), 4), "\n")
  cat("  t-statistic:", round(overall_excess_test$statistic, 3), "\n")
  cat("  p-value (one-tailed):", format.pval(overall_excess_test$p.value), "\n")
  cat("  Interpretation:", ifelse(overall_excess_test$p.value < 0.05,
                                  "Significant SELECTION FOR preferred codons",
                                  "No significant selection detected"), "\n")
  
  results$overall_excess_test <- overall_excess_test
  results$neutral_expectations <- neutral_expectations
  
  # 4. Polymorphism load: proportion of sites that are polymorphic
  cat("\n\n=== Polymorphism Load ===\n")
  
  poly_load <- full_dist %>%
    rowwise() %>%
    mutate(
      Total_Sites = sum(c_across(starts_with("Bin_")), na.rm = TRUE),
      Fixed_Sites = Bin_0_Fixed_NonPref + Bin_11_Fixed_Pref,
      Polymorphic_Sites = Total_Sites - Fixed_Sites,
      Prop_Polymorphic = Polymorphic_Sites / Total_Sites
    ) %>%
    ungroup() %>%
    dplyr::select(Amino_Acid, Family, Total_Sites, Fixed_Sites, 
           Polymorphic_Sites, Prop_Polymorphic)
  
  print(poly_load)
  results$polymorphism_load <- poly_load
  
  cat("\n\nMean polymorphism proportion:", 
      round(mean(poly_load$Prop_Polymorphic), 4), "\n")
  
  return(results)
}


#' Compare SFS between amino acid families
#' 
#' @param raw_freq Raw frequency data
#' 
#' @return List with comparison results
compare_cfs_by_family <- function(raw_freq) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("COMPARING CFS ACROSS DEGENERACY FAMILIES\n")
  cat(strrep("=", 70), "\n\n")
  
  results <- list()
  
  # Kruskal-Wallis test: do families differ in frequency distribution?
  kw_test <- kruskal.test(Preferred_Freq ~ Family, data = raw_freq)
  
  cat("Kruskal-Wallis test:\n")
  cat("  H0: All families have the same frequency distribution\n")
  cat("  Chi-squared:", round(kw_test$statistic, 3), "\n")
  cat("  df:", kw_test$parameter, "\n")
  cat("  p-value:", format.pval(kw_test$p.value), "\n\n")
  
  results$kruskal_wallis <- kw_test
  
  # Pairwise Wilcoxon tests
  if (kw_test$p.value < 0.05) {
    cat("Pairwise Wilcoxon tests (BH-adjusted):\n\n")
    
    pairwise_results <- pairwise.wilcox.test(
      raw_freq$Preferred_Freq,
      raw_freq$Family,
      p.adjust.method = "BH"
    )
    
    print(pairwise_results)
    results$pairwise_wilcox <- pairwise_results
  }
  
  # Effect sizes (median differences)
  cat("\n\nMedian frequencies by family:\n")
  median_by_family <- raw_freq %>%
    group_by(Family) %>%
    summarise(
      n = n(),
      median = median(Preferred_Freq),
      IQR = IQR(Preferred_Freq),
      .groups = "drop"
    ) %>%
    arrange(desc(median))
  
  print(median_by_family)
  results$median_by_family <- median_by_family
  
  return(results)
}


# =============================================================================
# SELECTION COEFFICIENT ESTIMATION
# =============================================================================

#' Estimate selection coefficient from SFS
#' 
#' Uses the ratio of fixed to polymorphic sites to estimate
#' the population-scaled selection coefficient (gamma = 2Ns)
#' 
#' @param full_dist Full distribution data
#' 
#' @return Data frame with selection estimates by amino acid
estimate_selection_from_cfs <- function(full_dist) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("ESTIMATING SELECTION FROM CODON FREQUENCY SPECTRUM\n")
  cat(strrep("=", 70), "\n\n")
  
  cat("Method: Comparing observed fixation bias to neutral expectation\n")
  cat("Under neutrality, fixed_pref / fixed_nonpref should equal 1/(n-1)\n")
  cat("where n = number of synonymous codons in the family\n")
  cat("Selection coefficient gamma = 2Ns can be estimated from deviation\n\n")
  
  selection_estimates <- full_dist %>%
    mutate(
      Fixed_Pref = as.numeric(Bin_11_Fixed_Pref),
      Fixed_NonPref = as.numeric(Bin_0_Fixed_NonPref)
    ) %>%
    rowwise() %>%
    mutate(
      # Total sites = sum of all bins (columns 3 onwards are the bin columns)
      Total = sum(c_across(3:14), na.rm = TRUE),
      # Polymorphic = Total - Fixed sites
      Polymorphic = Total - Fixed_Pref - Fixed_NonPref,
      
      # Fixation probability ratio
      Fix_Ratio = Fixed_Pref / Fixed_NonPref,
      Log_Fix_Ratio = log(Fix_Ratio),
      
      # Simple gamma estimate: gamma ~ ln(fixation_ratio)
      # This is a rough approximation; true estimate requires more complex modeling
      Gamma_Estimate = Log_Fix_Ratio,
      
      # Direction of selection
      Selection_Direction = case_when(
        Fix_Ratio > 1.5 ~ "Strong for Preferred",
        Fix_Ratio > 1.1 ~ "Weak for Preferred",
        Fix_Ratio > 0.9 ~ "Neutral",
        Fix_Ratio > 0.6 ~ "Weak against Preferred",
        TRUE ~ "Strong against Preferred"
      )
    ) %>%
    ungroup() %>%
    dplyr::select(Amino_Acid, Family, Fixed_Pref, Fixed_NonPref, Polymorphic,
           Fix_Ratio, Gamma_Estimate, Selection_Direction) %>%
    arrange(desc(Gamma_Estimate))
  
  print(selection_estimates)
  
  cat("\n\nInterpretation:\n")
  cat("  Gamma > 0: Selection favors preferred codons\n")
  cat("  Gamma < 0: Selection favors non-preferred codons\n")
  cat("  |Gamma| > 1: Strong selection (|2Ns| > 1)\n\n")
  
  return(selection_estimates)
}


# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

#' Run complete CFS analysis pipeline
#' 
#' @param data_dir Directory containing input files
#' @param output_dir Directory for output files
#' 
#' @return List with all analysis results
run_cfs_analysis <- function(
    data_dir = "./data",
    output_dir = "./results/cfs_analysis"
) {
  
  cat("\n")
  cat(strrep("#", 80), "\n")
  cat("CODON FREQUENCY SPECTRUM (CFS) ANALYSIS PIPELINE\n")
  cat(strrep("#", 80), "\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list()
  
  # 1. Load data
  data <- load_cfs_data(data_dir = data_dir)
  results$data <- data
  
  # 2. Create visualizations
  results$bridge_plots <- plot_suspension_bridge(
    data$binned_distribution,
    by_family = TRUE,
    output_dir = output_dir
  )
  
  results$polymorphic_plots <- plot_polymorphic_frequency_distributions(
    data$raw_frequencies,
    output_dir = output_dir
  )
  
  # 3. Statistical tests
  results$selection_tests <- test_selection_signature(
    data$raw_frequencies,
    data$binned_distribution
  )
  
  results$family_comparison <- compare_cfs_by_family(data$raw_frequencies)
  
  # 4. Selection coefficient estimation
  results$selection_estimates <- estimate_selection_from_cfs(data$binned_distribution)
  
  # 5. Save results
  write.csv(results$selection_tests$polymorphic_summary,
            file.path(output_dir, "polymorphic_site_summary.csv"),
            row.names = FALSE)
  
  write.csv(results$selection_tests$fixation_index,
            file.path(output_dir, "fixation_index.csv"),
            row.names = FALSE)
  
  write.csv(results$selection_estimates,
            file.path(output_dir, "selection_estimates.csv"),
            row.names = FALSE)
  
  # Summary
  cat("\n")
  cat(strrep("#", 80), "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(strrep("#", 80), "\n\n")
  
  cat("Key findings:\n")
  cat("  - Overall mean preferred frequency:", 
      round(mean(data$raw_frequencies$Preferred_Freq), 4), "\n")
  cat("  - Overall test for directional selection: p =",
      format.pval(results$selection_tests$overall_asymmetry$p.value), "\n")
  cat("  - Mean fixation index (Pref/NonPref):",
      round(mean(results$selection_tests$fixation_index$Fixation_Index), 3), "\n")
  
  cat("\nResults saved to:", output_dir, "\n\n")
  
  return(results)
}


# =============================================================================
# EXPRESSION-STRATIFIED SFS (if expression data available)
# =============================================================================

#' Analyze SFS stratified by expression level
#' 
#' @param sfs_data SFS data from load_sfs_data()
#' @param gene_expression Data frame with Gene_name and expression values
#' @param gene_codon_preference Data linking genes to codon preferences
#' 
#' @return List with stratified analysis results
analyze_sfs_by_expression <- function(
    sfs_data,
    gene_expression,
    gene_codon_preference = NULL
) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("EXPRESSION-STRATIFIED SFS ANALYSIS\n")
  cat(strrep("=", 70), "\n\n")
  
  cat("Note: This analysis requires per-gene codon frequency data\n")
  cat("that links polymorphic sites to specific genes.\n")
  cat("Currently using aggregated SFS data.\n\n")
  
  # If we had per-gene data, we could:
  # 1. Split genes into expression quantiles
  # 2. Compare SFS shape between high vs low expression genes
  # 3. Test if preferred codon frequency correlates with expression
  
  # For now, report that this analysis needs additional data
  cat("To perform expression-stratified SFS analysis,\n")
  cat("you would need per-gene codon frequency data from VCF processing.\n")
  
  return(NULL)
}
