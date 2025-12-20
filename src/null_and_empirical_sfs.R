project_sfs <- function(counts_df, target_n) {
  #' Projects SFS summary table down to a common sample size
  #' Input:
  #'   counts_df: Data frame with columns 'k', 'n', and 'count'
  #'   target_n: The sample size to project down to
  
  # Initialize empty SFS (indices 0 to target_n)
  proj_sfs <- numeric(target_n + 1) 
  
  # Filter sites that have enough samples to be projected
  valid_sites <- counts_df[counts_df$n >= target_n, ]
  
  for(i in 0:target_n) {
    # 1. Calculate prob of site mapping to 'i' (vectorized across all rows)
    probs <- dhyper(x = i, 
                    m = valid_sites$k, 
                    n = valid_sites$n - valid_sites$k, 
                    k = target_n)
    
    # 2. WEIGHT BY THE COUNT COLUMN (Crucial Fix)
    # This accounts for the fact that one row in your CSV represents multiple sites
    weighted_probs <- probs * valid_sites$count
    
    # Sum to get total sites in bin 'i'
    proj_sfs[i + 1] <- sum(weighted_probs)
  }
  
  # Return ALL bins (0 to target_n)
  return(proj_sfs) 
}

# --- 2. Corrected Expectation Function ---
generate_expected_counts <- function(neutral_param,
                                     observed_sfs_list, # New argument
                                     target_n = 90) {
  #' Generates expected neutral counts scaled to match the Observed SFS size
  #' 
  #' @param neutral_param List with alpha/beta parameters
  #' @param observed_sfs_list A list containing the projected OBSERVED SFS 
  #'        (e.g., list(G = vec_G, C = vec_C)). Used to get the total site count.
  #' @param target_n The projection size used
  
  expectations <- list()
  targets <- c("G", "C")
  
  for (t in targets) {
    # 1. Get Parameters
    alpha <- neutral_param[[paste0("alpha_", t)]]
    beta <- neutral_param[[paste0("beta_", t)]]
    
    # 2. Generate Probabilities (The Shape)
    k_values <- 0:(target_n)
    
    # Log-space calculation for stability
    log_probs <- lchoose(target_n, k_values) + 
      lbeta(k_values + alpha, target_n - k_values + beta) - 
      lbeta(alpha, beta)
    
    probs <- exp(log_probs)
    
    # Normalize probabilities to sum to 1 (conditional on segregating)
    probs_norm <- probs / sum(probs)
    
    # 3. SCALE BY OBSERVED TOTAL (Crucial Fix)
    # We multiply the neutral probability by the total number of variants 
    # in the observed data so the bars are comparable in the plot.
    total_observed_sites <- sum(observed_sfs_list[[t]])
    
    expectations[[t]] <- probs_norm * total_observed_sites
  }
  
  return(expectations)
}

validate_selection_signal <- function(plot_data, target_type) {
  
  # 1. Extract Data for the specific Target (C or G)
  # ------------------------------------------------
  target_data <- plot_data |>
    filter(grepl(target_type, Metric))
  
  obs_data <- target_data |> filter(grepl("Observed", Metric))
  exp_data <- target_data |> filter(grepl("Expected", Metric))
  
  # 2. Check Sample Size (Effective Number of Variants)
  # ---------------------------------------------------
  # The sum of the projected SFS is the total count of segregating sites
  # used in the analysis.
  n_variants <- sum(obs_data$Counts)
  
  message(paste0("\n=== Validation for ", target_type, "-ending Codons ==="))
  message(paste("Total Segregating Sites (Projected):", round(n_variants, 1)))
  
  # Heuristic Check: Do we have enough data?
  if(n_variants < 100) {
    warning("Sample size is low (< 100 sites). The Chi-Squared test may lack power.")
  } else {
    message("Sample size is sufficient for robust statistical testing.")
  }
  
  # 3. Chi-Squared Goodness of Fit Test
  # -----------------------------------
  # H0: The Observed SFS follows the Neutral Expectation
  # H1: The Observed SFS differs significantly (i.e., Selection)
  
  # We assume the 'Expected' counts are the theoretical values.
  # We compare 'Observed' against them.
  # Note: Since projection gives non-integers, we round for the strict chisq.test,
  # or calculate manually. Manual calculation is safer here.
  
  O <- obs_data$Counts
  E <- exp_data$Counts
  
  # Avoid division by zero in rare cases where E is tiny
  valid_bins <- E > 1e-6
  O <- O[valid_bins]
  E <- E[valid_bins]
  
  # Calculate Chi-Square Statistic
  chi_sq_stat <- sum((O - E)^2 / E)
  
  # Degrees of Freedom = (Number of Bins) - 1
  df <- length(O) - 1
  
  # P-value
  p_val <- pchisq(chi_sq_stat, df, lower.tail = FALSE)
  
  message(paste("Chi-Squared Statistic:", round(chi_sq_stat, 2)))
  message(paste("Degrees of Freedom:", df))
  message(paste("P-value:", format.pval(p_val, digits = 4)))
  
  if(p_val < 0.05) {
    message("RESULT: REJECT Null Hypothesis. Significant deviation from neutrality detected.")
  } else {
    message("RESULT: CANNOT Reject Null Hypothesis. Observed deviation may be noise.")
  }
  
  # 4. Quantify the "Shift" (Difference in Mean Allele Frequency)
  # -------------------------------------------------------------
  # Mean Frequency = Sum(Count * Freq) / Total_Count
  # Freq index is roughly 'Num_seq' / 90
  
  freqs <- obs_data$Num_seq / max(obs_data$Num_seq + 1) # approx frequency
  
  mean_freq_obs <- sum(obs_data$Counts * freqs) / sum(obs_data$Counts)
  mean_freq_exp <- sum(exp_data$Counts * freqs) / sum(exp_data$Counts)
  
  message(paste("Mean Allele Freq (Observed):", round(mean_freq_obs, 4)))
  message(paste("Mean Allele Freq (Neutral): ", round(mean_freq_exp, 4)))
  message(paste("Direction of Shift:", ifelse(mean_freq_obs > mean_freq_exp, 
                                              "Towards Fixation (Positive Selection)", 
                                              "Towards Loss (Negative Selection)")))
}

estimate_gamma <- function(observed_sfs, neutral_param, target_type, target_n = 90) {
  #' Estimate gamma (4Nes) from observed SFS using maximum likelihood
  #' 
  #' @param observed_sfs Numeric vector of observed SFS (length = target_n + 1)
  #' @param neutral_param List with alpha/beta parameters
  #' @param target_type "G" or "C" to select appropriate alpha/beta
  #' @param target_n Sample size used for projection
  #' @return List with gamma_mle, log_likelihood, and convergence info
  
  # Get neutral parameters
  alpha <- neutral_param[[paste0("alpha_", target_type)]]
  beta <- neutral_param[[paste0("beta_", target_type)]]
  
  # Define log-likelihood function (all calculations in log-space for numerical stability)
  log_likelihood <- function(gamma) {
    k_values <- 0:target_n  # All frequency bins including fixed sites
    
    # Log selection weight: gamma * k (in log space)
    log_selection_weight <- gamma * k_values
    
    # Neutral log-probability (beta-binomial)
    log_neutral_prob <- lchoose(target_n, k_values) + 
      lbeta(k_values + alpha, target_n - k_values + beta) - 
      lbeta(alpha, beta)
    
    # Combine in log space: log(P(k)) = gamma*k + log(BetaBin(k))
    log_prob_k_unnorm <- log_selection_weight + log_neutral_prob
    
    # Normalize using log-sum-exp trick for numerical stability
    max_log_prob <- max(log_prob_k_unnorm)
    log_prob_k <- log_prob_k_unnorm - max_log_prob - log(sum(exp(log_prob_k_unnorm - max_log_prob)))
    
    # Log-likelihood (observed counts are weights)
    ll <- sum(observed_sfs * log_prob_k)
    
    # Return -Inf if numerical issues
    if (!is.finite(ll)) {
      return(-Inf)
    }
    
    return(ll)
  }
  
  # Optimize to find MLE of gamma
  # Search range: gamma ∈ [-10, 50] (covers strong negative to strong positive selection)
  result <- tryCatch(
    optimize(
      f = log_likelihood,
      interval = c(-10, 50),
      maximum = TRUE  # We want to maximize log-likelihood
    ),
    error = function(e) {
      warning("Optimization failed: ", e$message)
      return(list(maximum = 0, objective = log_likelihood(0)))
    }
  )
  
  gamma_mle <- result$maximum
  max_ll <- result$objective
  
  # Also calculate likelihood at gamma = 0 (neutral)
  ll_neutral <- log_likelihood(0)
  
  # Likelihood ratio test statistic
  lr_statistic <- 2 * (max_ll - ll_neutral)
  
  return(list(
    gamma = gamma_mle,
    log_likelihood = max_ll,
    log_likelihood_neutral = ll_neutral,
    lr_statistic = lr_statistic,
    p_value = pchisq(lr_statistic, df = 1, lower.tail = FALSE)
  ))
}

plot_sfs_validation_4panel <- function(plot_data_top5, plot_data_bottom5, 
                                       gamma_top5_G, gamma_top5_C,
                                       gamma_bottom5_G, gamma_bottom5_C) {
  #' Create 4-panel validation plot with gamma estimates
  #' 
  #' @param plot_data_top5 Plot data for top 5% genes (contains G and C)
  #' @param plot_data_bottom5 Plot data for bottom 5% genes (contains G and C)
  #' @param gamma_top5_G Gamma estimate for top 5% G-ending codons
  #' @param gamma_top5_C Gamma estimate for top 5% C-ending codons
  #' @param gamma_bottom5_G Gamma estimate for bottom 5% G-ending codons
  #' @param gamma_bottom5_C Gamma estimate for bottom 5% C-ending codons
  
  require(ggplot2)
  require(gridExtra)
  
  # Helper function to create individual panel
  create_panel <- function(data, target_type, gamma_result, title_prefix) {
    target_data <- data |>
      filter(grepl(target_type, Metric))
    
    # Format gamma annotation (use 'gamma' instead of γ to avoid encoding issues)
    gamma_text <- sprintf(
      "gamma = %.3f (p = %.2e)\n%s",
      gamma_result$gamma,
      gamma_result$p_value,
      ifelse(gamma_result$p_value < 0.05, 
             ifelse(gamma_result$gamma > 0, "Positive Selection", "Negative Selection"),
             "Neutral")
    )
    
    # Separate observed and expected data for plotting
    obs_data <- target_data %>% filter(grepl("Observed", Metric))
    exp_data <- target_data %>% filter(grepl("Expected", Metric))
    
    p <- ggplot() +
      # Observed data as bars
      geom_bar(data = obs_data, aes(x = Num_seq, y = Counts, fill = "Observed"),
               stat = "identity", alpha = 0.7) +
      # Neutral expectation as a line overlay
      geom_line(data = exp_data, aes(x = Num_seq, y = Counts, color = "Neutral"),
                linewidth = 1.2) +
      geom_point(data = exp_data, aes(x = Num_seq, y = Counts, color = "Neutral"),
                 size = 1, alpha = 0.6) +
      scale_fill_manual(
        name = NULL,
        values = c("Observed" = "#E41A1C"),
        labels = c("Observed SFS")
      ) +
      scale_color_manual(
        name = NULL,
        values = c("Neutral" = "#377EB8"),
        labels = c("Neutral Expectation")
      ) +
      scale_y_log10() +  # Log scale to visualize both fixed (k=0, k=90) and polymorphic sites
      labs(
        title = paste0(title_prefix, ": ", target_type, "-ending Codons"),
        x = "Allele Count (k)",
        y = "Number of Sites (log scale)",
        fill = NULL
      ) +
      theme_custom() +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "top",
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(order = 1), 
             color = guide_legend(order = 2)) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = gamma_text,
        hjust = 1.1, vjust = 1.5,
        size = 3.5,
        fontface = "bold",
        color = ifelse(gamma_result$p_value < 0.05, "#E41A1C", "#666666")
      )
    
    return(p)
  }
  
  # Create 4 panels
  p1 <- create_panel(plot_data_top5, "G", gamma_top5_G, "Top 5% (High Expression)")
  p2 <- create_panel(plot_data_top5, "C", gamma_top5_C, "Top 5% (High Expression)")
  p3 <- create_panel(plot_data_bottom5, "G", gamma_bottom5_G, "Bottom 5% (Low Expression)")
  p4 <- create_panel(plot_data_bottom5, "C", gamma_bottom5_C, "Bottom 5% (Low Expression)")
  
  # Combine into 2x2 grid
  combined_plot <- gridExtra::grid.arrange(
    p1, p2, p3, p4,
    ncol = 2,
    nrow = 2,
    top = grid::textGrob("SFS Validation: Selection Signature in Preferred Codons",
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  return(combined_plot)
}

process_gene_set_sfs <- function(gene_list, gene_set_name, target_n = 90) {
  #' Process codon frequency data for a specific gene set
  #' 
  #' @param gene_list Character vector of gene names
  #' @param gene_set_name Description for logging (e.g., "Top 5%")
  #' @param target_n Sample size for SFS projection
  #' @return List with obs_sfs_G and obs_sfs_C
  
  cat(sprintf("\n=== Processing %s Genes ===\n", gene_set_name))
  cat(sprintf("%s genes: %d\n", gene_set_name, length(gene_list)))
  
  # Load and filter codon frequency data
  raw_variants <- read.delim("./data/all_chromosomes.codon_frequencies_preferred.txt", 
                             stringsAsFactors = FALSE) |>
    dplyr::mutate(Gene = paste0("MgIM767.", Gene)) |>
    dplyr::filter(Gene %in% gene_list)
  
  cat(sprintf("Codon positions in %s: %d\n", gene_set_name, nrow(raw_variants)))
  
  # Parse the 'Codon_Variants' column
  parsed_sfs_data <- raw_variants |>
    rowwise() |>
    dplyr::mutate(
      # A. Identify Target Base (G or C)
      Target_Base = substring(Preferred_Codon, 3, 3),
      
      # B. Parse Counts from string "AAA:100;AAG:2"
      entries = list(str_split(Codon_Variants, ";")[[1]]),
      
      # Extract counts for all variants at this site
      all_counts = list(as.numeric(sub(".*:", "", entries))),
      all_codons = list(sub(":.*", "", entries)),
      
      # Calculate Total Depth (n)
      n = sum(all_counts),
      
      # Calculate Count of Preferred Allele (k)
      k_index = match(Preferred_Codon, all_codons),
      k = ifelse(is.na(k_index), 0, all_counts[k_index])
    ) |>
    ungroup() |>
    # Filter for valid sites (n > 0) and only G/C targets
    dplyr::filter(n > 0, Target_Base %in% c("G", "C")) |>
    dplyr::select(Gene, Target_Base, n, k)
  
  # Project SFS for G-ending sites
  syn_G_summary <- parsed_sfs_data |>
    dplyr::filter(Target_Base == "G") |>
    dplyr::count(n, k, name = "count")
  
  obs_sfs_G <- project_sfs(syn_G_summary, target_n)
  
  # Project SFS for C-ending sites
  syn_C_summary <- parsed_sfs_data |>
    dplyr::filter(Target_Base == "C") |>
    dplyr::count(n, k, name = "count")
  
  obs_sfs_C <- project_sfs(syn_C_summary, target_n)
  
  return(list(obs_sfs_G = obs_sfs_G, obs_sfs_C = obs_sfs_C))
}