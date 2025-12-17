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
  
  # Return segregating sites only (1 to target_n-1)
  return(proj_sfs[2:target_n]) 
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
    k_values <- 1:(target_n - 1)
    
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