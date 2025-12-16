#' Auxiliary functions to estimate gamma from polymorphism data
#' 
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

# ******************************************************************************
# 1) Define the Probability Function (Wright's Distribution integrated) ----
# ______________________________________________________________________________

get_prob_k <- function(k, n, u, v, S, theta = 0.0312) {
  #' Calculates probability using internally scaled mutation rates
  #' theta: The observed population mutation rate (Pi from synonymous sites)
  
  # --- INTERNAL SCALING FIX ---
  # Force alpha + beta to equal the observed theta (0.0312)
  # Preserving the bias ratio (u vs v) from the Q-matrix
  rate_sum <- u + v
  alpha <- theta * (u / rate_sum)
  beta  <- theta * (v / rate_sum)
  
  # Define the unnormalized Wright density
  wright_density <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    # Using log scale for numerical stability
    val <- (alpha - 1)*log(x) + (beta - 1)*log(1-x) + S*x
    return(exp(val))
  }
  
  # Normalization Constant (Denominator)
  denom <- tryCatch(
    integrate(wright_density, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(denom) || denom == 0) return(0)
  
  # Numerator (Probability of sampling k alleles)
  numerator_func <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    # dbinom handles the combinatorial math
    return(dbinom(k, size=n, prob=x) * wright_density(x))
  }
  
  num <- tryCatch(
    integrate(numerator_func, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(num)) return(0)
  
  return(num / denom)
}

get_prob_k_analytical <- function(k, n, u, v, S, theta = 0.0312) {
  # u and v come from Q-matrix (per-generation rates)
  # theta = 4*N*mu_total where mu_total is the total mutation rate
  # We scale u and v by theta to get 4*N*u and 4*N*v
  
  rate_sum <- u + v
  alpha <- theta * (u / rate_sum)  # This is 4*N*u
  beta  <- theta * (v / rate_sum)  # This is 4*N*v
  
  A_prime <- k + alpha
  B_prime <- n - k + beta
  
  log_binom <- lchoose(n, k)
  log_beta_ratio <- lbeta(A_prime, B_prime) - lbeta(alpha, beta)
  log_hyperg_num <- log(gsl::hyperg_1F1(A_prime, A_prime + B_prime, S))
  log_hyperg_den <- log(gsl::hyperg_1F1(alpha, alpha + beta, S))
  
  log_prob <- log_binom + log_beta_ratio + (log_hyperg_num - log_hyperg_den)
  
  return(exp(log_prob))
}

get_prob_k_analytical_precomputed <- function(k, n, alpha, beta, S) {
  #' Calculate probability using pre-computed alpha and beta from introns
  #' This version directly uses empirically-derived 4N*u and 4N*v values
  #' 
  #' @param k Number of preferred alleles observed
  #' @param n Total sample size
  #' @param alpha Pre-computed 4N*u from intronic G or C sites
  #' @param beta Pre-computed 4N*v from intronic G or C sites
  #' @param S Selection coefficient (4Nes) to evaluate
  #' @return Probability of observing k preferred alleles
  
  A_prime <- k + alpha
  B_prime <- n - k + beta
  
  log_binom <- lchoose(n, k)
  log_beta_ratio <- lbeta(A_prime, B_prime) - lbeta(alpha, beta)
  log_hyperg_num <- log(gsl::hyperg_1F1(A_prime, A_prime + B_prime, S))
  log_hyperg_den <- log(gsl::hyperg_1F1(alpha, alpha + beta, S))
  
  log_prob <- log_binom + log_beta_ratio + (log_hyperg_num - log_hyperg_den)
  
  return(exp(log_prob))
}

solve_alpha_and_beta_from_introns <- function(sfs_file) {
  #' Estimate alpha and beta from neutral (intronic) site frequency spectra
  #' Uses Weighted Beta-Binomial Maximum Likelihood Estimation
  #' 
  #' @param sfs_file Path to CSV file with columns: n, k, count
  #' @return A list containing estimated alpha, beta, and convergence codes
  #' ___________________________________________________________________________
  
  require(stats)
  require(data.table)
  
  # Load SFS data
  sfs_data <- fread(sfs_file)
  
  # Filter out invalid rows if any
  dt <- sfs_data[n > 0]
  
  # Extract vectors for vectorized calculation
  k <- dt$k
  n <- dt$n
  w <- dt$count # The weights
  
  # Negative Log-Likelihood Function (Weighted)
  # L = Sum [ count * log( P(k|n, alpha, beta) ) ]
  nll <- function(params) {
    log_alpha <- params[1]
    log_beta <- params[2]
    
    alpha <- exp(log_alpha)
    beta <- exp(log_beta)
    
    # Beta-Binomial Log-Probability
    # const + lbeta(k+a, n-k+b) - lbeta(a, b)
    log_prob <- lchoose(n, k) + 
      lbeta(k + alpha, n - k + beta) - 
      lbeta(alpha, beta)
    
    # WEIGHTED sum (Critical optimization)
    return(-sum(w * log_prob))
  }
  
  # Initial guesses based on weighted mean frequency
  total_alleles <- sum(as.numeric(w) * n)
  total_k <- sum(as.numeric(w) * k)
  p_hat <- total_k / total_alleles
  
  start_alpha <- 0.01 * p_hat
  start_beta <- 0.01 * (1 - p_hat)
  
  # Optimize
  opt <- optim(par = c(log(start_alpha), log(start_beta)), 
               fn = nll, 
               method = "Nelder-Mead")
  
  results <- list(
    alpha = exp(opt$par[1]),
    beta = exp(opt$par[2]),
    convergence = opt$convergence,
    n_sites = sum(as.numeric(w)) # Total number of sites processed
  )
  
  return(results)
}

# ******************************************************************************
# 2) Likelihood Optimizer ----
# ______________________________________________________________________________

# Input: 
#   counts: Vector of observed counts of PREFERRED codons (k) across many sites
#   sample_sizes: Vector of sample sizes (n) for those sites (usually 187, because we have haplotypes)
#   u, v: Mutation rates for this specific Amino Acid family

estimate_gamma_for_AA <- function(counts, sample_sizes, alpha, beta, 
                                  S_interval = c(0, 50)) {
  #' Estimate gamma (4Nes) using pre-computed alpha and beta from introns
  #' Constrained to POSITIVE values only (gamma >= 0) to match AnaCoDa framework
  #' where preferred codon is optimal and gamma measures selection favorability
  #' 
  #' @param counts Vector of preferred codon counts (k)
  #' @param sample_sizes Vector of total sample sizes (n)
  #' @param alpha Pre-computed 4N*u (unpreferred -> preferred mutation rate)
  #' @param beta Pre-computed 4N*v (preferred -> unpreferred mutation rate)
  #' @param S_interval Search interval for gamma (default: [0, 50] for positive selection)
  #' @return Estimated gamma value (gamma >= 0)
  
  nll <- function(S) {
    log_probs <- mapply(function(k, n) {
      # Use analytical formula with pre-computed alpha, beta
      p <- get_prob_k_analytical_precomputed(k, n, alpha, beta, S)
      if(is.na(p) || p <= 0) return(-1e6) 
      return(log(p))
    }, counts, sample_sizes)
    
    return(-sum(log_probs))
  }
  
  opt <- optimize(nll, interval = S_interval)
  return(opt$minimum)
}

worker_gamma_corrected <- function(k, n, u, v, theta) {
  # FIX: Changed threshold from 5 to 1 to unlock 85% of data
  # Most genes have only ~2 sites per AA, so requiring 5 was too strict
  if (length(k) < 1) return(NA_real_)
  
  tryCatch(
    estimate_gamma_for_AA(
      counts = k, 
      sample_sizes = n,
      alpha = theta * (u / (u + v)),
      beta = theta * (v / (u + v)),
      S_interval = c(0, 50)  # Bound to positive selection
    ),
    error = function(e) NA_real_
  )
}

# ******************************************************************************
# 3) Calculate Expected Nucleotide Diversity (Pi) ----
# ______________________________________________________________________________

calculate_expected_pi <- function(alpha, beta, S) {
  #' Calculates the expected Nucleotide Diversity (Pi) for a population
  #' by integrating 2p(1-p) over the Wright distribution.
  #' 
  #' alpha: 4Neu (Rate Unpreferred -> Preferred)
  #' beta:  4Nev (Rate Preferred -> Unpreferred)
  #' S:     4Nes (Selection coefficient)
  #' ___________________________________________________________________________
  
  # 1. Define the Density Function (Same as your get_prob_k)
  wright_density <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    # Log-scale calculation for stability
    val <- (alpha - 1)*log(x) + (beta - 1)*log(1-x) + S*x
    return(exp(val))
  }
  
  # 2. Calculate Normalization Constant (The Denominator)
  # Integral of phi(x)
  denom <- tryCatch(
    integrate(wright_density, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(denom) || denom == 0) return(0)
  
  # 3. Calculate Numerator (Expected Heterozygosity)
  # Integral of 2*x*(1-x) * phi(x)
  numerator_func <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    heterozygosity <- 2 * x * (1 - x)
    return(heterozygosity * wright_density(x))
  }
  
  num <- tryCatch(
    integrate(numerator_func, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(num)) return(0)
  
  return(num / denom)
}

calculate_expected_pi_robust <- function(alpha, beta, S) {
  
  # Define the integrand: 2p(1-p) * exp(S*p) * p^(a-1) * (1-p)^(b-1)
  # This looks exactly like an unnormalized Beta distribution weighted by exp(S*p)
  
  # We use the fact that:
  # Integral[ x^(A-1) (1-x)^(B-1) ] is Beta(A,B)
  
  # We can pull out the '2' and integrating x(1-x) * density becomes:
  # Integral [ 2 * x^alpha * (1-x)^beta * exp(S*x) ]
  
  # Define the function for the Numerator (Heterozygosity)
  num_func <- function(x) {
    # Using log-exp to avoid underflow
    log_val <- log(2) + alpha*log(x) + beta*log(1-x) + S*x
    return(exp(log_val))
  }
  
  # Define the function for the Denominator (Normalization)
  den_func <- function(x) {
    log_val <- (alpha-1)*log(x) + (beta-1)*log(1-x) + S*x
    return(exp(log_val))
  }
  
  # Integrate with safer bounds (avoiding exact 0 and 1)
  # Using 1e-8 avoids the singularity at the edge
  num <- integrate(num_func, 1e-8, 1-1e-8)$value
  den <- integrate(den_func, 1e-8, 1-1e-8)$value
  
  return(num / den)
}

#' Calculate nucleotide diversity (pi) at synonymous sites
#' 
#' This function calculates TRUE nucleotide diversity by examining
#' nucleotide-level variation at each codon position, rather than
#' treating whole codons as alleles.
#' 
#' @author Luis Javier Madrigal-Roca
#' _____________________________________________________________________________

calculate_nucleotide_pi_at_synonymous_sites <- function(vcf_dt, genetic_code_df) {
  #' Calculate nucleotide-level pi at synonymous sites
  #' 
  #' For each codon position, we examine which nucleotide positions
  #' are synonymous (don't change the amino acid) and calculate
  #' nucleotide diversity at those positions.
  #' 
  #' @param vcf_dt data.table with columns: Gene, Codon_Pos, AA, 
  #'               Preferred_Codon, Codon_Variants
  #' @param genetic_code_df data.frame mapping codons to amino acids
  #' @return data.table with nucleotide-level pi at each synonymous position
  
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  
  # Create codon-to-AA lookup
  setkey(genetic_code_df, Codon)
  
  # Expand variants
  dt <- vcf_dt[, .(
    Variant_String = unlist(strsplit(Codon_Variants, ";", fixed = TRUE))
  ), by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  # Parse codon:count
  dt[, c("Variant_Codon", "Count") := tstrsplit(Variant_String, ":", fixed=TRUE)]
  dt[, Count := as.integer(Count)]
  dt[, Variant_String := NULL]
  
  # Get AA for each variant
  dt <- merge(dt, genetic_code_df, by.x = "Variant_Codon", by.y = "Codon", all.x = TRUE)
  setnames(dt, "AA.y", "Variant_AA")
  setnames(dt, "AA.x", "Site_AA")
  
  # Keep only synonymous variants
  dt_syn <- dt[Site_AA == Variant_AA]
  
  # For each site, calculate nucleotide diversity at each of the 3 positions
  result <- dt_syn[, {
    
    # Get all variant codons and their counts
    codons <- Variant_Codon
    counts <- Count
    total_n <- sum(counts)
    
    if (length(codons) < 2 || total_n < 2) {
      # No variation or insufficient sample
      list(
        Pos1_Pi = 0, Pos2_Pi = 0, Pos3_Pi = 0,
        Pos1_Syn = FALSE, Pos2_Syn = FALSE, Pos3_Syn = FALSE,
        n = total_n,
        Mean_Syn_Pi = 0,
        N_Syn_Positions = 0
      )
    } else {
      
      # Check each position
      pi_vals <- numeric(3)
      is_syn <- logical(3)
      
      for (pos in 1:3) {
        # Extract nucleotides at this position
        nucs <- substr(codons, pos, pos)
        
        # Calculate nucleotide frequencies
        nuc_counts <- tapply(counts, nucs, sum)
        nuc_freqs <- nuc_counts / total_n
        
        # Check if this position is synonymous (has variation)
        is_syn[pos] <- length(unique(nucs)) > 1
        
        # Calculate pi: (n/(n-1)) * (1 - sum(p_i^2))
        if (is_syn[pos]) {
          pi_vals[pos] <- (total_n / (total_n - 1)) * (1 - sum(nuc_freqs^2))
        } else {
          pi_vals[pos] <- 0
        }
      }
      
      # Calculate mean pi across synonymous positions only
      n_syn_pos <- sum(is_syn)
      mean_pi <- if (n_syn_pos > 0) mean(pi_vals[is_syn]) else 0
      
      list(
        Pos1_Pi = pi_vals[1],
        Pos2_Pi = pi_vals[2],
        Pos3_Pi = pi_vals[3],
        Pos1_Syn = is_syn[1],
        Pos2_Syn = is_syn[2],
        Pos3_Syn = is_syn[3],
        n = total_n,
        Mean_Syn_Pi = mean_pi,
        N_Syn_Positions = n_syn_pos
      )
    }
  }, by = .(Gene, Codon_Pos, Site_AA, Preferred_Codon)]
  
  return(result)
}


process_codon_vcf_with_nucleotide_pi <- function(vcf_dt, aa_mut_rates, genetic_code_df) {
  #' Combined function: biallelic codon analysis + nucleotide pi
  #' 
  #' Returns a data.table with:
  #' - k, n, p: biallelic codon frequencies (for SFS/gamma estimation)
  #' - Site_Pi_Codon: heterozygosity at codon level (biallelic)
  #' - Site_Pi_Nucleotide: true nucleotide diversity at synonymous positions
  #' - u, v: mutation rates
  
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  if (!is.data.table(aa_mut_rates)) setDT(aa_mut_rates)
  
  # --- PREPARE GENETIC CODE ---
  gen_code_mod <- copy(genetic_code_df)
  setnames(gen_code_mod, c("AA", "Codon"), c("Variant_AA", "Codon"))
  gen_code_mod[Codon %in% c("AGT", "AGC"), Variant_AA := "Z"]
  
  # --- PREPARE VCF ---
  dt <- vcf_dt[, .(
    Variant_String = unlist(strsplit(Codon_Variants, ";", fixed = TRUE))
  ), by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  dt[, c("Variant_Codon", "Count") := tstrsplit(Variant_String, ":", fixed=TRUE)]
  dt[, Count := as.integer(Count)]
  dt[, Variant_String := NULL]
  
  setnames(dt, "AA", "Site_AA")
  dt[Site_AA == "S" & (Preferred_Codon == "AGT" | Preferred_Codon == "AGC"), Site_AA := "Z"]
  
  # --- JOIN & FILTER ---
  setkey(gen_code_mod, Codon)
  setkey(dt, Variant_Codon)
  dt_merged <- gen_code_mod[dt, nomatch=0] 
  dt_syn <- dt_merged[Site_AA == Variant_AA]
  
  # --- BIALLELIC CODON ANALYSIS (for SFS) ---
  result_codon <- dt_syn[, .(
    k = sum(Count[Codon == Preferred_Codon]),
    n = sum(Count)
  ), by = .(Gene, Codon_Pos, Site_AA, Preferred_Codon)]
  
  result_codon[, p := ifelse(n > 0, k/n, 0)]
  result_codon[, Site_Pi_Codon := ifelse(n > 1, 2 * p * (1-p) * (n/(n-1)), 0)]
  
  # --- NUCLEOTIDE-LEVEL PI (BIALLELIC: Preferred vs Non-Preferred) ---
  # Calculate nucleotide diversity treating preferred/non-preferred as biallelic
  # This is consistent with the Wright-Fisher model used for gamma estimation
  result_nuc <- dt_syn[, {
    codons <- Codon
    counts <- Count
    pref_codon <- Preferred_Codon[1]
    total_n <- sum(counts)
    
    if (length(codons) < 2 || total_n < 2) {
      list(Site_Pi_Nucleotide = 0.0, N_Syn_Positions = 0L)
    } else {
      # Classify codons as preferred or non-preferred
      is_preferred <- codons == pref_codon
      n_pref <- sum(counts[is_preferred])
      n_nonpref <- sum(counts[!is_preferred])
      
      # Calculate biallelic frequency
      p_pref <- n_pref / total_n
      
      # Biallelic heterozygosity (same as Site_Pi_Codon but recalculated here for nucleotide interpretation)
      pi_codon_biallelic <- ifelse(total_n > 1, 
                                    2 * p_pref * (1 - p_pref) * (total_n / (total_n - 1)), 
                                    0.0)
      
      # For nucleotide-level pi: count how many nucleotide positions differ between
      # preferred and non-preferred codons on average
      if (n_pref > 0 && n_nonpref > 0) {
        # Get representative codons
        pref_codons_present <- unique(codons[is_preferred])
        nonpref_codons_present <- unique(codons[!is_preferred])
        
        # Count positions that differ
        n_diff_positions <- 0
        for (pos in 1:3) {
          nucs_pref <- unique(substr(pref_codons_present, pos, pos))
          nucs_nonpref <- unique(substr(nonpref_codons_present, pos, pos))
          # Position differs if there's no overlap between preferred and non-preferred nucleotides
          if (length(intersect(nucs_pref, nucs_nonpref)) == 0) {
            n_diff_positions <- n_diff_positions + 1
          } else if (length(nucs_pref) > 1 || length(nucs_nonpref) > 1) {
            # Position is polymorphic within groups
            n_diff_positions <- n_diff_positions + 1
          }
        }
        
        # Nucleotide pi is the codon-level heterozygosity scaled by proportion of differing positions
        # This gives expected nucleotide diversity per synonymous site
        pi_nucleotide <- pi_codon_biallelic * (n_diff_positions / 3)
        n_syn_pos <- as.integer(n_diff_positions)
      } else {
        pi_nucleotide <- 0.0
        n_syn_pos <- 0L
      }
      
      list(
        Site_Pi_Nucleotide = as.numeric(pi_nucleotide),
        N_Syn_Positions = n_syn_pos
      )
    }
  }, by = .(Gene, Codon_Pos, Site_AA, Preferred_Codon)]
  
  # --- MERGE RESULTS ---
  result <- merge(result_codon, result_nuc, 
                  by = c("Gene", "Codon_Pos", "Site_AA", "Preferred_Codon"))
  
  setnames(result, "Site_AA", "AA")
  
  # --- MERGE MUTATION RATES ---
  setkey(result, AA)
  setkey(aa_mut_rates, AA)
  final <- aa_mut_rates[result, nomatch=0]
  
  return(final)
}


calculate_pi_analytical <- function(alpha, beta, S) {
  require(gsl)
  
  # 1. Calculate the Denominator (Normalization)
  # Beta(a,b) * 1F1(a, a+b, S)
  # We use lbeta (log beta) to avoid underflow/overflow
  log_denom <- lbeta(alpha, beta) + log(hyperg_1F1(alpha, alpha + beta, S))
  
  # 2. Calculate the Numerator (Expected Heterozygosity)
  # We are integrating 2 * x * (1-x) * Density
  # This shifts the powers to alpha+1 and beta+1
  log_num_integral <- lbeta(alpha + 1, beta + 1) + log(hyperg_1F1(alpha + 1, alpha + beta + 2, S))
  
  # Add log(2) because of the "2pq" in heterozygosity
  log_numerator <- log(2) + log_num_integral
  
  # 3. Ratio
  return(exp(log_numerator - log_denom))
}

get_aa_mutation_rates <- function(Q, pref_codons_list, genetic_code) {
  
  results <- data.frame(AA = character(), u = numeric(), v = numeric(), 
                        stringsAsFactors = FALSE)
  
  for (aa in names(genetic_code)) {
    codons <- genetic_code[[aa]]
    if (length(codons) < 2) next # Skip Met, Trp
    
    # Identify Preferred and Unpreferred sets for this AA
    # You need to supply 'pref_codons_list' based on your AnaCoDa results
    pref <- intersect(codons, pref_codons_list[[aa]])
    unpref <- setdiff(codons, pref)
    
    if (length(pref) == 0 || length(unpref) == 0) next
    
    # --- Calculate u (Unpreferred -> Preferred) ---
    # We average the "flux" leaving each Unpreferred codon
    u_rates <- c()
    for (c_un in unpref) {
      rate_out <- 0
      for (c_p in pref) {
        # Check if 1-step mutation
        diffs <- 0
        nuc_from <- ""
        nuc_to <- ""
        for (i in 1:3) {
          if (substr(c_un, i, i) != substr(c_p, i, i)) {
            diffs <- diffs + 1
            nuc_from <- substr(c_un, i, i)
            nuc_to <- substr(c_p, i, i)
          }
        }
        
        # If it is a single step mutation, it is a valid synonymous path
        if (diffs == 1) {
          rate_out <- rate_out + Q[nuc_from, nuc_to]
        }
      }
      u_rates <- c(u_rates, rate_out)
    }
    u <- mean(u_rates) # Average rate per Unpreferred codon
    
    # --- Calculate v (Preferred -> Unpreferred) ---
    v_rates <- c()
    for (c_p in pref) {
      rate_out <- 0
      for (c_un in unpref) {
        # Check if 1-step mutation
        diffs <- 0
        nuc_from <- ""
        nuc_to <- ""
        for (i in 1:3) {
          if (substr(c_p, i, i) != substr(c_un, i, i)) {
            diffs <- diffs + 1
            nuc_from <- substr(c_p, i, i)
            nuc_to <- substr(c_un, i, i)
          }
        }
        
        if (diffs == 1) {
          rate_out <- rate_out + Q[nuc_from, nuc_to]
        }
      }
      v_rates <- c(v_rates, rate_out)
    }
    v <- mean(v_rates) # Average rate per Preferred codon
    
    results <- rbind(results, data.frame(AA = aa, u = u, v = v))
  }
  return(results)
}

calc_gamma_wrapper <- function(k_list, n_list, u, v) {
  # Don't run if too few sites
  if(length(k_list) < 2) return(NA) 
  
  # Run your optimization function
  gamma <- estimate_gamma_for_AA(counts = unlist(k_list), 
                                 sample_sizes = unlist(n_list), 
                                 u = u, v = v)
  return(gamma)
}

process_codon_vcf_corrected <- function(vcf_dt, aa_mut_rates, genetic_code_df) {
  
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  if (!is.data.table(aa_mut_rates)) setDT(aa_mut_rates)
  
  # Prepare genetic code with S/Z split
  gen_code_mod <- copy(genetic_code_df)
  setnames(gen_code_mod, c("AA", "Codon"), c("Variant_AA", "Codon"))
  gen_code_mod[Codon %in% c("AGT", "AGC"), Variant_AA := "Z"]
  
  # Expand VCF rows
  dt <- vcf_dt[, .(
    Variant_String = unlist(strsplit(Codon_Variants, ";", fixed = TRUE))
  ), by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  # Split "Codon:Count"
  dt[, c("Variant_Codon", "Count") := tstrsplit(Variant_String, ":", fixed=TRUE)]
  dt[, Count := as.integer(Count)]
  dt[, Variant_String := NULL]
  
  # Handle S/Z split for sites
  setnames(dt, "AA", "Site_AA")
  dt[Site_AA == "S" & (Preferred_Codon == "AGT" | Preferred_Codon == "AGC"), Site_AA := "Z"]
  
  # Join with genetic code
  setkey(gen_code_mod, Codon)
  setkey(dt, Variant_Codon)
  dt_merged <- gen_code_mod[dt, nomatch=0]
  
  # Filter for synonymous only
  dt_syn <- dt_merged[Site_AA == Variant_AA]
  
  # *** KEY FIX: Calculate n as total sample size, not sum of counts ***
  # The sample size should be consistent across all variants at a site
  result <- dt_syn[, .(
    k = sum(Count[Codon == Preferred_Codon]),  # Count of preferred allele
    n = sum(Count),  # Total alleles sampled at this site (should be ~187 based on your data)
    n_variants = .N  # Number of distinct variants (for QC)
  ), by = .(Gene, Codon_Pos, Site_AA, Preferred_Codon)]
  
  setnames(result, "Site_AA", "AA")
  
  # Calculate allele frequency and Pi
  result[, p := ifelse(n > 0, k/n, 0)]
  result[, Site_Pi := ifelse(n > 1, 2 * p * (1-p) * (n/(n-1)), 0)]
  
  # Merge mutation rates
  setkey(result, AA)
  setkey(aa_mut_rates, AA)
  final <- aa_mut_rates[result, nomatch=0]
  
  return(final)
}

calc_gamma_vectorized <- function(k_vec, n_vec, u, v) {
  # If gene has too few sites for this AA, return NA
  if(length(k_vec) == 0) return(NA)
  
  # Use the optimizer function defined previously
  tryCatch(
    estimate_gamma_for_AA(counts = k_vec, sample_sizes = n_vec, u = u, v = v),
    error = function(e) NA
  )
}

calc_gamma_wrapper <- function(k_list, n_list, u, v) {
  # Don't run if too few sites
  if(length(k_list) < 2) return(NA) 
  
  # Run your optimization function
  gamma <- estimate_gamma_for_AA(counts = unlist(k_list), 
                                 sample_sizes = unlist(n_list), 
                                 u = u, v = v)
  return(gamma)
}

process_codon_vcf_fast <- function(vcf_dt, aa_mut_rates, genetic_code_df) {
  
  # 0. Ensure inputs are data.tables
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  if (!is.data.table(aa_mut_rates)) setDT(aa_mut_rates)
  
  # --- PREPARE GENETIC CODE ---
  # Rename columns to avoid collision: AA -> Variant_AA
  gen_code_mod <- copy(genetic_code_df)
  setnames(gen_code_mod, c("AA", "Codon"), c("Variant_AA", "Codon"))
  
  # Teach genetic code that AGT/AGC are 'Z' (for the Variant check)
  gen_code_mod[Codon %in% c("AGT", "AGC"), Variant_AA := "Z"]
  
  # --- PREPARE VCF ---
  # 1. Expand Rows (Long Format)
  # Unlist splits the variants into rows
  dt <- vcf_dt[, .(
    Variant_String = unlist(strsplit(Codon_Variants, ";", fixed = TRUE))
  ), by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  # 2. Split "Codon:Count"
  dt[, c("Variant_Codon", "Count") := tstrsplit(Variant_String, ":", fixed=TRUE)]
  dt[, Count := as.integer(Count)]
  dt[, Variant_String := NULL]
  
  # 3. Handle S/Z Split for the SITES
  # Rename AA to Site_AA to be explicit
  setnames(dt, "AA", "Site_AA")
  dt[Site_AA == "S" & (Preferred_Codon == "AGT" | Preferred_Codon == "AGC"), Site_AA := "Z"]
  
  # --- JOIN & FILTER ---
  # Set keys for fast merge
  setkey(gen_code_mod, Codon)
  setkey(dt, Variant_Codon)
  
  # Merge on Codon == Variant_Codon
  # dt_merged will have: Site_AA, Preferred_Codon, Variant_AA, Count
  dt_merged <- gen_code_mod[dt, nomatch=0] 
  
  # Filter: The Variant's Amino Acid must match the Site's Amino Acid
  # (This excludes non-synonymous variants)
  dt_syn <- dt_merged[Site_AA == Variant_AA]
  
  # --- AGGREGATE ---
  # Now we have clear column names, no "i." needed
  # For SFS-based gamma estimation, we treat this as BIALLELIC:
  # Preferred codon(s) vs Non-preferred codon(s)
  result <- dt_syn[, .(
    k = sum(Count[Codon == Preferred_Codon]), # Count of preferred codon alleles
    n = sum(Count)                             # Total sample size at this site
  ), by = .(Gene, Codon_Pos, Site_AA, Preferred_Codon)]
  
  # Rename back to standard 'AA' for final merge
  setnames(result, "Site_AA", "AA")
  
  # --- CALCULATE PI (BIALLELIC HETEROZYGOSITY) ---
  # p = frequency of preferred codon
  # Site_Pi = expected heterozygosity assuming biallelic model
  result[, p := ifelse(n > 0, k/n, 0)]
  result[, Site_Pi := ifelse(n > 1, 2 * p * (1-p) * (n/(n-1)), 0)]
  
  # --- MERGE RATES ---
  setkey(result, AA)
  setkey(aa_mut_rates, AA)
  final <- aa_mut_rates[result, nomatch=0]
  
  return(final)
}

worker_gamma_est <- function(k, n, u, v) {
  # Quick check to skip useless rows
  if (length(k) < 5) return(NA_real_) 
  
  # Run your existing estimator
  tryCatch(
    estimate_gamma_for_AA(counts = k, sample_sizes = n, u = u, v = v),
    error = function(e) NA_real_
  )
}