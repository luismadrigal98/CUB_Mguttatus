#' Auxiliary functions to estimate gamma from polymorphism data
#' 
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

# ******************************************************************************
# 1) Define the Probability Function (Wright's Distribution integrated) ----
# ______________________________________________________________________________

get_prob_k <- function(k, n, alpha, beta, S) 
{
  #' Calculates the probability of observing k derived alleles in sample size n
  #' given mutation rates (alpha=4Neu, beta=4Nev) and selection (S=4Nes).
  #'
  #' u: rate from Unpreferred -> Preferred
  #' v: rate from Preferred -> Unpreferred
  #' alpha: Population level mutation rate towards preferred
  #' beta: Population level mutation rate towards unpreferred
  #' S: Selection coefficient favoring Preferred
  
  # Define the unnormalized Wright density
  # x is frequency of PREFERRED allele
  wright_density <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    # Using log scale for numerical stability then exp
    val <- (alpha - 1)*log(x) + (beta - 1)*log(1-x) + S*x
    return(exp(val))
  }
  
  # Calculate Normalization Constant (Denominator)
  # Integrate density from 0 to 1
  denom <- tryCatch(
    integrate(wright_density, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(denom) || denom == 0) return(0)
  
  # Calculate Numerator (Probability of sampling k alleles)
  # Integral of binom(n,k) * x^k * (1-x)^(n-k) * density
  numerator_func <- function(x) {
    if(x <= 0 | x >= 1) return(0)
    binom_prob <- dbinom(k, size=n, prob=x)
    return(binom_prob * wright_density(x))
  }
  
  num <- tryCatch(
    integrate(numerator_func, 0, 1)$value,
    error = function(e) NA
  )
  
  if(is.na(num)) return(0)
  
  return(num / denom)
}

# ******************************************************************************
# 2) Likelihood Optimizer ----
# ______________________________________________________________________________

# Input: 
#   counts: Vector of observed counts of PREFERRED codons (k) across many sites
#   sample_sizes: Vector of sample sizes (n) for those sites (usually 187, because we have haplotypes)
#   u, v: Mutation rates for this specific Amino Acid family

estimate_gamma_for_AA <- function(counts, sample_sizes, u, v,
                                  S_interval = c(-5, 20)) {
  
  # Negative Log Likelihood Function to minimize
  nll <- function(S) {
    # Calculate log-probability for each site
    log_probs <- mapply(function(k, n) {
      p <- get_prob_k(k, n, u, v, S)
      if(p <= 0) return(-1e6) # Penalty for impossible values
      return(log(p))
    }, counts, sample_sizes)
    
    return(-sum(log_probs))
  }
  
  # Optimize S (gamma)
  # Searching likely range -5 to 20
  opt <- optimize(nll, interval = S_interval)
  
  return(opt$minimum) # This is your Gamma (4Nes)
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

process_codon_vcf <- function(vcf_dt, aa_mut_rates, genetic_code_df) {
  
  # Ensure input is a data.table for speed
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  
  # A. Parse the "Codon_Variants" column
  # We select only necessary columns and expand the variants
  # "GCC:183;GCA:4" -> two rows
  long_vcf <- vcf_dt[, .(Gene, Codon_Pos, AA, Preferred_Codon, Codon_Variants)] %>%
    separate_rows(Codon_Variants, sep = ";") %>%
    separate(Codon_Variants, into = c("Variant_Codon", "Count"), sep = ":", convert = TRUE)
  
  # B. Handle the Serine S/Z Split Correction
  # AnaCoDa calls 2-fold Serines 'Z' (AGT, AGC).
  # The VCF likely labels them 'S'. We must re-label them based on the Preferred Codon.
  # If Preferred is AGT or AGC, change AA to 'Z'. Otherwise leave as 'S'.
  long_vcf <- long_vcf %>%
    mutate(
      AA = case_when(
        AA == "S" & (Preferred_Codon == "AGT" | Preferred_Codon == "AGC") ~ "Z",
        TRUE ~ AA
      )
    )
  
  # C. Filter for Synonymous Variants Only
  # Join variant codons with genetic code to check their AA identity
  clean_counts <- long_vcf %>%
    inner_join(genetic_code_df, by = c("Variant_Codon" = "Codon")) %>%
    # IMPORTANT: Filter condition
    # AA.x is the Site's AA (from VCF column, potentially corrected to Z)
    # AA.y is the Variant's AA (from genetic code lookup)
    filter(AA.x == AA.y) %>% 
    rename(AA_Site = AA.x) %>%
    dplyr::select(-AA.y) 
  
  # D. Aggregate to get k (Preferred) and Recalculated n (Total Synonymous)
  site_stats <- clean_counts %>%
    group_by(Gene, Codon_Pos, AA_Site, Preferred_Codon) %>%
    summarise(
      # k = Sum counts where Variant is the Preferred one
      k = sum(Count[Variant_Codon == Preferred_Codon]),
      
      # n = Sum counts of ALL synonymous variants (excluding non-syn noise)
      n = sum(Count),
      
      # Calculate Site Pi
      p = ifelse(n > 0, k/n, 0),
      # Note: n/(n-1) correction handles sample size bias
      Site_Pi = ifelse(n > 1, 2 * p * (1-p) * (n/(n-1)), 0),
      
      .groups = "drop"
    ) %>%
    rename(AA = AA_Site) # Rename back to AA for merging
  
  # E. Merge in the Mutation Rates (u, v)
  final_site_data <- site_stats %>%
    left_join(aa_mut_rates, by = "AA") %>%
    filter(!is.na(u)) # Removes Stop codons or AAs without defined rates
  
  return(final_site_data)
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
  
  # Ensure input is data.table
  if (!is.data.table(vcf_dt)) setDT(vcf_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  if (!is.data.table(aa_mut_rates)) setDT(aa_mut_rates)
  
  # 1. Select only needed columns to save memory
  # We subset immediately
  dt <- vcf_dt[, .(Gene, Codon_Pos, AA, Preferred_Codon, Codon_Variants)]
  
  # 2. Fast String Splitting (The heavy lifting)
  # This uses data.table's internal C-based splitter which is vastly faster than tidyr
  dt <- dt[, tstrsplit(Codon_Variants, ";", fixed=TRUE), by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  # Melt long to get single column of variants
  dt <- melt(dt, id.vars = c("Gene", "Codon_Pos", "AA", "Preferred_Codon"), 
             value.name = "Variant_String", na.rm = TRUE)
  
  # Split "Codon:Count" (e.g. "GCC:183")
  dt[, c("Variant_Codon", "Count") := tstrsplit(Variant_String, ":", fixed=TRUE)]
  dt[, Count := as.integer(Count)]
  dt[, Variant_String := NULL] # Clean up
  
  # 3. Handle S/Z Split (Vectorized)
  dt[AA == "S" & (Preferred_Codon == "AGT" | Preferred_Codon == "AGC"), AA := "Z"]
  
  # 4. Filter Synonymous (Fast Join)
  # Set keys for speed
  setkey(dt, Variant_Codon)
  setkey(genetic_code_df, Codon)
  
  # Inner Join
  dt <- genetic_code_df[dt, nomatch=0] # Match Variant_Codon to Codon
  
  # Filter where Site AA == Variant AA (renamed column handling)
  # data.table merge usually keeps the 'i' columns, check names
  # genetic_code_df has 'AA' and 'Codon'. dt has 'AA' (site). 
  # After merge, we likely have i.AA (site) and AA (variant)
  dt <- dt[AA == i.AA] 
  
  # 5. Aggregation
  result <- dt[, .(
    k = sum(Count[Variant_Codon == i.Preferred_Codon]),
    n = sum(Count)
  ), by = .(Gene, Codon_Pos, i.AA, i.Preferred_Codon)]
  
  # Rename for clarity
  setnames(result, c("i.AA", "i.Preferred_Codon"), c("AA", "Preferred_Codon"))
  
  # 6. Calculate Pi
  result[, p := ifelse(n > 0, k/n, 0)]
  result[, Site_Pi := ifelse(n > 1, 2 * p * (1-p) * (n/(n-1)), 0)]
  
  # 7. Merge Rates
  setkey(result, AA)
  setkey(aa_mut_rates, AA)
  final <- aa_mut_rates[result, nomatch=0]
  
  return(final)
}