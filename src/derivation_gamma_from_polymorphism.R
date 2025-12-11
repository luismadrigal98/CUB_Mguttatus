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