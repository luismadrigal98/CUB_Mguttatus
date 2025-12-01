#!/usr/bin/env Rscript
#
# ******************************************************************************
# Pure R Implementation of the ROC (Ribosome Overhead Cost) Model
# 
# This is a pure R implementation of the AnaCoDa ROC model for codon usage bias
# analysis. It avoids the C++ backend and associated bugs, at the cost of speed.
#
# @author Luis Javier Madrigal-Roca
# @date 12/01/2025
#
# REFERENCE:
# Gilchrist et al. (2015) "Estimating Gene Expression and Codon-Specific 
# Translational Efficiencies, Mutation Biases, and Selection Coefficients 
# from Genomic Data Alone"
#
# MODEL OVERVIEW:
# The ROC model assumes codon frequencies follow a multinomial distribution
# where the probability of each codon depends on:
#   - Delta M (mutation bias): reflects mutational pressure at wobble position
#   - Delta Eta (selection): reflects selection for translational efficiency
#   - Phi (expression level): genes with higher phi show stronger selection
#
# P(codon_i | phi, dM, dEta) = exp(-dM_i - dEta_i * phi) / Z
# where Z is the partition function (sum over all codons for that amino acid)
#
# PARALLELIZATION:
# Uses the 'parallel' package with mclapply for gene-level likelihood 
# calculations. Set the number of cores via --cores argument or n_cores parameter.
# On HPC systems with SLURM, will auto-detect from SLURM_CPUS_PER_TASK.
#
# ******************************************************************************

library(parallel)

# ==============================================================================
# GLOBAL PARALLEL CONFIGURATION
# ==============================================================================

# Detect available cores from environment or system
detect_cores <- function() {
  # First check SLURM environment
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (slurm_cpus != "") {
    return(as.integer(slurm_cpus))
  }
  
  # Then check PBS/Torque
  pbs_cpus <- Sys.getenv("PBS_NUM_PPN")
  if (pbs_cpus != "") {
    return(as.integer(pbs_cpus))
  }
  
  # Fall back to parallel detection (leave 1 core free)
  n <- parallel::detectCores()
  if (is.na(n)) n <- 1
  return(max(1, n - 1))
}

# Global variable for number of cores (can be overridden)
N_CORES <- 1  # Default to serial, will be set by run_roc_mcmc

# ==============================================================================
# GENETIC CODE AND CODON TABLES
# ==============================================================================

#' Get the standard genetic code
#' @return Named list mapping amino acids to their codons
get_genetic_code <- function() {
  # Standard genetic code (codon -> amino acid)
  codon_to_aa <- c(
    "TTT"="F", "TTC"="F",
    "TTA"="L", "TTG"="L", "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
    "ATT"="I", "ATC"="I", "ATA"="I",
    "ATG"="M",
    "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
    "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S", "AGT"="S", "AGC"="S",
    "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
    "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
    "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
    "TAT"="Y", "TAC"="Y",
    "TAA"="X", "TAG"="X", "TGA"="X",  # Stop codons
    "CAT"="H", "CAC"="H",
    "CAA"="Q", "CAG"="Q",
    "AAT"="N", "AAC"="N",
    "AAA"="K", "AAG"="K",
    "GAT"="D", "GAC"="D",
    "GAA"="E", "GAG"="E",
    "TGT"="C", "TGC"="C",
    "TGG"="W",
    "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R", "AGA"="R", "AGG"="R",
    "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
  )
  
  # Invert to get amino acid -> codons
  aa_to_codons <- list()
  for (codon in names(codon_to_aa)) {
    aa <- codon_to_aa[codon]
    if (is.null(aa_to_codons[[aa]])) {
      aa_to_codons[[aa]] <- codon
    } else {
      aa_to_codons[[aa]] <- c(aa_to_codons[[aa]], codon)
    }
  }
  
  # Sort codons alphabetically within each AA (last one is reference)
  for (aa in names(aa_to_codons)) {
    aa_to_codons[[aa]] <- sort(aa_to_codons[[aa]])
  }
  
  return(aa_to_codons)
}

#' Get amino acids that have synonymous codons (exclude M, W, stop codons)
#' @return Character vector of amino acid single-letter codes
get_synonymous_aas <- function() {
  aa_to_codons <- get_genetic_code()
  # Keep only AAs with 2+ codons, exclude stops
  syn_aas <- names(aa_to_codons)[sapply(aa_to_codons, length) >= 2]
  syn_aas <- syn_aas[!syn_aas %in% c("M", "W", "X")]
  return(sort(syn_aas))
}

# ==============================================================================
# SEQUENCE PROCESSING
# ==============================================================================

#' Read FASTA file and extract codon counts per gene
#' @param fasta_file Path to FASTA file
#' @return List with gene_ids and codon_counts matrix
read_fasta_codon_counts <- function(fasta_file) {
  lines <- readLines(fasta_file)
  
  gene_ids <- c()
  sequences <- c()
  current_seq <- ""
  current_id <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (current_id != "") {
        gene_ids <- c(gene_ids, current_id)
        sequences <- c(sequences, current_seq)
      }
      # Extract gene ID (first word after >)
      current_id <- sub("^>([^ ]+).*", "\\1", line)
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, toupper(gsub("\\s", "", line)))
    }
  }
  # Don't forget the last sequence
  if (current_id != "") {
    gene_ids <- c(gene_ids, current_id)
    sequences <- c(sequences, current_seq)
  }
  
  # Count codons for each gene
  aa_to_codons <- get_genetic_code()
  all_codons <- unlist(aa_to_codons)
  
  codon_counts <- matrix(0, nrow = length(sequences), ncol = length(all_codons))
  colnames(codon_counts) <- all_codons
  rownames(codon_counts) <- gene_ids
  
  for (i in seq_along(sequences)) {
    seq <- sequences[i]
    seq_len <- nchar(seq)
    
    # Extract codons (must be multiple of 3)
    if (seq_len %% 3 != 0) {
      warning(paste("Gene", gene_ids[i], "length not divisible by 3, trimming"))
      seq_len <- seq_len - (seq_len %% 3)
    }
    
    for (j in seq(1, seq_len - 2, by = 3)) {
      codon <- substr(seq, j, j + 2)
      if (codon %in% all_codons) {
        codon_counts[i, codon] <- codon_counts[i, codon] + 1
      }
    }
  }
  
  return(list(
    gene_ids = gene_ids,
    codon_counts = codon_counts,
    n_genes = length(gene_ids)
  ))
}

# ==============================================================================
# ROC MODEL LIKELIHOOD
# ==============================================================================

#' Calculate log probability of codon usage for one amino acid in one gene
#' @param codon_counts Named vector of codon counts for this AA
#' @param dM Vector of delta M values (mutation bias) for non-reference codons
#' @param dEta Vector of delta Eta values (selection) for non-reference codons
#' @param phi Gene expression level
#' @return Log likelihood contribution
calc_log_likelihood_aa <- function(codon_counts, dM, dEta, phi) {
  n_codons <- length(codon_counts)
  if (n_codons <= 1 || sum(codon_counts) == 0) {
    return(0)  # No information from single codon AAs or zero counts
  }
  
  # Reference codon is last (alphabetically), has dM=0, dEta=0
  # Calculate log probabilities
  # log P(codon_i) = -dM_i - dEta_i * phi - log(Z)
  
  # For numerical stability, shift by minimum value
  log_numerators <- c(-dM - dEta * phi, 0)  # Last is reference (0)
  max_log <- max(log_numerators)
  log_Z <- max_log + log(sum(exp(log_numerators - max_log)))
  
  log_probs <- log_numerators - log_Z
  
  # Multinomial log likelihood
  log_lik <- sum(codon_counts * log_probs)
  
  return(log_lik)
}

#' Calculate total log likelihood for one gene
#' @param gene_codon_counts Named vector of all codon counts for gene
#' @param dM_list List of dM vectors, one per amino acid
#' @param dEta_list List of dEta vectors, one per amino acid
#' @param phi Gene expression level
#' @param aa_to_codons Genetic code mapping
#' @return Total log likelihood for gene
calc_log_likelihood_gene <- function(gene_codon_counts, dM_list, dEta_list, 
                                     phi, aa_to_codons) {
  syn_aas <- get_synonymous_aas()
  log_lik <- 0
  
  for (aa in syn_aas) {
    codons <- aa_to_codons[[aa]]
    counts <- gene_codon_counts[codons]
    
    if (sum(counts) > 0) {
      n_non_ref <- length(codons) - 1
      dM <- dM_list[[aa]][1:n_non_ref]
      dEta <- dEta_list[[aa]][1:n_non_ref]
      
      log_lik <- log_lik + calc_log_likelihood_aa(counts, dM, dEta, phi)
    }
  }
  
  return(log_lik)
}

#' Calculate total log likelihood across all genes (PARALLEL VERSION)
#' @param codon_counts Matrix of codon counts (genes x codons)
#' @param dM_list List of dM vectors
#' @param dEta_list List of dEta vectors
#' @param phi Vector of expression levels
#' @param aa_to_codons Genetic code mapping
#' @param n_cores Number of cores for parallel computation
#' @return Total log likelihood
calc_total_log_likelihood <- function(codon_counts, dM_list, dEta_list, 
                                      phi, aa_to_codons, n_cores = 1) {
  n_genes <- nrow(codon_counts)
  
  if (n_cores > 1 && n_genes >= n_cores) {
    # Parallel computation across genes
    gene_logliks <- parallel::mclapply(
      1:n_genes,
      function(i) {
        calc_log_likelihood_gene(
          codon_counts[i, ], dM_list, dEta_list, phi[i], aa_to_codons
        )
      },
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )
    log_lik <- sum(unlist(gene_logliks))
  } else {
    # Serial computation
    log_lik <- 0
    for (i in 1:n_genes) {
      log_lik <- log_lik + calc_log_likelihood_gene(
        codon_counts[i, ], dM_list, dEta_list, phi[i], aa_to_codons
      )
    }
  }
  
  return(log_lik)
}

# ==============================================================================
# PRIOR DISTRIBUTIONS
# ==============================================================================

#' Log prior for phi (log-normal distribution)
#' @param phi Vector of expression levels
#' @param sphi Standard deviation of log(phi)
#' @return Log prior probability
log_prior_phi <- function(phi, sphi) {
  # phi ~ LogNormal(mu = -sphi^2/2, sigma = sphi)
  # This parameterization gives E[phi] = 1
  mu <- -sphi^2 / 2
  sum(dlnorm(phi, meanlog = mu, sdlog = sphi, log = TRUE))
}

#' Log prior for CSP (mutation/selection parameters)
#' @param params Vector of parameter values
#' @param prior_mean Prior mean (usually 0)
#' @param prior_sd Prior standard deviation
#' @return Log prior probability
log_prior_csp <- function(params, prior_mean = 0, prior_sd = 0.35) {
  sum(dnorm(params, mean = prior_mean, sd = prior_sd, log = TRUE))
}

#' Log prior for observation noise (when using empirical phi)
#' @param obs_phi Observed phi values
#' @param true_phi True phi values
#' @param sepsilon Observation noise SD
#' @return Log prior probability
log_prior_obs_phi <- function(obs_phi, true_phi, sepsilon) {
  # log(obs_phi) ~ Normal(log(true_phi), sepsilon)
  sum(dnorm(log(obs_phi), mean = log(true_phi), sd = sepsilon, log = TRUE))
}

# ==============================================================================
# MCMC PROPOSALS
# ==============================================================================

#' Propose new phi value (on log scale)
#' @param current_phi Current phi value
#' @param proposal_sd Proposal standard deviation
#' @return New proposed phi value
propose_phi <- function(current_phi, proposal_sd = 0.1) {
  log_phi_new <- rnorm(1, mean = log(current_phi), sd = proposal_sd)
  return(exp(log_phi_new))
}

#' Propose new CSP value
#' @param current_val Current parameter value
#' @param proposal_sd Proposal standard deviation
#' @return New proposed value
propose_csp <- function(current_val, proposal_sd = 0.1) {
  rnorm(1, mean = current_val, sd = proposal_sd)
}

# ==============================================================================
# MCMC SAMPLER
# ==============================================================================

#' Initialize parameters
#' @param n_genes Number of genes
#' @param aa_to_codons Genetic code
#' @param obs_phi Optional observed phi values for initialization
#' @return List of initial parameter values
initialize_parameters <- function(n_genes, aa_to_codons, obs_phi = NULL) {
  syn_aas <- get_synonymous_aas()
  
  # Initialize phi
  if (!is.null(obs_phi)) {
    phi <- obs_phi
    phi[phi <= 0] <- 0.01  # Ensure positive
  } else {
    phi <- rep(1, n_genes)
  }
  
  # Initialize dM and dEta (all zeros = no bias)
  dM_list <- list()
  dEta_list <- list()
  
  for (aa in syn_aas) {
    n_codons <- length(aa_to_codons[[aa]])
    n_params <- n_codons - 1  # Reference codon excluded
    dM_list[[aa]] <- rep(0, n_params)
    dEta_list[[aa]] <- rep(0, n_params)
  }
  
  # Hyperparameters
  sphi <- 1.0
  sepsilon <- 0.5  # If using observed phi
  
  return(list(
    phi = phi,
    dM = dM_list,
    dEta = dEta_list,
    sphi = sphi,
    sepsilon = sepsilon
  ))
}

#' Run MCMC for ROC model
#' @param codon_counts Matrix of codon counts
#' @param n_samples Number of MCMC samples
#' @param thin Thinning interval
#' @param obs_phi Optional observed phi values
#' @param fix_phi If TRUE, don't update phi
#' @param fix_dM If TRUE, don't update dM
#' @param fix_dEta If TRUE, don't update dEta
#' @param init_dM Optional initial dM values (list)
#' @param init_dEta Optional initial dEta values (list)
#' @param adapt_interval Interval for adapting proposal widths
#' @param n_cores Number of cores for parallel computation
#' @param verbose Print progress
#' @return List with MCMC samples and diagnostics
run_roc_mcmc <- function(codon_counts, 
                         n_samples = 1000,
                         thin = 10,
                         burnin_frac = 0.5,
                         obs_phi = NULL,
                         fix_phi = FALSE,
                         fix_dM = FALSE,
                         fix_dEta = FALSE,
                         init_dM = NULL,
                         init_dEta = NULL,
                         adapt_interval = 100,
                         n_cores = 1,
                         verbose = TRUE) {
  
  # Setup
  aa_to_codons <- get_genetic_code()
  syn_aas <- get_synonymous_aas()
  n_genes <- nrow(codon_counts)
  n_iter <- n_samples * thin
  
  # Configure parallelization
  if (n_cores <= 0) {
    n_cores <- detect_cores()
  }
  n_cores <- min(n_cores, n_genes)  # Don't use more cores than genes
  
  if (verbose) {
    message("============================================")
    message("     ROC Model MCMC (Pure R Implementation)")
    message("============================================")
    message(sprintf("Genes: %d", n_genes))
    message(sprintf("Samples: %d (thin=%d, total iter=%d)", n_samples, thin, n_iter))
    message(sprintf("Fix phi: %s, Fix dM: %s, Fix dEta: %s", fix_phi, fix_dM, fix_dEta))
    message(sprintf("Parallel cores: %d", n_cores))
    message("============================================")
  }
  
  # Initialize
  params <- initialize_parameters(n_genes, aa_to_codons, obs_phi)
  
  # Override with provided initial values
  if (!is.null(init_dM)) {
    params$dM <- init_dM
  }
  if (!is.null(init_dEta)) {
    params$dEta <- init_dEta
  }
  
  # If phi is fixed, use observed values
  if (fix_phi && !is.null(obs_phi)) {
    params$phi <- obs_phi
  }
  
  # Storage for samples
  n_store <- n_samples
  phi_samples <- matrix(NA, nrow = n_store, ncol = n_genes)
  sphi_samples <- numeric(n_store)
  log_posterior_trace <- numeric(n_store)
  
  # Store CSP samples (flattened)
  dM_samples <- list()
  dEta_samples <- list()
  for (aa in syn_aas) {
    n_params <- length(params$dM[[aa]])
    dM_samples[[aa]] <- matrix(NA, nrow = n_store, ncol = n_params)
    dEta_samples[[aa]] <- matrix(NA, nrow = n_store, ncol = n_params)
  }
  
  # Proposal widths (will be adapted)
  phi_prop_sd <- rep(0.3, n_genes)
  dM_prop_sd <- 0.1
  dEta_prop_sd <- 0.1
  sphi_prop_sd <- 0.1
  
  # Acceptance counters
  phi_accept <- rep(0, n_genes)
  dM_accept <- list()
  dEta_accept <- list()
  for (aa in syn_aas) {
    n_params <- length(params$dM[[aa]])
    dM_accept[[aa]] <- rep(0, n_params)
    dEta_accept[[aa]] <- rep(0, n_params)
  }
  sphi_accept <- 0
  
  # Current log posterior
  current_log_lik <- calc_total_log_likelihood(
    codon_counts, params$dM, params$dEta, params$phi, aa_to_codons, n_cores
  )
  current_log_prior_phi <- log_prior_phi(params$phi, params$sphi)
  
  current_log_prior_dM <- 0
  current_log_prior_dEta <- 0
  for (aa in syn_aas) {
    current_log_prior_dM <- current_log_prior_dM + log_prior_csp(params$dM[[aa]])
    current_log_prior_dEta <- current_log_prior_dEta + log_prior_csp(params$dEta[[aa]])
  }
  
  current_log_posterior <- current_log_lik + current_log_prior_phi + 
    current_log_prior_dM + current_log_prior_dEta
  
  sample_idx <- 0
  
  # MCMC loop
  for (iter in 1:n_iter) {
    
    # --- Update phi (gene by gene) ---
    if (!fix_phi) {
      for (g in 1:n_genes) {
        phi_new <- propose_phi(params$phi[g], phi_prop_sd[g])
        
        if (phi_new > 0) {
          # Calculate new likelihood for this gene only
          old_lik_g <- calc_log_likelihood_gene(
            codon_counts[g, ], params$dM, params$dEta, params$phi[g], aa_to_codons
          )
          new_lik_g <- calc_log_likelihood_gene(
            codon_counts[g, ], params$dM, params$dEta, phi_new, aa_to_codons
          )
          
          # Prior ratio (on log scale)
          mu <- -params$sphi^2 / 2
          old_prior <- dlnorm(params$phi[g], meanlog = mu, sdlog = params$sphi, log = TRUE)
          new_prior <- dlnorm(phi_new, meanlog = mu, sdlog = params$sphi, log = TRUE)
          
          # Jacobian for log-scale proposal
          log_jacobian <- log(phi_new) - log(params$phi[g])
          
          log_accept_ratio <- (new_lik_g - old_lik_g) + (new_prior - old_prior) + log_jacobian
          
          if (log(runif(1)) < log_accept_ratio) {
            current_log_lik <- current_log_lik + (new_lik_g - old_lik_g)
            current_log_prior_phi <- current_log_prior_phi + (new_prior - old_prior)
            params$phi[g] <- phi_new
            phi_accept[g] <- phi_accept[g] + 1
          }
        }
      }
    }
    
    # --- Update dEta (selection coefficients) ---
    if (!fix_dEta) {
      for (aa in syn_aas) {
        for (p in seq_along(params$dEta[[aa]])) {
          old_val <- params$dEta[[aa]][p]
          new_val <- propose_csp(old_val, dEta_prop_sd)
          
          # Store old, try new
          params$dEta[[aa]][p] <- new_val
          
          new_log_lik <- calc_total_log_likelihood(
            codon_counts, params$dM, params$dEta, params$phi, aa_to_codons, n_cores
          )
          new_prior <- log_prior_csp(new_val)
          old_prior <- log_prior_csp(old_val)
          
          log_accept_ratio <- (new_log_lik - current_log_lik) + (new_prior - old_prior)
          
          if (log(runif(1)) < log_accept_ratio) {
            current_log_lik <- new_log_lik
            dEta_accept[[aa]][p] <- dEta_accept[[aa]][p] + 1
          } else {
            params$dEta[[aa]][p] <- old_val  # Revert
          }
        }
      }
    }
    
    # --- Update dM (mutation coefficients) ---
    if (!fix_dM) {
      for (aa in syn_aas) {
        for (p in seq_along(params$dM[[aa]])) {
          old_val <- params$dM[[aa]][p]
          new_val <- propose_csp(old_val, dM_prop_sd)
          
          params$dM[[aa]][p] <- new_val
          
          new_log_lik <- calc_total_log_likelihood(
            codon_counts, params$dM, params$dEta, params$phi, aa_to_codons, n_cores
          )
          new_prior <- log_prior_csp(new_val)
          old_prior <- log_prior_csp(old_val)
          
          log_accept_ratio <- (new_log_lik - current_log_lik) + (new_prior - old_prior)
          
          if (log(runif(1)) < log_accept_ratio) {
            current_log_lik <- new_log_lik
            dM_accept[[aa]][p] <- dM_accept[[aa]][p] + 1
          } else {
            params$dM[[aa]][p] <- old_val
          }
        }
      }
    }
    
    # --- Update sphi (hyperparameter) ---
    if (!fix_phi) {
      sphi_new <- exp(rnorm(1, log(params$sphi), sphi_prop_sd))
      
      if (sphi_new > 0.01 && sphi_new < 10) {
        new_prior_phi <- log_prior_phi(params$phi, sphi_new)
        # Jacobian for log-scale proposal
        log_jacobian <- log(sphi_new) - log(params$sphi)
        
        log_accept_ratio <- (new_prior_phi - current_log_prior_phi) + log_jacobian
        
        if (log(runif(1)) < log_accept_ratio) {
          params$sphi <- sphi_new
          current_log_prior_phi <- new_prior_phi
          sphi_accept <- sphi_accept + 1
        }
      }
    }
    
    # Update total log posterior
    current_log_posterior <- current_log_lik + current_log_prior_phi + 
      current_log_prior_dM + current_log_prior_dEta
    
    # --- Store samples ---
    if (iter %% thin == 0) {
      sample_idx <- sample_idx + 1
      phi_samples[sample_idx, ] <- params$phi
      sphi_samples[sample_idx] <- params$sphi
      log_posterior_trace[sample_idx] <- current_log_posterior
      
      for (aa in syn_aas) {
        dM_samples[[aa]][sample_idx, ] <- params$dM[[aa]]
        dEta_samples[[aa]][sample_idx, ] <- params$dEta[[aa]]
      }
    }
    
    # --- Adapt proposal widths ---
    if (iter %% adapt_interval == 0 && iter < n_iter * burnin_frac) {
      # Adapt phi proposals
      for (g in 1:n_genes) {
        accept_rate <- phi_accept[g] / adapt_interval
        if (accept_rate < 0.2) {
          phi_prop_sd[g] <- phi_prop_sd[g] * 0.8
        } else if (accept_rate > 0.4) {
          phi_prop_sd[g] <- phi_prop_sd[g] * 1.2
        }
        phi_accept[g] <- 0
      }
      
      # Could adapt CSP proposals similarly...
    }
    
    # Progress
    if (verbose && iter %% (n_iter / 10) == 0) {
      message(sprintf("Iteration %d/%d (%.0f%%) - Log posterior: %.2f", 
                      iter, n_iter, 100 * iter / n_iter, current_log_posterior))
    }
  }
  
  # Calculate posterior means (excluding burn-in)
  burnin_samples <- floor(n_samples * burnin_frac)
  post_samples <- (burnin_samples + 1):n_samples
  
  phi_mean <- colMeans(phi_samples[post_samples, , drop = FALSE])
  phi_sd <- apply(phi_samples[post_samples, , drop = FALSE], 2, sd)
  sphi_mean <- mean(sphi_samples[post_samples])
  
  dM_mean <- list()
  dEta_mean <- list()
  for (aa in syn_aas) {
    dM_mean[[aa]] <- colMeans(dM_samples[[aa]][post_samples, , drop = FALSE])
    dEta_mean[[aa]] <- colMeans(dEta_samples[[aa]][post_samples, , drop = FALSE])
  }
  
  if (verbose) {
    message("============================================")
    message("MCMC Complete!")
    message(sprintf("Final log posterior: %.2f", current_log_posterior))
    message(sprintf("Posterior mean sphi: %.3f", sphi_mean))
    message("============================================")
  }
  
  return(list(
    # Traces
    phi_samples = phi_samples,
    sphi_samples = sphi_samples,
    dM_samples = dM_samples,
    dEta_samples = dEta_samples,
    log_posterior = log_posterior_trace,
    
    # Posterior summaries
    phi_mean = phi_mean,
    phi_sd = phi_sd,
    sphi_mean = sphi_mean,
    dM_mean = dM_mean,
    dEta_mean = dEta_mean,
    
    # Settings
    n_samples = n_samples,
    thin = thin,
    burnin_frac = burnin_frac,
    gene_ids = rownames(codon_counts)
  ))
}

# ==============================================================================
# OUTPUT FUNCTIONS
# ==============================================================================

#' Convert dM/dEta list to data frame
#' @param param_list List of parameter vectors by amino acid
#' @param aa_to_codons Genetic code
#' @param param_name Name of parameter (dM or dEta)
#' @return Data frame with columns AA, Codon, Value
params_to_dataframe <- function(param_list, aa_to_codons, param_name = "Value") {
  syn_aas <- get_synonymous_aas()
  
  result <- data.frame(
    AA = character(),
    Codon = character(),
    Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (aa in syn_aas) {
    codons <- aa_to_codons[[aa]]
    n_codons <- length(codons)
    values <- c(param_list[[aa]], 0)  # Add 0 for reference codon
    
    for (i in 1:n_codons) {
      result <- rbind(result, data.frame(
        AA = aa,
        Codon = codons[i],
        Value = values[i],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  names(result)[3] <- param_name
  return(result)
}

#' Save MCMC results to files
#' @param results MCMC results from run_roc_mcmc
#' @param output_dir Output directory
save_results <- function(results, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  aa_to_codons <- get_genetic_code()
  
  # Save phi estimates
  phi_df <- data.frame(
    GeneID = results$gene_ids,
    Phi_Mean = results$phi_mean,
    Phi_SD = results$phi_sd
  )
  write.csv(phi_df, file.path(output_dir, "phi_estimates.csv"), row.names = FALSE)
  
  # Save dM estimates
  dM_df <- params_to_dataframe(results$dM_mean, aa_to_codons, "DeltaM")
  write.csv(dM_df, file.path(output_dir, "dM_estimates.csv"), row.names = FALSE)
  
  # Save dEta estimates
  dEta_df <- params_to_dataframe(results$dEta_mean, aa_to_codons, "DeltaEta")
  write.csv(dEta_df, file.path(output_dir, "dEta_estimates.csv"), row.names = FALSE)
  
  # Save hyperparameters
  hyper_df <- data.frame(
    Parameter = c("sphi"),
    Mean = c(results$sphi_mean)
  )
  write.csv(hyper_df, file.path(output_dir, "hyperparameters.csv"), row.names = FALSE)
  
  # Save log posterior trace
  write.csv(
    data.frame(Sample = seq_along(results$log_posterior), LogPosterior = results$log_posterior),
    file.path(output_dir, "log_posterior_trace.csv"), 
    row.names = FALSE
  )
  
  message(paste("Results saved to:", output_dir))
}

# ==============================================================================
# COMMAND LINE INTERFACE
# ==============================================================================

# Only run CLI when executed directly (not when sourced)
# Check if this script is the main script being run
.is_main_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(FALSE)
  script_name <- basename(sub("^--file=", "", file_arg))
  return(grepl("ROC_model_pure_R\\.R$", script_name))
}

if (!interactive() && .is_main_script()) {
  library(argparse)
  
  parser <- ArgumentParser(description = "Pure R ROC Model MCMC")
  
  parser$add_argument("-i", "--input", required = TRUE,
                      help = "Input FASTA file with CDS sequences")
  parser$add_argument("-o", "--output", default = "./roc_results",
                      help = "Output directory")
  parser$add_argument("-s", "--samples", type = "integer", default = 1000,
                      help = "Number of MCMC samples")
  parser$add_argument("-t", "--thin", type = "integer", default = 10,
                      help = "Thinning interval")
  parser$add_argument("-c", "--cores", type = "integer", default = 0,
                      help = "Number of parallel cores (0 = auto-detect)")
  parser$add_argument("--phi", default = NULL,
                      help = "CSV file with observed phi values (GeneID, expression columns)")
  parser$add_argument("--phi_col", default = NULL,
                      help = "Column name for phi values (default: use column 2)")
  parser$add_argument("--fix_phi", action = "store_true",
                      help = "Fix phi at observed values (requires --phi)")
  parser$add_argument("--dM", default = NULL,
                      help = "CSV file with initial dM values (AA, Codon, dM columns)")
  parser$add_argument("--fix_dM", action = "store_true",
                      help = "Fix dM at initial values (requires --dM)")
  
  args <- parser$parse_args()
  
  # Set up parallel cores
  n_cores <- args$cores
  if (n_cores <= 0) {
    n_cores <- detect_cores()
    message(sprintf("Auto-detected %d cores", n_cores))
  }
  
  # Read data
  message("Reading FASTA file...")
  data <- read_fasta_codon_counts(args$input)
  
  # Read observed phi if provided
  obs_phi <- NULL
  if (!is.null(args$phi)) {
    message("Reading observed phi values...")
    phi_df <- read.csv(args$phi)
    
    # Determine which column to use for phi
    if (!is.null(args$phi_col) && args$phi_col %in% colnames(phi_df)) {
      phi_col <- args$phi_col
    } else {
      phi_col <- colnames(phi_df)[2]  # Default: second column
    }
    message(sprintf("Using column '%s' for phi values", phi_col))
    
    # Match to gene order
    gene_col <- colnames(phi_df)[1]  # First column is gene ID
    obs_phi <- phi_df[match(data$gene_ids, phi_df[[gene_col]]), phi_col]
    obs_phi[is.na(obs_phi)] <- 1  # Default for missing
  }
  
  # Read initial dM if provided
  init_dM <- NULL
  if (!is.null(args$dM)) {
    message("Reading initial dM values...")
    dM_df <- read.csv(args$dM)
    aa_to_codons <- get_genetic_code()
    syn_aas <- get_synonymous_aas()
    
    # Handle column name variations (dM, DM, DeltaM)
    dM_col <- NULL
    for (col_name in c("dM", "DM", "DeltaM", "deltaM")) {
      if (col_name %in% colnames(dM_df)) {
        dM_col <- col_name
        break
      }
    }
    if (is.null(dM_col)) {
      stop("Could not find dM column in input file. Expected one of: dM, DM, DeltaM, deltaM")
    }
    message(sprintf("Using column '%s' for dM values", dM_col))
    
    init_dM <- list()
    for (aa in syn_aas) {
      codons <- aa_to_codons[[aa]]
      non_ref_codons <- codons[-length(codons)]  # Exclude reference
      
      values <- numeric(length(non_ref_codons))
      for (i in seq_along(non_ref_codons)) {
        row_idx <- which(dM_df$Codon == non_ref_codons[i])
        if (length(row_idx) > 0) {
          values[i] <- dM_df[[dM_col]][row_idx[1]]
        }
      }
      init_dM[[aa]] <- values
    }
  }
  
  # Run MCMC
  results <- run_roc_mcmc(
    codon_counts = data$codon_counts,
    n_samples = args$samples,
    thin = args$thin,
    obs_phi = obs_phi,
    fix_phi = args$fix_phi,
    fix_dM = args$fix_dM,
    init_dM = init_dM,
    n_cores = n_cores,
    verbose = TRUE
  )
  
  # Save results
  save_results(results, args$output)
}
