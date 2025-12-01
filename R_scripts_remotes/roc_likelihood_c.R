# ==============================================================================
# ROC Likelihood - C-accelerated Functions
#
# This file provides R wrappers for the C implementations of the ROC model
# likelihood calculations. Falls back to pure R if C library not available.
#
# Usage:
#   source("roc_likelihood_c.R")
#   # Functions like calc_total_log_likelihood_c() are now available
#
# Compile the C code first with:
#   system("R CMD SHLIB src/roc_likelihood.c -o src/roc_likelihood.so")
#
# Author: Luis Javier Madrigal-Roca
# Date: 2025-12-01
# ==============================================================================

# Try to load the shared library
.load_roc_c <- function() {
  # Find the shared library
  script_dir <- dirname(sys.frame(1)$ofile)
  if (is.null(script_dir) || script_dir == "") {
    script_dir <- getwd()
  }
  
  # Look for .so (Linux/Mac) or .dll (Windows)
  so_paths <- c(
    file.path(script_dir, "src", "roc_likelihood.so"),
    file.path(script_dir, "roc_likelihood.so"),
    file.path("R_scripts_remotes", "src", "roc_likelihood.so"),
    file.path("src", "roc_likelihood.so"),
    "roc_likelihood.so"
  )
  
  dll_paths <- gsub("\\.so$", ".dll", so_paths)
  all_paths <- c(so_paths, dll_paths)
  
  for (path in all_paths) {
    if (file.exists(path)) {
      tryCatch({
        dyn.load(path)
        message(sprintf("Loaded C library: %s", path))
        return(TRUE)
      }, error = function(e) {
        message(sprintf("Failed to load %s: %s", path, e$message))
      })
    }
  }
  
  message("C library not found. Using pure R (slower).")
  message("To compile: R CMD SHLIB src/roc_likelihood.c -o src/roc_likelihood.so")
  return(FALSE)
}

# Global flag for whether C code is available
USE_C_CODE <- FALSE

# ==============================================================================
# Amino Acid Info Structure
# 
# Creates a matrix that tells C code about the codon structure:
# - Which column indices correspond to which amino acid
# - How many codons per amino acid
# - Where in the parameter vectors to find dM/dEta values
# ==============================================================================

#' Build amino acid info matrix for C code
#' @param codon_counts Matrix with codon names as column names
#' @param aa_to_codons Genetic code mapping
#' @return Matrix with columns [aa_idx, start_codon, n_codons, start_param]
build_aa_codon_info <- function(codon_counts, aa_to_codons) {
  syn_aas <- get_synonymous_aas()
  all_codons <- colnames(codon_counts)
  
  n_aa <- length(syn_aas)
  info <- matrix(0L, nrow = n_aa, ncol = 4)
  colnames(info) <- c("aa_idx", "start_codon", "n_codons", "start_param")
  
  param_offset <- 0L
  
  for (i in seq_along(syn_aas)) {
    aa <- syn_aas[i]
    codons <- aa_to_codons[[aa]]
    
    # Find column indices for these codons (0-indexed for C)
    codon_indices <- match(codons, all_codons) - 1L
    
    info[i, "aa_idx"] <- i - 1L
    info[i, "start_codon"] <- codon_indices[1]
    info[i, "n_codons"] <- length(codons)
    info[i, "start_param"] <- param_offset
    
    # Parameters are for non-reference codons (n_codons - 1)
    param_offset <- param_offset + length(codons) - 1L
  }
  
  return(info)
}

#' Flatten dM/dEta list to vector (for C code)
#' @param param_list List of parameter vectors by amino acid
#' @return Numeric vector
flatten_params <- function(param_list) {
  syn_aas <- get_synonymous_aas()
  unlist(param_list[syn_aas])
}

#' Unflatten parameter vector back to list
#' @param param_vec Flattened parameter vector
#' @param aa_to_codons Genetic code
#' @return List of parameter vectors by amino acid
unflatten_params <- function(param_vec, aa_to_codons) {
  syn_aas <- get_synonymous_aas()
  result <- list()
  offset <- 1
  
  for (aa in syn_aas) {
    n_params <- length(aa_to_codons[[aa]]) - 1
    result[[aa]] <- param_vec[offset:(offset + n_params - 1)]
    offset <- offset + n_params
  }
  
  return(result)
}

# ==============================================================================
# C-Accelerated Likelihood Functions
# ==============================================================================

#' Calculate log likelihood for all genes using C code
#' 
#' @param codon_counts Integer matrix (n_genes x n_codons)
#' @param dM_list List of dM vectors
#' @param dEta_list List of dEta vectors  
#' @param phi Numeric vector of expression levels
#' @param aa_to_codons Genetic code mapping
#' @param n_cores Number of cores (ignored - C code is single-threaded but fast)
#' @return Total log likelihood
calc_total_log_likelihood_c <- function(codon_counts, dM_list, dEta_list,
                                         phi, aa_to_codons, n_cores = 1) {
  if (!USE_C_CODE) {
    # Fall back to pure R
    return(calc_total_log_likelihood(codon_counts, dM_list, dEta_list,
                                     phi, aa_to_codons, n_cores))
  }
  
  # Prepare data for C
  aa_info <- build_aa_codon_info(codon_counts, aa_to_codons)
  dM_vec <- flatten_params(dM_list)
  dEta_vec <- flatten_params(dEta_list)
  n_aa <- nrow(aa_info)
  
  # Ensure integer matrix
  if (!is.integer(codon_counts)) {
    storage.mode(codon_counts) <- "integer"
  }
  
  # Call C function
  gene_logliks <- .Call("C_calc_log_lik_all_genes",
                        codon_counts, dM_vec, dEta_vec, phi,
                        aa_info, as.integer(n_aa))
  
  return(sum(gene_logliks))
}

#' Calculate log likelihood for one gene using C code
#' 
#' @param gene_codon_counts Named vector of codon counts
#' @param dM_list List of dM vectors
#' @param dEta_list List of dEta vectors
#' @param phi Expression level for this gene
#' @param aa_to_codons Genetic code mapping
#' @return Log likelihood for gene
calc_log_likelihood_gene_c <- function(gene_codon_counts, dM_list, dEta_list,
                                        phi, aa_to_codons) {
  if (!USE_C_CODE) {
    return(calc_log_likelihood_gene(gene_codon_counts, dM_list, dEta_list,
                                    phi, aa_to_codons))
  }
  
  # Build aa_info for single gene (need codon ordering)
  syn_aas <- get_synonymous_aas()
  all_codons <- names(gene_codon_counts)
  
  n_aa <- length(syn_aas)
  aa_info <- matrix(0L, nrow = n_aa, ncol = 4)
  param_offset <- 0L
  
  for (i in seq_along(syn_aas)) {
    aa <- syn_aas[i]
    codons <- aa_to_codons[[aa]]
    codon_indices <- match(codons, all_codons) - 1L
    
    aa_info[i, 1] <- i - 1L
    aa_info[i, 2] <- codon_indices[1]
    aa_info[i, 3] <- length(codons)
    aa_info[i, 4] <- param_offset
    param_offset <- param_offset + length(codons) - 1L
  }
  
  dM_vec <- flatten_params(dM_list)
  dEta_vec <- flatten_params(dEta_list)
  
  counts <- as.integer(gene_codon_counts)
  
  return(.Call("C_calc_log_lik_one_gene",
               counts, dM_vec, dEta_vec, phi, aa_info, as.integer(n_aa)))
}

#' Batch update phi values using C code
#' 
#' This is the main speedup - does all gene phi updates in C without
#' R interpreter overhead.
#'
#' @param codon_counts Integer matrix
#' @param dM_list List of dM vectors
#' @param dEta_list List of dEta vectors
#' @param phi Current phi values
#' @param obs_phi_matrix Observed phi matrix (or NULL)
#' @param sepsilon Observation noise vector
#' @param sphi Phi prior SD
#' @param prop_sd Proposal SD vector
#' @param aa_to_codons Genetic code
#' @param with_phi Include observation likelihood
#' @param gene_indices Which genes to update (1-indexed)
#' @return List with updated phi and acceptance indicators
batch_update_phi_c <- function(codon_counts, dM_list, dEta_list, phi,
                                obs_phi_matrix, sepsilon, sphi, prop_sd,
                                aa_to_codons, with_phi, gene_indices = NULL) {
  
  if (is.null(gene_indices)) {
    gene_indices <- seq_len(nrow(codon_counts))
  }
  
  if (!USE_C_CODE) {
    # Fall back to pure R (slow)
    return(.batch_update_phi_r(codon_counts, dM_list, dEta_list, phi,
                               obs_phi_matrix, sepsilon, sphi, prop_sd,
                               aa_to_codons, with_phi, gene_indices))
  }
  
  aa_info <- build_aa_codon_info(codon_counts, aa_to_codons)
  dM_vec <- flatten_params(dM_list)
  dEta_vec <- flatten_params(dEta_list)
  n_aa <- nrow(aa_info)
  
  if (!is.integer(codon_counts)) {
    storage.mode(codon_counts) <- "integer"
  }
  
  # Handle NULL obs_phi_matrix
  if (is.null(obs_phi_matrix)) {
    obs_phi_matrix <- matrix(NA_real_, nrow = length(phi), ncol = 1)
    sepsilon <- 1.0
    with_phi <- FALSE
  }
  
  result <- .Call("C_batch_update_phi",
                  codon_counts, dM_vec, dEta_vec, phi,
                  obs_phi_matrix, as.numeric(sepsilon),
                  as.numeric(sphi), as.numeric(prop_sd),
                  aa_info, as.integer(n_aa),
                  as.logical(with_phi), as.integer(gene_indices))
  
  return(result)
}

#' Pure R fallback for batch_update_phi
.batch_update_phi_r <- function(codon_counts, dM_list, dEta_list, phi,
                                 obs_phi_matrix, sepsilon, sphi, prop_sd,
                                 aa_to_codons, with_phi, gene_indices) {
  n_tissues <- if (!is.null(obs_phi_matrix)) ncol(obs_phi_matrix) else 0
  accept <- integer(length(gene_indices))
  phi_new <- phi
  mu <- -sphi^2 / 2
  
  for (i in seq_along(gene_indices)) {
    g <- gene_indices[i]
    
    phi_prop <- exp(rnorm(1, log(phi[g]), prop_sd[g]))
    if (phi_prop <= 0) next
    
    # Codon likelihood
    old_lik <- calc_log_likelihood_gene(codon_counts[g, ], dM_list, dEta_list,
                                        phi[g], aa_to_codons)
    new_lik <- calc_log_likelihood_gene(codon_counts[g, ], dM_list, dEta_list,
                                        phi_prop, aa_to_codons)
    
    # Observation likelihood
    old_obs <- 0
    new_obs <- 0
    if (with_phi && n_tissues > 0) {
      for (k in seq_len(n_tissues)) {
        obs_k <- obs_phi_matrix[g, k]
        if (!is.na(obs_k) && obs_k > 0) {
          old_obs <- old_obs + dnorm(log(obs_k), log(phi[g]), sepsilon[k], log = TRUE)
          new_obs <- new_obs + dnorm(log(obs_k), log(phi_prop), sepsilon[k], log = TRUE)
        }
      }
    }
    
    # Prior
    old_prior <- dlnorm(phi[g], mu, sphi, log = TRUE)
    new_prior <- dlnorm(phi_prop, mu, sphi, log = TRUE)
    
    # Jacobian
    log_jacobian <- log(phi_prop) - log(phi[g])
    
    log_alpha <- (new_lik - old_lik) + (new_obs - old_obs) + 
                 (new_prior - old_prior) + log_jacobian
    
    if (log(runif(1)) < log_alpha) {
      phi_new[g] <- phi_prop
      accept[i] <- 1L
    }
  }
  
  return(list(phi = phi_new, accept = accept))
}

#' Calculate observation likelihood for multi-tissue phi using C
calc_obs_phi_likelihood_c <- function(obs_phi_matrix, true_phi, sepsilon) {
  if (!USE_C_CODE) {
    return(log_lik_obs_phi_multitissue(obs_phi_matrix, true_phi, sepsilon))
  }
  
  return(.Call("C_log_lik_obs_phi", obs_phi_matrix, true_phi, as.numeric(sepsilon)))
}

# ==============================================================================
# Initialization
# ==============================================================================

# Try to load C code when this file is sourced
USE_C_CODE <- .load_roc_c()

#' Check if C acceleration is available
is_c_available <- function() {
  USE_C_CODE
}

#' Compile the C code
compile_roc_c <- function(src_dir = NULL) {
  if (is.null(src_dir)) {
    src_dir <- file.path(dirname(sys.frame(1)$ofile), "src")
    if (!dir.exists(src_dir)) {
      src_dir <- "R_scripts_remotes/src"
    }
  }
  
  c_file <- file.path(src_dir, "roc_likelihood.c")
  if (!file.exists(c_file)) {
    stop(sprintf("C source file not found: %s", c_file))
  }
  
  old_wd <- getwd()
  setwd(src_dir)
  
  tryCatch({
    system2("R", c("CMD", "SHLIB", "roc_likelihood.c"), 
            stdout = TRUE, stderr = TRUE)
    message("Compilation successful!")
    
    # Try to load
    USE_C_CODE <<- .load_roc_c()
  }, finally = {
    setwd(old_wd)
  })
  
  return(USE_C_CODE)
}

# Print status
if (USE_C_CODE) {
  message("ROC likelihood: C acceleration ENABLED")
} else {
  message("ROC likelihood: Using pure R (run compile_roc_c() to enable C acceleration)")
}
