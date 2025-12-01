/*
 * roc_likelihood.c - Optimized C implementations for ROC model likelihood
 *
 * These functions are the bottleneck operations from the pure R ROC MCMC.
 * Provides ~10-50x speedup over pure R for likelihood calculations.
 *
 * Compile with: R CMD SHLIB roc_likelihood.c
 * Or in R: system("R CMD SHLIB R_scripts_remotes/src/roc_likelihood.c")
 *
 * Author: Luis Javier Madrigal-Roca
 * Date: 2025-12-01
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* 
 * Log-sum-exp trick for numerical stability
 * Computes log(sum(exp(x))) without overflow
 */
static inline double log_sum_exp(double *x, int n) {
    if (n == 0) return R_NegInf;
    
    // Find maximum
    double max_val = x[0];
    for (int i = 1; i < n; i++) {
        if (x[i] > max_val) max_val = x[i];
    }
    
    if (!R_FINITE(max_val)) return max_val;
    
    // Compute sum of exp(x - max)
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += exp(x[i] - max_val);
    }
    
    return max_val + log(sum);
}

/*
 * Calculate log likelihood for one amino acid in one gene
 * 
 * codon_counts: counts for each codon (length n_codons)
 * dM: mutation bias for non-reference codons (length n_codons - 1)
 * dEta: selection coefficients for non-reference codons (length n_codons - 1)
 * phi: expression level
 * n_codons: number of codons for this amino acid
 * 
 * Returns: log likelihood contribution
 */
static double calc_log_lik_aa_c(int *codon_counts, double *dM, double *dEta,
                                 double phi, int n_codons) {
    if (n_codons <= 1) return 0.0;
    
    // Check if total count is zero
    int total = 0;
    for (int i = 0; i < n_codons; i++) {
        total += codon_counts[i];
    }
    if (total == 0) return 0.0;
    
    // Calculate log numerators for each codon
    // log P(codon_i) = -dM_i - dEta_i * phi - log(Z)
    // Reference codon (last) has dM=0, dEta=0
    double *log_numer = (double *)R_alloc(n_codons, sizeof(double));
    
    for (int i = 0; i < n_codons - 1; i++) {
        log_numer[i] = -dM[i] - dEta[i] * phi;
    }
    log_numer[n_codons - 1] = 0.0;  // Reference codon
    
    // Partition function (normalizer)
    double log_Z = log_sum_exp(log_numer, n_codons);
    
    // Log probabilities
    double log_lik = 0.0;
    for (int i = 0; i < n_codons; i++) {
        if (codon_counts[i] > 0) {
            log_lik += codon_counts[i] * (log_numer[i] - log_Z);
        }
    }
    
    return log_lik;
}

/*
 * SEXP wrapper: Calculate log likelihood for all genes
 *
 * Arguments:
 *   codon_counts_matrix: integer matrix (n_genes x n_codons)
 *   dM_vec: numeric vector of all dM values (concatenated across AAs)
 *   dEta_vec: numeric vector of all dEta values (concatenated across AAs)
 *   phi_vec: numeric vector of expression levels (length n_genes)
 *   aa_codon_info: integer matrix with columns [aa_idx, start_codon, n_codons, start_param]
 *   n_aa: number of amino acids
 *
 * Returns: numeric vector of log likelihoods per gene
 */
SEXP C_calc_log_lik_all_genes(SEXP codon_counts_matrix, SEXP dM_vec, SEXP dEta_vec,
                               SEXP phi_vec, SEXP aa_codon_info, SEXP n_aa_sexp) {
    // Get dimensions
    int n_genes = Rf_nrows(codon_counts_matrix);
    int n_codons_total = Rf_ncols(codon_counts_matrix);
    int n_aa = INTEGER(n_aa_sexp)[0];
    
    // Get pointers
    int *counts = INTEGER(codon_counts_matrix);
    double *dM = REAL(dM_vec);
    double *dEta = REAL(dEta_vec);
    double *phi = REAL(phi_vec);
    int *aa_info = INTEGER(aa_codon_info);  // 4 columns: aa_idx, start_codon, n_codons, start_param
    
    // Allocate output
    SEXP result = PROTECT(Rf_allocVector(REALSXP, n_genes));
    double *log_liks = REAL(result);
    
    // Process each gene
    for (int g = 0; g < n_genes; g++) {
        double gene_log_lik = 0.0;
        
        // Process each amino acid
        for (int aa = 0; aa < n_aa; aa++) {
            int start_codon = aa_info[aa + n_aa];      // Column 1: start codon index
            int n_codons = aa_info[aa + 2 * n_aa];     // Column 2: number of codons
            int start_param = aa_info[aa + 3 * n_aa];  // Column 3: start param index
            
            if (n_codons <= 1) continue;
            
            // Get codon counts for this gene and AA
            int *gene_counts = (int *)R_alloc(n_codons, sizeof(int));
            for (int c = 0; c < n_codons; c++) {
                // R matrices are column-major
                gene_counts[c] = counts[(start_codon + c) * n_genes + g];
            }
            
            // Get dM and dEta for this AA (n_codons - 1 parameters)
            double *aa_dM = dM + start_param;
            double *aa_dEta = dEta + start_param;
            
            gene_log_lik += calc_log_lik_aa_c(gene_counts, aa_dM, aa_dEta,
                                               phi[g], n_codons);
        }
        
        log_liks[g] = gene_log_lik;
    }
    
    UNPROTECT(1);
    return result;
}

/*
 * Calculate log likelihood for a SINGLE gene (for Metropolis-Hastings updates)
 * 
 * This is used in the inner loop when updating phi gene by gene.
 */
SEXP C_calc_log_lik_one_gene(SEXP codon_counts_vec, SEXP dM_vec, SEXP dEta_vec,
                              SEXP phi_sexp, SEXP aa_codon_info, SEXP n_aa_sexp) {
    int n_aa = INTEGER(n_aa_sexp)[0];
    
    int *counts = INTEGER(codon_counts_vec);
    double *dM = REAL(dM_vec);
    double *dEta = REAL(dEta_vec);
    double phi = REAL(phi_sexp)[0];
    int *aa_info = INTEGER(aa_codon_info);
    
    double log_lik = 0.0;
    
    for (int aa = 0; aa < n_aa; aa++) {
        int start_codon = aa_info[aa + n_aa];
        int n_codons = aa_info[aa + 2 * n_aa];
        int start_param = aa_info[aa + 3 * n_aa];
        
        if (n_codons <= 1) continue;
        
        int *aa_counts = (int *)R_alloc(n_codons, sizeof(int));
        for (int c = 0; c < n_codons; c++) {
            aa_counts[c] = counts[start_codon + c];
        }
        
        log_lik += calc_log_lik_aa_c(aa_counts, dM + start_param,
                                      dEta + start_param, phi, n_codons);
    }
    
    return Rf_ScalarReal(log_lik);
}

/*
 * Fast batch update of phi for multiple genes
 * 
 * Performs Metropolis-Hastings updates for a batch of genes at once,
 * which is more cache-friendly than gene-by-gene in R.
 *
 * Arguments:
 *   codon_counts_matrix: integer matrix (n_genes x n_codons)
 *   dM_vec, dEta_vec: parameter vectors
 *   phi_vec: current phi values (modified in place)
 *   obs_phi_matrix: observed phi (n_genes x n_tissues) or NULL
 *   sepsilon_vec: observation noise per tissue
 *   sphi: std dev of phi prior
 *   prop_sd_vec: proposal std dev per gene
 *   aa_codon_info, n_aa_sexp: amino acid structure info
 *   with_phi_sexp: include observation likelihood
 *   gene_indices: which genes to update (1-indexed)
 *
 * Returns: list(phi = updated phi, accept = acceptance vector)
 */
SEXP C_batch_update_phi(SEXP codon_counts_matrix, SEXP dM_vec, SEXP dEta_vec,
                         SEXP phi_vec, SEXP obs_phi_matrix, SEXP sepsilon_vec,
                         SEXP sphi_sexp, SEXP prop_sd_vec,
                         SEXP aa_codon_info, SEXP n_aa_sexp,
                         SEXP with_phi_sexp, SEXP gene_indices) {
    
    GetRNGstate();
    
    int n_genes = Rf_nrows(codon_counts_matrix);
    int n_codons_total = Rf_ncols(codon_counts_matrix);
    int n_aa = INTEGER(n_aa_sexp)[0];
    int n_update = Rf_length(gene_indices);
    int with_phi = Rf_asLogical(with_phi_sexp);
    
    int n_tissues = 0;
    if (!Rf_isNull(obs_phi_matrix) && with_phi) {
        n_tissues = Rf_ncols(obs_phi_matrix);
    }
    
    int *counts = INTEGER(codon_counts_matrix);
    double *dM = REAL(dM_vec);
    double *dEta = REAL(dEta_vec);
    double *phi = REAL(phi_vec);
    double *prop_sd = REAL(prop_sd_vec);
    double *sepsilon = (n_tissues > 0) ? REAL(sepsilon_vec) : NULL;
    double *obs_phi = (n_tissues > 0) ? REAL(obs_phi_matrix) : NULL;
    int *aa_info = INTEGER(aa_codon_info);
    int *idx = INTEGER(gene_indices);
    double sphi = REAL(sphi_sexp)[0];
    
    // Output: updated phi and acceptance indicators
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP phi_out = PROTECT(Rf_duplicate(phi_vec));
    SEXP accept_out = PROTECT(Rf_allocVector(INTSXP, n_update));
    double *phi_new = REAL(phi_out);
    int *accept = INTEGER(accept_out);
    
    // Prior mean for phi
    double mu_phi = -sphi * sphi / 2.0;
    
    for (int i = 0; i < n_update; i++) {
        int g = idx[i] - 1;  // Convert to 0-indexed
        accept[i] = 0;
        
        // Propose new phi (log-scale random walk)
        double log_phi_old = log(phi[g]);
        double log_phi_prop = log_phi_old + norm_rand() * prop_sd[g];
        double phi_prop = exp(log_phi_prop);
        
        if (phi_prop <= 0 || !R_FINITE(phi_prop)) continue;
        
        // Calculate old and new codon log likelihood for this gene
        double old_lik = 0.0, new_lik = 0.0;
        
        for (int aa = 0; aa < n_aa; aa++) {
            int start_codon = aa_info[aa + n_aa];
            int n_codons = aa_info[aa + 2 * n_aa];
            int start_param = aa_info[aa + 3 * n_aa];
            
            if (n_codons <= 1) continue;
            
            int *aa_counts = (int *)R_alloc(n_codons, sizeof(int));
            for (int c = 0; c < n_codons; c++) {
                aa_counts[c] = counts[(start_codon + c) * n_genes + g];
            }
            
            old_lik += calc_log_lik_aa_c(aa_counts, dM + start_param,
                                          dEta + start_param, phi[g], n_codons);
            new_lik += calc_log_lik_aa_c(aa_counts, dM + start_param,
                                          dEta + start_param, phi_prop, n_codons);
        }
        
        // Observation likelihood (multi-tissue)
        double old_obs_lik = 0.0, new_obs_lik = 0.0;
        if (with_phi && n_tissues > 0) {
            for (int k = 0; k < n_tissues; k++) {
                double obs_k = obs_phi[k * n_genes + g];
                if (!ISNA(obs_k) && obs_k > 0) {
                    old_obs_lik += dnorm(log(obs_k), log(phi[g]), sepsilon[k], 1);
                    new_obs_lik += dnorm(log(obs_k), log(phi_prop), sepsilon[k], 1);
                }
            }
        }
        
        // Prior
        double old_prior = dlnorm(phi[g], mu_phi, sphi, 1);
        double new_prior = dlnorm(phi_prop, mu_phi, sphi, 1);
        
        // Jacobian for log-scale proposal
        double log_jacobian = log_phi_prop - log_phi_old;
        
        // Acceptance ratio
        double log_alpha = (new_lik - old_lik) + 
                           (new_obs_lik - old_obs_lik) +
                           (new_prior - old_prior) + log_jacobian;
        
        if (log(unif_rand()) < log_alpha) {
            phi_new[g] = phi_prop;
            accept[i] = 1;
        }
    }
    
    SET_VECTOR_ELT(result, 0, phi_out);
    SET_VECTOR_ELT(result, 1, accept_out);
    
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("phi"));
    SET_STRING_ELT(names, 1, Rf_mkChar("accept"));
    Rf_setAttrib(result, R_NamesSymbol, names);
    
    PutRNGstate();
    
    UNPROTECT(4);
    return result;
}

/*
 * Calculate observation likelihood for multi-tissue phi
 */
SEXP C_log_lik_obs_phi(SEXP obs_phi_matrix, SEXP true_phi_vec, SEXP sepsilon_vec) {
    int n_genes = Rf_length(true_phi_vec);
    int n_tissues = Rf_ncols(obs_phi_matrix);
    
    double *obs_phi = REAL(obs_phi_matrix);
    double *true_phi = REAL(true_phi_vec);
    double *sepsilon = REAL(sepsilon_vec);
    
    double total_log_lik = 0.0;
    
    for (int k = 0; k < n_tissues; k++) {
        for (int g = 0; g < n_genes; g++) {
            double obs = obs_phi[k * n_genes + g];
            if (!ISNA(obs) && obs > 0 && true_phi[g] > 0) {
                total_log_lik += dnorm(log(obs), log(true_phi[g]), sepsilon[k], 1);
            }
        }
    }
    
    return Rf_ScalarReal(total_log_lik);
}

/*
 * Registration table for .Call
 */
static const R_CallMethodDef CallEntries[] = {
    {"C_calc_log_lik_all_genes", (DL_FUNC) &C_calc_log_lik_all_genes, 6},
    {"C_calc_log_lik_one_gene", (DL_FUNC) &C_calc_log_lik_one_gene, 6},
    {"C_batch_update_phi", (DL_FUNC) &C_batch_update_phi, 12},
    {"C_log_lik_obs_phi", (DL_FUNC) &C_log_lik_obs_phi, 3},
    {NULL, NULL, 0}
};

void R_init_roclikelihood(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
