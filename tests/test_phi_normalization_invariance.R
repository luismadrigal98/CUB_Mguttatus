#!/usr/bin/env Rscript
#
# Mathematical proof + numerical verification that ROC_eff is invariant to
# multiplicative scaling of phi (observed expression).
#
# Background
# ----------
# The AnaCoDa ROC model pipeline divides all CPM values by their grand mean
# (~60.85) before passing them to the MCMC.  This raises the question: does
# scaling phi change the estimated selection coefficients and, ultimately, the
# ROC_eff values used downstream?
#
# Mathematical proof
# ------------------
# The ROC model codon probability for codon i of synonymous family a, given
# synthesis rate phi, mutation mu, and selection eta, is:
#
#   P(codon_i | phi, mu, eta)  prop  exp( -mu_i  -  eta_i * phi )
#
# Claim: the likelihood is EXACTLY invariant under the joint substitution
#   phi'  = phi / k     (scaling observed AND estimated phi by 1/k)
#   eta'  = eta * k
# for any k > 0.
#
# Proof (one line):
#   -eta'_i * phi'  =  -(k * eta_i) * (phi / k)  =  -eta_i * phi   QED.
#
# Consequence for ROC_eff
# ---------------------
# AnaCoDa's calculateSelectionCoefficients computes (for each codon i):
#
#   S_i  =  phi_estimated * ( eta_i - min_a(eta) )
#
# Under the scaling above:
#
#   S'_i  =  phi' * ( eta'_i - min(eta') )
#          =  (phi / k) * k * ( eta_i - min(eta) )
#          =  phi * ( eta_i - min(eta) )
#          =  S_i
#
# ROC_eff is therefore EXACTLY invariant to multiplicative scaling of phi.
#
# Observation noise model
# -----------------------
# The model is:  log(obsPhi) ~ Normal( log(phi_est) + A_phi, sEpsilon )
#
# When BOTH obsPhi and phi_est are scaled by 1/k:
#
#   log(obsPhi/k) - log(phi_est/k) - A_phi
#   = ( log(obsPhi) - log(k) ) - ( log(phi_est) - log(k) ) - A_phi
#   = log(obsPhi) - log(phi_est) - A_phi        (identical to original)
#
# So the observation likelihood is EXACTLY invariant with noiseOffset (A_phi)
# UNCHANGED -- both sides shift by -log(k) and cancel.
#
# Usage
# -----
#   Rscript tests/test_phi_normalization_invariance.R
# ____________________________________________________________________________

cat("=================================================================\n")
cat("Test: ROC_eff invariance to phi normalization\n")
cat("=================================================================\n\n")

n_pass <- 0L
n_fail <- 0L

assert_near <- function(label, a, b, tol = 1e-10) {
  diff <- max(abs(a - b))
  if (diff < tol) {
    cat(sprintf("  PASS  %s  (max|diff| = %.2e)\n", label, diff))
    n_pass <<- n_pass + 1L
  } else {
    cat(sprintf("  FAIL  %s  (max|diff| = %.2e, tol = %.2e)\n",
                label, diff, tol))
    n_fail <<- n_fail + 1L
  }
}

# ── ROC log-likelihood for one amino-acid family, one gene ──────────────────
#
# counts : integer vector length n_codons (last = reference)
# phi    : scalar synthesis rate
# mu     : numeric vector length (n_codons - 1)  [reference mu = 0]
# eta    : numeric vector length (n_codons - 1)  [reference eta = 0]
#
roc_log_lik <- function(counts, phi, mu, eta) {
  log_p_unnorm <- c(-mu - eta * phi, 0)         # reference codon = 0
  log_z        <- log(sum(exp(log_p_unnorm - max(log_p_unnorm)))) +
                  max(log_p_unnorm)
  sum(counts * (log_p_unnorm - log_z))
}

# ── ROC_eff per-codon for one gene ─────────────────────────────────────────────
#
# phi_est : estimated synthesis rate (scalar)
# eta     : full selection vector length n_codons (reference appended as 0)
#
roc_eff_vector <- function(phi_est, eta) {
  phi_est * (eta - min(eta))
}


# ── Test parameters ──────────────────────────────────────────────────────────
set.seed(2024)

n_codons <- 4L  # synonymous family size

# True parameters on CPM scale (realistic Mimulus range)
mu_true  <- c(-0.6,  0.3, -0.1)          # mutation  (n_codons - 1)
eta_true <- c(-0.02, -0.008, -0.001)     # selection (n_codons - 1)

phi_true <- c(0.5, 1, 5, 20, 60, 200, 500, 1000, 5000, 37000)

# k values: actual pipeline value plus edge cases
k_values <- c(60.85, 1000, mean(phi_true), 1)

# Simulate codon counts under the true model
sim_counts <- vapply(phi_true, function(phi) {
  lp <- c(-mu_true - eta_true * phi, 0)
  lp <- lp - max(lp)
  p  <- exp(lp)
  p  <- pmax(p / sum(p), .Machine$double.eps)
  as.integer(rmultinom(1L, size = 300L, prob = p))
}, integer(n_codons))   # n_codons x n_genes


# ── Section 1: Codon likelihood is exactly invariant ────────────────────────
cat("--- Section 1: Codon likelihood invariance ---\n")

for (k in k_values) {
  llik_orig   <- mapply(roc_log_lik,
                        as.data.frame(sim_counts), phi_true,
                        MoreArgs = list(mu = mu_true, eta = eta_true))
  llik_scaled <- mapply(roc_log_lik,
                        as.data.frame(sim_counts), phi_true / k,
                        MoreArgs = list(mu = mu_true, eta = eta_true * k))
  assert_near(
    sprintf("log-likelihood invariant  (k = %8.2f)", k),
    llik_orig, llik_scaled
  )
}


# ── Section 2: ROC_eff = phi * eta is exactly invariant ───────────────────────
cat("\n--- Section 2: ROC_eff invariance ---\n")

for (k in k_values) {
  eta_full    <- c(eta_true, 0)
  eta_scaled  <- c(eta_true * k, 0)

  roc_eff_orig   <- vapply(phi_true,
                        roc_eff_vector, numeric(n_codons), eta = eta_full)
  roc_eff_scaled <- vapply(phi_true / k,
                        roc_eff_vector, numeric(n_codons), eta = eta_scaled)
  assert_near(
    sprintf("ROC_eff invariant           (k = %8.2f)", k),
    roc_eff_orig, roc_eff_scaled
  )
}


# ── Section 3: Observation noise likelihood invariant (same noiseOffset) ────
cat("\n--- Section 3: Observation noise model invariance ---\n")
#
# The key identity: both obsPhi and phi_est are scaled by 1/k, so the
# deviation log(obsPhi) - log(phi_est) is unchanged.  noiseOffset does NOT
# need to change.
#
s_epsilon    <- 0.5
noise_offset <- 0.3

# Fixed observed phi (same for both versions of the test)
set.seed(42)
obs_phi_vals <- phi_true * exp(rnorm(length(phi_true), noise_offset, s_epsilon))

for (k in k_values) {
  obs_log_lik_orig <- dnorm(
    log(obs_phi_vals),
    log(phi_true)      + noise_offset,
    s_epsilon, log = TRUE
  )
  obs_log_lik_scaled <- dnorm(
    log(obs_phi_vals / k),          # obsPhi scaled by 1/k
    log(phi_true / k) + noise_offset,  # phi_est also scaled by 1/k; offset unchanged
    s_epsilon, log = TRUE
  )
  assert_near(
    sprintf("Obs. lik. invariant       (k = %8.2f)", k),
    obs_log_lik_orig, obs_log_lik_scaled
  )
}


# ── Section 4: Analytical gradient rescaling (chain rule) ───────────────────
cat("\n--- Section 4: Analytical gradient rescaling ---\n")
#
# Closed-form score d/d(phi) log L for the ROC model (one AA family):
#
#   score = -sum_i n_i * eta_i  +  N * E[eta | phi, mu, eta]
#
# where the expectation is over the model probabilities P(codon_i | phi).
#
# Under the substitution phi' = phi/k, eta' = k*eta, the model probabilities
# are UNCHANGED (Section 1 proof), so E[eta' | phi'] = k * E[eta | phi].
# Therefore:
#
#   score(phi', eta')
#   = -sum_i n_i * (k*eta_i) + N * E[k*eta | phi']
#   = k * ( -sum_i n_i * eta_i + N * E[eta | phi] )
#   = k * score(phi, eta)       QED.
#
# This test verifies the identity exactly (no finite differences involved).

roc_score <- function(counts, phi, mu, eta) {
  eta_full    <- c(eta, 0)          # append 0 for reference codon
  lp          <- c(-mu - eta * phi, 0)
  lp          <- lp - max(lp)
  p           <- exp(lp) / sum(exp(lp))
  n_total     <- sum(counts)
  -sum(counts * eta_full) + n_total * sum(eta_full * p)
}

for (k in k_values) {
  scores_orig   <- mapply(roc_score,
                          as.data.frame(sim_counts), phi_true,
                          MoreArgs = list(mu = mu_true, eta = eta_true))
  scores_scaled <- mapply(roc_score,
                          as.data.frame(sim_counts), phi_true / k,
                          MoreArgs = list(mu = mu_true, eta = eta_true * k))

  # score(phi', k*eta) must equal k * score(phi, eta)
  assert_near(
    sprintf("score(phi',k*eta) = k*score(phi,eta)  (k = %8.2f)", k),
    scores_scaled, k * scores_orig
  )
}


# ── Summary ──────────────────────────────────────────────────────────────────
cat(sprintf("\n=================================================================\n"))
cat(sprintf("Results: %d passed, %d failed\n", n_pass, n_fail))
cat(sprintf("=================================================================\n"))

if (n_fail > 0L) {
  quit(status = 1L)
} else {
  cat("\nConclusion: phi normalization is mathematically safe.\n")
  cat("ROC_eff = phi_estimated * eta is invariant to multiplicative\n")
  cat("scaling of phi.  The ROC_eff > 1 threshold retains its\n")
  cat("population-genetic meaning regardless of the CPM scale.\n")
}
