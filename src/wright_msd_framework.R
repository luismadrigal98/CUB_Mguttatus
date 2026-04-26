##' @title  Wright's mutation-selection-drift framework for CUS
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' @date   2026-04-25
##'
##' Theoretical machinery for predicting equilibrium preferred-codon
##' frequency Q and per-site nucleotide diversity pi as functions of the
##' scaled selection coefficient S = 2 N s (haploid convention) under a
##' two-allele Wright model with biased mutation.
##'
##' The stationary density (Wright 1937, "Evolution and the genetics of
##' populations") is
##'
##'     phi(p) ∝ p^(V-1) * (1-p)^(U-1) * exp(S * p)
##'
##' with V = 2*N*v (scaled mutation rate TO the preferred allele) and
##' U = 2*N*u (scaled mutation rate FROM the preferred allele).  The
##' normalization and moments are available in closed form via Kummer's
##' confluent hypergeometric function M(a, b, z) = 1F1(a; b; z), so we
##' avoid the numerical pitfalls of integrating across the boundary
##' singularities at p = 0 and p = 1.
##'
##' Identities used (gsl::hyperg_1F1):
##'   Z(S, U, V)         = B(V, U)         * M(V,   U+V,   S)
##'   E[p]               = (V/(U+V))       * M(V+1, U+V+1, S) / M(V, U+V, S)
##'   E[p (1-p)]         = (V*U / ((U+V)*(U+V+1))) *
##'                        M(V+1, U+V+2, S) / M(V, U+V, S)
##'
##' At S = 0 these reduce to the analytic neutral results
##'   Q_neutral  = V / (U + V)
##'   pi_neutral = 2 * V * U / ((U + V) * (U + V + 1))
##' which we use to solve for U and V given empirical Q and pi at putatively
##' neutral genes (lowest-S_ROC bin).

## ---------------------------------------------------------------------------
## Core moments of the Wright stationary distribution
## ---------------------------------------------------------------------------

#' Q(S) = E[p] under Wright's stationary distribution
#'
#' @param S Numeric vector of scaled selection coefficients (2 N s).
#' @param U Scaled mutation rate FROM preferred (2 N u). Must be > 0.
#' @param V Scaled mutation rate TO preferred (2 N v). Must be > 0.
#' @return Numeric vector of equilibrium preferred-codon frequencies.
wright_Q <- function(S, U, V) {
  stopifnot(U > 0, V > 0, all(is.finite(S)))
  num <- vapply(S, function(s) gsl::hyperg_1F1(V + 1, U + V + 1, s), numeric(1))
  den <- vapply(S, function(s) gsl::hyperg_1F1(V,     U + V,     s), numeric(1))
  (V / (U + V)) * num / den
}

#' pi(S) = 2 * E[p (1 - p)] under Wright's stationary distribution
#'
#' This is the per-site heterozygosity in a two-allele model.
wright_pi <- function(S, U, V) {
  stopifnot(U > 0, V > 0, all(is.finite(S)))
  num <- vapply(S, function(s) gsl::hyperg_1F1(V + 1, U + V + 2, s), numeric(1))
  den <- vapply(S, function(s) gsl::hyperg_1F1(V,     U + V,     s), numeric(1))
  prefactor <- 2 * V * U / ((U + V) * (U + V + 1))
  prefactor * num / den
}

## ---------------------------------------------------------------------------
## Closed-form neutral solver: recover (U, V) from observed Q and pi at S = 0
## ---------------------------------------------------------------------------

#' Recover scaled mutation parameters from observed neutral Q and pi
#'
#' Given Q_neutral = V/(U+V) and pi_neutral = 2*V*U/((U+V)*(U+V+1)),
#' the system has a unique solution in U, V > 0 whenever
#' 0 < Q_neutral < 1 and pi_neutral < 2 * Q_neutral * (1 - Q_neutral).
#'
#' @param Q_neutral Observed preferred-codon frequency in the (near-)neutral pool.
#' @param pi_neutral Observed nucleotide diversity in the same pool.
#' @return Named numeric vector with elements `U`, `V`, `W` (= U+V, total scaled mutation).
wright_solve_UV <- function(Q_neutral, pi_neutral) {
  stopifnot(Q_neutral > 0, Q_neutral < 1, pi_neutral > 0)
  hardy_max <- 2 * Q_neutral * (1 - Q_neutral)
  if (pi_neutral >= hardy_max) {
    stop(sprintf(
      "pi_neutral (%.5f) must be < 2*Q*(1-Q) = %.5f for a finite-mutation solution.",
      pi_neutral, hardy_max
    ))
  }
  r <- pi_neutral / hardy_max          # in (0, 1)
  W <- r / (1 - r)                     # W = U + V
  V <- Q_neutral * W
  U <- (1 - Q_neutral) * W
  c(U = U, V = V, W = W)
}

## ---------------------------------------------------------------------------
## Calibrate S_ROC -> S_Wright via a single proportional scale factor alpha
## ---------------------------------------------------------------------------

#' Fit alpha such that Q(alpha * S_ROC; U, V) approximates observed Q
#'
#' Uses bin-aggregated empirical Q (more stable than per-gene fits, since
#' per-gene Q has high sampling noise from small 4-fold-site counts).
#'
#' @param S_ROC_bin   Numeric vector of mean S_ROC per bin.
#' @param Q_obs_bin   Numeric vector of mean preferred-codon frequency per bin.
#' @param weights     Optional bin weights (e.g. site counts) for weighted LS.
#' @param U,V         Wright mutation parameters.
#' @param alpha_init  Starting value for alpha.
#' @return List with fitted `alpha`, residual sd, and the nls fit.
wright_calibrate_alpha <- function(S_ROC_bin, Q_obs_bin,
                                   weights = NULL, U, V, alpha_init = 1) {
  stopifnot(length(S_ROC_bin) == length(Q_obs_bin))
  if (is.null(weights)) weights <- rep(1, length(S_ROC_bin))
  df <- data.frame(s = S_ROC_bin, q = Q_obs_bin, w = weights)
  fit <- tryCatch(
    nls(q ~ wright_Q(alpha * s, U = U, V = V),
        data = df, start = list(alpha = alpha_init), weights = w,
        control = nls.control(warnOnly = TRUE, maxiter = 200)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(alpha = NA_real_, residual_sd = NA_real_, fit = NULL))
  }
  alpha_hat <- unname(coef(fit)["alpha"])
  resid_sd  <- sqrt(mean(residuals(fit)^2))
  list(alpha = alpha_hat, residual_sd = resid_sd, fit = fit)
}

## ---------------------------------------------------------------------------
## Threshold derivation on the Wright Q(S) curve
## ---------------------------------------------------------------------------

#' Solve S such that Q(S; U, V) = Q_target
#'
#' Used to map a chosen Q-criterion (e.g. midpoint, half-fixation) to the
#' corresponding S in the Wright scale.  Once alpha is estimated, the
#' equivalent S_ROC threshold is S_target / alpha.
wright_invert_Q <- function(Q_target, U, V, S_lower = 0, S_upper = 50,
                            S_upper_max = 400) {
  q_low  <- wright_Q(S_lower, U, V)
  if (Q_target <= q_low) return(NA_real_)
  # If Q_target sits above the current upper bound, expand S_upper adaptively
  # so high-Q outlier genes can still be inverted (slowly-saturating regime).
  q_high <- wright_Q(S_upper, U, V)
  while (Q_target >= q_high && S_upper < S_upper_max) {
    S_upper <- min(S_upper * 2, S_upper_max)
    q_high  <- wright_Q(S_upper, U, V)
  }
  if (Q_target >= q_high) return(NA_real_)
  uniroot(function(s) wright_Q(s, U, V) - Q_target,
          interval = c(S_lower, S_upper), tol = 1e-6)$root
}

#' Three candidate Wright-S thresholds with biological rationale
#'
#' @return data.frame with columns name, Q_target, S_Wright, rationale.
wright_threshold_candidates <- function(U, V) {
  Q_neutral <- V / (U + V)
  Q_mid     <- (Q_neutral + 1) / 2          # halfway from neutral to fixation
  Q_half    <- 0.5                           # majority-preferred crossover
  Q_2x      <- min(0.95, 2 * Q_neutral)      # twice neutral expectation, capped
  cands <- data.frame(
    name      = c("Q_midpoint", "Q_majority", "Q_2x_neutral"),
    Q_target  = c(Q_mid, Q_half, Q_2x),
    rationale = c(
      "Halfway from neutral to fixation; balances drift and selection regimes",
      "Q crosses 0.5 (preferred codon now wins more often than not)",
      "Q reaches twice the neutral expectation; departure clearly above mutation"
    ),
    stringsAsFactors = FALSE
  )
  cands$S_Wright <- vapply(cands$Q_target,
                           wright_invert_Q, numeric(1), U = U, V = V)
  cands
}
