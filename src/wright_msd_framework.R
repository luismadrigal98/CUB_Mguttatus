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
##' neutral genes (lowest-L_ROC bin).

## ***************************************************************************
## Core moments of the Wright stationary distribution ----
## ___________________________________________________________________________

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

## ***************************************************************************
## Closed-form neutral solver: recover (U, V) from observed Q and pi at S = 0
## ___________________________________________________________________________

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

## ***************************************************************************
## Calibrate L_ROC -> S_Wright via a single proportional scale factor alpha
## ___________________________________________________________________________

#' Fit alpha such that Q(alpha * L_ROC; U, V) approximates observed Q
#'
#' Uses bin-aggregated empirical Q (more stable than per-gene fits, since
#' per-gene Q has high sampling noise from small 4-fold-site counts).
#'
#' @param S_ROC_bin   Numeric vector of mean L_ROC per bin.
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

## ***************************************************************************
## Threshold derivation on the Wright Q(S) curve
## ___________________________________________________________________________

#' Solve S such that Q(S; U, V) = Q_target
#'
#' Q(S; U, V) is strictly monotone-increasing on the whole real line:
#' Q(-Inf) = 0, Q(0) = V/(U+V) (the neutral mean), Q(+Inf) = 1.  So *any*
#' Q_target in (0, 1) has a unique S in R.  This solver searches a signed
#' bracket and adaptively expands it on either side, so genes with
#' Q_obs < Q_neutral correctly return a *negative* S_Wright (selection
#' against the preferred allele) instead of being floored to zero or NA.
#'
#' If `floor_at_zero = TRUE`, the solver reports S = 0 for any Q_target
#' below the neutral mean Q(0) = V/(U+V).  This is the operational
#' "drift cap" used in main.R §8.3.4: per-gene Q below Q_neutral is
#' interpreted as drift acting (no detectable selection FOR the preferred
#' allele) rather than as an inferred selection coefficient against it.
#' The signed value remains the right thing to compute when probing the
#' empirical (U, V) calibration (e.g. for `bin_eta$S_Wright_bin`).
#'
#' @param Q_target  Target preferred-base frequency, in (0, 1).
#' @param U,V       Wright mutation parameters (must be > 0).
#' @param S_lower,S_upper  Initial bracket; expanded outward as needed.
#' @param S_lower_min,S_upper_max  Hard limits on the expansion.
#' @param floor_at_zero  If TRUE, return 0 (not a negative S) when
#'        Q_target falls below the neutral mean V/(U+V).  Default FALSE.
#' @return Numeric S, or NA_real_ if Q_target falls outside (0, 1) or
#'         the expansion cannot bracket the root within the hard limits.
wright_invert_Q <- function(Q_target, U, V,
                            S_lower = -50, S_upper = 50,
                            S_lower_min = -400, S_upper_max = 400,
                            floor_at_zero = FALSE) {
  if (!is.finite(Q_target) || Q_target <= 0 || Q_target >= 1) return(NA_real_)
  if (floor_at_zero && Q_target <= V / (U + V)) return(0)
  # Expand the lower bracket outward until Q(S_lower) < Q_target
  q_low <- wright_Q(S_lower, U, V)
  while (Q_target <= q_low && S_lower > S_lower_min) {
    S_lower <- max(S_lower * 2, S_lower_min)
    q_low   <- wright_Q(S_lower, U, V)
  }
  if (Q_target <= q_low) return(NA_real_)
  # Expand the upper bracket outward until Q(S_upper) > Q_target
  q_high <- wright_Q(S_upper, U, V)
  while (Q_target >= q_high && S_upper < S_upper_max) {
    S_upper <- min(S_upper * 2, S_upper_max)
    q_high  <- wright_Q(S_upper, U, V)
  }
  if (Q_target >= q_high) return(NA_real_)
  uniroot(function(s) wright_Q(s, U, V) - Q_target,
          interval = c(S_lower, S_upper), tol = 1e-6)$root
}

#' Solve S such that pi(S; U, V) = pi_target
#'
#' pi(S; U, V) is non-monotone: it rises from pi(0) to a maximum (typically
#' near S ~ ln(U/V)) and then decays as selection drives the preferred
#' allele toward fixation and heterozygosity collapses. This solver finds
#' the *rising-flank* crossing (low S, selection-favoring regime), which
#' is appropriate for per-gene π values observed in the data.
#'
#' The solver uses a log-spaced grid to locate the first sign change on the
#' rising flank, then refines with uniroot. This is analogous to
#' derive_S_barrier() but inverts individual per-gene π values rather than
#' computing a population-level threshold.
#'
#' If `floor_at_zero = TRUE`, the solver reports S = 0 for any pi_target
#' below the neutral π(0) = 2*V*U/((U+V)*(U+V+1)), matching the operational
#' drift-cap logic of wright_invert_Q(). Otherwise, returns NA if pi_target
#' falls outside the feasible range or no rising-flank crossing is found.
#'
#' @param pi_target  Target nucleotide diversity, pi_target > 0.
#' @param U,V        Wright mutation parameters (must be > 0).
#' @param S_grid_max Upper end of the log-spaced search grid (default 50).
#' @param S_grid_n   Resolution of the search grid (default 1000).
#' @param floor_at_zero  If TRUE, return 0 (not NA) when pi_target falls
#'        below the neutral π(0). Default FALSE.
#' @param tol        uniroot tolerance.
#' @return Numeric S at which pi(S; U, V) first crosses pi_target on the
#'         rising flank, or NA_real_ if pi_target is outside the feasible
#'         range or no rising-flank crossing exists.
wright_invert_pi <- function(pi_target, U, V,
                             S_grid_max = 50, S_grid_n = 1000,
                             floor_at_zero = FALSE,
                             tol = 1e-6) {
  if (!is.finite(pi_target) || pi_target <= 0) return(NA_real_)
  pi_neutral <- wright_pi(0, U, V)
  if (floor_at_zero && pi_target <= pi_neutral) return(0)
  # Log-spaced grid to find rising-flank crossing
  s_grid    <- exp(seq(log(max(tol, 1e-5)), log(S_grid_max),
                       length.out = S_grid_n))
  pi_grid   <- wright_pi(s_grid, U, V)
  diff_grid <- pi_grid - pi_target
  # First grid index where pi_grid > pi_target -> rising-flank crossing
  first_pos <- which(diff_grid > 0)[1]
  if (is.na(first_pos)) return(NA_real_)
  s_lo <- if (first_pos == 1) 0 else s_grid[first_pos - 1]
  s_hi <- s_grid[first_pos]
  tryCatch(
    uniroot(function(s) wright_pi(s, U, V) - pi_target,
            interval = c(s_lo, s_hi), tol = tol)$root,
    error = function(e) NA_real_
  )
}

## ***************************************************************************
## Dynamic drift-barrier threshold on the Wright pi(S) curve
## ___________________________________________________________________________

#' Solve S such that pi(S; U, V) = (1 + fraction) * pi_neutral
#' on the rising flank of Wright's pi(S).
#'
#' pi(S) is non-monotone: it rises from pi(0) to a hump apex (typically
#' near S ~ ln(U/V)) and then collapses as selection drives the preferred
#' allele toward fixation.  We want the *first* (low-S, rising-flank)
#' crossing of (1 + fraction) * pi_neutral, NOT the descending-flank
#' crossing on the far side of the hump.  We locate the first sign change
#' on a log-spaced grid, then refine with uniroot inside that bracket.
#'
#' This replaces the hand-set S_BARRIER = 0.1 from JK's Mathematica
#' pi-rise simulation with a value derived directly from the empirical
#' (U, V).  For the M. guttatus calibration (U_emp = 0.0698, V_emp =
#' 0.0245), a 1% rise corresponds to S = 0.044, while JK's prior 0.10
#' value corresponds to roughly a 2.3% rise on the same curve.
#'
#' @param U,V         Wright mutation parameters (>0).
#' @param pi_neutral  Reference pi at S = 0 (typically wright_pi(0, U, V)).
#' @param fraction    Relative rise over pi_neutral defining the barrier.
#'                    Default 0.01 (1% rise).
#' @param S_grid_max  Upper end of the log-spaced search grid (default 50).
#' @param S_grid_n    Resolution of the search grid (default 1000).
#' @param tol         uniroot tolerance.
#' @return Numeric S at which pi crosses (1 + fraction) * pi_neutral on
#'         its rising flank, or NA_real_ if no rising-flank crossing
#'         exists within the grid.
derive_S_barrier <- function(U, V, pi_neutral,
                             fraction = 0.01,
                             S_grid_max = 50, S_grid_n = 1000,
                             tol = 1e-6) {
  stopifnot(U > 0, V > 0, pi_neutral > 0, fraction > 0)
  target <- (1 + fraction) * pi_neutral
  f <- function(s) wright_pi(s, U, V) - target
  if (f(0) >= 0) {
    stop(sprintf(
      "pi(0) = %.5g already at/above target (%.5g); check pi_neutral input.",
      wright_pi(0, U, V), target
    ))
  }
  # Log-spaced grid (rising flank is concentrated at small S; log-spacing
  # gives ~3-decade dynamic range without wasting points at large S).
  s_grid    <- exp(seq(log(max(tol, 1e-5)), log(S_grid_max),
                       length.out = S_grid_n))
  pi_grid   <- wright_pi(s_grid, U, V)
  diff_grid <- pi_grid - target
  # First grid index where pi_grid > target -> rising-flank crossing
  # lies between (s_grid[first_pos - 1] or 0, s_grid[first_pos]).
  first_pos <- which(diff_grid > 0)[1]
  if (is.na(first_pos)) return(NA_real_)
  s_lo <- if (first_pos == 1) 0 else s_grid[first_pos - 1]
  s_hi <- s_grid[first_pos]
  uniroot(f, interval = c(s_lo, s_hi), tol = tol)$root
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
