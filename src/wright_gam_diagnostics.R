##' @title  Diagnostic battery for GAM-RE-shrunk Wright S inversion
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' @date   2026-04-28
##'
##' Tier 1-3 diagnostics for deciding whether the per-gene GAM-RE-shrunk
##' inversion can replace the raw per-gene wright_invert_Q (and
##' wright_invert_pi) pipeline.
##'
##' Tier 1 (adoption-blocking):
##'   1.1 Simulation recovery: known-S simulation drawn from a parametric
##'       prior anchored to the bin-level S_Wright (low-noise reference);
##'       compare RMSE of raw vs shrunk recovery stratified by N_4fold
##'       quartile.
##'   1.2 Bin-level consistency: bin-mean(per-gene shrunk S) vs the
##'       bin-pooled S_Wright_bin (gold standard within data).
##'   1.3 Threshold robustness: bootstrap thr_eta both ways, compare CIs.
##'
##' Tier 2 (model fit):
##'   2.4 RE BLUP QQ-plot and REML variance.
##'   2.5 Shrinkage pattern: predicted vs observed Q (or π), colored by N.
##'   2.6 Dispersion check (binomial scale ≈ 1; gaussian σ²).
##'
##' Tier 3 (signal preservation):
##'   3.7 Spearman(S_eta, S_Wright) raw vs shrunk vs π-branch.
##'   3.8 Rank-displacement audit: which genes shift most, and is the
##'       shift concentrated at low N (expected) vs spread across N
##'       (red flag)?

# --------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------
# --------------------------------------------------------------------------

.invert_Q_safe <- function(q_vec, U, V) {
  vapply(q_vec, function(q) {
    if (!is.finite(q)) return(NA_real_)
    tryCatch(wright_invert_Q(q, U = U, V = V),
             error = function(e) NA_real_)
  }, numeric(1))
}

.invert_pi_safe <- function(pi_vec, U, V) {
  vapply(pi_vec, function(p) {
    if (!is.finite(p)) return(NA_real_)
    tryCatch(wright_invert_pi(p, U = U, V = V),
             error = function(e) NA_real_)
  }, numeric(1))
}

# --------------------------------------------------------------------------
# Tier 1.1 — Simulation recovery -------------------------------------------
# --------------------------------------------------------------------------
#
# Truth model: S = f(S_eta) + Normal(0, sigma), where f is a smooth on
# bin-pooled S_Wright_bin (low-noise reference, NOT the GAM-RE fit being
# tested). This gives a non-circular ground truth that respects the broad
# covariate-S relationship without conflating the simulation with the model
# under test.

.build_truth_model <- function(bin_eta, msd_data) {
  ok <- is.finite(bin_eta$mean_S_eta) & is.finite(bin_eta$S_Wright_bin)
  if (sum(ok) < 5) {
    stop("Not enough finite (mean_S_eta, S_Wright_bin) bins to build truth model.")
  }
  f_true <- approxfun(bin_eta$mean_S_eta[ok], bin_eta$S_Wright_bin[ok],
                      rule = 2)
  S_resid <- msd_data$S_Wright_signed - f_true(msd_data$S_eta)
  sigma_S <- stats::sd(S_resid, na.rm = TRUE)
  list(f = f_true, sigma = sigma_S)
}

.simulate_recovery_Q <- function(gam_Q_pool, truth, U_emp, V_emp, seed) {
  set.seed(seed)
  d <- gam_Q_pool
  d$S_true <- truth$f(d$S_eta) + stats::rnorm(nrow(d), 0, truth$sigma)
  d$Q_true <- wright_Q(d$S_true, U = U_emp, V = V_emp)
  d$N_pref_sim <- stats::rbinom(nrow(d), d$N_4fold_sites, d$Q_true)
  d$N_other_sim <- d$N_4fold_sites - d$N_pref_sim
  d$Q_obs_sim <- d$N_pref_sim / d$N_4fold_sites
  d$S_raw_recov <- .invert_Q_safe(d$Q_obs_sim, U_emp, V_emp)

  gam_sim <- tryCatch(
    mgcv::bam(
      cbind(N_pref_sim, N_other_sim) ~
        s(Max_Log10_Exp, k = 8) + s(Exp_breadth, k = 8) +
        s(Log_CDS_length_nt, k = 8) + s(Gene_name_f, bs = "re"),
      data = d, family = stats::binomial(link = "logit"),
      method = "fREML", discrete = TRUE
    ),
    error = function(e) { message("[sim Q] bam failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(gam_sim)) {
    d$Q_GAM_sim     <- NA_real_
    d$S_shrunk_recov <- NA_real_
  } else {
    d$Q_GAM_sim     <- as.numeric(predict(gam_sim, newdata = d, type = "response"))
    d$S_shrunk_recov <- .invert_Q_safe(d$Q_GAM_sim, U_emp, V_emp)
  }
  d
}

.simulate_recovery_pi <- function(gam_pi_pool, truth, U_emp, V_emp, seed) {
  set.seed(seed)
  d <- gam_pi_pool
  # π is non-monotone in S; restrict simulated truth to S >= 0 (rising
  # flank) so wright_invert_pi has a well-defined inverse.
  d$S_true <- pmax(truth$f(d$S_eta) +
                     stats::rnorm(nrow(d), 0, truth$sigma), 0)
  d$pi_true <- wright_pi(d$S_true, U = U_emp, V = V_emp)
  # Sample π_obs from a Beta surrogate matched to π_true with variance
  # π(1 - π/2) / N (binomial heterozygosity SE), avoiding a per-site loop.
  d$pi_obs_sim <- mapply(function(pi_t, N) {
    if (!is.finite(pi_t) || pi_t <= 0 || N < 1) return(NA_real_)
    var_pi <- pi_t * (1 - pi_t / 2) / N
    var_pi <- max(var_pi, 1e-12)
    p_norm   <- pi_t / 0.5
    var_norm <- var_pi / 0.25
    if (p_norm * (1 - p_norm) <= var_norm) return(pi_t)
    nu <- p_norm * (1 - p_norm) / var_norm - 1
    a  <- p_norm * nu; b <- (1 - p_norm) * nu
    0.5 * stats::rbeta(1, a, b)
  }, d$pi_true, d$N_4fold_sites)
  d$S_raw_recov <- .invert_pi_safe(d$pi_obs_sim, U_emp, V_emp)

  gam_sim <- tryCatch(
    mgcv::bam(
      pi_obs_sim ~
        s(Max_Log10_Exp, k = 8) + s(Exp_breadth, k = 8) +
        s(Log_CDS_length_nt, k = 8) + s(Gene_name_f, bs = "re"),
      data = d, family = stats::gaussian(link = "identity"),
      weights = N_4fold_sites,
      method = "fREML", discrete = TRUE
    ),
    error = function(e) { message("[sim π] bam failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(gam_sim)) {
    d$pi_GAM_sim    <- NA_real_
    d$S_shrunk_recov <- NA_real_
  } else {
    d$pi_GAM_sim    <- as.numeric(predict(gam_sim, newdata = d, type = "response"))
    d$S_shrunk_recov <- .invert_pi_safe(d$pi_GAM_sim, U_emp, V_emp)
  }
  d
}

.summarize_recovery <- function(sim_data, branch_label) {
  n_q <- stats::quantile(sim_data$N_4fold_sites,
                         probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sim_data$N_strata <- cut(sim_data$N_4fold_sites,
                           breaks  = c(-Inf, n_q, Inf),
                           labels  = c("Q1 (smallest N)", "Q2", "Q3",
                                       "Q4 (largest N)"),
                           include.lowest = TRUE)

  long <- dplyr::bind_rows(
    data.frame(S_true = sim_data$S_true, S_recov = sim_data$S_raw_recov,
               estimator = "raw", N_strata = sim_data$N_strata),
    data.frame(S_true = sim_data$S_true, S_recov = sim_data$S_shrunk_recov,
               estimator = "shrunk (GAM-RE)", N_strata = sim_data$N_strata)
  ) |> dplyr::filter(is.finite(S_true), is.finite(S_recov))

  rmse_table <- long |>
    dplyr::group_by(N_strata, estimator) |>
    dplyr::summarize(
      rmse     = sqrt(mean((S_recov - S_true)^2, na.rm = TRUE)),
      bias     = mean(S_recov - S_true, na.rm = TRUE),
      spearman = stats::cor(S_recov, S_true, method = "spearman",
                            use = "complete.obs"),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = S_true, y = S_recov,
                                          color = estimator)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(size = 0.4, alpha = 0.25) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.6) +
    ggplot2::facet_wrap(~ N_strata, ncol = 2) +
    ggplot2::scale_color_manual(values = c("raw" = "#d95f02",
                                           "shrunk (GAM-RE)" = "#1b9e77")) +
    ggplot2::labs(x = "True S (simulated)", y = "Recovered S",
                  title = sprintf("Tier 1.1 — Simulation recovery: %s branch",
                                  branch_label),
                  color = NULL) +
    theme_custom() +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, rmse_table = rmse_table)
}

# --------------------------------------------------------------------------
# Tier 1.2 — Bin-level consistency -----------------------------------------
# --------------------------------------------------------------------------

.plot_bin_consistency <- function(msd_data, bin_eta, shrunk_col,
                                  branch_label) {
  d <- msd_data |>
    dplyr::filter(!is.na(S_eta), is.finite(.data[[shrunk_col]])) |>
    dplyr::arrange(S_eta) |>
    dplyr::mutate(Seta_bin = dplyr::ntile(S_eta, 30)) |>
    dplyr::group_by(Seta_bin) |>
    dplyr::summarize(mean_S_eta  = mean(S_eta),
                     mean_shrunk = mean(.data[[shrunk_col]], na.rm = TRUE),
                     n           = dplyr::n(),
                     .groups = "drop")

  cmp <- d |>
    dplyr::inner_join(bin_eta |>
                        dplyr::select(Seta_bin, S_Wright_bin),
                      by = "Seta_bin") |>
    dplyr::filter(is.finite(mean_shrunk), is.finite(S_Wright_bin))

  if (nrow(cmp) < 5) {
    return(list(plot = NULL, slope = NA_real_, r2 = NA_real_,
                pearson = NA_real_))
  }

  fit       <- stats::lm(mean_shrunk ~ S_Wright_bin, data = cmp)
  slope     <- unname(stats::coef(fit)[2])
  r2        <- summary(fit)$r.squared
  pearson_r <- stats::cor(cmp$mean_shrunk, cmp$S_Wright_bin)

  p <- ggplot2::ggplot(cmp, ggplot2::aes(x = S_Wright_bin, y = mean_shrunk)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(size = 2, color = "#1b9e77") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#1b9e77",
                         fill = "#1b9e7733", linewidth = 0.6) +
    ggplot2::labs(
      x = "S_Wright_bin (bin-pooled, low-noise reference)",
      y = "mean(per-gene shrunk S) within bin",
      title = sprintf(
        "Tier 1.2 — Bin-level consistency: %s\nslope = %.3f, R² = %.3f, r = %.3f",
        branch_label, slope, r2, pearson_r)
    ) +
    theme_custom()

  list(plot = p, slope = slope, r2 = r2, pearson = pearson_r)
}

# --------------------------------------------------------------------------
# Tier 1.3 — Threshold robustness (Q-pipeline only) ------------------------
# --------------------------------------------------------------------------

.bootstrap_thr_eta <- function(data, sw_col, S_BARRIER, n_boot = 200) {
  d <- data |>
    dplyr::filter(is.finite(.data[[sw_col]]), is.finite(S_eta))
  if (nrow(d) < 100) return(rep(NA_real_, n_boot))
  vapply(seq_len(n_boot), function(b) {
    idx <- sample(seq_len(nrow(d)), replace = TRUE)
    bd  <- d[idx, ]
    fit <- tryCatch(
      mgcv::bam(stats::as.formula(paste0(sw_col, " ~ s(S_eta, k = 8)")),
                data = bd, method = "fREML", discrete = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    seta_grid <- seq(min(bd$S_eta, na.rm = TRUE),
                     max(bd$S_eta, na.rm = TRUE), length.out = 401)
    sw_pred <- as.numeric(stats::predict(fit,
                                         newdata = data.frame(S_eta = seta_grid)))
    crossings <- which(sw_pred >= S_BARRIER)
    if (length(crossings) == 0) return(NA_real_)
    i <- crossings[1]
    if (i == 1) return(seta_grid[i])
    sw_lo <- sw_pred[i - 1]; sw_hi <- sw_pred[i]
    e_lo  <- seta_grid[i - 1]; e_hi <- seta_grid[i]
    e_lo + (e_hi - e_lo) * (S_BARRIER - sw_lo) / (sw_hi - sw_lo)
  }, numeric(1))
}

.plot_threshold_robustness <- function(thr_raw, thr_shr, S_BARRIER) {
  long <- dplyr::bind_rows(
    data.frame(thr = thr_raw, source = "raw (per-gene Q)"),
    data.frame(thr = thr_shr, source = "shrunk (GAM-RE)")
  ) |> dplyr::filter(is.finite(thr))

  if (nrow(long) < 20) {
    return(list(plot = NULL, ci_table = NULL))
  }

  ci_table <- long |>
    dplyr::group_by(source) |>
    dplyr::summarize(median = stats::median(thr, na.rm = TRUE),
                     lo     = stats::quantile(thr, 0.025, na.rm = TRUE),
                     hi     = stats::quantile(thr, 0.975, na.rm = TRUE),
                     .groups = "drop")

  p <- ggplot2::ggplot(long, ggplot2::aes(x = thr, fill = source)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(data = ci_table,
                        ggplot2::aes(xintercept = median, color = source),
                        linewidth = 0.6, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("raw (per-gene Q)" = "#d95f02",
                                          "shrunk (GAM-RE)" = "#1b9e77")) +
    ggplot2::scale_color_manual(values = c("raw (per-gene Q)" = "#d95f02",
                                           "shrunk (GAM-RE)" = "#1b9e77")) +
    ggplot2::labs(x = "thr_eta (bootstrap)", y = "density",
                  title = sprintf(
                    "Tier 1.3 — Threshold robustness (S_BARRIER = %.4f, B = %d)",
                    S_BARRIER, length(thr_raw)),
                  fill = NULL, color = NULL) +
    theme_custom() +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, ci_table = ci_table)
}

# --------------------------------------------------------------------------
# Tier 1.4 — Bulmer/Li closed-form approximation vs exact Wright -----------
# --------------------------------------------------------------------------
#
# When U, V << 1 (M. guttatus: U_emp ≈ 0.07, V_emp ≈ 0.025), Wright's
# stationary distribution concentrates at p = 0 and p = 1, and the
# equilibrium preferred-base frequency reduces to the Bulmer/Li form:
#
#     Q ≈ (V/U) * exp(S) / (1 + (V/U) * exp(S))
#
# which inverts in closed form:
#
#     S ≈ logit(Q) - log(V/U)
#
# If S_Bulmer agrees with the exact Wright inversion across the bulk of the
# data, we can adopt the closed-form expression as the canonical operational
# mapping, citing Bulmer (1991), with this panel as the validation against
# the exact 1F1 machinery.

.compute_S_bulmer <- function(Q, U, V) {
  Q_clip <- pmin(pmax(Q, 1e-9), 1 - 1e-9)
  log(Q_clip / (1 - Q_clip)) - log(V / U)
}

.plot_bulmer_check <- function(msd_data, U_emp, V_emp) {
  d <- msd_data |>
    dplyr::filter(is.finite(Q_GAM), is.finite(S_Wright_GAM_signed))
  if (nrow(d) < 30) {
    return(list(plot = NULL, cor_bulmer_gam = NA_real_,
                slope = NA_real_, intercept = NA_real_,
                max_abs_diff = NA_real_, mean_abs_diff = NA_real_))
  }
  d$S_Bulmer <- .compute_S_bulmer(d$Q_GAM, U_emp, V_emp)
  d$abs_diff <- abs(d$S_Bulmer - d$S_Wright_GAM_signed)

  cor_bg    <- stats::cor(d$S_Bulmer, d$S_Wright_GAM_signed,
                          method = "pearson", use = "complete.obs")
  fit       <- stats::lm(S_Bulmer ~ S_Wright_GAM_signed, data = d)
  slope     <- unname(stats::coef(fit)[2])
  intercept <- unname(stats::coef(fit)[1])

  p <- ggplot2::ggplot(d, ggplot2::aes(x = S_Wright_GAM_signed,
                                       y = S_Bulmer,
                                       color = abs_diff)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::scale_color_viridis_c(name = "|Δ|", option = "plasma") +
    ggplot2::labs(
      x = "S_Wright_GAM (exact 1F1 inversion)",
      y = "S_Bulmer = logit(Q_GAM) − log(V/U)",
      title = sprintf(
        "Tier 1.4 — Bulmer/Li closed form vs exact Wright\nr = %.4f, slope = %.3f, intercept = %.3f, max |Δ| = %.3f, mean |Δ| = %.3f",
        cor_bg, slope, intercept,
        max(d$abs_diff, na.rm = TRUE), mean(d$abs_diff, na.rm = TRUE))
    ) +
    theme_custom()

  list(plot = p,
       cor_bulmer_gam = cor_bg,
       slope = slope, intercept = intercept,
       max_abs_diff  = max(d$abs_diff, na.rm = TRUE),
       mean_abs_diff = mean(d$abs_diff, na.rm = TRUE))
}

# --------------------------------------------------------------------------
# Tier 2 — Model-fit sanity ------------------------------------------------
# --------------------------------------------------------------------------

.plot_re_blup_qq <- function(gam_fit, branch_label) {
  cf <- stats::coef(gam_fit)
  blup_idx <- grep("Gene_name_f", names(cf), fixed = TRUE)
  if (length(blup_idx) == 0) {
    return(list(plot = NULL, re_sd = NA_real_, blups = numeric(0)))
  }
  blups <- unname(cf[blup_idx])
  vc <- tryCatch(mgcv::gam.vcomp(gam_fit, rescale = FALSE, conf.lev = 0.95),
                 error = function(e) NULL)
  re_sd <- if (is.null(vc)) NA_real_ else {
    re_row <- grep("Gene_name_f", rownames(vc), fixed = TRUE)
    if (length(re_row) == 0) NA_real_ else as.numeric(vc[re_row, "std.dev"])
  }

  qq_df <- data.frame(theoretical = stats::qnorm(stats::ppoints(length(blups))),
                      sample = sort(blups))
  p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_abline(slope = stats::sd(blups),
                         intercept = mean(blups),
                         color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(size = 0.5, alpha = 0.4, color = "#1b9e77") +
    ggplot2::labs(
      x = "Theoretical quantiles (Normal)",
      y = "Sample BLUP",
      title = sprintf("Tier 2.4 — RE BLUP QQ: %s\nRE sd (REML) = %s",
                      branch_label,
                      if (is.na(re_sd)) "NA" else sprintf("%.3f", re_sd))
    ) +
    theme_custom()
  list(plot = p, re_sd = re_sd, blups = blups)
}

.plot_shrinkage_pattern <- function(pool, q_obs_col, q_gam_col,
                                    branch_label) {
  d <- pool |>
    dplyr::filter(is.finite(.data[[q_obs_col]]),
                  is.finite(.data[[q_gam_col]]))
  d$N_log10 <- log10(d$N_4fold_sites)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[q_obs_col]],
                                       y = .data[[q_gam_col]],
                                       color = N_log10)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::scale_color_viridis_c(name = "log10(N_4fold)") +
    ggplot2::labs(
      x = sprintf("%s (observed)", q_obs_col),
      y = sprintf("%s (shrunk GAM-RE)", q_gam_col),
      title = sprintf("Tier 2.5 — Shrinkage pattern: %s", branch_label)
    ) +
    theme_custom()
  list(plot = p)
}

# --------------------------------------------------------------------------
# Tier 3 — Signal preservation ---------------------------------------------
# --------------------------------------------------------------------------

.plot_rank_displacement <- function(msd_data, raw_col, shrunk_col,
                                    branch_label) {
  d <- msd_data |>
    dplyr::filter(is.finite(.data[[raw_col]]),
                  is.finite(.data[[shrunk_col]]),
                  !is.na(N_4fold_sites))
  d$rank_raw    <- rank(d[[raw_col]])
  d$rank_shrunk <- rank(d[[shrunk_col]])
  d$delta       <- abs(d$rank_raw - d$rank_shrunk)
  d$mover       <- d$delta >= 1000

  p <- ggplot2::ggplot(d, ggplot2::aes(x = N_4fold_sites, fill = mover)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey60", "TRUE" = "#d95f02"),
      labels = c("FALSE" = "|Δrank| < 1000",
                 "TRUE"  = "|Δrank| ≥ 1000")
    ) +
    ggplot2::labs(
      x = "N_4fold_sites (log10)", y = "density",
      title = sprintf("Tier 3.8 — Rank displacement vs N: %s",
                      branch_label),
      fill = NULL
    ) +
    theme_custom() +
    ggplot2::theme(legend.position = "bottom")

  pct_movers_low_half <- if (any(d$mover, na.rm = TRUE)) {
    median_N <- stats::quantile(d$N_4fold_sites, 0.5, na.rm = TRUE)
    mean(d$N_4fold_sites[d$mover] < median_N, na.rm = TRUE)
  } else NA_real_

  list(plot = p,
       n_movers = sum(d$mover, na.rm = TRUE),
       pct_movers_low_half = pct_movers_low_half)
}

# --------------------------------------------------------------------------
# Main entry point ---------------------------------------------------------
# --------------------------------------------------------------------------

run_wright_gam_diagnostics <- function(
    msd_data,
    gam_Q_pool,  gam_Q_wright,
    gam_pi_pool, gam_pi_wright,
    bin_eta,
    U_emp, V_emp, S_BARRIER,
    output_pdf  = "results/Wright_GAM_diagnostics.pdf",
    n_sim_reps  = 3,
    n_boot      = 200,
    seed        = 42L) {

  set.seed(seed)
  truth <- .build_truth_model(bin_eta, msd_data)
  cat(sprintf("[GAM diagnostics] Truth model: residual sigma = %.3f\n",
              truth$sigma))

  cat("[GAM diagnostics] Tier 1.1 — simulation recovery (Q)...\n")
  sim_runs_Q <- lapply(seq_len(n_sim_reps), function(r) {
    .simulate_recovery_Q(gam_Q_pool, truth, U_emp, V_emp,
                         seed = seed + r)
  })
  sim_Q  <- dplyr::bind_rows(sim_runs_Q, .id = "rep")
  sim_Qx <- .summarize_recovery(sim_Q, "Q")

  cat("[GAM diagnostics] Tier 1.1 — simulation recovery (π)...\n")
  sim_runs_pi <- lapply(seq_len(n_sim_reps), function(r) {
    .simulate_recovery_pi(gam_pi_pool, truth, U_emp, V_emp,
                          seed = seed + r + 1000L)
  })
  sim_pi  <- dplyr::bind_rows(sim_runs_pi, .id = "rep")
  sim_pix <- .summarize_recovery(sim_pi, "π")

  cat("[GAM diagnostics] Tier 1.2 — bin-level consistency...\n")
  bin_Q  <- .plot_bin_consistency(msd_data, bin_eta,
                                  "S_Wright_GAM_signed", "Q")
  bin_pi <- .plot_bin_consistency(msd_data, bin_eta,
                                  "S_Wright_pi_signed", "π")

  cat(sprintf("[GAM diagnostics] Tier 1.3 — threshold bootstrap (B = %d)...\n",
              n_boot))
  thr_raw <- .bootstrap_thr_eta(msd_data, "S_Wright_signed",
                                S_BARRIER, n_boot)
  thr_shr <- .bootstrap_thr_eta(msd_data, "S_Wright_GAM_signed",
                                S_BARRIER, n_boot)
  thr_x   <- .plot_threshold_robustness(thr_raw, thr_shr, S_BARRIER)

  cat("[GAM diagnostics] Tier 1.4 — Bulmer closed form vs exact Wright...\n")
  bulmer_x <- .plot_bulmer_check(msd_data, U_emp, V_emp)

  cat("[GAM diagnostics] Tier 2 — model-fit panels...\n")
  blup_Q   <- .plot_re_blup_qq(gam_Q_wright, "Q")
  blup_pi  <- .plot_re_blup_qq(gam_pi_wright, "π")
  shrink_Q <- .plot_shrinkage_pattern(gam_Q_pool,  "Q_pref_base", "Q_GAM", "Q")
  shrink_pi <- .plot_shrinkage_pattern(gam_pi_pool, "pi_2allele",  "pi_GAM", "π")

  disp_Q  <- gam_Q_wright$deviance / gam_Q_wright$df.residual
  disp_pi <- summary(gam_pi_wright)$scale

  cat("[GAM diagnostics] Tier 3 — signal-preservation panels...\n")
  spear_eta_raw <- stats::cor(msd_data$S_eta, msd_data$S_Wright_signed,
                              method = "spearman", use = "complete.obs")
  spear_eta_gam <- stats::cor(msd_data$S_eta, msd_data$S_Wright_GAM_signed,
                              method = "spearman", use = "complete.obs")
  spear_eta_pi  <- stats::cor(msd_data$S_eta, msd_data$S_Wright_pi_signed,
                              method = "spearman", use = "complete.obs")
  rank_Q  <- .plot_rank_displacement(msd_data, "S_Wright_signed",
                                     "S_Wright_GAM_signed", "Q")
  rank_pi <- .plot_rank_displacement(msd_data, "S_Wright_signed",
                                     "S_Wright_pi_signed",
                                     "π (vs raw Q)")

  # Pass/fail evaluation -----------------------------------------------------
  shr_col <- "shrunk (GAM-RE)"

  rmse_Q <- sim_Qx$rmse_table |>
    tidyr::pivot_wider(id_cols = N_strata, names_from = estimator,
                       values_from = rmse)
  rmse_pi <- sim_pix$rmse_table |>
    tidyr::pivot_wider(id_cols = N_strata, names_from = estimator,
                       values_from = rmse)

  pass_sim_Q  <- isTRUE(rmse_Q[[shr_col]][1]  <= rmse_Q[["raw"]][1])
  pass_sim_pi <- isTRUE(rmse_pi[[shr_col]][1] <= rmse_pi[["raw"]][1])

  pass_bin_Q  <- !is.null(bin_Q$plot)  && isTRUE(bin_Q$pearson  > 0.95) &&
                 isTRUE(bin_Q$slope  >= 0.9 && bin_Q$slope  <= 1.1)
  pass_bin_pi <- !is.null(bin_pi$plot) && isTRUE(bin_pi$pearson > 0.95) &&
                 isTRUE(bin_pi$slope >= 0.9 && bin_pi$slope <= 1.1)

  pass_thr <- !is.null(thr_x$ci_table) && {
    raw_ci <- thr_x$ci_table |> dplyr::filter(grepl("raw",    source))
    shr_ci <- thr_x$ci_table |> dplyr::filter(grepl("shrunk", source))
    isTRUE(is.finite(shr_ci$median) && is.finite(raw_ci$lo) &&
             is.finite(raw_ci$hi) &&
             shr_ci$median >= raw_ci$lo && shr_ci$median <= raw_ci$hi)
  }

  pass_bulmer <- !is.null(bulmer_x$plot) &&
                 isTRUE(bulmer_x$cor_bulmer_gam > 0.99) &&
                 isTRUE(bulmer_x$slope >= 0.95 && bulmer_x$slope <= 1.05)

  cat("\n=========================================================\n")
  cat("[GAM diagnostics] Tier 1 pass/fail summary\n")
  cat("=========================================================\n")
  cat(sprintf("1.1 Q  simulation (low-N RMSE shrunk ≤ raw):    %s\n",
              if (pass_sim_Q)  "PASS" else "FAIL"))
  cat(sprintf("1.1 π  simulation (low-N RMSE shrunk ≤ raw):    %s\n",
              if (pass_sim_pi) "PASS" else "FAIL"))
  cat(sprintf("1.2 Q  bin consistency (r>0.95, slope∈[0.9,1.1]): %s  (r=%.3f, slope=%.3f)\n",
              if (pass_bin_Q) "PASS" else "FAIL",
              bin_Q$pearson, bin_Q$slope))
  cat(sprintf("1.2 π  bin consistency (r>0.95, slope∈[0.9,1.1]): %s  (r=%.3f, slope=%.3f)\n",
              if (pass_bin_pi) "PASS" else "FAIL",
              bin_pi$pearson, bin_pi$slope))
  cat(sprintf("1.3 thr_eta robustness (shrunk median in raw 95%% CI): %s\n",
              if (pass_thr) "PASS" else "FAIL"))
  cat(sprintf("1.4 Bulmer vs exact Wright (r>0.99, slope∈[0.95,1.05]): %s  (r=%.4f, slope=%.3f, max|Δ|=%.3f)\n",
              if (pass_bulmer) "PASS" else "FAIL",
              bulmer_x$cor_bulmer_gam, bulmer_x$slope,
              bulmer_x$max_abs_diff))

  cat(sprintf("\n[Tier 2] Q-GAM dispersion (binomial, target ≈ 1):    %.3f\n",
              disp_Q))
  cat(sprintf("[Tier 2] π-GAM scale (gaussian σ², weighted):        %.5g\n",
              disp_pi))
  cat(sprintf("[Tier 2] Q  RE sd (REML, logit scale):  %s\n",
              if (is.na(blup_Q$re_sd)) "NA" else sprintf("%.3f", blup_Q$re_sd)))
  cat(sprintf("[Tier 2] π  RE sd (REML, π scale):      %s\n",
              if (is.na(blup_pi$re_sd)) "NA" else sprintf("%.3f", blup_pi$re_sd)))

  cat(sprintf("\n[Tier 3] Spearman(S_eta, S_Wright_raw):     %.3f\n",
              spear_eta_raw))
  cat(sprintf("[Tier 3] Spearman(S_eta, S_Wright_GAM):     %.3f  (Δ = %+.3f)\n",
              spear_eta_gam, spear_eta_gam - spear_eta_raw))
  cat(sprintf("[Tier 3] Spearman(S_eta, S_Wright_pi):      %.3f  (Δ = %+.3f)\n",
              spear_eta_pi,  spear_eta_pi  - spear_eta_raw))

  if (!is.na(rank_Q$pct_movers_low_half)) {
    cat(sprintf("[Tier 3] Q  movers (|Δrank|≥1000): %d, %.0f%% in low-N half\n",
                rank_Q$n_movers, 100 * rank_Q$pct_movers_low_half))
  }
  cat("=========================================================\n\n")

  # Save multi-panel PDF -----------------------------------------------------
  panels <- list(sim_Qx$plot,    sim_pix$plot,
                 bin_Q$plot,     bin_pi$plot,
                 thr_x$plot,     bulmer_x$plot,
                 blup_Q$plot,    blup_pi$plot,
                 shrink_Q$plot,  shrink_pi$plot,
                 rank_Q$plot,    rank_pi$plot)
  panels <- panels[!vapply(panels, is.null, logical(1))]

  if (length(panels) > 0) {
    combined <- patchwork::wrap_plots(panels, ncol = 2)
    ggplot2::ggsave(output_pdf, combined,
                    width = 14, height = 4 * ceiling(length(panels) / 2),
                    limitsize = FALSE, device = cairo_pdf)
    cat(sprintf("[GAM diagnostics] Diagnostic PDF written to: %s\n",
                output_pdf))
  }

  invisible(list(
    sim_Q = sim_Qx, sim_pi = sim_pix,
    bin_Q = bin_Q,  bin_pi = bin_pi,
    thr_raw = thr_raw, thr_shr = thr_shr, thr = thr_x,
    bulmer = bulmer_x,
    blup_Q = blup_Q, blup_pi = blup_pi,
    shrink_Q = shrink_Q, shrink_pi = shrink_pi,
    rank_Q = rank_Q, rank_pi = rank_pi,
    disp_Q = disp_Q, disp_pi = disp_pi,
    pass = list(sim_Q = pass_sim_Q, sim_pi = pass_sim_pi,
                bin_Q = pass_bin_Q, bin_pi = pass_bin_pi,
                thr   = pass_thr,
                bulmer = pass_bulmer),
    spear = list(eta_raw = spear_eta_raw,
                 eta_gam = spear_eta_gam,
                 eta_pi  = spear_eta_pi),
    truth = truth
  ))
}
