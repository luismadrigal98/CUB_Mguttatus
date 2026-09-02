## ============================================================================
## detect_preferred_codons.R
##
## ROC-free identification of candidate optimal (preferred) codons.
##
## Replaces the AnaCoDa/ROC-SEMPPR call `which.min(Delta_eta)` with an
## estimator fitted directly to the data of this study.
##
## ---------------------------------------------------------------------------
## The two-regime problem (this is the whole design rationale)
## ---------------------------------------------------------------------------
## Codon share is NOT a monotone function of expression in M. guttatus.  Pooled
## over all genes, GC3-ending codons *decline* steadily with expression across
## roughly the bottom 95% of genes, and then *reverse sharply* in the extreme
## upper tail.  Alanine is typical: the GCC share falls 0.262 -> 0.185 from the
## first decile to the 90-92nd percentile, then climbs back to 0.283 in the top
## 0.5% of genes.  Phenylalanine TTC does the same (0.564 -> 0.451 -> 0.633).
##
## The two regimes have different causes and only one of them is selection:
##   * bulk regime  - drift/mutation dominated.  ~99.8% of genes sit below the
##     drift barrier, so their codon composition tracks AT-biased mutation, not
##     translational optimality.
##   * tail regime  - the small "elite" set where Ns is large enough for
##     selection on codon usage to be effective.  This is where the direction
##     of selection is identifiable.
##
## A single linear expression slope averages the two and returns the *bulk*
## answer (T-ending codons), which is a statement about mutation bias, not
## about optimality.  So the preferred codon is defined here by the REGIME
## CONTRAST: the derivative of a smooth fit at a high expression quantile MINUS
## the derivative at the median.  Ranking on the tail derivative alone is not
## enough -- a codon can be rising at high expression simply because it is also
## rising through the drift-dominated bulk (Arg_4 CGT and Gly GGT both do).
## What identifies an optimal codon is that its trajectory bends upward
## specifically where selection becomes effective.
##
## This is also why the ROC-SEMPPR trajectory figure drew the referee's
## objection that "most of the curve dynamics are beyond the extent of the data
## points": ROC asks which codon wins as phi -> infinity, i.e. it reports the
## tail regime, but it draws that answer as a sigmoid over a phi range the data
## barely populate.  The estimator below reports the same quantity on the scale
## of the data, with the regime structure made explicit instead of extrapolated.
##
## Why GC12 and not GC3s as the composition control: GC3s is largely a
## restatement of synonymous codon choice, so conditioning on it would remove
## the signal being estimated.  GC12 sits at non-degenerate positions and
## indexes the regional base composition a gene sits in (mutation bias, and any
## GC-biased gene conversion acting on its neighbourhood) without being a
## direct function of the response.  Same control as
## `fit_pairwise_binomial_models.R`.
## ============================================================================

detect_preferred_codons <- function(codon_counts,
                                    gene_meta,
                                    genetic_code,
                                    expression_var = "Mean_Log10_Exp",
                                    composition_var = "GC12",
                                    length_var = "CDS_length_nt",
                                    tail_quantile = 0.99,
                                    bulk_quantile = 0.50,
                                    min_family_count = 5L,
                                    k_smooth = 8L,
                                    fdr_level = 0.05,
                                    control_composition = TRUE,
                                    verbose = TRUE) {
  #' Identify candidate optimal codons from the high-expression regime
  #'
  #' For every codon, fits a quasi-binomial GAM of that codon's share within its
  #' synonymous family against a smooth of expression, controlling for
  #' background composition and CDS length.  The preferred codon of a family is
  #' the one with the largest, FDR-significant REGIME CONTRAST: the derivative
  #' at `tail_quantile` minus the derivative at `bulk_quantile`.  Note this is
  #' not the same as the largest tail derivative, and the two can disagree.
  #'
  #' @param codon_counts data.table/data.frame with `Gene_name` plus one integer
  #'   column per sense codon (`data/codon_counts_from_fasta.rds`).
  #' @param gene_meta data.table/data.frame with `Gene_name`, the expression
  #'   variable, the composition variable and the CDS length variable.
  #' @param genetic_code named character vector, codon -> family label.  Use
  #'   `genetic_code_dna_long` for split 2/4-fold families (Leu_2/Leu_4,
  #'   Ser_2/Ser_4, Arg_2/Arg_4).
  #' @param expression_var log10 expression column.  Project convention:
  #'   `Mean_Log10_Exp` standalone, `Max_Log10_Exp` only alongside breadth.
  #' @param composition_var background composition control (default `GC12`).
  #' @param length_var CDS length in nucleotides; entered as log10.
  #' @param tail_quantile expression quantile at which the selection-regime
  #'   derivative is evaluated (default 0.99).
  #' @param bulk_quantile expression quantile for the drift-regime derivative,
  #'   reported for contrast (default 0.50).
  #' @param min_family_count minimum codons of the family a gene must carry to
  #'   contribute (default 5; guards against shares estimated from 1-2 codons).
  #' @param k_smooth basis dimension for the expression smooth.
  #' @param fdr_level Benjamini-Hochberg level applied across all codons.
  #' @param control_composition if FALSE, drops `composition_var` (sensitivity:
  #'   is the call driven by the GC control?).
  #' @param verbose print per-family progress.
  #' @return list with `codon_table` (per codon: tail and bulk derivatives with
  #'   SEs, p, FDR, fitted shares), `preferred` (per family: the winning codon,
  #'   separation from the runner-up, and the call), `regime` (model-free
  #'   elite-vs-bulk contrast) and `settings`.

  suppressPackageStartupMessages({
    require(data.table)
    require(mgcv)
  })

  codon_counts <- data.table::as.data.table(codon_counts)
  gene_meta    <- data.table::as.data.table(gene_meta)

  needed <- c("Gene_name", expression_var, length_var,
              if (control_composition) composition_var)
  missing_cols <- setdiff(needed, names(gene_meta))
  if (length(missing_cols)) {
    stop("gene_meta is missing: ", paste(missing_cols, collapse = ", "))
  }

  meta <- gene_meta[, ..needed]
  data.table::setnames(meta, expression_var, "Expr")
  data.table::setnames(meta, length_var, "CDS_len")
  if (control_composition) data.table::setnames(meta, composition_var, "Comp")

  dat <- merge(codon_counts, meta, by = "Gene_name")
  dat[, LogLen := log10(CDS_len)]
  dat <- dat[stats::complete.cases(
    dat[, c("Expr", "LogLen", if (control_composition) "Comp"), with = FALSE]
  )]

  expr_tail <- unname(stats::quantile(dat$Expr, tail_quantile, na.rm = TRUE))
  expr_bulk <- unname(stats::quantile(dat$Expr, bulk_quantile, na.rm = TRUE))
  ref_len   <- stats::median(dat$LogLen, na.rm = TRUE)
  ref_comp  <- if (control_composition) stats::median(dat$Comp, na.rm = TRUE) else NA_real_

  code <- genetic_code[genetic_code != "STOP"]
  code <- code[names(code) %in% names(dat)]
  families <- split(names(code), unname(code))
  families <- families[vapply(families, length, 1L) >= 2L]

  rhs <- paste(c(sprintf("s(Expr, k = %d)", k_smooth),
                 if (control_composition) "Comp", "LogLen"), collapse = " + ")

  # Finite-difference derivative rows of the linear predictor (standard mgcv
  # lpmatrix approach).  Returned as the contrast row itself so that the tail
  # derivative, the bulk derivative and their *difference* can all be tested
  # against the same model covariance -- the two derivatives come from one fit
  # and are correlated, so the difference needs the joint covariance, not
  # sqrt(se1^2 + se2^2).
  deriv_row <- function(fit, e) {
    eps <- 1e-5
    nd0 <- data.frame(Expr = e,       LogLen = ref_len)
    nd1 <- data.frame(Expr = e + eps, LogLen = ref_len)
    if (control_composition) { nd0$Comp <- ref_comp; nd1$Comp <- ref_comp }
    X0 <- stats::predict(fit, newdata = nd0, type = "lpmatrix")
    X1 <- stats::predict(fit, newdata = nd1, type = "lpmatrix")
    (X1 - X0) / eps
  }

  est_se <- function(fit, Xrow) {
    est <- as.numeric(Xrow %*% stats::coef(fit))
    se  <- sqrt(as.numeric(rowSums((Xrow %*% stats::vcov(fit)) * Xrow)))
    c(est = est, se = se)
  }

  fitted_share_at <- function(fit, e) {
    nd <- data.frame(Expr = e, LogLen = ref_len)
    if (control_composition) nd$Comp <- ref_comp
    as.numeric(stats::predict(fit, newdata = nd, type = "response"))
  }

  rows <- list()
  for (fam in names(families)) {
    codons <- families[[fam]]
    fam_dat <- dat[, c("Expr", "LogLen", if (control_composition) "Comp", codons),
                   with = FALSE]
    fam_dat[, N_family := rowSums(.SD), .SDcols = codons]
    fam_dat <- fam_dat[N_family >= min_family_count]

    if (nrow(fam_dat) < 200L) {
      if (verbose) message(sprintf("  [skip] %s: only %d usable genes", fam, nrow(fam_dat)))
      next
    }

    for (cod in codons) {
      fam_dat[, Y := get(cod)]
      fit <- try(
        suppressWarnings(mgcv::gam(
          stats::as.formula(paste("cbind(Y, N_family - Y) ~", rhs)),
          family = stats::quasibinomial(), data = fam_dat, method = "REML"
        )),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        rows[[length(rows) + 1L]] <- data.table::data.table(
          Family = fam, Codon = cod, Third_base = substr(cod, 3, 3),
          N_genes = nrow(fam_dat), Mean_share = NA_real_,
          Slope_tail = NA_real_, SE_tail = NA_real_, CI_low = NA_real_,
          CI_high = NA_real_, p_value = NA_real_, Slope_bulk = NA_real_,
          SE_bulk = NA_real_, Share_bulk = NA_real_, Share_tail = NA_real_,
          Converged = FALSE
        )
        next
      }

      Xt <- deriv_row(fit, expr_tail)
      Xb <- deriv_row(fit, expr_bulk)
      dt <- est_se(fit, Xt)
      db <- est_se(fit, Xb)
      dc <- est_se(fit, Xt - Xb)          # regime contrast, joint covariance
      z  <- dt[["est"]] / dt[["se"]]
      zc <- dc[["est"]] / dc[["se"]]

      rows[[length(rows) + 1L]] <- data.table::data.table(
        Family = fam, Codon = cod, Third_base = substr(cod, 3, 3),
        N_genes = nrow(fam_dat),
        Mean_share = sum(fam_dat[[cod]]) / sum(fam_dat$N_family),
        Slope_tail = dt[["est"]], SE_tail = dt[["se"]],
        CI_low  = dt[["est"]] - 1.96 * dt[["se"]],
        CI_high = dt[["est"]] + 1.96 * dt[["se"]],
        p_value = 2 * stats::pnorm(-abs(z)),
        Slope_bulk = db[["est"]], SE_bulk = db[["se"]],
        Slope_contrast = dc[["est"]], SE_contrast = dc[["se"]],
        Contrast_CI_low  = dc[["est"]] - 1.96 * dc[["se"]],
        Contrast_CI_high = dc[["est"]] + 1.96 * dc[["se"]],
        p_contrast = 2 * stats::pnorm(-abs(zc)),
        Share_bulk = fitted_share_at(fit, expr_bulk),
        Share_tail = fitted_share_at(fit, expr_tail),
        Converged = TRUE
      )
    }
    if (verbose) message(sprintf("  [ok] %s (%d codons, %d genes)",
                                 fam, length(codons), nrow(fam_dat)))
  }

  codon_table <- data.table::rbindlist(rows)
  codon_table[, FDR := stats::p.adjust(p_value, method = "BH")]
  codon_table[, FDR_contrast := stats::p.adjust(p_contrast, method = "BH")]
  codon_table[, Significant := !is.na(FDR) & FDR < fdr_level]
  # Does this codon reverse direction between the two regimes?
  codon_table[, Regime_switch := !is.na(Slope_bulk) & !is.na(Slope_tail) &
                sign(Slope_bulk) != sign(Slope_tail)]
  data.table::setorder(codon_table, Family, -Slope_contrast)

  # Ranking criterion: the *regime contrast* (tail derivative minus bulk
  # derivative), not the tail derivative alone.  A codon can have a positive
  # tail slope simply because it is also rising through the drift-dominated
  # bulk -- Arg_4 CGT and Gly GGT both do.  What identifies an optimal codon is
  # that its trajectory bends upward specifically where selection becomes
  # effective, i.e. a negative-to-positive reversal.  Ranking on the contrast
  # separates Gly (GGC vs GGT, tail slopes 0.150 vs 0.147 -- effectively tied)
  # and Arg_4 (CGT leads on tail slope but rises in both regimes, while CGC
  # shows the reversal).
  preferred <- codon_table[Converged == TRUE][
    order(Family, -Slope_contrast),
    {
      top <- .SD[1L]
      runner <- if (.N >= 2L) .SD[2L] else NULL
      list(
        Preferred_Codon = top$Codon,
        Third_base      = top$Third_base,
        Slope_contrast  = top$Slope_contrast,
        Contrast_CI_low = top$Contrast_CI_low,
        Contrast_CI_high = top$Contrast_CI_high,
        FDR_contrast    = top$FDR_contrast,
        Slope_tail      = top$Slope_tail,
        Slope_bulk      = top$Slope_bulk,
        Regime_switch   = top$Regime_switch,
        Share_bulk      = top$Share_bulk,
        Share_tail      = top$Share_tail,
        Runner_up       = if (is.null(runner)) NA_character_ else runner$Codon,
        Runner_up_contrast = if (is.null(runner)) NA_real_ else runner$Slope_contrast,
        Separated       = if (is.null(runner)) NA else
                            top$Contrast_CI_low > runner$Contrast_CI_high,
        Call            = if (top$Slope_contrast <= 0 || top$Slope_tail <= 0) {
                            "none (no codon favoured at high expression)"
                          } else if (!isTRUE(top$FDR_contrast < fdr_level)) {
                            "unresolved (regime contrast not significant)"
                          } else {
                            "preferred"
                          }
      )
    },
    by = Family
  ]

  list(
    codon_table = codon_table[],
    preferred   = preferred[],
    regime      = elite_vs_bulk_contrast(codon_counts, gene_meta, genetic_code,
                                         expression_var = expression_var,
                                         min_family_count = min_family_count),
    settings    = list(
      expression_var = expression_var,
      composition_var = if (control_composition) composition_var else NA_character_,
      control_composition = control_composition,
      tail_quantile = tail_quantile, bulk_quantile = bulk_quantile,
      expr_at_tail = expr_tail, expr_at_bulk = expr_bulk,
      min_family_count = min_family_count,
      fdr_level = fdr_level, n_genes = nrow(dat)
    )
  )
}


elite_vs_bulk_contrast <- function(codon_counts, gene_meta, genetic_code,
                                   expression_var = "Mean_Log10_Exp",
                                   elite_quantile = 0.99,
                                   bulk_range = c(0.50, 0.95),
                                   min_family_count = 5L) {
  #' Model-free elite-vs-bulk codon share contrast
  #'
  #' Companion to `detect_preferred_codons()` that needs no link function and no
  #' smooth: compares each codon's share in the top `1 - elite_quantile` of
  #' genes against its share in the mid-to-upper bulk.  Reports both the pooled
  #' share and the per-gene median (the pooled figure alone would be vulnerable
  #' to a handful of codon-rich genes), plus a Wilcoxon test on per-gene shares.
  #' Written for the Methods: it is the whole result in one table.

  suppressPackageStartupMessages(require(data.table))

  codon_counts <- data.table::as.data.table(codon_counts)
  meta <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(meta, expression_var, "Expr")
  dat <- merge(codon_counts, meta, by = "Gene_name")[!is.na(Expr)]

  q_elite <- stats::quantile(dat$Expr, elite_quantile)
  q_lo    <- stats::quantile(dat$Expr, bulk_range[1])
  q_hi    <- stats::quantile(dat$Expr, bulk_range[2])

  code <- genetic_code[genetic_code != "STOP"]
  code <- code[names(code) %in% names(dat)]
  families <- split(names(code), unname(code))
  families <- families[vapply(families, length, 1L) >= 2L]

  out <- lapply(names(families), function(fam) {
    codons <- families[[fam]]
    x <- dat[, c("Expr", codons), with = FALSE]
    x[, N_family := rowSums(.SD), .SDcols = codons]
    x <- x[N_family >= min_family_count]
    elite <- x[Expr >= q_elite]
    bulk  <- x[Expr >= q_lo & Expr < q_hi]
    if (nrow(elite) < 30L || nrow(bulk) < 30L) return(NULL)

    data.table::rbindlist(lapply(codons, function(cod) {
      se <- elite[[cod]] / elite$N_family
      sb <- bulk[[cod]]  / bulk$N_family
      data.table::data.table(
        Family = fam, Codon = cod, Third_base = substr(cod, 3, 3),
        N_elite = nrow(elite), N_bulk = nrow(bulk),
        Share_elite_pooled  = sum(elite[[cod]]) / sum(elite$N_family),
        Share_bulk_pooled   = sum(bulk[[cod]])  / sum(bulk$N_family),
        Share_elite_median  = stats::median(se),
        Share_bulk_median   = stats::median(sb),
        Delta_pooled = sum(elite[[cod]]) / sum(elite$N_family) -
                       sum(bulk[[cod]])  / sum(bulk$N_family),
        p_wilcox = tryCatch(stats::wilcox.test(se, sb)$p.value,
                            error = function(e) NA_real_)
      )
    }))
  })

  res <- data.table::rbindlist(out)
  if (nrow(res)) {
    res[, FDR := stats::p.adjust(p_wilcox, method = "BH")]
    data.table::setorder(res, Family, -Delta_pooled)
  }
  res[]
}


codon_share_by_expression_bin <- function(codon_counts, gene_meta, genetic_code,
                                          expression_var = "Mean_Log10_Exp",
                                          probs = c(seq(0, 0.90, by = 0.10),
                                                    0.92, 0.94, 0.96, 0.98,
                                                    0.99, 0.995, 1),
                                          min_family_count = 5L) {
  #' Model-free codon share across expression bins
  #'
  #' The raw pattern the GAM summarises, with the default bin edges deliberately
  #' dense in the upper tail — a uniform decile grid hides the regime reversal
  #' inside its top bin.  Returns a long data.table
  #' (Family, Codon, Third_base, Bin, Share, Mean_expression, N_genes).

  suppressPackageStartupMessages(require(data.table))

  codon_counts <- data.table::as.data.table(codon_counts)
  meta <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(meta, expression_var, "Expr")
  dat <- merge(codon_counts, meta, by = "Gene_name")[!is.na(Expr)]

  code <- genetic_code[genetic_code != "STOP"]
  code <- code[names(code) %in% names(dat)]
  families <- split(names(code), unname(code))
  families <- families[vapply(families, length, 1L) >= 2L]

  out <- lapply(names(families), function(fam) {
    codons <- families[[fam]]
    x <- dat[, c("Expr", codons), with = FALSE]
    x[, N_family := rowSums(.SD), .SDcols = codons]
    x <- x[N_family >= min_family_count]
    x[, Bin := cut(Expr, stats::quantile(Expr, probs),
                   include.lowest = TRUE, labels = FALSE)]
    agg <- x[, c(lapply(.SD, sum),
                 list(N_family = sum(N_family),
                      Mean_expression = mean(Expr), N_genes = .N)),
             .SDcols = codons, by = Bin]
    long <- data.table::melt(agg,
                             id.vars = c("Bin", "N_family", "Mean_expression", "N_genes"),
                             variable.name = "Codon", value.name = "Count")
    long[, `:=`(Family = fam,
                Third_base = substr(as.character(Codon), 3, 3),
                Share = Count / N_family)]
    long[]
  })

  data.table::rbindlist(out)[order(Family, Codon, Bin)]
}


compare_preferred_codon_sets <- function(preferred_new, roc_codons, genetic_code) {
  #' Concordance between the expression-regime call and the AnaCoDa/ROC list
  #'
  #' @param preferred_new the `preferred` element of `detect_preferred_codons()`
  #' @param roc_codons character vector of ROC-SEMPPR preferred codons
  #'   (`results/preferred_codons.txt`)
  #' @return data.table, one row per family, with both calls and whether they agree.

  suppressPackageStartupMessages(require(data.table))

  roc_dt <- data.table::data.table(
    Family = unname(genetic_code[roc_codons]),
    ROC_Codon = roc_codons
  )[!is.na(Family)]

  out <- merge(data.table::as.data.table(preferred_new), roc_dt,
               by = "Family", all = TRUE)
  out[, Agrees := !is.na(ROC_Codon) & !is.na(Preferred_Codon) &
        ROC_Codon == Preferred_Codon]
  out[, `:=`(Third_base_new = substr(Preferred_Codon, 3, 3),
             Third_base_ROC = substr(ROC_Codon, 3, 3))]
  out[]
}


plot_codon_regimes <- function(bin_table, families = NULL, output_file = NULL,
                               title = "Codon share across the expression range") {
  #' Two-regime figure: codon share vs expression, on the scale of the data
  #'
  #' This is the honest replacement for the ROC-SEMPPR trajectory panel the
  #' referee objected to.  Every point is an observed pooled share in an
  #' expression bin; nothing is extrapolated, and the drift/selection regime
  #' reversal in the upper tail is visible rather than implied.

  suppressPackageStartupMessages({ require(ggplot2); require(data.table) })

  d <- data.table::as.data.table(bin_table)
  if (!is.null(families)) d <- d[Family %in% families]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = Mean_expression, y = Share,
                                       colour = Third_base, group = Codon)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::facet_wrap(~ Family, scales = "free_y", ncol = 5) +
    ggplot2::scale_colour_manual(values = c("A" = "#3B7DD8", "C" = "#E08214",
                                            "G" = "#1B7837", "T" = "#B2182B"),
                                 name = "3rd base") +
    ggplot2::labs(title = title,
                  x = expression("Mean log"[10] * " expression"),
                  y = "Share within synonymous family") +
    theme_custom() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 14, height = 11, dpi = 300)
  }
  p
}


plot_preferred_codon_slopes <- function(codon_table, preferred,
                                        output_file = NULL,
                                        title = "Codon-share derivative at high expression, by synonymous family") {
  #' Forest plot of per-codon tail derivatives, faceted by family
  #'
  #' The winning codon of each family is highlighted.  Bulk-regime derivatives
  #' are overlaid as hollow points so the regime reversal is legible per codon.

  suppressPackageStartupMessages({ require(ggplot2); require(data.table) })

  d <- data.table::as.data.table(codon_table)[Converged == TRUE]
  win <- data.table::as.data.table(preferred)[, .(Family, Codon = Preferred_Codon)]
  win[, Is_preferred := TRUE]
  d <- merge(d, win, by = c("Family", "Codon"), all.x = TRUE)
  d[is.na(Is_preferred), Is_preferred := FALSE]
  d[, Status := data.table::fifelse(
    Is_preferred, "Preferred (largest drift-to-selection shift)",
    data.table::fifelse(Significant & Slope_tail > 0, "Rises at high expression",
                        data.table::fifelse(Significant, "Falls at high expression",
                                            "Not significant")))]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(Codon, Slope_tail),
                                       y = Slope_tail, colour = Status)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = CI_low, ymax = CI_high), size = 0.35) +
    ggplot2::geom_point(ggplot2::aes(y = Slope_bulk), shape = 21, fill = NA,
                        colour = "grey35", size = 1.6) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ Family, scales = "free_y", ncol = 5) +
    ggplot2::scale_colour_manual(values = c(
      "Preferred (largest drift-to-selection shift)" = "#1B7837",
      "Rises at high expression"                     = "#7FBC41",
      "Falls at high expression"                     = "#D6604D",
      "Not significant"                              = "grey65"
    )) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Filled = derivative at the 99th expression percentile (selection regime);",
                       "hollow = at the median (drift regime).\nThe call is made on the GAP between them,",
                       "so the preferred codon is not always the highest filled point."),
      x = NULL,
      y = "d(logit codon share) / d(log10 expression)",
      colour = NULL) +
    theme_custom() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 14, height = 11, dpi = 300)
  }
  p
}
