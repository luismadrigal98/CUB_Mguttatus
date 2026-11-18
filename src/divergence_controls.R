#' Polymorphism over divergence (pi / d), with M. jungermannioides (JUNG1) as the
#' outgroup and M. glabratus (KONZA) as a check
#'
#' Why pi / d. Under neutrality both pi and d are proportional to the mutational
#' input of a site class; for a two-allele site at compositional equilibrium both
#' are proportional to 2uv / (u + v). So pi / d does not move with the mutation rate
#' OR the mutational bias, at any gene, in any sequence. That makes it the direct
#' answer to Referee 2's proposal of a different mutational input in highly
#' expressed genes, without the exon-specific mutation model we constructed for
#' the exon / intron ratio:
#'   mutation (rate or spectrum, shared or exon-specific)  pi / d flat
#'   weak selection (Wright; McVean & Charlesworth 1999)   pi rises with S while the
#'                                                         substitution rate falls,
#'                                                         so pi / d rises with S
#'   linked selection                                      pi falls, d does not, so
#'                                                         pi / d falls, introns too
#' Introns of the same genes are the class that should stay flat with S.
#'
#' Divergence is IM767 reference vs outgroup, from the AnchorWave proAli
#' alignments (miscellanea_code/divergence_by_gene.py). pi and d are always taken
#' over the same sites: the 4-fold sites joined by coordinate, and the intronic
#' sites through the aligned-position mask given to calculate_intronic_twostatepi.py.
#'
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

wright_pi_2allele <- function(S, U, V) {
  #' Expected two-allele heterozygosity under Wright's distribution (as in
  #' exon_intron_pi_ratio_by_S()).
  vapply(S, function(s) 2 * V * U / ((U + V) * (U + V + 1)) *
           exp(log_kummer_1F1(V + 1, U + V + 2, s) - log_kummer_1F1(V, U + V, s)),
         numeric(1))
}

wright_substitution_rate <- function(S, U, V) {
  #' Substitution rate of a two-allele site at stationarity, weak-mutation limit:
  #' K(S) = 2 Q U f(-S), with Q = V e^S / (V e^S + U) and f(x) = x / (1 - e^-x),
  #' the fixation probability of a new mutant relative to a neutral one. K(0) is
  #' the neutral 2UV / (U + V). Same S as Wright's density exp(S p).
  f <- function(x) ifelse(abs(x) < 1e-8, 1, x / (1 - exp(-x)))
  Q <- V * exp(S) / (V * exp(S) + U)
  2 * Q * U * f(-S)
}

# ---------------------------------------------------------------------------
# Per-gene pi and d on shared sites
# ---------------------------------------------------------------------------

fourfold_pi_d_by_gene <- function(fourfold_sites, div_sites, taxon, fixed = TRUE) {
  #' Per gene: pi and d at the 4-fold sites aligned to `taxon`.
  #'
  #' @param fourfold_sites extract_fourfold_sites() output (needs pos3, ref3).
  #' @param div_sites path or table: divergence_by_gene.py fourfold_sites.tsv.gz.
  #' @param fixed TRUE (default): d counts FIXED differences, the outgroup base
  #'   absent from every line (C system: outgroup C and no line C, or outgroup
  #'   non-C and every line C). FALSE: IM767 reference vs outgroup, which picks up
  #'   polymorphism (37% of JUNG1 and 72% of KONZA reference differences fall on
  #'   segregating sites) and is kept only for comparison.
  #' @return data.table Gene, pi_n, pi_C, pi_tot, d_n, d_C, d_tot (sums over
  #'   sites; pi_n = d_n), with the join diagnostics as attributes.
  b  <- paste0(taxon, "_base")
  dv <- if (is.character(div_sites))
    data.table::fread(div_sites, select = c("Gene", "Pos", "Ref", b)) else
    data.table::as.data.table(div_sites)[, c("Gene", "Pos", "Ref", b), with = FALSE]
  data.table::setnames(dv, b, "Out")
  x <- merge(fourfold_sites[, .(Gene, pos3, ref3, n, A, C, G, T)], dv,
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))
  ref_match <- x[, mean(ref3 == Ref)]
  a <- x[Out %chin% c("A", "C", "G", "T") & ref3 == Ref]
  a[, `:=`(
    pi_tot = n / (n - 1) * (1 - ((A / n)^2 + (C / n)^2 + (G / n)^2 + (T / n)^2)),
    pi_C   = 2 * C * (n - C) / (n * (n - 1)))]
  if (fixed) {
    a[, oc := data.table::fcase(Out == "A", A, Out == "C", C, Out == "G", G, Out == "T", T)]
    a[, `:=`(d_tot = as.numeric(oc == 0L),
             d_C   = as.numeric((Out == "C" & C == 0L) | (Out != "C" & C == n)))]
  } else {
    a[, `:=`(d_tot = as.numeric(Ref != Out),
             d_C   = as.numeric((Ref == "C") != (Out == "C")))]
  }
  g <- a[, .(pi_n = .N, pi_C = sum(pi_C), pi_tot = sum(pi_tot),
             d_n = .N, d_C = sum(d_C), d_tot = sum(d_tot)), by = Gene]
  attr(g, "join") <- c(codon_sites = nrow(fourfold_sites), joined = nrow(x),
                       ref_match = ref_match, aligned = nrow(a))
  g[]
}

intron_pi_d_by_gene <- function(intron_div, intron_masked, taxon, fixed = TRUE) {
  #' Per gene: intronic pi and d over the sites aligned to `taxon`. Same columns
  #' as the 4-fold table.
  #' With fixed = TRUE (default) d is the FIXED-difference count that
  #' calculate_intronic_twostatepi.py writes when run with --outgroup (columns
  #' n_out, fixed_sub, fixed_subC); an older masked table without them falls
  #' back to the reference-based counts of divergence_by_gene.py, with a warning.
  pm <- data.table::fread(intron_masked)
  if (!"pi_total" %in% names(pm))
    stop("[intron pi/d] ", intron_masked, " lacks pi_total: rerun calculate_intronic_twostatepi.py")
  has_fixed <- all(c("n_out", "fixed_sub", "fixed_subC") %in% names(pm))
  if (fixed && has_fixed) {
    out <- pm[, .(Gene, d_n = n_out, d_tot = fixed_sub, d_C = fixed_subC,
                  pi_n = n_sites, pi_C = pi_2allele_C * n_sites, pi_tot = pi_total * n_sites)]
    return(out[d_n > 0 & pi_n > 0])
  }
  if (fixed)
    warning("[intron pi/d] ", intron_masked, " has no fixed-difference columns; intronic d is ",
            "reference-based (rerun the intronic script with OUTGROUPS=...)", call. = FALSE)
  dv <- data.table::fread(intron_div)
  dv <- dv[, .(Gene, d_n = get(paste0(taxon, "_n_aln")), d_tot = get(paste0(taxon, "_sub")),
               d_C = get(paste0(taxon, "_subC")))]
  pm <- pm[, .(Gene, pi_n = n_sites, pi_C = pi_2allele_C * n_sites, pi_tot = pi_total * n_sites)]
  merge(dv, pm, by = "Gene")[d_n > 0 & pi_n > 0]
}

# ---------------------------------------------------------------------------
# pi / d by S_Wright bin and by expression category
# ---------------------------------------------------------------------------

.pi_d_stats <- function(g, i, m) {
  pi <- sum(g[[paste0("pi_", m)]][i]) / sum(g$pi_n[i])
  d  <- sum(g[[paste0("d_", m)]][i]) / sum(g$d_n[i])
  c(pi = pi, d = d, pd = pi / d)
}

pi_d_by_S <- function(gene_tab, meta, U, V, measure = "C",
                      S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                      min_d_n = 20L, n_boot = 2000L, seed = 1L) {
  #' pi, d and pi / d per S_Wright bin, each relative to S = 0 (linear
  #' interpolation between the [-Inf,0) and [0,0.25) bins, as in
  #' exon_intron_pi_ratio_by_S()), with the Wright predictions.
  #' Genes are resampled within bins; the S = 0 reference is formed inside each
  #' replicate, so intervals carry its uncertainty.
  #' @param measure "C" (two-allele C vs other, Wright's system) or "tot".
  set.seed(seed)
  g <- merge(data.table::as.data.table(gene_tab), meta[, .(Gene, S)], by = "Gene")
  g <- g[is.finite(S) & d_n >= min_d_n]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  bins <- levels(droplevels(g$Sbin))
  idx  <- lapply(bins, function(b) which(g$Sbin == b)); names(idx) <- bins
  mS   <- vapply(idx, function(i) mean(g$S[i]), numeric(1))
  stopifnot(all(c("[-Inf,0)", "[0,0.25)") %in% bins))
  w0   <- -mS["[-Inf,0)"] / (mS["[0,0.25)"] - mS["[-Inf,0)"])
  est  <- sapply(idx, function(i) .pi_d_stats(g, i, measure))          # 3 x bins
  ref0 <- (1 - w0) * est[, "[-Inf,0)"] + w0 * est[, "[0,0.25)"]
  bt   <- replicate(n_boot, sapply(idx, function(i)
    .pi_d_stats(g, i[sample.int(length(i), length(i), TRUE)], measure)))   # 3 x bins x B
  ref0_b <- (1 - w0) * bt[, "[-Inf,0)", ] + w0 * bt[, "[0,0.25)", ]
  rel_b  <- sweep(bt, c(1, 3), ref0_b, "/")
  q <- function(k, p) apply(rel_b[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  pred_pi <- wright_pi_2allele(mS, U, V) / wright_pi_2allele(0, U, V)
  pred_d  <- wright_substitution_rate(mS, U, V) / wright_substitution_rate(0, U, V)
  data.table::data.table(
    Sbin = bins, n_genes = lengths(idx), mean_S = mS,
    pi = est["pi", ], d = est["d", ], pi_d = est["pd", ],
    pi_rel = est["pi", ] / ref0["pi"], pi_lo = q("pi", 0.025), pi_hi = q("pi", 0.975),
    d_rel = est["d", ] / ref0["d"], d_lo = q("d", 0.025), d_hi = q("d", 0.975),
    pd_rel = est["pd", ] / ref0["pd"], pd_lo = q("pd", 0.025), pd_hi = q("pd", 0.975),
    pred_pi = pred_pi, pred_d = pred_d, pred_pd = pred_pi / pred_d)
}

pi_d_by_expression <- function(gene_tab, meta, measure = "tot", bin_width = 0.2,
                               bulk_max = 1.2, min_genes = 40L, min_d_n = 20L,
                               n_boot = 2000L, seed = 1L) {
  #' pi, d and pi / d per expression category, relative to the drift-dominated
  #' bulk (categories <= bulk_max), bootstrapped jointly. Separates a mutational
  #' rise in the elite genes (pi and d up together, pi / d flat) from relief of
  #' linked selection (pi up, d flat, pi / d up), and gives divergence against
  #' expression directly (Referee 2's long-term question).
  set.seed(seed)
  g <- merge(data.table::as.data.table(gene_tab), meta[, .(Gene, Mean_Log10_Exp)], by = "Gene")
  g <- g[is.finite(Mean_Log10_Exp) & d_n >= min_d_n]
  g[, Exp_cat := round(round(Mean_Log10_Exp / bin_width) * bin_width, 1)]
  cats <- sort(unique(g$Exp_cat))
  idx  <- lapply(cats, function(cc) which(g$Exp_cat == cc))
  keep <- lengths(idx) >= min_genes
  cats <- cats[keep]; idx <- idx[keep]
  bulk <- which(g$Exp_cat <= bulk_max)
  est  <- sapply(idx, function(i) .pi_d_stats(g, i, measure))
  ref  <- .pi_d_stats(g, bulk, measure)
  bt   <- replicate(n_boot, {
    r <- .pi_d_stats(g, bulk[sample.int(length(bulk), length(bulk), TRUE)], measure)
    sapply(idx, function(i) .pi_d_stats(g, i[sample.int(length(i), length(i), TRUE)], measure)) / r
  })
  q <- function(k, p) apply(bt[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  data.table::data.table(
    Exp_cat = cats, n_genes = lengths(idx),
    pi = est["pi", ], d = est["d", ], pi_d = est["pd", ],
    pi_rel = est["pi", ] / ref["pi"], pi_lo = q("pi", 0.025), pi_hi = q("pi", 0.975),
    d_rel = est["d", ] / ref["d"], d_lo = q("d", 0.025), d_hi = q("d", 0.975),
    pd_rel = est["pd", ] / ref["pd"], pd_lo = q("pd", 0.025), pd_hi = q("pd", 0.975))
}

# ---------------------------------------------------------------------------
# Plots (no titles: the legend text lives in the manuscript)
# ---------------------------------------------------------------------------

.seq_pal <- c("4-fold sites" = "#3B7DD8", "Introns, same genes" = "grey35")

plot_pi_over_d_by_S <- function(tab, check = NULL, title = NULL) {
  #' pi / d relative to S = 0 by S_Wright bin, 4-fold sites and introns, with
  #' Wright's selection prediction for 4-fold sites (solid) and the prediction of
  #' any mutational explanation (1, dashed). `check` (same layout, the KONZA
  #' outgroup) is drawn as hollow symbols.
  #' @param tab data.table with Sequence, mean_S, pd_rel, pd_lo, pd_hi, pred_pd.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  d[, Outgroup := "JUNG1"]
  if (!is.null(check)) d <- rbind(d, data.table::as.data.table(check)[, Outgroup := "KONZA"], fill = TRUE)
  dodge <- ifelse(d$Sequence == "4-fold sites", -0.02, 0.02) + ifelse(d$Outgroup == "KONZA", 0.04, 0)
  d[, x := mean_S + dodge]
  pred <- data.table::as.data.table(tab)[Sequence == "4-fold sites", .(mean_S, pred_pd)]
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey65", linewidth = 0.4) +
    geom_line(data = pred, aes(mean_S, pred_pd, linetype = "Selection (Wright), 4-fold"),
              colour = "grey15", linewidth = 0.8) +
    geom_errorbar(data = d, aes(x = x, ymin = pd_lo, ymax = pd_hi, colour = Sequence),
                  width = 0.04, linewidth = 0.5) +
    geom_point(data = d, aes(x = x, y = pd_rel, colour = Sequence, shape = Sequence,
                             fill = interaction(Sequence, Outgroup)), size = 2.4, stroke = 0.9) +
    scale_colour_manual(values = .seq_pal, name = NULL) +
    scale_shape_manual(values = c("4-fold sites" = 21, "Introns, same genes" = 24), name = NULL) +
    scale_fill_manual(values = c("4-fold sites.JUNG1" = "#3B7DD8", "Introns, same genes.JUNG1" = "grey35",
                                 "4-fold sites.KONZA" = "white", "Introns, same genes.KONZA" = "white"),
                      guide = "none") +
    scale_linetype_manual(values = c("Selection (Wright), 4-fold" = "solid"), name = NULL) +
    labs(x = expression(S[Wright] * " (bin mean)"),
         y = expression(pi / d * ", relative to " * S == 0), title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
}

plot_pi_d_by_expression <- function(tab, what = c("pd", "d", "pi"), title = NULL) {
  #' pi / d (or d, or pi) by expression category relative to the bulk, 4-fold
  #' sites and introns, gene-bootstrap bands.
  what <- match.arg(what)
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  y <- paste0(what, "_rel"); lo <- paste0(what, "_lo"); hi <- paste0(what, "_hi")
  ylab <- switch(what, pd = expression(pi / d * ", relative to the bulk"),
                 d = "Divergence, relative to the bulk", pi = expression(pi * ", relative to the bulk"))
  ggplot(d, aes(Exp_cat, .data[[y]], colour = Sequence, fill = Sequence, shape = Sequence)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_colour_manual(values = .seq_pal, name = NULL) +
    scale_fill_manual(values = .seq_pal, name = NULL) +
    scale_shape_manual(values = c("4-fold sites" = 16, "Introns, same genes" = 17), name = NULL) +
    scale_x_continuous(breaks = seq(0, 2.6, 0.4), labels = scales::number_format(accuracy = 0.1)) +
    labs(x = expression("Expression level category (log"[10] * ")"), y = ylab, title = title) +
    theme_custom() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
}

# ---------------------------------------------------------------------------
# The selection test: 4-fold pi/d over intronic pi/d of the same genes, by S
# ---------------------------------------------------------------------------

pi_d_double_ratio_by_S <- function(g4, gI, meta, U, V,
                                   S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, Inf),
                                   expr_range = c(-Inf, Inf), bin_width = 0.2,
                                   min_d_n = 20L, n_boot = 2000L, seed = 1L) {
  #' (pi4 / d4) / (piI / dI), C-vs-other system, per S_Wright bin, relative to
  #' S = 0 (interpolated between the bins flanking zero), genes resampled with
  #' both sequence classes together.
  #'
  #' Why this and not 4-fold pi / d alone. Linked selection lowers pi but not d,
  #' and it is stronger in high-S genes: intronic pi / d falls to ~0.92 at
  #' S ~ 1.2. Dividing by the intron of the same gene removes it, together with
  #' anything mutational. Mutation (rate or spectrum, shared or exon-specific)
  #' and linked selection all predict 1; Wright's selection model predicts
  #' pi(S)/pi(0) over K(S)/K(0).
  #' @param expr_range c(lo, hi) on the expression category: S is not monotone
  #'   with expression, so the test is repeated within expression strata.
  set.seed(seed)
  g <- merge(merge(data.table::as.data.table(g4), data.table::as.data.table(gI),
                   by = "Gene", suffixes = c("4", "I")),
             data.table::as.data.table(meta)[, .(Gene, S, Exp = Mean_Log10_Exp)], by = "Gene")
  g <- g[d_n4 >= min_d_n & d_nI >= min_d_n & is.finite(S) & is.finite(Exp)]
  g[, Exp_cat := round(round(Exp / bin_width) * bin_width, 1)]
  g <- g[Exp_cat >= expr_range[1] & Exp_cat <= expr_range[2]]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  rat <- function(i, cl) (sum(g[[paste0("pi_C", cl)]][i]) / sum(g[[paste0("pi_n", cl)]][i])) /
                         (sum(g[[paste0("d_C", cl)]][i]) / sum(g[[paste0("d_n", cl)]][i]))
  st <- function(i) { a <- rat(i, "4"); b <- rat(i, "I"); c(r4 = a, rI = b, dbl = a / b) }
  bins <- levels(droplevels(g$Sbin))
  idx  <- lapply(bins, function(b) which(g$Sbin == b)); names(idx) <- bins
  stopifnot(all(c("[-Inf,0)", "[0,0.25)") %in% bins))
  mS  <- vapply(idx, function(i) mean(g$S[i]), numeric(1))
  w0  <- -mS["[-Inf,0)"] / (mS["[0,0.25)"] - mS["[-Inf,0)"])
  est <- sapply(idx, st)
  ref <- (1 - w0) * est[, "[-Inf,0)"] + w0 * est[, "[0,0.25)"]
  bt  <- replicate(n_boot, sapply(idx, function(i) st(i[sample.int(length(i), length(i), TRUE)])))
  ref_b <- (1 - w0) * bt[, "[-Inf,0)", ] + w0 * bt[, "[0,0.25)", ]
  rel <- sweep(bt, c(1, 3), ref_b, "/")
  q  <- function(k, p) apply(rel[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  pv <- apply(rel["dbl", , , drop = FALSE], 2, function(z) min(1, 2 * min(mean(z <= 1), mean(z >= 1))))
  pred <- (wright_pi_2allele(mS, U, V) / wright_pi_2allele(0, U, V)) /
          (wright_substitution_rate(mS, U, V) / wright_substitution_rate(0, U, V))
  data.table::data.table(
    Sbin = bins, n_genes = lengths(idx), mean_S = mS,
    pd4_rel = est["r4", ] / ref["r4"], pdI_rel = est["rI", ] / ref["rI"],
    double = est["dbl", ] / ref["dbl"], lo = q("dbl", 0.025), hi = q("dbl", 0.975),
    P = pv, pred_pd = pred)
}

plot_pi_d_double_ratio <- function(tab, check = NULL, title = NULL) {
  #' The selection test: 4-fold over intronic pi / d by S_Wright bin, relative
  #' to S = 0, with gene-bootstrap intervals; filled = M. jungermannioides,
  #' open = M. glabratus (`check`). Solid line: Wright's selection prediction.
  #' Dashed line at 1: any mutational explanation, and linked selection.
  #' With a `stratum` column the panel is faceted by expression stratum.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::copy(data.table::as.data.table(tab))[, Outgroup := "JUNG1"]
  if (!is.null(check))
    d <- rbind(d, data.table::copy(data.table::as.data.table(check))[, Outgroup := "KONZA"], fill = TRUE)
  d[, x := mean_S + ifelse(Outgroup == "KONZA", 0.03, 0)]
  d[, Outgroup := factor(Outgroup, levels = c("JUNG1", "KONZA"))]
  keep <- intersect(c("stratum", "mean_S", "pred_pd"), names(tab))
  pred <- data.table::as.data.table(tab)[, ..keep]
  p <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey65", linewidth = 0.4) +
    geom_line(data = pred, aes(mean_S, pred_pd, linetype = "Selection (Wright)"),
              colour = "grey15", linewidth = 0.8) +
    geom_errorbar(data = d, aes(x = x, ymin = lo, ymax = hi), colour = "#3B7DD8",
                  width = 0.04, linewidth = 0.5) +
    geom_point(data = d, aes(x = x, y = double, shape = Outgroup), colour = "#3B7DD8",
               size = 2.4, stroke = 0.9) +
    scale_shape_manual(values = c(JUNG1 = 16, KONZA = 1), drop = TRUE,
                       labels = c(JUNG1 = expression(italic("M. jungermannioides")),
                                  KONZA = expression(italic("M. glabratus"))),
                       name = "Outgroup") +
    scale_linetype_manual(values = c("Selection (Wright)" = "solid"), name = NULL) +
    labs(x = expression(S[Wright] * " (bin mean)"),
         y = expression((pi / d)["4-fold"] / (pi / d)["intron"] * ", relative to " * S == 0),
         title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
  if ("stratum" %in% names(d) && data.table::uniqueN(d$stratum) > 1) {
    lv <- unique(d$stratum)
    p$layers <- lapply(p$layers, function(l) {
      if (is.data.frame(l$data) && "stratum" %in% names(l$data))
        l$data$stratum <- factor(l$data$stratum, levels = lv)
      l
    })
    p <- p + facet_wrap(~ stratum, nrow = 1)
  }
  p
}

# ---------------------------------------------------------------------------
# Divergence on its own: raw means by expression, and the turnover test by S
# ---------------------------------------------------------------------------
# Neutral divergence measures the mutation rate (plus a share of ancestral
# polymorphism, ~pi/d of it: ~0.17 at 4-fold sites for JUNG1, ~0.45 for KONZA).
# The neutral rate is not constant across genes: it falls ~20% with expression
# and rises ~30-50% with S. So the selection test needs a yardstick inside each
# gene. Under Wright's model, C/non-C turnover at 4-fold sites per opportunity
# is slower than the neutral rate by
#   f(-S) (V e^S + U) / (U + V)        (~0.70 at S = 1.2 for U/V = 5.7)
# because selection slows the loss of C more than it speeds its gain. Two
# yardsticks: changes among A, G and T at the same 4-fold sites (neutral under
# the C-preference model), and the introns of the same genes.

divergence_per_gene <- function(fourfold_sites, zerofold_div, intron_div, intron_comp, taxon,
                                pop_sites = NULL) {
  #' Per gene, over every site aligned to `taxon`: aligned sites and
  #' substitutions at 0-fold, 4-fold and intronic sites, plus the counts of the
  #' turnover test (C vs other; and, among sites non-C in both, the A/G/T
  #' composition of each side and the A/G/T changes).
  #'
  #' Two ways to count a 4-fold difference:
  #'   pop_sites = NULL  reference: IM767's base vs the outgroup's. Picks up
  #'                     polymorphism: 37% of JUNG1 and 72% of KONZA differences
  #'                     fall on sites segregating among the 187 lines.
  #'   pop_sites given   fixed differences: the outgroup base is absent from every
  #'                     IM line; IM composition is the population frequency. This
  #'                     is divergence that does not share sites with pi.
  #' @param fourfold_sites divergence_by_gene.py fourfold_sites.tsv.gz, as a table.
  #' @param zerofold_div,intron_div paths to its per-gene 0-fold and intronic
  #'   tables; zerofold_div may be absent (runs before 2026-11-17).
  #' @param intron_comp calculate_intronic_twostatepi.py output masked to
  #'   `taxon` (q_pref_C: intronic C over the same aligned sites).
  #' @param pop_sites extract_fourfold_sites() output (Gene, pos3, ref3, n, A..T).
  #' @return per-gene table; with pop_sites it also carries piC and piT (sums of
  #'   two-allele and total pi over the same sites), and the attribute "overlap".
  b  <- paste0(taxon, "_base")
  ff <- data.table::as.data.table(fourfold_sites)[, .(Gene, Pos, r = Ref, o = get(b))][
    r %chin% c("A", "C", "G", "T") & o %chin% c("A", "C", "G", "T")]
  if (is.null(pop_sites)) {
    g <- ff[, .(n4 = .N, d4 = sum(r != o), C4_ref = sum(r == "C"), C4_out = sum(o == "C"),
                d4C = sum((r == "C") != (o == "C")),
                nN = sum(r != "C" & o != "C"), dN = sum(r != "C" & o != "C" & r != o),
                rA = sum(r == "A" & o != "C"), rG = sum(r == "G" & o != "C"), rT = sum(r == "T" & o != "C"),
                oA = sum(o == "A" & r != "C"), oG = sum(o == "G" & r != "C"), oT = sum(o == "T" & r != "C")),
            by = Gene]
    ov <- NULL
  } else {
    x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T)],
               ff, by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r]
    x[, oc := data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T)]
    x[, `:=`(seg = (A > 0) + (C > 0) + (G > 0) + (T > 0) > 1L, noC = C == 0L,
             piT = n / (n - 1) * (1 - ((A / n)^2 + (C / n)^2 + (G / n)^2 + (T / n)^2)),
             piC = 2 * C * (n - C) / (n * (n - 1)))]
    ov <- x[, c(sites = .N, segregating = mean(seg), differ_from_ref = mean(r != o),
                fixed = mean(oc == 0L), ref_differences_on_segregating_sites = mean(seg[r != o]),
                expected_if_independent = mean(seg),
                share_of_pi_on_ref_differences = sum(piT[r != o]) / sum(piT))]
    g <- x[, .(n4 = .N, d4 = sum(oc == 0L), C4_ref = sum(C / n), C4_out = sum(o == "C"),
               d4C = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
               nN = sum(noC & o != "C"), dN = sum(noC & o != "C" & oc == 0L),
               rA = sum((A / n)[noC & o != "C"]), rG = sum((G / n)[noC & o != "C"]),
               rT = sum((T / n)[noC & o != "C"]),
               oA = sum(noC & o == "A"), oG = sum(noC & o == "G"), oT = sum(noC & o == "T"),
               piC = sum(piC), piT = sum(piT)), by = Gene]
    rm(x)
  }
  # Introns: fixed differences when the masked table carries them (intronic
  # script run with --outgroup) and fixed counting was asked for; otherwise the
  # reference-based counts of divergence_by_gene.py. QIo = outgroup intronic C.
  pm <- data.table::fread(intron_comp)
  has_fixed <- all(c("n_out", "out_C", "fixed_sub", "fixed_subC") %in% names(pm))
  if (!is.null(pop_sites) && has_fixed) {
    dI <- pm[, .(Gene, nI = n_out, dI = fixed_sub, dIC = fixed_subC, QI = q_pref_C,
                 QIo = out_C / n_out)]
  } else {
    if (!is.null(pop_sites))
      warning("[divergence] ", intron_comp, " has no fixed-difference columns; intronic d is ",
              "reference-based (rerun the intronic script with OUTGROUPS=...)", call. = FALSE)
    dI <- data.table::fread(intron_div)[, .(Gene, nI = get(paste0(taxon, "_n_aln")),
                                           dI = get(paste0(taxon, "_sub")),
                                           dIC = get(paste0(taxon, "_subC")))]
    dI <- merge(dI, pm[, .(Gene, QI = q_pref_C,
                           QIo = if ("out_C" %in% names(pm)) out_C / n_out else NA_real_)],
                by = "Gene", all.x = TRUE)
  }
  g <- merge(g, dI, by = "Gene", all = TRUE)
  if (!is.null(zerofold_div) && file.exists(zerofold_div)) {
    z <- data.table::fread(zerofold_div)[, .(Gene, n0 = get(paste0(taxon, "_n_aln")),
                                            d0 = get(paste0(taxon, "_sub")))]
    g <- merge(g, z, by = "Gene", all = TRUE)
  }
  if (!is.null(ov)) attr(g, "overlap") <- ov
  attr(g, "intron_counting") <- if (!is.null(pop_sites) && has_fixed) "fixed" else "reference"
  g[]
}

divergence_raw_by_expression <- function(gd, meta, bin_width = 0.2, bulk_max = 1.2,
                                         min_sites = 20L, min_genes = 40L) {
  #' Raw divergence (p, substitutions per aligned site) and its Jukes-Cantor
  #' correction at 0-fold, 4-fold and intronic sites, per expression category,
  #' on one gene set (every class aligned at >= min_sites), with each class
  #' relative to the drift-dominated bulk (categories <= bulk_max) and the
  #' 0-fold / 4-fold and 4-fold / intron ratios.
  g <- merge(data.table::as.data.table(gd), meta[, .(Gene, Mean_Log10_Exp)], by = "Gene")
  has0 <- "n0" %in% names(g)
  g <- g[is.finite(Mean_Log10_Exp) & n4 >= min_sites & nI >= min_sites &
           (if (has0) n0 >= min_sites else TRUE)]
  g[, Exp_cat := round(round(Mean_Log10_Exp / bin_width) * bin_width, 1)]
  jc <- function(p) -0.75 * log(1 - 4 * p / 3)
  cls <- c(`0-fold` = "0", `4-fold` = "4", intron = "I")
  if (!has0) cls <- cls[-1]
  summ <- function(x) {
    out <- list(genes = nrow(x))
    for (k in names(cls)) {
      p <- sum(x[[paste0("d", cls[k])]]) / sum(x[[paste0("n", cls[k])]])
      out[[paste0("p_", k)]] <- p
      out[[paste0("JC_", k)]] <- jc(p)
    }
    data.table::as.data.table(out)
  }
  tab <- rbind(g[, summ(.SD), by = Exp_cat][genes >= min_genes][order(Exp_cat)][
                 , Exp_cat := as.character(Exp_cat)],
               data.table::data.table(Exp_cat = sprintf("bulk (<= %.1f)", bulk_max),
                                      summ(g[Exp_cat <= bulk_max])),
               data.table::data.table(Exp_cat = "all genes", summ(g)))
  bulk <- tab[startsWith(Exp_cat, "bulk")]
  for (k in names(cls))
    tab[, (paste0("rel_", k)) := get(paste0("p_", k)) / bulk[[paste0("p_", k)]]]
  if (has0) tab[, p0_over_p4 := `p_0-fold` / `p_4-fold`]
  tab[, p4_over_pI := `p_4-fold` / p_intron]
  tab[]
}

divergence_turnover_by_S <- function(gd, meta, U, V,
                                     S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                     min_sites = 20L, n_boot = 2000L, seed = 1L) {
  #' The divergence-only selection test, per S_Wright bin, relative to S = 0
  #' (interpolated between the bins flanking zero), genes resampled.
  #' Each turnover is divergence over the mismatch expected between two
  #' sequences drawn independently with the observed compositions, so it does
  #' not move with composition:
  #'   C/non-C at 4-fold sites   d_C / [Q_ref (1 - Q_out) + Q_out (1 - Q_ref)]
  #'   A/G/T at 4-fold sites     (sites non-C in both) d / [1 - sum_b f_ref(b) f_out(b)]
  #'   introns, C/non-C          d_C / [Q_I (1 - Q_Io) + Q_Io (1 - Q_I)], with
  #'                             Q_Io the outgroup's intronic C when known, else Q_I
  #' ratio_nonC = C/non-C over A/G/T at the same 4-fold sites (all genes);
  #' ratio_intron = C/non-C at 4-fold sites over introns (genes with introns).
  #' Mutation, in any form shared with the yardstick: 1. Wright:
  #' f(-S) (V e^S + U) / (U + V). Raw C/non-C divergence relative to S = 0, and
  #' Wright's K(S) / K(0), are returned alongside.
  set.seed(seed)
  g <- merge(data.table::as.data.table(gd), meta[, .(Gene, S)], by = "Gene")
  g <- g[is.finite(S) & n4 >= min_sites]
  g[, hasI := !is.na(nI) & nI >= min_sites & is.finite(QI)]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  sm <- function(v, i) sum(g[[v]][i])
  st <- function(i) {
    Qr <- sm("C4_ref", i) / sm("n4", i); Qo <- sm("C4_out", i) / sm("n4", i)
    tC <- (sm("d4C", i) / sm("n4", i)) / (Qr * (1 - Qo) + Qo * (1 - Qr))
    fr <- c(sm("rA", i), sm("rG", i), sm("rT", i)) / sm("nN", i)
    fo <- c(sm("oA", i), sm("oG", i), sm("oT", i)) / sm("nN", i)
    tN <- (sm("dN", i) / sm("nN", i)) / (1 - sum(fr * fo))
    j  <- i[g$hasI[i]]
    QI <- sum(g$QI[j] * g$nI[j]) / sm("nI", j)
    QIo <- if (all(is.finite(g$QIo[j]))) sum(g$QIo[j] * g$nI[j]) / sm("nI", j) else QI
    tI <- (sm("dIC", j) / sm("nI", j)) / (QI * (1 - QIo) + QIo * (1 - QI))
    c(d4C = sm("d4C", i) / sm("n4", i), dIC = sm("dIC", j) / sm("nI", j),
      tC = tC, tN = tN, tI = tI, rN = tC / tN, rI = tC / tI)
  }
  bins <- levels(droplevels(g$Sbin))
  idx  <- lapply(bins, function(b) which(g$Sbin == b)); names(idx) <- bins
  stopifnot(all(c("[-Inf,0)", "[0,0.25)") %in% bins))
  mS  <- vapply(idx, function(i) mean(g$S[i]), numeric(1))
  w0  <- -mS["[-Inf,0)"] / (mS["[0,0.25)"] - mS["[-Inf,0)"])
  est <- sapply(idx, st)
  ref <- (1 - w0) * est[, "[-Inf,0)"] + w0 * est[, "[0,0.25)"]
  bt  <- replicate(n_boot, sapply(idx, function(i) st(i[sample.int(length(i), length(i), TRUE)])))
  rel <- sweep(bt, c(1, 3), (1 - w0) * bt[, "[-Inf,0)", ] + w0 * bt[, "[0,0.25)", ], "/")
  q  <- function(k, p) apply(rel[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  pv <- function(k) apply(rel[k, , , drop = FALSE], 2,
                          function(z) min(1, 2 * min(mean(z <= 1), mean(z >= 1))))
  f  <- function(x) ifelse(abs(x) < 1e-8, 1, x / (1 - exp(-x)))
  comp <- function(a, b) vapply(idx, function(i) sm(a, i) / sm(b, i), numeric(1))
  data.table::data.table(
    Sbin = bins, n_genes = lengths(idx), n_genes_intron = vapply(idx, function(i) sum(g$hasI[i]), 1L),
    mean_S = mS, Q4_ref = comp("C4_ref", "n4"), Q4_out = comp("C4_out", "n4"),
    G_share_nonC = comp("rG", "nN"),
    d4C_rel = est["d4C", ] / ref["d4C"], dIC_rel = est["dIC", ] / ref["dIC"],
    Wright_K = wright_substitution_rate(mS, U, V) / wright_substitution_rate(0, U, V),
    turnover_C = est["tC", ] / ref["tC"], turnover_nonC = est["tN", ] / ref["tN"],
    turnover_intron = est["tI", ] / ref["tI"],
    ratio_nonC = est["rN", ] / ref["rN"], lo_nonC = q("rN", 0.025), hi_nonC = q("rN", 0.975),
    P_nonC = pv("rN"),
    ratio_intron = est["rI", ] / ref["rI"], lo_intron = q("rI", 0.025), hi_intron = q("rI", 0.975),
    P_intron = pv("rI"),
    Wright_turnover = f(-mS) * (V * exp(mS) + U) / (U + V))
}

plot_divergence_turnover <- function(tab, check = NULL, yardstick = c("nonC", "intron"),
                                     title = NULL) {
  #' The divergence-only selection test (Fig 7D): C/non-C substitution rate at
  #' 4-fold sites over the neutral yardstick of the same genes, by S_Wright bin,
  #' relative to S = 0, with gene-bootstrap intervals; filled = M.
  #' jungermannioides, open = M. glabratus (`check`). Solid line: Wright's
  #' selection prediction. Dashed line at 1: any mutational explanation.
  #' @param yardstick "nonC" (A/G/T changes at the same 4-fold sites, Fig 7D) or
  #'   "intron" (introns of the same genes, supplement).
  yardstick <- match.arg(yardstick)
  suppressPackageStartupMessages(require(ggplot2))
  cols <- if (yardstick == "nonC") c("ratio_nonC", "lo_nonC", "hi_nonC") else
    c("ratio_intron", "lo_intron", "hi_intron")
  pick <- function(x, lab) {
    x <- data.table::as.data.table(x)
    data.table::data.table(Outgroup = lab, mean_S = x$mean_S, y = x[[cols[1]]],
                           lo = x[[cols[2]]], hi = x[[cols[3]]])
  }
  d <- pick(tab, "JUNG1")
  if (!is.null(check) && nrow(check)) d <- rbind(d, pick(check, "KONZA"))
  d[, x := mean_S + ifelse(Outgroup == "KONZA", 0.03, 0)]
  d[, Outgroup := factor(Outgroup, levels = c("JUNG1", "KONZA"))]
  pred <- data.table::as.data.table(tab)[, .(mean_S, Wright_turnover)]
  # K = substitution rate per opportunity (divergence over the mismatch expected
  # from composition alone); the legend text defines it.
  ylab <- if (yardstick == "nonC")
    expression((K["C/non-C"] / K["A/G/T"])["4-fold"] * ", relative to " * S == 0) else
    expression(K["C/non-C, 4-fold"] / K["C/non-C, intron"] * ", relative to " * S == 0)
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey65", linewidth = 0.4) +
    geom_line(data = pred, aes(mean_S, Wright_turnover, linetype = "Selection (Wright)"),
              colour = "grey15", linewidth = 0.8) +
    geom_errorbar(data = d, aes(x = x, ymin = lo, ymax = hi), colour = "#3B7DD8",
                  width = 0.04, linewidth = 0.5) +
    geom_point(data = d, aes(x = x, y = y, shape = Outgroup), colour = "#3B7DD8",
               size = 2.4, stroke = 0.9) +
    scale_shape_manual(values = c(JUNG1 = 16, KONZA = 1), drop = TRUE,
                       labels = c(JUNG1 = expression(italic("M. jungermannioides")),
                                  KONZA = expression(italic("M. glabratus"))),
                       name = "Outgroup") +
    scale_linetype_manual(values = c("Selection (Wright)" = "solid"), name = NULL) +
    labs(x = expression(S[Wright] * " (bin mean)"), y = ylab, title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
}

# ---------------------------------------------------------------------------
# Numerators and denominators: pi and fixed divergence by S, raw and per unit
# of neutral rate, against Wright (John Kelly's request, 2026-11-18)
# ---------------------------------------------------------------------------

.neutral_rate <- function(g, i) {
  # A/G/T turnover at the same 4-fold sites, per opportunity (see the turnover test)
  s  <- function(v) sum(g[[v]][i])
  fr <- c(s("rA"), s("rG"), s("rT")) / s("nN")
  fo <- c(s("oA"), s("oG"), s("oT")) / s("nN")
  (s("dN") / s("nN")) / (1 - sum(fr * fo))
}

pi_d_components_by_S <- function(gf, meta, U, V,
                                 S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                 min_sites = 20L, n_boot = 2000L, seed = 1L) {
  #' C-vs-other pi and fixed C/non-C divergence at 4-fold sites per S_Wright bin,
  #' relative to S = 0, each raw and divided by the neutral rate of the same genes
  #' (A/G/T changes at the same sites), with Wright's pi(S) and K(S).
  #' Absolute predictions are returned too: pi from U and V directly, and
  #' divergence as t K(S) with t set so that t K(0) equals the observed fixed
  #' divergence at S = 0 (the time cancels in every relative quantity).
  #' @param gf divergence_per_gene(..., pop_sites = ) output (fixed differences).
  set.seed(seed)
  g <- merge(data.table::as.data.table(gf), meta[, .(Gene, S)], by = "Gene")
  g <- g[is.finite(S) & n4 >= min_sites]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  st <- function(i) {
    piC <- sum(g$piC[i]) / sum(g$n4[i]); dC <- sum(g$d4C[i]) / sum(g$n4[i])
    nr  <- .neutral_rate(g, i)
    c(pi = piC, d = dC, neutral = nr, pi_sc = piC / nr, d_sc = dC / nr, pd = piC / dC)
  }
  bins <- levels(droplevels(g$Sbin))
  idx  <- lapply(bins, function(b) which(g$Sbin == b)); names(idx) <- bins
  mS   <- vapply(idx, function(i) mean(g$S[i]), numeric(1))
  w0   <- -mS["[-Inf,0)"] / (mS["[0,0.25)"] - mS["[-Inf,0)"])
  est  <- sapply(idx, st)
  ref  <- (1 - w0) * est[, "[-Inf,0)"] + w0 * est[, "[0,0.25)"]
  bt   <- replicate(n_boot, sapply(idx, function(i) st(i[sample.int(length(i), length(i), TRUE)])))
  rel  <- sweep(bt, c(1, 3), (1 - w0) * bt[, "[-Inf,0)", ] + w0 * bt[, "[0,0.25)", ], "/")
  q    <- function(k, p) apply(rel[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  t_sc <- ref["d"] / wright_substitution_rate(0, U, V)
  out <- data.table::data.table(Sbin = bins, n_genes = lengths(idx), mean_S = mS,
    pi = est["pi", ], pi_Wright = wright_pi_2allele(mS, U, V),
    d = est["d", ], d_Wright = t_sc * wright_substitution_rate(mS, U, V),
    Wright_pi_rel = wright_pi_2allele(mS, U, V) / wright_pi_2allele(0, U, V),
    Wright_d_rel = wright_substitution_rate(mS, U, V) / wright_substitution_rate(0, U, V))
  out[, Wright_pd_rel := Wright_pi_rel / Wright_d_rel]
  for (k in rownames(est)) {
    out[, (paste0(k, "_rel")) := est[k, ] / ref[k]]
    out[, (paste0(k, "_lo")) := q(k, 0.025)]
    out[, (paste0(k, "_hi")) := q(k, 0.975)]
  }
  attr(out, "time") <- unname(t_sc)
  out[]
}

plot_pi_d_components <- function(tab, title = NULL) {
  #' Two panels, pi and fixed divergence (C vs other, 4-fold sites), relative to
  #' S = 0: observed (filled), observed per unit of the neutral rate (open), the
  #' neutral rate itself (grey), and Wright's prediction (line).
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  mk <- function(q, lab) data.table::rbindlist(list(
    d[, .(panel = lab, series = "Observed", mean_S, y = get(paste0(q, "_rel")),
          lo = get(paste0(q, "_lo")), hi = get(paste0(q, "_hi")))],
    d[, .(panel = lab, series = "Observed / neutral rate", mean_S, y = get(paste0(q, "_sc_rel")),
          lo = get(paste0(q, "_sc_lo")), hi = get(paste0(q, "_sc_hi")))],
    d[, .(panel = lab, series = "Neutral rate (A/G/T changes)", mean_S, y = neutral_rel,
          lo = neutral_lo, hi = neutral_hi)]))
  labs_p <- c("pi[C]", "d[C]~'(fixed differences)'")
  pts <- rbind(mk("pi", labs_p[1]), mk("d", labs_p[2]))
  pts[, panel := factor(panel, levels = labs_p)]
  pts[, x := mean_S + data.table::fcase(series == "Observed", -0.025,
                                        series == "Observed / neutral rate", 0.025, default = 0)]
  wr <- rbind(d[, .(panel = labs_p[1], mean_S, y = Wright_pi_rel)],
              d[, .(panel = labs_p[2], mean_S, y = Wright_d_rel)])
  wr[, panel := factor(panel, levels = labs_p)]
  cols <- c("Observed" = "#3B7DD8", "Observed / neutral rate" = "#3B7DD8",
            "Neutral rate (A/G/T changes)" = "grey50")
  shp  <- c("Observed" = 16, "Observed / neutral rate" = 1, "Neutral rate (A/G/T changes)" = 17)
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_line(data = pts[series == "Neutral rate (A/G/T changes)"], aes(x, y),
              colour = "grey60", linewidth = 0.5) +
    geom_line(data = wr, aes(mean_S, y, linetype = "Selection (Wright)"),
              colour = "grey15", linewidth = 0.8) +
    geom_errorbar(data = pts, aes(x = x, ymin = lo, ymax = hi, colour = series),
                  width = 0.04, linewidth = 0.45) +
    geom_point(data = pts, aes(x = x, y = y, colour = series, shape = series),
               size = 2.2, stroke = 0.9) +
    facet_wrap(~ panel, nrow = 1, labeller = label_parsed) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_linetype_manual(values = c("Selection (Wright)" = "solid"), name = NULL) +
    labs(x = expression(S[Wright] * " (bin mean)"), y = expression("Relative to " * S == 0),
         title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"),
          legend.background = element_blank(), legend.key = element_blank())
}

# ---------------------------------------------------------------------------
# Referee 2's partition figure, per unit of neutral rate
# ---------------------------------------------------------------------------

partition_scaled_by_neutral <- function(partition_genes, gf, bulk_max = 1.2, min_genes = 40L,
                                        n_boot = 2000L, seed = 1L) {
  #' Both polymorphism classes (pi per 4-fold site; partition_genes from
  #' gene_level_partition_table()) per expression category relative to the bulk,
  #' raw and divided by the neutral rate of the same genes (A/G/T fixed
  #' differences at the same 4-fold sites, gf from divergence_per_gene(...,
  #' pop_sites = )). The neutral rate follows the same U shape as divergence in
  #' introns: it falls ~20% to categories 2.0-2.2 and rises again above, which
  #' is what lifts both classes at the top.
  set.seed(seed)
  g <- merge(data.table::as.data.table(partition_genes), data.table::as.data.table(gf),
             by = "Gene")[n4 >= 20]
  st <- function(i) {
    nr <- .neutral_rate(g, i)
    c(pref = sum(g$pi_sum_pref[i]) / sum(g$tot_sites[i]),
      ctrl = sum(g$pi_sum_ctrl[i]) / sum(g$tot_sites[i]), neutral = nr)
  }
  cats <- sort(unique(g$Exp_cat))
  idx  <- lapply(cats, function(cc) which(g$Exp_cat == cc))
  keep <- lengths(idx) >= min_genes; cats <- cats[keep]; idx <- idx[keep]
  bulk <- which(g$Exp_cat <= bulk_max)
  sc   <- function(v) c(v, pref_sc = v[["pref"]] / v[["neutral"]], ctrl_sc = v[["ctrl"]] / v[["neutral"]])
  est  <- sapply(idx, function(i) sc(st(i))); ref <- sc(st(bulk))
  bt   <- replicate(n_boot, {
    r <- sc(st(bulk[sample.int(length(bulk), length(bulk), TRUE)]))
    sapply(idx, function(i) sc(st(i[sample.int(length(i), length(i), TRUE)])) / r)
  })
  q <- function(k, p) apply(bt[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  out <- data.table::data.table(Exp_cat = cats, n_genes = lengths(idx))
  for (k in rownames(est)) {
    out[, (paste0(k, "_rel")) := est[k, ] / ref[k]]
    out[, (paste0(k, "_lo")) := q(k, 0.025)]
    out[, (paste0(k, "_hi")) := q(k, 0.975)]
  }
  out[]
}

plot_partition_scaled <- function(tab, title = NULL) {
  #' The two classes relative to the bulk: raw (thin dashed) and per unit of
  #' neutral rate (points, ribbons); grey, the neutral rate itself.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  lv <- c("preferred / unpreferred", "unpreferred / unpreferred")
  long <- rbind(
    d[, .(Exp_cat, Class = lv[1], raw = pref_rel, y = pref_sc_rel, lo = pref_sc_lo, hi = pref_sc_hi)],
    d[, .(Exp_cat, Class = lv[2], raw = ctrl_rel, y = ctrl_sc_rel, lo = ctrl_sc_lo, hi = ctrl_sc_hi)])
  long[, Class := factor(Class, levels = lv)]
  pal <- c("preferred / unpreferred" = "#B2182B", "unpreferred / unpreferred" = "#2166AC")
  ggplot(long, aes(Exp_cat, colour = Class, fill = Class)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_line(data = d, aes(Exp_cat, neutral_rel, linetype = "Neutral rate (A/G/T changes)"),
              inherit.aes = FALSE, colour = "grey50", linewidth = 0.7) +
    geom_line(aes(y = raw, linetype = "Observed"), linewidth = 0.45, alpha = 0.7) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = y, linetype = "Observed / neutral rate"), linewidth = 0.7) +
    geom_point(aes(y = y, shape = Class), size = 2.1) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_linetype_manual(values = c("Observed" = "22", "Observed / neutral rate" = "solid",
                                     "Neutral rate (A/G/T changes)" = "solid"), name = NULL) +
    guides(linetype = guide_legend(override.aes = list(colour = c("grey50", "grey20", "grey20")))) +
    scale_x_continuous(breaks = seq(0, max(d$Exp_cat), 0.4)) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "Relative to the bulk", title = title) +
    theme_custom() +
    theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
}

# ---------------------------------------------------------------------------
# Robustness of the divergence-only test (Fig 7D)
# ---------------------------------------------------------------------------
# Five versions of the same test, each closing one route to a false positive:
#   baseline          A/G/T yardstick, S from the whole gene, global U/V
#   A<->T yardstick   G may itself be mildly favoured, so drop it from the yardstick
#   intron yardstick  the gene's own introns, fixed differences
#   split half        S from half of each gene's 4-fold sites, everything measured
#                     on the other half (no shared sampling noise); Wright is read
#                     at the composition of the half measured, since S from half the
#                     sites is noisy and bins regress toward the mean
#   local bias        S from each gene's own intronic C content as its mutation bias
#                     instead of the global U/V (Referee 2's alternative)
# The mutation RATE barely enters S (exactly not in the weak-mutation limit; with
# the full Wright inversion a 0.8x-1.4x rate moves S by <= 0.03, 2x by 0.07); the
# mutation BIAS does (~0.19 per 20%), which is what "local bias" tests. Proofs:
# formal/wright_checks; code-vs-proof check: formal/numeric_crosscheck.R.
# Mutation bias alone, producing the same C content, predicts turnover
# (U e^-S + V) / (U + V) (less mutation away from C) to (U + V e^S) / (U + V)
# (more toward C), and pi/d = 1 for any mixture; Wright predicts f(-S)(V e^S + U)/(U + V).

.wright_Q_kummer <- function(S, U, V)
  vapply(S, function(s) (V / (U + V)) *
           exp(log_kummer_1F1(V + 1, U + V + 1, s) - log_kummer_1F1(V, U + V, s)), numeric(1))

divergence_turnover_robustness <- function(pop_sites, div_sites, meta, intron_comp, U, V, taxon,
                                           S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                           min_sites = 20L, shrink = 300, n_boot = 1000L,
                                           seed = 1L) {
  #' The turnover test (C/non-C over the neutral yardstick, per opportunity,
  #' fixed differences, relative to S = 0) in the five versions above.
  #' @param pop_sites extract_fourfold_sites() output.
  #' @param div_sites divergence_by_gene.py fourfold_sites.tsv.gz, as a table.
  #' @param meta Gene, S (S_Wright_signed).
  #' @param intron_comp the intronic table masked to `taxon` (with --outgroup
  #'   columns for the intron yardstick; q_pref_C for the local bias).
  #' @return data.table variant, Sbin, n_genes, S (where Wright is read), turnover,
  #'   lo, hi, Wright, and the two mutation-bias predictions.
  set.seed(seed)
  f    <- function(z) ifelse(abs(z) < 1e-8, 1, z / (1 - exp(-z)))
  W    <- function(S) f(-S) * (V * exp(S) + U) / (U + V)
  grid <- seq(-6, 9, by = 0.002); qg <- .wright_Q_kummer(grid, U, V)
  invS <- function(q) stats::approx(qg, grid, xout = q, ties = mean)$y
  b  <- paste0(taxon, "_base")
  x  <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T)],
              data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref, o = get(b))],
              by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[
                ref3 == r & o %chin% c("A", "C", "G", "T")]
  data.table::setorder(x, Gene, pos3)
  x[, half := ifelse(seq_len(.N) %% 2L == 1L, "A", "B"), by = Gene]
  x[, oc := data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T)]
  x[, `:=`(noC = C == 0L, W = C == 0L & G == 0L)]
  agg <- function(by) x[, .(
    n4 = .N, Cr = sum(r == "C"), QIM = sum(C / n), Qo = sum(o == "C"),
    dC = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
    nN = sum(noC & o != "C"), dN = sum(noC & o != "C" & oc == 0L),
    fA = sum((A / n)[noC & o != "C"]), fG = sum((G / n)[noC & o != "C"]),
    fT = sum((T / n)[noC & o != "C"]),
    oA = sum(noC & o == "A"), oG = sum(noC & o == "G"), oT = sum(noC & o == "T"),
    nW = sum(W & o %chin% c("A", "T")), dW = sum(W & o %chin% c("A", "T") & oc == 0L),
    wA = sum((A / n)[W & o %chin% c("A", "T")]), wT = sum((T / n)[W & o %chin% c("A", "T")]),
    owA = sum(W & o == "A"), owT = sum(W & o == "T")), by = by]
  stat <- function(g, i, yard) {
    s  <- function(v) sum(g[[v]][i])
    QI <- s("QIM") / s("n4"); Qo <- s("Qo") / s("n4")
    tC <- (s("dC") / s("n4")) / (QI * (1 - Qo) + Qo * (1 - QI))
    nr <- switch(yard,
      AGT = { fr <- c(s("fA"), s("fG"), s("fT")) / s("nN"); fo <- c(s("oA"), s("oG"), s("oT")) / s("nN")
              (s("dN") / s("nN")) / (1 - sum(fr * fo)) },
      AT  = { fr <- c(s("wA"), s("wT")) / s("nW"); fo <- c(s("owA"), s("owT")) / s("nW")
              (s("dW") / s("nW")) / (1 - sum(fr * fo)) },
      intron = { qi <- sum(g$QI[i] * g$nI[i]) / s("nI"); qo <- s("outCI") / s("nI")
                 (s("dIC") / s("nI")) / (qi * (1 - qo) + qo * (1 - qi)) })
    tC / nr
  }
  run <- function(g, Scol, yard, label, at_measured = FALSE) {
    g <- g[is.finite(get(Scol)) & n4 >= 10]
    g[, Sbin := cut(get(Scol), S_breaks, right = FALSE)]
    bins <- levels(droplevels(g$Sbin))
    idx  <- lapply(bins, function(bb) which(g$Sbin == bb)); names(idx) <- bins
    if (!all(c("[-Inf,0)", "[0,0.25)") %in% bins)) return(NULL)
    gi   <- split(seq_len(nrow(g)), g$Gene); genes <- names(gi)
    mS   <- vapply(idx, function(i) mean(g[[Scol]][i]), numeric(1))
    xS   <- if (at_measured) vapply(idx, function(i) invS(sum(g$Cr[i]) / sum(g$n4[i])), numeric(1)) else mS
    w0   <- -mS[1] / (mS[2] - mS[1])
    est  <- vapply(idx, function(i) stat(g, i, yard), numeric(1))
    rel  <- est / ((1 - w0) * est[1] + w0 * est[2])
    bt <- replicate(n_boot, {
      rows <- unlist(gi[sample(genes, length(genes), TRUE)], use.names = FALSE)
      bb <- g$Sbin[rows]
      e <- vapply(bins, function(k) stat(g, rows[bb == k], yard), numeric(1))
      e / ((1 - w0) * e[1] + w0 * e[2])
    })
    # Predictions are scaled to the same reference as the data: the value
    # interpolated at binning S = 0, read at the S where each bin truly sits
    # (with split-half S that reference sits near S = 0.2, not 0).
    sc <- function(v) v / ((1 - w0) * v[1] + w0 * v[2])
    data.table::data.table(variant = label, Sbin = bins, n_genes = lengths(idx), S = xS,
      turnover = rel, lo = apply(bt, 1, stats::quantile, 0.025, na.rm = TRUE),
      hi = apply(bt, 1, stats::quantile, 0.975, na.rm = TRUE), Wright = sc(W(xS)),
      bias_less_C_loss = sc((U * exp(-xS) + V) / (U + V)),
      bias_more_C_gain = sc((U + V * exp(xS)) / (U + V)))
  }
  m  <- data.table::as.data.table(meta)[, .(Gene, S)]
  gw <- merge(agg("Gene"), m, by = "Gene")[n4 >= min_sites]
  pm <- data.table::fread(intron_comp)
  out <- list(run(gw, "S", "AGT", "Baseline"), run(gw, "S", "AT", "A/T changes as yardstick"))
  if (all(c("n_out", "out_C", "fixed_subC") %in% names(pm))) {
    gI <- merge(gw, pm[n_out >= min_sites, .(Gene, nI = n_out, dIC = fixed_subC, QI = q_pref_C,
                                             outCI = out_C)], by = "Gene")
    out[[length(out) + 1]] <- run(gI, "S", "intron", "Introns as yardstick")
  } else {
    message("[robustness] ", intron_comp, ": no fixed-difference columns; intron yardstick skipped")
  }
  qh <- x[, .(Qh = mean(r == "C"), nh = .N), by = .(Gene, half)][nh >= min_sites %/% 2]
  qh[, `:=`(S_other = invS(Qh), half = ifelse(half == "A", "B", "A"))]
  gh <- merge(agg(c("Gene", "half")), qh[, .(Gene, half, S_other)], by = c("Gene", "half"))
  out[[length(out) + 1]] <- run(gh, "S_other", "AGT", "S from the other half of the sites",
                                at_measured = TRUE)
  ib <- pm[n_sites >= 200, .(Gene, nIs = n_sites, QIg = q_pref_C)]
  QIbar <- ib[, sum(QIg * nIs) / sum(nIs)]
  gl <- merge(gw, ib, by = "Gene")
  gl[, S_local := S - (stats::qlogis((QIg * nIs + shrink * QIbar) / (nIs + shrink)) -
                         stats::qlogis(V / (U + V)))]
  out[[length(out) + 1]] <- run(gl, "S_local", "AGT", "S with each gene's intronic mutation bias")
  res <- data.table::rbindlist(out)
  res[, variant := factor(variant, levels = unique(variant))]
  res[]
}

plot_turnover_robustness <- function(tab, title = NULL) {
  #' One panel per version of the test; filled = M. jungermannioides, open =
  #' M. glabratus; solid line Wright; dotted lines mutation bias alone (the two
  #' extremes); dashed line at 1 a change in mutation rate.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  d[, x := S + ifelse(Outgroup == "KONZA", 0.03, 0)]
  d[, Outgroup := factor(Outgroup, levels = c("JUNG1", "KONZA"))]
  ln <- d[Outgroup == "JUNG1"]
  alt <- data.table::rbindlist(list(
    ln[, .(variant, S, y = bias_less_C_loss, k = "a")], ln[, .(variant, S, y = bias_more_C_gain, k = "b")]))
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_line(data = alt, aes(S, y, group = k, linetype = "Mutation bias alone"),
              colour = "grey45", linewidth = 0.5) +
    geom_line(data = ln, aes(S, Wright, linetype = "Selection (Wright)"), colour = "grey15",
              linewidth = 0.8) +
    geom_errorbar(data = d, aes(x = x, ymin = lo, ymax = hi), colour = "#3B7DD8",
                  width = 0.04, linewidth = 0.45) +
    geom_point(data = d, aes(x = x, y = turnover, shape = Outgroup), colour = "#3B7DD8",
               size = 2, stroke = 0.85) +
    facet_wrap(~ variant, ncol = 3) +
    coord_cartesian(ylim = c(0.3, 1.25)) +
    scale_shape_manual(values = c(JUNG1 = 16, KONZA = 1), drop = TRUE,
                       labels = c(JUNG1 = expression(italic("M. jungermannioides")),
                                  KONZA = expression(italic("M. glabratus"))), name = "Outgroup") +
    scale_linetype_manual(values = c("Selection (Wright)" = "solid", "Mutation bias alone" = "13"),
                          name = NULL) +
    labs(x = expression(S[Wright]),
         y = expression((K["C/non-C"] / K["neutral"])["4-fold"] * ", relative to " * S == 0),
         title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
}

# ---------------------------------------------------------------------------
# The high-S divergence excess (S > ~1.5, ~5% of genes)
# ---------------------------------------------------------------------------
# Above S ~ 1.5, fixed divergence exceeds Wright's K(S) while pi matches pi(S).
# Two cheap explanations, each with a direct test:
#  (1) IM gained C recently (composition not at equilibrium on the IM lineage).
#      Then IM carries more C than its sister M. glabratus (KONZA) at the same
#      sites. Topology (((guttatus, nasutus), glabratus), jungermannioides).
#  (2) Ranking genes on IM's own composition picks genes whose S rose on the IM
#      side, which have extra recent C gains (regression to the mean across
#      lineages). Then ranking on the outgroup's composition, or on the mean of
#      both, removes the excess, and the IM-vs-outgroup asymmetry with it.
# Every ranking is split half (rank on one half of a gene's sites, measure on the
# other), so per-site sampling noise is not shared between ranking and measure.
# Wright is read at the composition actually measured (mean of the IM reference
# and the outgroup on the measured half), scaled to the same S = 0 reference.
# Result (2026-11-18, local data): IM - M. glabratus <= 0.006 at every S; ranking
# on either lineage or both leaves the top-bin turnover at 0.79-0.81 (Wright 0.65)
# while the asymmetry falls from 0.24 to 0.06. Neither explanation holds.

divergence_high_S_check <- function(pop_sites, div_sites, U, V, outgroup = "JUNG1",
                                    sister = "KONZA",
                                    S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                    min_sites = 20L, n_boot = 1000L, seed = 1L) {
  #' @param pop_sites extract_fourfold_sites() output.
  #' @param div_sites divergence_by_gene.py fourfold_sites.tsv.gz, as a table, with
  #'   `<outgroup>_base` and `<sister>_base`.
  #' @return data.table ranking, Sbin, n_genes, S (measured), turnover / d and pi
  #'   per neutral rate relative to S = 0 with CIs and Wright, asym (IM-side C
  #'   minus outgroup-side C among fixed C/non-C differences, over their sum), and
  #'   C-content differences IM - outgroup, sister - outgroup, IM - sister (sites
  #'   aligned in all three; IM = mean over the 187 lines) with CIs.
  set.seed(seed)
  grid <- seq(-6, 9, by = 0.002); qg <- .wright_Q_kummer(grid, U, V)
  invS <- function(q) stats::approx(qg, grid, xout = q, ties = mean)$y
  f    <- function(z) ifelse(abs(z) < 1e-8, 1, z / (1 - exp(-z)))
  Wt   <- function(S) f(-S) * (V * exp(S) + U) / (U + V)
  acgt <- c("A", "C", "G", "T")
  x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T)],
             data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref,
                                                      o = get(paste0(outgroup, "_base")),
                                                      s = get(paste0(sister, "_base")))],
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r & o %chin% acgt]
  data.table::setorder(x, Gene, pos3)
  x[, half := ifelse(seq_len(.N) %% 2L == 1L, "A", "B"), by = Gene]
  x[, `:=`(oc = data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T),
           noC = C == 0L, s3 = s %chin% acgt)]
  g <- x[, .(
    n4 = .N, Cr = sum(r == "C"), QIM = sum(C / n), Qo = sum(o == "C"),
    piC = sum(2 * C * (n - C) / (n * (n - 1))),
    dC = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
    aIM = sum(o != "C" & C == n), bO = sum(o == "C" & C == 0L),
    nN = sum(noC & o != "C"), dN = sum(noC & o != "C" & oc == 0L),
    fA = sum((A / n)[noC & o != "C"]), fG = sum((G / n)[noC & o != "C"]),
    fT = sum((T / n)[noC & o != "C"]),
    oA = sum(noC & o == "A"), oG = sum(noC & o == "G"), oT = sum(noC & o == "T"),
    n3 = sum(s3), q3IM = sum((C / n)[s3]), q3s = sum(s3 & s == "C"), q3o = sum(s3 & o == "C")),
    by = .(Gene, half)]
  # ranking scores from the OTHER half: IM reference, outgroup, and their mean
  # (single genomes each, so the three are on the same footing)
  oth <- g[, .(Gene, half = ifelse(half == "A", "B", "A"), m = n4,
               S_IM = invS(Cr / n4), S_out = invS(Qo / n4), S_both = invS((Cr + Qo) / (2 * n4)))]
  g <- merge(g, oth[m >= min_sites %/% 2], by = c("Gene", "half"))[n4 >= min_sites %/% 2]
  stat <- function(i) {
    s  <- function(v) sum(g[[v]][i])
    QI <- s("QIM") / s("n4"); Qo <- s("Qo") / s("n4")
    fr <- c(s("fA"), s("fG"), s("fT")) / s("nN"); fo <- c(s("oA"), s("oG"), s("oT")) / s("nN")
    nr <- (s("dN") / s("nN")) / (1 - sum(fr * fo))
    c(turn = (s("dC") / s("n4")) / (QI * (1 - Qo) + Qo * (1 - QI)) / nr,
      d_sc = (s("dC") / s("n4")) / nr, pi_sc = (s("piC") / s("n4")) / nr,
      asym = (s("aIM") - s("bO")) / (s("aIM") + s("bO")),
      IM_minus_out = (s("q3IM") - s("q3o")) / s("n3"),
      sister_minus_out = (s("q3s") - s("q3o")) / s("n3"),
      IM_minus_sister = (s("q3IM") - s("q3s")) / s("n3"))
  }
  rel_k <- c("turn", "d_sc", "pi_sc")
  genes <- unique(g$Gene); gi <- split(seq_len(nrow(g)), g$Gene)
  run <- function(Scol, label) {
    b <- cut(g[[Scol]], S_breaks, right = FALSE)
    bins <- levels(droplevels(b[!is.na(b)]))
    idx  <- lapply(bins, function(z) which(b == z)); names(idx) <- bins
    mS   <- vapply(idx, function(i) mean(g[[Scol]][i]), numeric(1))
    w0   <- -mS[1] / (mS[2] - mS[1])
    xS   <- vapply(idx, function(i) invS((sum(g$Cr[i]) + sum(g$Qo[i])) / (2 * sum(g$n4[i]))), numeric(1))
    sc   <- function(v) v / ((1 - w0) * v[1] + w0 * v[2])
    norm <- function(e) { e[rel_k, ] <- e[rel_k, ] / ((1 - w0) * e[rel_k, 1] + w0 * e[rel_k, 2]); e }
    est  <- norm(sapply(idx, stat))
    bt <- replicate(n_boot, {
      rows <- unlist(gi[sample(genes, length(genes), TRUE)], use.names = FALSE)
      bb <- b[rows]
      norm(sapply(bins, function(z) stat(rows[!is.na(bb) & bb == z])))
    })
    q <- function(k, p) apply(bt[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
    out <- data.table::data.table(ranking = label, Sbin = bins,
      n_genes = vapply(idx, function(i) data.table::uniqueN(g$Gene[i]), integer(1)), S = xS,
      Wright_turnover = sc(Wt(xS)), Wright_K = sc(wright_substitution_rate(xS, U, V)),
      Wright_pi = sc(wright_pi_2allele(xS, U, V)))
    for (k in rownames(est)) {
      out[, (k) := est[k, ]]
      out[, (paste0(k, "_lo")) := q(k, 0.025)]
      out[, (paste0(k, "_hi")) := q(k, 0.975)]
    }
    data.table::setnames(out, c("turn", "d_sc", "pi_sc"), c("turnover", "d_per_neutral", "pi_per_neutral"))
    data.table::setnames(out, function(z) sub("^turn_", "turnover_", sub("^d_sc_", "d_per_neutral_",
                                              sub("^pi_sc_", "pi_per_neutral_", z))))
    out
  }
  res <- data.table::rbindlist(list(run("S_IM", "IM"), run("S_out", outgroup),
                                    run("S_both", paste("IM and", outgroup))))
  res[, ranking := factor(ranking, levels = unique(ranking))]
  attr(res, "taxa") <- c(outgroup = outgroup, sister = sister)
  res[]
}

plot_high_S_check <- function(tab, title = NULL) {
  #' Two panels. composition: C content of IM (mean of 187 lines) and of the
  #' sister M. glabratus, each minus M. jungermannioides, at sites aligned in all
  #' three, genes ranked on IM (split half). ranking: Fig 7D's turnover with genes
  #' ranked on IM, on the outgroup, or on both; line = Wright at the measured S.
  suppressPackageStartupMessages(require(ggplot2))
  d  <- data.table::as.data.table(tab)
  pal <- RColorBrewer::brewer.pal(3, "Dark2")
  cm <- d[ranking == levels(ranking)[1]]
  cm <- data.table::rbindlist(list(
    cm[, .(S, who = "IM (187 lines)", y = IM_minus_out, lo = IM_minus_out_lo, hi = IM_minus_out_hi)],
    cm[, .(S = S + 0.03, who = "M. glabratus", y = sister_minus_out, lo = sister_minus_out_lo,
           hi = sister_minus_out_hi)]))
  p1 <- ggplot(cm, aes(S, y, colour = who)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.04, linewidth = 0.45) +
    geom_line(linewidth = 0.5) + geom_point(size = 2) +
    scale_colour_manual(values = c("IM (187 lines)" = pal[1], "M. glabratus" = pal[2]),
                        labels = c(expression("IM (187 lines)"), expression(italic("M. glabratus"))),
                        name = NULL) +
    labs(x = expression(S[Wright]),
         y = expression("C content minus " * italic("M. jungermannioides")), title = title) +
    theme_custom() + theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 0.5))
  d[, x := S + (as.integer(ranking) - 2) * 0.03]
  ln <- d[ranking == levels(ranking)[3]]
  p2 <- ggplot(d) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_line(data = ln, aes(S, Wright_turnover), colour = "grey15", linewidth = 0.8) +
    geom_errorbar(aes(x = x, ymin = turnover_lo, ymax = turnover_hi, colour = ranking),
                  width = 0.04, linewidth = 0.45) +
    geom_point(aes(x = x, y = turnover, colour = ranking), size = 2) +
    scale_colour_manual(values = stats::setNames(pal, levels(d$ranking)), name = "Genes ranked on",
                        labels = c(expression("IM"), expression(italic("M. jungermannioides")),
                                   expression("Both"))) +
    labs(x = expression(S[Wright]),
         y = expression((K["C/non-C"] / K["neutral"])["4-fold"] * ", relative to " * S == 0)) +
    theme_custom() + theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 0.5))
  list(composition = p1, ranking = p2)
}

# ---------------------------------------------------------------------------
# Fixed 0-fold differences (needs zerofold_sites.tsv.gz, divergence_by_gene.py
# from 2026-11-18 on)
# ---------------------------------------------------------------------------

.ZEROFOLD_OFFSETS <- local({
  # Same rule as divergence_by_gene.py and describe_gene_positions_by_degeneracy.py:
  # a codon position is 0-fold when every other base changes the amino acid
  # (stop counts as an amino acid). Value: the 0-fold positions (1-3) as a string.
  b   <- c("T", "C", "A", "G")
  cod <- with(expand.grid(k = b, j = b, i = b, stringsAsFactors = FALSE), paste0(i, j, k))
  aa  <- stats::setNames(strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]], cod)
  vapply(cod, function(cd) {
    ks <- Filter(function(k) {
      alt <- setdiff(c("A", "C", "G", "T"), substr(cd, k, k))
      all(aa[vapply(alt, function(x) { y <- cd; substr(y, k, k) <- x; y }, "")] != aa[[cd]])
    }, 1:3)
    paste(ks, collapse = "")
  }, "")
})

extract_zerofold_sites <- function(codon_freq_file, min_called = 50L, genes = NULL) {
  #' One row per 0-fold site with the 187 lines' base counts, from the codon file.
  #' @return data.table Gene, pos (1-based genomic coordinate, the key of
  #'   zerofold_sites.tsv.gz), ref (reference base, coding strand), A, C, G, T, n.
  d <- data.table::fread(codon_freq_file, select = c("Gene", "Ref_Codon", "Genomic_Positions",
                                                    "Codon_Variants"), showProgress = FALSE)
  if (!is.null(genes)) d <- d[Gene %chin% genes]
  d <- d[Ref_Codon %chin% names(.ZEROFOLD_OFFSETS)]
  d[, `:=`(row = .I, off = .ZEROFOLD_OFFSETS[Ref_Codon])]
  gp <- data.table::tstrsplit(d$Genomic_Positions, ",", fixed = TRUE, type.convert = TRUE)
  d[, Genomic_Positions := NULL]
  parts <- strsplit(d$Codon_Variants, ";", fixed = TRUE)
  d[, Codon_Variants := NULL]
  long <- data.table::data.table(row = rep(d$row, lengths(parts)), kv = unlist(parts, use.names = FALSE))
  rm(parts); gc()
  long[, `:=`(codon = substr(kv, 1L, 3L), cnt = as.integer(sub("^[^:]+:", "", kv)))]
  long[, kv := NULL]
  long <- long[codon != "NNN" & !is.na(cnt) & cnt > 0L]
  data.table::setkey(long, row)
  out <- vector("list", 3L)
  for (k in 1:3) {
    rk <- d$row[grepl(as.character(k), d$off, fixed = TRUE)]
    lk <- long[data.table::J(rk), nomatch = 0L]
    lk[, base := substr(codon, k, k)]
    lk <- lk[base %chin% c("A", "C", "G", "T"), .(cnt = sum(cnt)), by = .(row, base)]
    bc <- data.table::dcast(lk, row ~ base, value.var = "cnt", fill = 0L)
    for (b in c("A", "C", "G", "T")) if (!b %in% names(bc)) data.table::set(bc, j = b, value = 0L)
    out[[k]] <- bc[, .(Gene = d$Gene[row], pos = gp[[k]][row], ref = substr(d$Ref_Codon[row], k, k),
                       A, C, G, T)]
    rm(lk, bc)
  }
  rm(long, d, gp); gc()
  res <- data.table::rbindlist(out)
  res[, n := A + C + G + T]
  res[n >= min_called]
}

zerofold_divergence_by_gene <- function(zerofold_sites, zf_pop, taxon) {
  #' Per gene, over 0-fold sites aligned to `taxon` and called in >= min_called
  #' lines: n0 sites, d0_all (IM767 reference differs from the outgroup) and
  #' d0_fix (the outgroup base is carried by none of the lines). Same sites for both.
  b <- paste0(taxon, "_base")
  z <- data.table::as.data.table(zerofold_sites)[, .(Gene, Pos, r = Ref, o = get(b))][
    o %chin% c("A", "C", "G", "T") & r %chin% c("A", "C", "G", "T")]
  x <- merge(data.table::as.data.table(zf_pop), z, by.x = c("Gene", "pos"), by.y = c("Gene", "Pos"))
  ref_match <- x[, mean(ref == r)]
  x <- x[ref == r]
  x[, oc := data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T)]
  g <- x[, .(n0 = .N, d0_all = sum(r != o), d0_fix = sum(oc == 0L)), by = Gene]
  attr(g, "join") <- c(aligned = nrow(z), joined = nrow(x), ref_match = ref_match)
  g[]
}

# ---------------------------------------------------------------------------
# Raw divergence by expression category: all aligned sites vs fixed differences
# (John Kelly, 2026-11-18)
# ---------------------------------------------------------------------------

divergence_counting_by_expression <- function(pop_sites, div_sites, meta, taxon, intron_div,
                                              intron_comp, zerofold = NULL, zerofold_div = NULL,
                                              bin_width = 0.2, bulk_max = 1.2, min_sites = 20L,
                                              min_genes = 40L) {
  #' Raw divergence (differences per aligned site, no correction) per expression
  #' category, for one outgroup, counted two ways on the SAME sites:
  #'   "All aligned sites"  the IM767 reference differs from the outgroup
  #'   "Fixed differences"  the outgroup base is carried by none of the 187 lines
  #'                        (the counting every test in the paper uses)
  #' Classes: 0-fold, 4-fold (any change), 4-fold C/non-C (the selected system),
  #' 4-fold A/G/T (the neutral yardstick; sites non-C on both sides) and introns.
  #' One gene set for every class and both countings. Columns for classes whose
  #' input is missing are NA (0-fold fixed needs zerofold_sites.tsv.gz; intronic
  #' fixed needs the intronic script run with OUTGROUPS).
  acgt <- c("A", "C", "G", "T")
  x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T)],
             data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref, o = get(paste0(taxon, "_base")))],
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r & o %chin% acgt]
  x[, `:=`(oc = data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T), noC = C == 0L)]
  g <- x[, .(n4 = .N, d4_all = sum(r != o), d4_fix = sum(oc == 0L),
             d4C_all = sum((r == "C") != (o == "C")),
             d4C_fix = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
             nN_all = sum(r != "C" & o != "C"), dN_all = sum(r != "C" & o != "C" & r != o),
             nN_fix = sum(noC & o != "C"), dN_fix = sum(noC & o != "C" & oc == 0L)), by = Gene]
  rm(x)
  pm <- data.table::fread(intron_comp)
  if (all(c("n_out", "sub_ref", "fixed_sub") %in% names(pm))) {
    gI <- pm[, .(Gene, nI = n_out, dI_all = sub_ref, dI_fix = fixed_sub)]
  } else {
    gI <- data.table::fread(intron_div)[, .(Gene, nI = get(paste0(taxon, "_n_aln")),
                                           dI_all = get(paste0(taxon, "_sub")), dI_fix = NA_real_)]
  }
  g <- merge(g, gI, by = "Gene")
  has0 <- FALSE
  if (!is.null(zerofold)) {
    g <- merge(g, data.table::as.data.table(zerofold)[, .(Gene, n0, d0_all, d0_fix)], by = "Gene")
    has0 <- TRUE
  } else if (!is.null(zerofold_div) && file.exists(zerofold_div)) {
    z <- data.table::fread(zerofold_div)[, .(Gene, n0 = get(paste0(taxon, "_n_aln")),
                                            d0_all = get(paste0(taxon, "_sub")), d0_fix = NA_real_)]
    g <- merge(g, z, by = "Gene"); has0 <- TRUE
  }
  g <- merge(g, data.table::as.data.table(meta)[, .(Gene, Mean_Log10_Exp)], by = "Gene")
  g <- g[is.finite(Mean_Log10_Exp) & n4 >= min_sites & nI >= min_sites &
           (if (has0) n0 >= min_sites else TRUE)]
  g[, Exp_cat := round(round(Mean_Log10_Exp / bin_width) * bin_width, 1)]
  cls <- list(`0-fold` = c("n0", "d0"), `4-fold` = c("n4", "d4"), `4-fold C/non-C` = c("n4", "d4C"),
              `4-fold A/G/T` = c("nN", "dN"), intron = c("nI", "dI"))
  if (!has0) cls[["0-fold"]] <- NULL
  summ <- function(s, cnt) {
    out <- list(genes = nrow(s))
    for (k in names(cls)) {
      nn <- cls[[k]][1]; dd <- paste0(cls[[k]][2], "_", cnt)
      if (nn == "nN") nn <- paste0("nN_", cnt)
      out[[paste0("sites_", k)]] <- sum(s[[nn]])
      out[[paste0("p_", k)]] <- if (all(is.na(s[[dd]]))) NA_real_ else sum(s[[dd]]) / sum(s[[nn]])
    }
    data.table::as.data.table(out)
  }
  rows <- list()
  for (cnt in c("all", "fix")) {
    lab <- if (cnt == "all") "All aligned sites" else "Fixed differences"
    byc <- g[, summ(.SD, cnt), by = Exp_cat][genes >= min_genes][order(Exp_cat)]
    byc[, Exp_cat := sprintf("%.1f", Exp_cat)]
    rows[[length(rows) + 1]] <- data.table::data.table(Counting = lab, rbind(
      byc,
      data.table::data.table(Exp_cat = sprintf("bulk (<= %.1f)", bulk_max), summ(g[Exp_cat <= bulk_max], cnt)),
      data.table::data.table(Exp_cat = "all genes", summ(g, cnt))))
  }
  res <- data.table::data.table(Outgroup = taxon, data.table::rbindlist(rows))
  ord <- c(sprintf("%.1f", sort(unique(g$Exp_cat))), sprintf("bulk (<= %.1f)", bulk_max), "all genes")
  res[, Exp_cat := factor(Exp_cat, levels = ord)]
  data.table::setorder(res, Exp_cat, Counting)
  res[, Exp_cat := as.character(Exp_cat)]
  data.table::setcolorder(res, c("Outgroup", "Exp_cat", "Counting"))
  res[]
}

# ---------------------------------------------------------------------------
# Fig 7C/7D without the ratios: absolute pi and divergence with a time
# (John Kelly, 2026-11-18)
# ---------------------------------------------------------------------------

divergence_absolute_by_S <- function(gf, meta, U, V,
                                     S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                     min_sites = 20L, n_boot = 1000L, seed = 1L) {
  #' Absolute pi and fixed divergence per 4-fold site in the C/non-C system, by
  #' S_Wright bin, with Wright's predictions in absolute units.
  #'  pi:  Wright's two-allele heterozygosity pi(S) from U and V; nothing fitted.
  #'  d:   needs the divergence time t (4N generations, summed over both
  #'       lineages). Exact two-state solution: per C/non-C opportunity the chance
  #'       of a difference is 1 - exp(-lambda(S) t), lambda(S) = f(-S)(V e^S + U)
  #'       (= U + V at S = 0), so per site d(S) = 2Q(1-Q)(1 - exp(-lambda t)),
  #'       which saturates; t K(S) is its linear (small-t) form.
  #'  t:   (i) fitted at S = 0 from the 4-fold data (interpolated between the bins
  #'       flanking zero): one parameter, which the ratios cancel; (ii) from the
  #'       introns alone, when gf carries fixed intronic differences: introns are
  #'       neutral C/non-C sites under the same U and V, so this t is independent
  #'       of the 4-fold data and the prediction then has no free parameter.
  #' @param gf divergence_per_gene(..., pop_sites = ) output (fixed differences).
  #' @return per-bin table (observed with 95% CIs from resampling genes within
  #'   bins; predictions); attributes t (with CI) and t_intron.
  set.seed(seed)
  g <- merge(data.table::as.data.table(gf), data.table::as.data.table(meta)[, .(Gene, S)], by = "Gene")
  g <- g[is.finite(S) & n4 >= min_sites]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  f   <- function(z) ifelse(abs(z) < 1e-8, 1, z / (1 - exp(-z)))
  lam <- function(S) f(-S) * (V * exp(S) + U)
  st <- function(i) {
    s  <- function(v) sum(g[[v]][i])
    n  <- s("n4"); Q1 <- s("C4_ref") / n; Q2 <- s("C4_out") / n
    opp <- Q1 * (1 - Q2) + Q2 * (1 - Q1)
    fr <- c(s("rA"), s("rG"), s("rT")) / s("nN"); fo <- c(s("oA"), s("oG"), s("oT")) / s("nN")
    c(Q_IM = Q1, Q_out = Q2, pi = s("piC") / n, d = s("d4C") / n, opp = opp,
      d_per_opp = s("d4C") / n / opp, dAGT = s("dN") / s("nN"),
      dAGT_per_opp = (s("dN") / s("nN")) / (1 - sum(fr * fo)))
  }
  bins <- levels(droplevels(g$Sbin))
  idx  <- lapply(bins, function(b) which(g$Sbin == b)); names(idx) <- bins
  mS   <- vapply(idx, function(i) mean(g$S[i]), numeric(1))
  w0   <- -mS[1] / (mS[2] - mS[1])
  tfit <- function(e) -log(1 - ((1 - w0) * e["d_per_opp", 1] + w0 * e["d_per_opp", 2])) / (U + V)
  est  <- sapply(idx, st)
  t0   <- tfit(est)
  bt   <- replicate(n_boot, sapply(idx, function(i) st(i[sample.int(length(i), length(i), TRUE)])))
  tb   <- apply(bt, 3, tfit)
  q    <- function(k, p) apply(bt[k, , , drop = FALSE], 2, stats::quantile, p, na.rm = TRUE)
  QW   <- .wright_Q_kummer(mS, U, V)
  t_I  <- NA_real_
  if (isTRUE(attr(gf, "intron_counting") == "fixed")) {
    h <- g[is.finite(QI) & is.finite(QIo) & nI > 0]
    qi <- sum(h$QI * h$nI) / sum(h$nI); qo <- sum(h$QIo * h$nI) / sum(h$nI)
    t_I <- -log(1 - (sum(h$dIC) / sum(h$nI)) / (qi * (1 - qo) + qo * (1 - qi))) / (U + V)
  }
  pred <- function(t) 2 * QW * (1 - QW) * (1 - exp(-lam(mS) * t))
  out <- data.table::data.table(Sbin = bins, n_genes = lengths(idx),
    n_sites = vapply(idx, function(i) sum(g$n4[i]), numeric(1)), mean_S = mS,
    Q_IM = est["Q_IM", ], Q_out = est["Q_out", ], Q_Wright = QW,
    pi = est["pi", ], pi_lo = q("pi", 0.025), pi_hi = q("pi", 0.975),
    pi_Wright = wright_pi_2allele(mS, U, V),
    d = est["d", ], d_lo = q("d", 0.025), d_hi = q("d", 0.975),
    d_Wright = pred(t0), d_Wright_linear = t0 * wright_substitution_rate(mS, U, V),
    d_Wright_intron_t = if (is.finite(t_I)) pred(t_I) else NA_real_,
    opportunity = est["opp", ],
    d_per_opp = est["d_per_opp", ], d_per_opp_lo = q("d_per_opp", 0.025),
    d_per_opp_hi = q("d_per_opp", 0.975), d_per_opp_Wright = 1 - exp(-lam(mS) * t0),
    dAGT = est["dAGT", ], dAGT_per_opp = est["dAGT_per_opp", ],
    dAGT_per_opp_lo = q("dAGT_per_opp", 0.025), dAGT_per_opp_hi = q("dAGT_per_opp", 0.975),
    pi_over_d = est["pi", ] / est["d", ], pi_over_d_Wright = wright_pi_2allele(mS, U, V) / pred(t0))
  attr(out, "t") <- c(t = unname(t0), lo = unname(stats::quantile(tb, 0.025)),
                      hi = unname(stats::quantile(tb, 0.975)))
  attr(out, "t_intron") <- t_I
  out[]
}

plot_divergence_absolute <- function(tab, title = NULL) {
  #' Three panels, absolute units, by S_Wright: pi_C per site with Wright's pi(S)
  #' (no fitted parameter); fixed C/non-C divergence per site with Wright's
  #' saturating prediction for the fitted t (dashed: the linear t K(S)); and per
  #' opportunity, C/non-C next to the neutral A/G/T rate of the same sites.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  pal <- c(obs = "#3B7DD8", neu = "grey45")
  pan <- c("pi[C]~per~site", "d[C]~per~site~(fixed)", "Divergence~per~opportunity")
  pts <- data.table::rbindlist(list(
    d[, .(S = mean_S, y = pi, lo = pi_lo, hi = pi_hi, k = "C/non-C", p = pan[1])],
    d[, .(S = mean_S, y = d, lo = d_lo, hi = d_hi, k = "C/non-C", p = pan[2])],
    d[, .(S = mean_S, y = d_per_opp, lo = d_per_opp_lo, hi = d_per_opp_hi, k = "C/non-C", p = pan[3])],
    d[, .(S = mean_S + 0.03, y = dAGT_per_opp, lo = dAGT_per_opp_lo, hi = dAGT_per_opp_hi,
          k = "A/G/T (neutral)", p = pan[3])]))
  lns <- data.table::rbindlist(list(
    d[, .(S = mean_S, y = pi_Wright, l = "Wright", p = pan[1])],
    d[, .(S = mean_S, y = d_Wright, l = "Wright", p = pan[2])],
    d[, .(S = mean_S, y = d_Wright_linear, l = "Wright, linear t K(S)", p = pan[2])],
    d[, .(S = mean_S, y = d_per_opp_Wright, l = "Wright", p = pan[3])]))
  pts[, p := factor(p, levels = pan)]; lns[, p := factor(p, levels = pan)]
  ggplot() +
    geom_line(data = lns, aes(S, y, linetype = l), colour = "grey15", linewidth = 0.7) +
    geom_errorbar(data = pts, aes(S, ymin = lo, ymax = hi, colour = k), width = 0.04, linewidth = 0.45) +
    geom_point(data = pts, aes(S, y, colour = k, shape = k), size = 2) +
    facet_wrap(~ p, scales = "free_y", labeller = label_parsed) +
    scale_colour_manual(values = c(`C/non-C` = pal[["obs"]], `A/G/T (neutral)` = pal[["neu"]]), name = NULL) +
    scale_shape_manual(values = c(`C/non-C` = 16, `A/G/T (neutral)` = 17), name = NULL) +
    scale_linetype_manual(values = c(Wright = "solid", `Wright, linear t K(S)` = "22"), name = NULL) +
    labs(x = expression(S[Wright] ~ "(bin mean)"), y = NULL, title = title) +
    theme_custom() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 0.5))
}
