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

fourfold_pi_d_by_gene <- function(fourfold_sites, div_sites, taxon) {
  #' Per gene: pi and d at the 4-fold sites aligned to `taxon`.
  #'
  #' @param fourfold_sites extract_fourfold_sites() output (needs pos3, ref3).
  #' @param div_sites path or table: divergence_by_gene.py fourfold_sites.tsv.gz.
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
    pi_C   = 2 * C * (n - C) / (n * (n - 1)),
    d_tot  = as.numeric(Ref != Out),
    d_C    = as.numeric((Ref == "C") != (Out == "C")))]
  g <- a[, .(pi_n = .N, pi_C = sum(pi_C), pi_tot = sum(pi_tot),
             d_n = .N, d_C = sum(d_C), d_tot = sum(d_tot)), by = Gene]
  attr(g, "join") <- c(codon_sites = nrow(fourfold_sites), joined = nrow(x),
                       ref_match = ref_match, aligned = nrow(a))
  g[]
}

intron_pi_d_by_gene <- function(intron_div, intron_masked, taxon) {
  #' Per gene: intronic d over the sites aligned to `taxon` (divergence_by_gene.py)
  #' and intronic pi over the same sites (calculate_intronic_twostatepi.py run with
  #' --mask <taxon>=intron_aligned_<taxon>.bed). Same columns as the 4-fold table.
  dv <- data.table::fread(intron_div)
  dv <- dv[, .(Gene, d_n = get(paste0(taxon, "_n_aln")), d_tot = get(paste0(taxon, "_sub")),
               d_C = get(paste0(taxon, "_subC")))]
  pm <- data.table::fread(intron_masked)
  if (!"pi_total" %in% names(pm))
    stop("[intron pi/d] ", intron_masked, " lacks pi_total: rerun calculate_intronic_twostatepi.py")
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
