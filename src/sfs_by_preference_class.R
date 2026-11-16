#' Site frequency spectrum of 4-fold polymorphisms, polarised by codon preference
#'
#' Referee 2 (2026-11 decision) asked whether elevated 4-fold diversity in highly
#' expressed genes could reflect a different mutation rate or spectrum rather than
#' selection, and proposed the SFS as the discriminating test.
#'
#' The test rests on one property of Wright's stationary distribution. For a
#' two-allele site, the expected number of segregating sites at which the
#' preferred allele is carried by k of n chromosomes is
#'
#'     C(n,k) B(k+V, n-k+U) / B(V,U)  ~=  [UV/(U+V)] * (1/k + 1/(n-k))
#'
#' when U and V are small (here U = 0.029, V = 0.005). Mutation sets the AMOUNT
#' of polymorphism; it does not set which allele is the rare one. Under any
#' mutation rate or bias at compositional equilibrium the preference-polarised
#' spectrum is symmetric, and so is it under any demography, because new
#' preferred and new unpreferred mutations enter at equal rates once composition
#' is at equilibrium. Selection (or gBGC) tilts it by exp(S p).
#'
#' Because the polarisation is by preference rather than by ancestry, no
#' outgroup is needed. Site classes, with what each was DESIGNED to isolate
#' under a single-preferred-base model:
#'   C/W  preferred vs weak         codon selection + gBGC
#'   C/G  preferred vs strong       codon selection, invisible to gBGC
#'   G/W  unpreferred strong/weak   gBGC only
#'   A/T  weak/weak                 neither
#' Two of those assignments failed on the data (2026-11-13): introns are
#' symmetric, so the G/W skew is not gBGC; and A/T trends with expression in
#' parallel with C/G, so it is not a flat control. Read the classes as base
#' pairs. The robust statements are (i) asymmetry at 4-fold sites and none in
#' introns, and (ii) the paired C/W - G/W contrast, which shares every
#' gene-level factor.
#'
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

FOURFOLD_PREFIX <- c(GC = "Ala", GG = "Gly", CC = "Pro", AC = "Thr",
                     GT = "Val", CT = "Leu_4", TC = "Ser_4", CG = "Arg_4")

# ---------------------------------------------------------------------------
# Wright two-allele SFS with binomial sampling
# ---------------------------------------------------------------------------

log_kummer_1F1 <- function(a, b, z, tol = 1e-15, max_terms = 5000L) {
  #' log 1F1(a; b; z), vectorised over `a`, scalar `b` and `z`.
  #'
  #' Direct series for z >= 0. For z < 0 the Kummer transformation
  #' 1F1(a;b;z) = e^z 1F1(b-a;b;-z) keeps every term positive, so there is no
  #' cancellation. Avoids a dependency on the gsl system library.
  if (z < 0) return(z + log_kummer_1F1(b - a, b, -z, tol, max_terms))
  term <- rep(1, length(a)); s <- term
  for (j in 0:(max_terms - 1L)) {
    term <- term * (a + j) / (b + j) * z / (j + 1)
    s <- s + term
    if (all(term <= tol * s)) break
  }
  log(s)
}

wright_sfs_logprob <- function(n, U, V, S) {
  #' log P(k | n) for k = 0..n: preferred-allele count in a sample of n under
  #' phi(p) ~ p^(V-1) (1-p)^(U-1) exp(S p), with binomial sampling.
  #' S = 4Ns, U = 4N(mutation away from preferred), V = 4N(mutation to preferred).
  k  <- 0:n
  lp <- lchoose(n, k) + lbeta(k + V, n - k + U) - lbeta(V, U) +
        log_kummer_1F1(k + V, n + U + V, S) - log_kummer_1F1(V, U + V, S)
  lp - (max(lp) + log(sum(exp(lp - max(lp)))))          # renormalise
}

wright_segregating_mean_freq <- function(n, U, V, S) {
  #' Expected preferred-allele frequency among segregating sites (0 < k < n).
  #' 0.5 exactly under neutrality in the small-theta limit, whatever U and V.
  p <- exp(wright_sfs_logprob(n, U, V, S))[2:n]
  sum(p * (1:(n - 1)) / n) / sum(p)
}

fit_wright_segregating_S <- function(sfs, U, V, interval = c(-15, 15)) {
  #' ML estimate of S from the SEGREGATING part of a preference-polarised SFS.
  #'
  #' Conditioning on 0 < k < n drops the monomorphic classes, which carry the
  #' information about U and V (composition), so S is estimated from shape
  #' alone. `sfs` is a vector of counts for k = 0..n.
  n   <- length(sfs) - 1L
  seg <- 2:n
  nll <- function(S) {
    lp <- wright_sfs_logprob(n, U, V, S)[seg]
    lp <- lp - (max(lp) + log(sum(exp(lp - max(lp)))))
    -sum(sfs[seg] * lp)
  }
  o <- stats::optimize(nll, interval)
  c(S = o$minimum, loglik = -o$objective, loglik0 = -nll(0))
}

# ---------------------------------------------------------------------------
# Per-site extraction
# ---------------------------------------------------------------------------

extract_fourfold_sites <- function(codon_freq_file, pref_base = "C",
                                   min_called = 50L) {
  #' One row per 4-fold site with base counts and local context.
  #'
  #' @param codon_freq_file `all_chromosomes.codon_frequencies.txt`.
  #' @param pref_base preferred 3rd-position base. All eight four-fold families
  #'   prefer C in this study (results/preferred_codon_contrast_vs_ROC.csv).
  #' @param min_called minimum called (non-NNN) chromosomes, as in the
  #'   polymorphism partition.
  #' @return data.table: Gene, Family, n, A, C, G, T, k_pref, n_alleles, Type,
  #'   next_base (first base of the following codon, CDS orientation), pos2, and
  #'   pos3 (1-based genomic coordinate of the 3rd codon position, the key that
  #'   joins these sites to the outgroup alignment in divergence_by_gene.py) and
  #'   ref3 (the reference base there, which the join checks).
  suppressPackageStartupMessages(require(data.table))
  d <- data.table::fread(codon_freq_file,
                         select = c("Gene", "Codon_Pos", "Ref_Codon", "Genomic_Positions",
                                    "Codon_Variants"),
                         showProgress = FALSE)
  data.table::setorder(d, Gene, Codon_Pos)
  # Context must be taken before the 4-fold filter removes the neighbours.
  d[, next_base := data.table::shift(substr(Ref_Codon, 1L, 1L), type = "lead"), by = Gene]
  d[, Family := FOURFOLD_PREFIX[substr(Ref_Codon, 1L, 2L)]]
  d <- d[!is.na(Family)]
  d[, `:=`(site = .I, pos2 = substr(Ref_Codon, 2L, 2L),
           pos3 = as.integer(sub("^.*,", "", Genomic_Positions)),
           ref3 = substr(Ref_Codon, 3L, 3L))]
  d[, Genomic_Positions := NULL]

  parts <- strsplit(d$Codon_Variants, ";", fixed = TRUE)
  long  <- data.table::data.table(site = rep(d$site, lengths(parts)),
                                  kv   = unlist(parts, use.names = FALSE))
  rm(parts); gc()
  long[, `:=`(codon = substr(kv, 1L, 3L),
              cnt   = as.integer(sub("^[^:]+:", "", kv)))]
  long <- long[codon != "NNN" & !is.na(cnt) & cnt > 0L]
  long[, base := substr(codon, 3L, 3L)]
  bc <- data.table::dcast(long[, .(cnt = sum(cnt)), by = .(site, base)],
                          site ~ base, value.var = "cnt", fill = 0L)
  rm(long); gc()
  for (b in c("A", "C", "G", "T"))
    if (!b %in% names(bc)) data.table::set(bc, j = b, value = 0L)

  out <- merge(d[, .(site, Gene, Family, next_base, pos2, pos3, ref3)],
               bc[, .(site, A, C, G, T)], by = "site")
  rm(d, bc); gc()
  out[, n := A + C + G + T]
  out <- out[n >= min_called]
  out[, n_alleles := (A > 0) + (C > 0) + (G > 0) + (T > 0)]
  out[, k_pref := get(pref_base)]
  # Polymorphism type for biallelic sites; everything else is mono or multi.
  out[, Type := data.table::fcase(
    n_alleles == 1L, "mono",
    n_alleles > 2L,  "multi",
    C > 0 & (A > 0 | T > 0), "C/W",
    C > 0 & G > 0,           "C/G",
    G > 0 & (A > 0 | T > 0), "G/W",
    A > 0 & T > 0,           "A/T")]
  out[, site := NULL]
  out[]
}

# ---------------------------------------------------------------------------
# Shape statistics by expression category, gene-clustered bootstrap
# ---------------------------------------------------------------------------

polarised_site_table <- function(sites, types = c("C/W", "C/G", "G/W", "A/T")) {
  #' Biallelic segregating sites with the focal allele's frequency.
  #' Focal allele: C for C/W and C/G (preferred), G for G/W (the strong allele,
  #' i.e. the one gBGC favours), A for A/T (arbitrary; neutral reference).
  s <- sites[Type %in% types]
  s[, k := data.table::fcase(Type %in% c("C/W", "C/G"), C,
                             Type == "G/W", G,
                             Type == "A/T", A)]
  s[, `:=`(f = k / n,
           # Tajima's estimators per site, variable n
           theta_pi = 2 * k * (n - k) / (n * (n - 1)),
           theta_w  = 1 / (digamma(n) - digamma(1)))]   # 1 / H_{n-1}
  s[]
}

gene_shape_table <- function(polarised, gene_meta, expression_var = "Mean_Log10_Exp",
                             bin_width = 0.2) {
  #' Per gene x type sums: n_seg, sum of focal frequency, theta_pi, theta_w.
  g <- polarised[, .(n_seg = .N, f_sum = sum(f),
                     theta_pi = sum(theta_pi), theta_w = sum(theta_w)),
                 by = .(Gene, Type)]
  meta <- data.table::as.data.table(gene_meta)[, .(Gene, Expr = get(expression_var))]
  g <- merge(g, meta, by = "Gene")
  g[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  g[]
}

bootstrap_shape_by_category <- function(gene_shape, n_boot = 2000L, seed = 1L,
                                        conf = 0.95) {
  #' Mean focal-allele frequency and theta_pi/theta_w per type x category,
  #' with intervals from resampling genes within each category.
  set.seed(seed)
  a <- (1 - conf) / 2
  gene_shape[, {
    fs <- f_sum; ns <- n_seg; tp <- theta_pi; tw <- theta_w; G <- .N
    bf <- numeric(n_boot); br <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      i <- sample.int(G, G, replace = TRUE)
      bf[b] <- sum(fs[i]) / sum(ns[i]); br[b] <- sum(tp[i]) / sum(tw[i])
    }
    .(n_genes = G, n_seg = sum(ns),
      mean_f = sum(fs) / sum(ns),
      mean_f_lo = stats::quantile(bf, a, na.rm = TRUE),
      mean_f_hi = stats::quantile(bf, 1 - a, na.rm = TRUE),
      pi_over_w = sum(tp) / sum(tw),
      pi_over_w_lo = stats::quantile(br, a, na.rm = TRUE),
      pi_over_w_hi = stats::quantile(br, 1 - a, na.rm = TRUE))
  }, by = .(Type, Exp_cat)][order(Type, Exp_cat)]
}

bootstrap_shape_regions <- function(gene_shape, regions, n_boot = 10000L,
                                    seed = 1L, conf = 0.95) {
  #' Pooled regions, and the change from the first region to each later one.
  #'
  #' The estimand mirrors the polymorphism partition: not the level of the
  #' asymmetry in any one region, but how it changes with expression, which a
  #' process indifferent to expression (gBGC, demography) cannot produce.
  #' @param regions named list of c(lo, hi) on Exp_cat, first = reference.
  set.seed(seed)
  a <- (1 - conf) / 2
  reg <- names(regions)
  gene_shape[, Region := NA_character_]
  for (r in reg) gene_shape[Exp_cat >= regions[[r]][1] & Exp_cat <= regions[[r]][2],
                            Region := r]
  gs <- gene_shape[!is.na(Region)]
  res <- list()
  for (ty in unique(gs$Type)) {
    x <- gs[Type == ty]
    idx <- lapply(reg, function(r) which(x$Region == r)); names(idx) <- reg
    stat <- function(ii) sum(x$f_sum[ii]) / sum(x$n_seg[ii])
    pt <- vapply(idx, stat, numeric(1))
    bt <- replicate(n_boot, vapply(idx, function(ii)
      stat(ii[sample.int(length(ii), length(ii), replace = TRUE)]), numeric(1)))
    for (j in seq_along(reg)) {
      d  <- bt[j, ] - bt[1, ]
      res[[length(res) + 1L]] <- data.table::data.table(
        Type = ty, Region = reg[j], n_genes = length(idx[[j]]),
        n_seg = sum(x$n_seg[idx[[j]]]),
        mean_f = pt[j],
        mean_f_lo = stats::quantile(bt[j, ], a), mean_f_hi = stats::quantile(bt[j, ], 1 - a),
        shift_vs_ref = pt[j] - pt[1],
        shift_lo = if (j > 1) stats::quantile(d, a) else NA_real_,
        shift_hi = if (j > 1) stats::quantile(d, 1 - a) else NA_real_,
        P_shift = if (j > 1) min(1, 2 * min(mean(d <= 0), mean(d >= 0))) else NA_real_)
    }
  }
  data.table::rbindlist(res)
}

# ---------------------------------------------------------------------------
# Projected preference-polarised SFS and shape-based S
# ---------------------------------------------------------------------------

project_gene_sfs <- function(sites, target_n = 100L, types = c("C/W", "C/G")) {
  #' Per-gene preference-polarised SFS projected to `target_n` chromosomes
  #' (hypergeometric), segregating sites of the given types plus monomorphic
  #' sites. Returns a gene x (k = 0..target_n) matrix.
  s <- sites[n >= target_n & (Type == "mono" | Type %in% types)]
  # Monomorphic sites project to the boundary without computation.
  s[, k_pref := data.table::fifelse(Type == "mono", data.table::fifelse(C > 0, n, 0L), C)]
  agg <- s[, .(count = .N), by = .(Gene, n, k_pref)]
  genes <- sort(unique(agg$Gene))
  M <- matrix(0, nrow = length(genes), ncol = target_n + 1L,
              dimnames = list(genes, 0:target_n))
  gi <- match(agg$Gene, genes)
  for (r in seq_len(nrow(agg))) {
    kk <- agg$k_pref[r]; nn <- agg$n[r]
    if (kk == 0L)      M[gi[r], 1L] <- M[gi[r], 1L] + agg$count[r]
    else if (kk == nn) M[gi[r], target_n + 1L] <- M[gi[r], target_n + 1L] + agg$count[r]
    else M[gi[r], ] <- M[gi[r], ] + agg$count[r] *
      stats::dhyper(0:target_n, kk, nn - kk, target_n)
  }
  M
}

fit_shape_S_by_region <- function(gene_sfs, gene_meta, regions, U, V,
                                  expression_var = "Mean_Log10_Exp",
                                  bin_width = 0.2, n_boot = 500L, seed = 1L) {
  #' Shape-based S per region, with a gene-clustered bootstrap interval.
  #' U and V only enter through the (negligible) small-theta correction; the
  #' sensitivity of S to them is reported separately.
  set.seed(seed)
  meta <- data.table::as.data.table(gene_meta)[, .(Gene, Expr = get(expression_var))]
  meta[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  ec <- meta$Exp_cat[match(rownames(gene_sfs), meta$Gene)]
  out <- list()
  for (r in names(regions)) {
    rows <- which(ec >= regions[[r]][1] & ec <= regions[[r]][2])
    if (!length(rows)) next
    sfs <- colSums(gene_sfs[rows, , drop = FALSE])
    f   <- fit_wright_segregating_S(sfs, U, V)
    bs  <- replicate(n_boot, fit_wright_segregating_S(
      colSums(gene_sfs[sample(rows, length(rows), replace = TRUE), , drop = FALSE]),
      U, V)[["S"]])
    out[[r]] <- data.table::data.table(
      Region = r, n_genes = length(rows),
      n_seg = sum(sfs[-c(1, length(sfs))]),
      S_shape = f[["S"]], S_lo = stats::quantile(bs, 0.025), S_hi = stats::quantile(bs, 0.975),
      LR_vs_neutral = 2 * (f[["loglik"]] - f[["loglik0"]]))
  }
  data.table::rbindlist(out)
}

bootstrap_type_contrast <- function(gene_shape, type_a, type_b, regions,
                                    n_boot = 10000L, seed = 1L, conf = 0.95) {
  #' Paired contrast of focal-allele frequency between two site types, by
  #' region, and its change from the reference (first) region.
  #'
  #' Genes are resampled jointly for both types, so a gene's C/W and G/W sites
  #' stay together and the shared gene-level factors (demography, linked
  #' selection, local gBGC) cancel in the difference.
  #' @param type_a,type_b site types, e.g. "C/W" and "G/W"; or type_b = NULL to
  #'   contrast type_a against the neutral value 0.5.
  set.seed(seed)
  a <- (1 - conf) / 2
  w <- data.table::dcast(gene_shape[Type %in% c(type_a, type_b)],
                         Gene + Exp_cat ~ Type, value.var = c("f_sum", "n_seg"), fill = 0)
  fa <- paste0("f_sum_", type_a); na <- paste0("n_seg_", type_a)
  stat <- if (is.null(type_b)) {
    function(ii) sum(w[[fa]][ii]) / sum(w[[na]][ii]) - 0.5
  } else {
    fb <- paste0("f_sum_", type_b); nb <- paste0("n_seg_", type_b)
    function(ii) sum(w[[fa]][ii]) / sum(w[[na]][ii]) - sum(w[[fb]][ii]) / sum(w[[nb]][ii])
  }
  idx <- lapply(regions, function(r) which(w$Exp_cat >= r[1] & w$Exp_cat <= r[2]))
  pt  <- vapply(idx, stat, numeric(1))
  bt  <- replicate(n_boot, vapply(idx, function(ii)
    stat(ii[sample.int(length(ii), length(ii), replace = TRUE)]), numeric(1)))
  data.table::rbindlist(lapply(seq_along(regions), function(j) {
    d <- bt[j, ] - bt[1, ]
    data.table::data.table(
      Contrast = if (is.null(type_b)) paste(type_a, "- 0.5") else paste(type_a, "-", type_b),
      Region = names(regions)[j], n_genes = length(idx[[j]]),
      estimate = pt[j], lo = stats::quantile(bt[j, ], a), hi = stats::quantile(bt[j, ], 1 - a),
      change_vs_ref = pt[j] - pt[1],
      change_lo = if (j > 1) stats::quantile(d, a) else NA_real_,
      change_hi = if (j > 1) stats::quantile(d, 1 - a) else NA_real_,
      P_change = if (j > 1) min(1, 2 * min(mean(d <= 0), mean(d >= 0))) else NA_real_)
  }))
}

# ---------------------------------------------------------------------------
# Plots (no titles: the legend text lives in the manuscript)
# ---------------------------------------------------------------------------

pool_shape_types <- function(gene_shape, types = c("C/W", "C/G"), label = "C/any") {
  #' Per-gene sums over several site types, as one new type. C/W + C/G pooled
  #' is the 4-fold counterpart of the intronic C vs non-C spectrum, so the two
  #' can be contrasted gene by gene with bootstrap_type_contrast().
  g <- gene_shape[Type %in% types,
                  .(n_seg = sum(n_seg), f_sum = sum(f_sum),
                    theta_pi = sum(theta_pi), theta_w = sum(theta_w)),
                  by = .(Gene, Expr, Exp_cat)]
  g[, Type := label]
  data.table::setcolorder(g, names(gene_shape))
  rbind(gene_shape, g)
}

plot_preference_asymmetry <- function(category_table, intron_table = NULL,
                                      min_genes = 40L, neutral = 0.5, title = NULL) {
  #' Focal-allele frequency among segregating sites, by site type and
  #' expression category, with gene-bootstrap bands.
  #' @param category_table output of `bootstrap_shape_by_category()`.
  #' @param intron_table the same for the introns of the same genes
  #'   (intronic_gene_shape() then bootstrap_shape_by_category()); drawn in grey
  #'   as the no-selection reference.
  #' @param min_genes categories with fewer genes for a type are dropped.
  suppressPackageStartupMessages(require(ggplot2))
  # Labelled by base pair, not by the process each class was meant to isolate:
  # introns show no strong/weak skew, so G/W does not measure gBGC, and A/T is
  # not flat, so "neither" would overstate the control.
  lab <- c("C/G" = "C vs G (focal C)", "C/W" = "C vs A/T (focal C)",
           "G/W" = "G vs A/T (focal G)", "A/T" = "A vs T (focal A)")
  intron_lab <- "Introns, C vs other (focal C)"
  pal <- c(stats::setNames(RColorBrewer::brewer.pal(4, "Dark2"), lab), grey35 = "grey35")
  names(pal)[5] <- intron_lab
  d <- data.table::as.data.table(category_table)[Type %in% names(lab)]
  d[, Type := lab[Type]]
  if (!is.null(intron_table)) {
    it <- data.table::as.data.table(intron_table)
    it[, Type := intron_lab]
    d <- rbind(d, it, fill = TRUE)
  }
  d <- d[n_genes >= min_genes]
  d[, Type := factor(Type, levels = intersect(names(pal), unique(Type)))]
  ggplot(d, aes(Exp_cat, mean_f, colour = Type, fill = Type)) +
    geom_hline(yintercept = neutral, linetype = "dashed", colour = "grey40") +
    geom_ribbon(aes(ymin = mean_f_lo, ymax = mean_f_hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "Focal-allele frequency\namong segregating sites", title = title) +
    theme_custom() + theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 3), fill = guide_legend(nrow = 3))
}

composition_exon_vs_intron <- function(gene_q, intron_comp, bin_width = 0.2,
                                       min_genes = 40L, n_boot = 2000L, seed = 1L) {
  #' C content at 4-fold sites (population Q) and in the introns of the same
  #' genes (reference sequence), per expression category, relative to the
  #' drift-dominated bulk (category <= 1.2). Gene-bootstrap intervals.
  #' @param gene_q data.table: Gene, Q4, n4, Exp.
  #' @param intron_comp data.table: Gene, C, n (intronic base counts).
  set.seed(seed)
  d <- merge(gene_q, intron_comp[, .(Gene, C, n)], by = "Gene")
  d[, Exp_cat := round(round(Exp / bin_width) * bin_width, 1)]
  b4 <- d[Exp_cat <= 1.2, sum(Q4 * n4) / sum(n4)]
  bI <- d[Exp_cat <= 1.2, sum(C) / sum(n)]
  out <- d[, {
    f <- function(i) c(sum(Q4[i] * n4[i]) / sum(n4[i]) / b4, sum(C[i]) / sum(n[i]) / bI)
    pt <- f(seq_len(.N)); bs <- replicate(n_boot, f(sample.int(.N, .N, TRUE)))
    .(Sequence = c("4-fold sites", "Introns, same genes"), rel = pt,
      lo = apply(bs, 1, stats::quantile, 0.025), hi = apply(bs, 1, stats::quantile, 0.975),
      n_genes = .N)
  }, by = Exp_cat]
  out[n_genes >= min_genes][order(Sequence, Exp_cat)]
}

plot_composition_exon_vs_intron <- function(tab, title = NULL) {
  #' C content of 4-fold sites and of the introns of the same genes, each
  #' relative to its own drift-dominated bulk, with gene-bootstrap bands.
  #' 4-fold sites in the blue of Fig 7A; introns in the grey used for the intronic
  #' reference elsewhere (plot_preference_asymmetry()).
  suppressPackageStartupMessages(require(ggplot2))
  pal <- c("4-fold sites" = "#3B7DD8", "Introns, same genes" = "grey35")
  ggplot(tab, aes(Exp_cat, rel, colour = Sequence, fill = Sequence, shape = Sequence)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_shape_manual(values = c("4-fold sites" = 16, "Introns, same genes" = 17), name = NULL) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "C content relative to the bulk", title = title) +
    theme_custom() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_blank())
}
