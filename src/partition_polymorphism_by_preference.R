## ============================================================================
## partition_polymorphism_by_preference.R
##
## Referee 1's line-535 test.
##
##   "The real test would be to partition polymorphic sites by those with
##    preferred/non-preferred polymorphisms vs those with non-preferred/
##    non-preferred polymorphisms. The prediction is that elevated polymorphism
##    would be present only for sites in the first scenario. If both scenarios
##    show elevated polymorphism (i.e. due to skew toward intermediate
##    frequencies), then it would point to some other explanation."
##
## The logic: codon-usage selection can only act on a segregating site if the
## alleles differ in preference. Where both segregating alleles are
## unpreferred, selection under a single-preferred-base model is indifferent
## between them, so that class is an internal negative control. If the
## diversity elevation at high expression is caused by selection opposing
## mutation, it must be confined to the preferred/non-preferred class.
## ============================================================================

classify_fourfold_polymorphism <- function(codon_freq_file,
                                           preferred_codons,
                                           genetic_code,
                                           min_called = 50L) {
  #' Partition segregating 4-fold sites by whether the preferred base is present
  #'
  #' @param codon_freq_file per-site table with Gene, AA, Ref_Codon,
  #'   Codon_Variants and Frequencies (`all_chromosomes.codon_frequencies.txt`).
  #' @param preferred_codons data.frame with Family and Preferred_Codon, as
  #'   returned by `detect_preferred_codons()$preferred`.
  #' @param genetic_code named vector codon -> family label.
  #' @param min_called minimum called (non-NNN) chromosomes for a site to count.
  #' @return data.table, one row per segregating 4-fold site: Gene, Family,
  #'   Class ("pref/non-pref" or "non-pref/non-pref"), pi_site, n_called.

  suppressPackageStartupMessages(require(data.table))

  fourfold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Leu_4", "Ser_4", "Arg_4")
  pref <- data.table::as.data.table(preferred_codons)[Family %in% fourfold]
  pref_base <- stats::setNames(substr(pref$Preferred_Codon, 3, 3), pref$Family)

  d <- data.table::fread(codon_freq_file,
                         select = c("Gene", "AA", "Ref_Codon", "Codon_Variants"),
                         showProgress = FALSE)
  # Family from the codon itself (handles the split six-fold families).
  d[, Family := unname(genetic_code[Ref_Codon])]
  d <- d[Family %in% fourfold]
  # Denominator for the density measure: EVERY 4-fold site of the gene,
  # monomorphic included. Must be taken before the segregating filter below.
  gene_totals <- d[, .(tot_sites = .N), by = Gene]
  d <- d[grepl(";", Codon_Variants)]

  # Parse "COD:count;COD:count", dropping the NNN missing-data entry.
  parts <- strsplit(d$Codon_Variants, ";", fixed = TRUE)
  res <- vapply(seq_along(parts), function(i) {
    kv <- parts[[i]]
    cod <- substr(kv, 1L, 3L)
    cnt <- as.numeric(sub("^[^:]+:", "", kv))
    keep <- cod != "NNN" & !is.na(cnt) & cnt > 0
    cod <- cod[keep]; cnt <- cnt[keep]
    if (!length(cnt)) return(c(NA_real_, NA_real_, NA_real_))
    n <- sum(cnt)
    b <- substr(cod, 3L, 3L)                      # 3rd position is the degenerate one
    agg <- tapply(cnt, b, sum)
    f <- agg / n
    # unbiased per-site heterozygosity
    pi_site <- if (n > 1) (n / (n - 1)) * (1 - sum(f^2)) else NA_real_
    has_pref <- as.numeric(pref_base[[d$Family[i]]] %in% names(agg))
    c(pi_site, n, has_pref)
  }, numeric(3))

  out <- data.table::data.table(
    Gene     = d$Gene,
    Family   = d$Family,
    pi_site  = res[1, ],
    n_called = res[2, ],
    has_pref = res[3, ]
  )
  out <- out[!is.na(pi_site) & n_called >= min_called & pi_site > 0]
  out[, Class := data.table::fifelse(has_pref == 1,
                                     "pref/non-pref", "non-pref/non-pref")]
  out[, has_pref := NULL]

  # Attached rather than returned separately so the two can never drift apart.
  data.table::setattr(out, "gene_totals", gene_totals)
  out[]
}


summarise_polymorphism_partition_density <- function(site_table, gene_meta,
                                                     expression_var = "Mean_Log10_Exp",
                                                     bin_width = 0.2,
                                                     min_genes = 30L) {
  #' pi CONTRIBUTED PER 4-FOLD SITE in each preference class
  #'
  #' The quantity that carries the signal. Heterozygosity conditional on a site
  #' already being polymorphic is flat-to-declining in both classes and would
  #' read as a null; the elevation at high expression lives in the DENSITY of
  #' segregating sites. Both classes are divided by the same denominator (all
  #' 4-fold sites of the genes in the bin), so the two are directly comparable.

  suppressPackageStartupMessages(require(data.table))
  totals <- attr(site_table, "gene_totals")
  if (is.null(totals)) stop("site_table carries no 'gene_totals' attribute")

  meta <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(meta, expression_var, "Expr")
  tot <- data.table::copy(data.table::as.data.table(totals))
  tot[, Gene_name := paste0("MgIM767.", Gene)]
  tot <- merge(tot, meta, by = "Gene_name")[!is.na(Expr)]
  tot[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  denom <- tot[, .(tot_sites = sum(tot_sites), n_genes = .N), by = Exp_cat]

  st <- merge(data.table::as.data.table(site_table),
              tot[, .(Gene, Exp_cat)], by = "Gene")
  num <- st[, .(n_seg = .N, pi_sum = sum(pi_site)), by = .(Exp_cat, Class)]

  agg <- merge(num, denom, by = "Exp_cat")[n_genes >= min_genes]
  agg[, pi_per_site := pi_sum / tot_sites]
  data.table::setorder(agg, Class, Exp_cat)
  agg[]
}



bootstrap_partition_density_bands <- function(gene_table,
                                              n_boot    = 2000L,
                                              conf      = 0.95,
                                              min_genes = 30L,
                                              seed      = 1L) {
  #' Per-category confidence bands for the two preference classes
  #'
  #' Genes are resampled with replacement WITHIN each expression category, for
  #' the same reason the region contrasts are bootstrapped that way: linked
  #' sites inside a gene are not independent, so a binomial or Poisson interval
  #' on site counts would be far too narrow. Categories holding fewer than
  #' `min_genes` genes are dropped, matching the density table.
  #'
  #' @return data.table Exp_cat, Class, pi_per_site, ci_low, ci_high, n_genes

  suppressPackageStartupMessages(require(data.table))
  g <- data.table::as.data.table(gene_table)[, .(Gene, Exp_cat, tot_sites,
                                                 pi_sum_pref, pi_sum_ctrl)]
  a <- (1 - conf) / 2
  set.seed(seed)

  out <- lapply(split(g, by = "Exp_cat", keep.by = TRUE), function(d) {
    n <- nrow(d)
    if (n < min_genes) return(NULL)
    tot <- d$tot_sites; pp <- d$pi_sum_pref; pc <- d$pi_sum_ctrl
    bs <- vapply(seq_len(n_boot), function(i) {
      k <- sample.int(n, n, replace = TRUE)
      s <- sum(tot[k])
      c(sum(pp[k]) / s, sum(pc[k]) / s)
    }, numeric(2))
    # The class contrast is taken WITHIN each replicate, so it keeps the
    # correlation between the two classes: they are measured on the same genes
    # and share a denominator, and the paired difference is far better resolved
    # than the overlap of the two marginal bands suggests.
    ct <- bs[1, ] / bs[2, ] - 1
    list(
      classes = data.table::data.table(
        Exp_cat     = d$Exp_cat[1],
        Class       = c("pref/non-pref", "non-pref/non-pref"),
        pi_per_site = c(sum(pp) / sum(tot), sum(pc) / sum(tot)),
        ci_low      = c(stats::quantile(bs[1, ], a), stats::quantile(bs[2, ], a)),
        ci_high     = c(stats::quantile(bs[1, ], 1 - a), stats::quantile(bs[2, ], 1 - a)),
        n_genes     = n),
      contrast = data.table::data.table(
        Exp_cat  = d$Exp_cat[1],
        contrast = sum(pp) / sum(pc) - 1,
        ci_low   = stats::quantile(ct, a),
        ci_high  = stats::quantile(ct, 1 - a),
        n_genes  = n))
  })
  out <- Filter(Negate(is.null), out)
  res <- data.table::rbindlist(lapply(out, `[[`, "classes"))
  con <- data.table::rbindlist(lapply(out, `[[`, "contrast"))
  data.table::setorder(res, Class, Exp_cat)
  data.table::setorder(con, Exp_cat)
  data.table::setattr(res, "contrast", con[])
  res[]
}



pool_partition_contrast <- function(gene_table, regions,
                                    n_boot = 50000L, conf = 0.95, seed = 1L) {
  #' Class contrast pooled over named regions of the expression range
  #'
  #' Pooling IS the multiplicity correction. Testing all fourteen expression
  #' categories separately answers a question nobody asked and spends the power
  #' on bins of 83-184 genes; the analysis is specified on the two regimes the
  #' paper is about, so the regions are the unit of inference. Per-category
  #' p-values are reported alongside (Benjamini-Hochberg) purely to show what a
  #' bin-by-bin reading would and would not support.
  #'
  #' @param regions named list of c(min_category, max_category); Inf allowed.
  #' @return data.table region, lo_cat, hi_cat, n_genes, estimate, ci_low,
  #'   ci_high, p_two_sided

  suppressPackageStartupMessages(require(data.table))
  g <- data.table::as.data.table(gene_table)
  a <- (1 - conf) / 2
  set.seed(seed)

  out <- lapply(names(regions), function(nm) {
    r <- regions[[nm]]
    d <- g[Exp_cat >= r[1] & Exp_cat <= r[2]]
    if (!nrow(d)) return(NULL)
    pp <- d$pi_sum_pref; pc <- d$pi_sum_ctrl; n <- nrow(d)
    bs <- vapply(seq_len(n_boot), function(i) {
      k <- sample.int(n, n, replace = TRUE)
      sum(pp[k]) / sum(pc[k]) - 1
    }, numeric(1))
    data.table::data.table(
      region = nm, lo_cat = r[1], hi_cat = r[2], n_genes = n,
      estimate = sum(pp) / sum(pc) - 1,
      ci_low = stats::quantile(bs, a), ci_high = stats::quantile(bs, 1 - a),
      p_two_sided = min(1, 2 * min(mean(bs <= 0), mean(bs >= 0))))
  })
  data.table::rbindlist(out)[]
}


plot_partition_contrast <- function(contrast_table, pooled = NULL, title = NULL) {
  #' Class contrast per expression category, with paired bootstrap band
  #'
  #' The quantity the analysis actually tests. Zero is the line at which the two
  #' classes carry equal diversity per site; the reversal is the crossing.
  #'
  #' `pooled`, as returned by `pool_partition_contrast()`, overlays the regions
  #' that are actually tested: a horizontal segment at the pooled estimate with
  #' its interval as a box. Single categories in the tail hold 83-184 genes and
  #' are individually underpowered, so drawing the pooled regions keeps the eye
  #' on the inference rather than on the noisiest points.

  suppressPackageStartupMessages({require(ggplot2); require(data.table)})
  d <- data.table::as.data.table(contrast_table)
  p <- ggplot(d, aes(x = Exp_cat, y = 100 * contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
               colour = "grey40")

  if (!is.null(pooled) && nrow(pooled)) {
    pl <- data.table::copy(data.table::as.data.table(pooled))
    rng <- range(d$Exp_cat)
    pl[, `:=`(x1 = pmax(lo_cat, rng[1]), x2 = pmin(hi_cat, rng[2]))]
    p <- p +
      geom_rect(data = pl, inherit.aes = FALSE,
                aes(xmin = x1, xmax = x2, ymin = 100 * ci_low, ymax = 100 * ci_high),
                fill = "grey35", alpha = 0.16) +
      geom_segment(data = pl, inherit.aes = FALSE,
                   aes(x = x1, xend = x2, y = 100 * estimate, yend = 100 * estimate),
                   colour = "grey20", linewidth = 0.7)
  }

  p +
    geom_ribbon(aes(ymin = 100 * ci_low, ymax = 100 * ci_high),
                fill = "#B2182B", alpha = 0.18) +
    geom_line(linewidth = 0.6, colour = "#B2182B") +
    geom_point(size = 2.2, colour = "#B2182B") +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "Excess diversity of the\npreferred / unpreferred class (%)",
         title = title) +
    theme_custom()
}



diagnose_partition_overlap <- function(gene_table, n_boot = 4000L,
                                       conf = 0.95, min_genes = 30L, seed = 1L) {
  #' Why the two bands in panel A overlap while the contrast in panel B does not
  #'
  #' Supplementary Figure 8's caption asserts that the class estimates are
  #' strongly correlated and that the paired difference therefore carries about
  #' half the sampling variance the marginal bands imply. This computes those
  #' numbers so the claim is reproducible rather than asserted:
  #'
  #'   Var(A - B) = Var(A) + Var(B) - 2 Cov(A, B)
  #'
  #' `width_ratio` is sd(paired difference) / sd(difference under independence).
  #' `bands_overlap` says whether the two marginal intervals touch, and
  #' `contrast_excludes_zero` whether the paired contrast does not - the rows
  #' where these disagree are the point.
  #'
  #' @return data.table Exp_cat, n_genes, r, bands_overlap, width_ratio,
  #'   contrast, ci_low, ci_high, contrast_excludes_zero

  suppressPackageStartupMessages(require(data.table))
  g <- data.table::as.data.table(gene_table)
  a <- (1 - conf) / 2
  set.seed(seed)

  res <- data.table::rbindlist(lapply(split(g, by = "Exp_cat"), function(d) {
    n <- nrow(d)
    if (n < min_genes) return(NULL)
    tot <- d$tot_sites; pp <- d$pi_sum_pref; pc <- d$pi_sum_ctrl
    bs <- vapply(seq_len(n_boot), function(i) {
      k <- sample.int(n, n, replace = TRUE); s <- sum(tot[k])
      c(sum(pp[k]) / s, sum(pc[k]) / s)
    }, numeric(2))
    ct <- bs[1, ] / bs[2, ] - 1
    qa <- stats::quantile(bs[1, ], c(a, 1 - a))
    qb <- stats::quantile(bs[2, ], c(a, 1 - a))
    lo <- stats::quantile(ct, a); hi <- stats::quantile(ct, 1 - a)
    data.table::data.table(
      Exp_cat = d$Exp_cat[1], n_genes = n,
      r = stats::cor(bs[1, ], bs[2, ]),
      bands_overlap = !(qa[2] < qb[1] || qb[2] < qa[1]),
      width_ratio = stats::sd(bs[1, ] - bs[2, ]) /
                    sqrt(stats::var(bs[1, ]) + stats::var(bs[2, ])),
      contrast = sum(pp) / sum(pc) - 1, ci_low = lo, ci_high = hi,
      contrast_excludes_zero = lo > 0 | hi < 0)
  }))
  data.table::setorder(res, Exp_cat)
  res[]
}


plot_polymorphism_partition <- function(density_table, title = NULL) {
  #' Diversity per 4-fold site in each preference class across the expression range
  #'
  #' If `density_table` carries `ci_low`/`ci_high` — as returned by
  #' `bootstrap_partition_density_bands()` — the intervals are drawn as ribbons.
  #' Without them the function falls back to lines and points, so the older
  #' density table still plots.

  suppressPackageStartupMessages({require(ggplot2); require(data.table)})
  d <- data.table::copy(data.table::as.data.table(density_table))
  d[, Class := factor(Class,
        levels = c("pref/non-pref", "non-pref/non-pref"),
        labels = c("preferred / unpreferred", "unpreferred / unpreferred"))]
  pal <- c("preferred / unpreferred"   = "#B2182B",
           "unpreferred / unpreferred" = "#2166AC")

  p <- ggplot(d, aes(x = Exp_cat, y = pi_per_site,
                     colour = Class, fill = Class, shape = Class))
  if (all(c("ci_low", "ci_high") %in% names(d)))
    p <- p + geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
                         alpha = 0.18, colour = NA)
  p +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = expression(pi ~ "per 4-fold site"),
         colour = NULL, fill = NULL, shape = NULL, title = title) +
    theme_custom() +
    theme(legend.position = "top")
}


summarise_polymorphism_partition <- function(site_table, gene_meta,
                                             expression_var = "Mean_Log10_Exp",
                                             bin_width = 0.2,
                                             min_sites = 200L) {
  #' Mean per-site diversity in each class, across the expression range
  #'
  #' The discriminating comparison: if codon-usage selection drives the
  #' diversity elevation at high expression, it appears in the
  #' preferred/non-preferred class only.

  suppressPackageStartupMessages(require(data.table))
  meta <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(meta, expression_var, "Expr")
  # site table carries the bare gene id; metadata carries the MgIM767 prefix
  st <- data.table::copy(data.table::as.data.table(site_table))
  st[, Gene_name := paste0("MgIM767.", Gene)]
  st <- merge(st, meta, by = "Gene_name")[!is.na(Expr)]

  st[, Exp_cat := round(Expr / bin_width) * bin_width]
  agg <- st[, .(n_sites = .N, pi_mean = mean(pi_site),
                pi_se = stats::sd(pi_site) / sqrt(.N)),
            by = .(Class, Exp_cat)][n_sites >= min_sites]
  data.table::setorder(agg, Class, Exp_cat)
  agg[]
}


gene_level_partition_table <- function(site_table, gene_meta,
                                       expression_var = "Mean_Log10_Exp",
                                       bin_width = 0.2) {
  #' One row per gene: pi contributed by each preference class, and the shared
  #' denominator.
  #'
  #' The unit of resampling for every inference below. Segregating sites within
  #' a gene are linked and cannot be treated as independent draws; the gene is
  #' the smallest unit that is approximately exchangeable.

  suppressPackageStartupMessages(require(data.table))
  totals <- attr(site_table, "gene_totals")
  if (is.null(totals)) stop("site_table carries no 'gene_totals' attribute")

  meta <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(meta, expression_var, "Expr")
  tot <- data.table::copy(data.table::as.data.table(totals))
  tot[, Gene_name := paste0("MgIM767.", Gene)]
  tot <- merge(tot, meta, by = "Gene_name")[!is.na(Expr)]
  tot[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]

  st <- data.table::as.data.table(site_table)[, .(pi_sum = sum(pi_site), n_seg = .N),
                                              by = .(Gene, Class)]
  st[, Cls := data.table::fifelse(Class == "pref/non-pref", "pref", "ctrl")]
  w <- data.table::dcast(st, Gene ~ Cls,
                         value.var = c("pi_sum", "n_seg"), fill = 0)
  # dcast omits a column entirely if a class is absent from the whole table.
  for (j in c("pi_sum_pref", "pi_sum_ctrl", "n_seg_pref", "n_seg_ctrl"))
    if (!j %in% names(w)) data.table::set(w, j = j, value = 0)

  g <- merge(tot[, .(Gene, Exp_cat, tot_sites)], w, by = "Gene", all.x = TRUE)
  for (j in setdiff(names(g), c("Gene", "Exp_cat", "tot_sites")))
    data.table::set(g, which(is.na(g[[j]])), j, 0)
  g[]
}


bootstrap_partition_contrast <- function(gene_table,
                                         bulk_max  = 1.2,
                                         elite_min = 2.4,
                                         elite_max = Inf,
                                         n_boot    = 10000L,
                                         conf      = 0.95,
                                         seed      = 1L) {
  #' Gene-clustered bootstrap of the class contrast, and of the reversal itself
  #'
  #' Three things the chi-square on the top category could not do:
  #'
  #'  1. It tests the INTERACTION.  The claim is not that the selected class is
  #'     more diverse in the elite genes, it is that the two classes CHANGE
  #'     ORDER.  The statistic is therefore `reversal = contrast(elite) -
  #'     contrast(bulk)`, and the CI on that is the test.
  #'  2. It resamples GENES, not sites.  Linked sites within a gene are not
  #'     independent, so a site-level test is anticonservative.
  #'  3. It puts an interval on the elite contrast, which rests on few genes.
  #'
  #' Because both classes share the same per-gene denominator, the denominator
  #' cancels in the class ratio and the contrast is simply
  #' `sum(pi_pref) / sum(pi_ctrl) - 1` over the genes of a region.
  #'
  #' `elite_max` caps the top region. The figure plots only categories holding
  #' at least 30 genes, so capping at the last plotted category makes the
  #' interval describe exactly what the reader can see; leaving it at Inf uses
  #' the whole tail, which is the more powerful test.
  #'
  #' @return list(estimate, ci, boot) — `estimate` a named numeric of the three
  #'   statistics, `ci` a data.frame with percentile limits and a two-sided
  #'   bootstrap p-value for the reversal.

  suppressPackageStartupMessages(require(data.table))
  g <- data.table::as.data.table(gene_table)
  bulk  <- g[Exp_cat <= bulk_max]
  elite <- g[Exp_cat >= elite_min & Exp_cat <= elite_max]
  if (!nrow(bulk) || !nrow(elite))
    stop("bulk or elite region is empty; check bulk_max / elite_min")

  contrast <- function(d) sum(d$pi_sum_pref) / sum(d$pi_sum_ctrl) - 1

  est <- c(bulk     = contrast(bulk),
           elite    = contrast(elite),
           reversal = contrast(elite) - contrast(bulk))

  set.seed(seed)
  nb <- nrow(bulk); ne <- nrow(elite)
  bp <- bulk$pi_sum_pref;  bc <- bulk$pi_sum_ctrl
  ep <- elite$pi_sum_pref; ec <- elite$pi_sum_ctrl

  boot <- vapply(seq_len(n_boot), function(i) {
    ib <- sample.int(nb, nb, replace = TRUE)
    ie <- sample.int(ne, ne, replace = TRUE)
    cb <- sum(bp[ib]) / sum(bc[ib]) - 1
    ce <- sum(ep[ie]) / sum(ec[ie]) - 1
    c(cb, ce, ce - cb)
  }, numeric(3))
  rownames(boot) <- c("bulk", "elite", "reversal")

  a  <- (1 - conf) / 2
  ci <- data.frame(
    statistic = rownames(boot),
    estimate  = est,
    ci_low    = apply(boot, 1, stats::quantile, probs = a),
    ci_high   = apply(boot, 1, stats::quantile, probs = 1 - a),
    row.names = NULL
  )
  # Two-sided bootstrap p: how often does a replicate fall on the null side of 0?
  ci$p_two_sided <- vapply(rownames(boot), function(s) {
    b <- boot[s, ]
    min(1, 2 * min(mean(b <= 0), mean(b >= 0)))
  }, numeric(1))
  ci$n_genes <- c(nb, ne, NA_integer_)

  list(estimate = est, ci = ci, boot = boot,
       regions = list(bulk_max = bulk_max, elite_min = elite_min,
                      elite_max = elite_max, n_bulk = nb, n_elite = ne,
                      n_boot = n_boot))
}
