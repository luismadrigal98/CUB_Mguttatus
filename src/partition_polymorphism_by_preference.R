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


plot_polymorphism_partition <- function(density_table, title = NULL) {
  #' Diversity per 4-fold site in each preference class across the expression range
  suppressPackageStartupMessages({require(ggplot2); require(data.table)})
  d <- data.table::copy(data.table::as.data.table(density_table))
  d[, Class := factor(Class,
        levels = c("pref/non-pref", "non-pref/non-pref"),
        labels = c("preferred / unpreferred", "unpreferred / unpreferred"))]
  ggplot(d, aes(x = Exp_cat, y = pi_per_site, colour = Class, shape = Class)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c("preferred / unpreferred"   = "#B2182B",
                                   "unpreferred / unpreferred" = "#2166AC")) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = expression(pi ~ "per 4-fold site"),
         colour = NULL, shape = NULL, title = title) +
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
