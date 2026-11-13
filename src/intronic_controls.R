#' Introns of the same genes as the no-selection control for 4-fold sites
#'
#' An intron shares its gene's chromatin, transcription, recombination
#' environment and linked selection, but carries no codons. Any explanation of
#' the 4-fold pattern that runs through those shared factors (a mutation rate or
#' spectrum keyed to open chromatin, background selection, demography) predicts
#' the same pattern in the introns. Codon selection predicts it at 4-fold sites
#' only. Three tests use that contrast:
#'
#'   * composition: C content of 4-fold sites vs the introns of the same genes,
#'     from the reference sequence (no population data needed);
#'   * SFS shape: the C-polarised spectrum of segregating intronic sites, by
#'     expression region, against the 4-fold one (Referee 2);
#'   * diversity: 4-fold two-allele pi over intronic two-allele pi, by S_Wright,
#'     against Wright's prediction (Referee 1's counterfactual).
#'
#' Population data come from miscellanea_code/calculate_intronic_twostatepi.py,
#' which must be run on the ALL-SITES VCF. Its intron definition (primary
#' transcript, exonic features merged, 30 bp trimmed from each end, trimmed
#' width >= 86 bp) is reproduced exactly by intronic_reference_composition().
#'
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

# ---------------------------------------------------------------------------
# Reference-sequence composition of introns
# ---------------------------------------------------------------------------

intronic_reference_composition <- function(gff_file, genome_fasta,
                                           trim_bp = 30L, min_width = 86L) {
  #' Base counts of each gene's trimmed introns, in transcript orientation.
  #'
  #' @return data.table: Gene, A, C, G, T, n, fC, fG, fGC
  suppressPackageStartupMessages({ require(data.table); require(Biostrings) })

  gff <- data.table::fread(cmd = sprintf("grep -v '^#' %s", shQuote(gff_file)),
                           header = FALSE, sep = "\t",
                           col.names = c("chr", "src", "type", "start", "end",
                                         "score", "strand", "phase", "attr"))
  tx <- gff[type == "mRNA"]
  tx[, `:=`(tx_id   = sub(".*ID=([^;]+).*", "\\1", attr),
            Gene    = sub("^MgIM767\\.", "",
                          sub("\\.v[0-9]+\\.[0-9]+$", "",
                              sub(".*Parent=([^;]+).*", "\\1", attr))),
            primary = grepl("longest=1", attr, fixed = TRUE))]
  if (any(tx$primary)) tx <- tx[primary == TRUE]

  ex <- gff[type %in% c("exon", "CDS", "five_prime_UTR", "three_prime_UTR")]
  ex[, tx_id := sub(".*Parent=([^;,]+).*", "\\1", attr)]
  ex <- merge(ex[, .(tx_id, chr, start, end, strand)],
              tx[, .(tx_id, Gene)], by = "tx_id")
  rm(gff); invisible(gc())

  # CDS and UTR pieces of one exon abut: merge into exons, then take the gaps.
  data.table::setorder(ex, tx_id, start)
  ex[, grp := cumsum(c(1L, start[-1] > cummax(end)[-.N] + 1L)), by = tx_id]
  exm <- ex[, .(start = min(start), end = max(end)), by = .(tx_id, Gene, chr, strand, grp)]
  data.table::setorder(exm, tx_id, start)
  intr <- exm[, .(istart = end[-.N] + 1L, iend = start[-1] - 1L),
              by = .(tx_id, Gene, chr, strand)]
  intr <- intr[(iend - istart + 1L) >= 2L * trim_bp + min_width]
  intr[, `:=`(istart = istart + trim_bp, iend = iend - trim_bp)]
  rm(ex, exm); invisible(gc())

  g <- Biostrings::readDNAStringSet(genome_fasta)
  names(g) <- sub(" .*", "", names(g))
  intr <- intr[chr %in% names(g)]
  seqs <- Biostrings::subseq(g[intr$chr], start = intr$istart, end = intr$iend)
  rm(g); invisible(gc())
  minus <- intr$strand == "-"
  seqs[minus] <- Biostrings::reverseComplement(seqs[minus])
  bf <- Biostrings::alphabetFrequency(seqs, baseOnly = TRUE)[, c("A", "C", "G", "T")]
  rm(seqs); invisible(gc())

  intr <- cbind(intr[, .(Gene)], data.table::as.data.table(bf))
  comp <- intr[, .(A = sum(A), C = sum(C), G = sum(G), T = sum(T)), by = Gene]
  comp[, n := A + C + G + T]
  comp[, `:=`(fC = C / n, fG = G / n, fGC = (C + G) / n)]
  comp[]
}

composition_change_between <- function(gene_q, intron_comp, from, to,
                                       bin_width = 0.2, n_boot = 10000L,
                                       seed = 1L, conf = 0.95) {
  #' Relative change in C content between two expression regions, at 4-fold
  #' sites and in the introns of the same genes, and the difference between
  #' the two changes. Genes are resampled within each region.
  #'
  #' @param gene_q data.table: Gene, Q4, n4, Exp.
  #' @param intron_comp data.table: Gene, C, n.
  #' @param from,to c(lo, hi) on the expression category.
  set.seed(seed)
  a <- (1 - conf) / 2
  d <- merge(gene_q, intron_comp[, .(Gene, C, n)], by = "Gene")
  d[, Exp_cat := round(round(Exp / bin_width) * bin_width, 1)]
  i1 <- which(d$Exp_cat >= from[1] & d$Exp_cat <= from[2])
  i2 <- which(d$Exp_cat >= to[1]   & d$Exp_cat <= to[2])
  lvl <- function(i) c(sum(d$Q4[i] * d$n4[i]) / sum(d$n4[i]), sum(d$C[i]) / sum(d$n[i]))
  chg <- function(j1, j2) { x <- lvl(j1); y <- lvl(j2); c(y / x - 1, x, y) }
  pt <- chg(i1, i2)
  bs <- replicate(n_boot, chg(sample(i1, length(i1), TRUE),
                              sample(i2, length(i2), TRUE))[1:2])
  diff_bs <- bs[1, ] - bs[2, ]
  data.table::data.table(
    from = paste(from, collapse = "-"), to = paste(to, collapse = "-"),
    n_genes_from = length(i1), n_genes_to = length(i2),
    Q4_from = pt[3], Q4_to = pt[5], intronC_from = pt[4], intronC_to = pt[6],
    change_4fold  = pt[1], change_4fold_lo  = stats::quantile(bs[1, ], a),
    change_4fold_hi = stats::quantile(bs[1, ], 1 - a),
    change_intron = pt[2], change_intron_lo = stats::quantile(bs[2, ], a),
    change_intron_hi = stats::quantile(bs[2, ], 1 - a),
    difference = pt[1] - pt[2],
    difference_lo = stats::quantile(diff_bs, a),
    difference_hi = stats::quantile(diff_bs, 1 - a),
    P_difference = min(1, 2 * min(mean(diff_bs <= 0), mean(diff_bs >= 0))))
}

# ---------------------------------------------------------------------------
# Population data from calculate_intronic_twostatepi.py
# ---------------------------------------------------------------------------

check_intron_two_allele <- function(tab, Q_expected, tol = 0.05) {
  #' TRUE when the per-gene intronic table came from an all-sites VCF.
  #'
  #' The May 2026 file was built from the variants-only VCF: its C frequency
  #' (~0.32) was the composition of variable sites, its denominators counted
  #' only those sites, and every per-site quantity in it was wrong by a factor
  #' that varied among genes. Two signatures catch it: the columns added with
  #' the fix are absent, and the C frequency is far from V/(U+V).
  tab <- data.table::as.data.table(tab)
  if (!all(c("n_seg_C", "n_fixed_C") %in% names(tab))) {
    message("[intron check] intron_2_allele.csv predates the all-sites fix (no n_seg_C column). ",
            "Regenerate with Bash_scripts/run_intronic_twostate.sh on Lines187.to.IM767v2.1.vcf.gz.")
    return(FALSE)
  }
  q <- tab[, sum(q_pref_C * n_sites) / sum(n_sites)]
  ok <- abs(q - Q_expected) <= tol
  message(sprintf("[intron check] intronic C frequency %.4f vs V/(U+V) = %.4f: %s",
                  q, Q_expected, if (ok) "consistent" else "INCONSISTENT, not used"))
  ok
}

intronic_gene_shape <- function(sfs_file, gene_meta, focal = "C",
                                expression_var = "Mean_Log10_Exp", bin_width = 0.2,
                                min_n = 50L, type_label = "Intron C/non-C") {
  #' Per-gene shape sums for segregating intronic sites, in the same layout as
  #' gene_shape_table(), so the two can be bound and bootstrapped together.
  #' Focal allele = `focal` (C by default), matching C/W and C/G.
  #' `min_n` matches the 4-fold extraction (50 called lines).
  s <- data.table::fread(sfs_file)
  s <- s[s$System == focal & s$n >= min_n]
  g <- s[, .(n_seg    = sum(count),
             f_sum    = sum(count * k / n),
             theta_pi = sum(count * 2 * k * (n - k) / (n * (n - 1))),
             theta_w  = sum(count / (digamma(n) - digamma(1)))), by = Gene]
  meta <- data.table::as.data.table(gene_meta)[, .(Gene, Expr = get(expression_var))]
  g <- merge(g, meta, by = "Gene")
  g[, `:=`(Type = type_label, Exp_cat = round(round(Expr / bin_width) * bin_width, 1))]
  data.table::setcolorder(g, c("Gene", "Type"))
  g[]
}

project_intronic_sfs <- function(sfs_file, gene_table, target_n = 100L, focal = "C") {
  #' Per-gene C-polarised intronic SFS projected to `target_n` chromosomes,
  #' laid out like project_gene_sfs() so fit_shape_S_by_region() applies.
  #' Monomorphic sites go straight to the boundary classes.
  s <- data.table::fread(sfs_file)
  s <- s[s$System == focal & s$n >= target_n]
  gt <- data.table::as.data.table(gene_table)
  seg_col <- paste0("n_seg_", focal); fix_col <- paste0("n_fixed_", focal)
  genes <- sort(intersect(unique(s$Gene), gt$Gene))
  M <- matrix(0, nrow = length(genes), ncol = target_n + 1L,
              dimnames = list(genes, 0:target_n))
  s <- s[Gene %in% genes]
  gi <- match(s$Gene, genes)
  for (r in seq_len(nrow(s)))
    M[gi[r], ] <- M[gi[r], ] + s$count[r] *
      stats::dhyper(0:target_n, s$k[r], s$n[r] - s$k[r], target_n)
  gm <- gt[match(genes, Gene)]
  M[, target_n + 1L] <- M[, target_n + 1L] + gm[[fix_col]]
  M[, 1L] <- M[, 1L] + (gm$n_sites - gm[[seg_col]] - gm[[fix_col]])
  M
}

# ---------------------------------------------------------------------------
# Referee 1: diversity against the no-selection counterfactual
# ---------------------------------------------------------------------------

exon_intron_pi_ratio_by_S <- function(wright_genes, intron_two, U, V,
                                      S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                      ref_bin = "[0,0.25)", min_4fold = 20L, min_intron = 200L,
                                      n_boot = 2000L, seed = 1L) {
  #' Ratio of 4-fold to intronic two-allele (C vs non-C) diversity by S_Wright
  #' bin, relative to a reference bin, against Wright's prediction.
  #'
  #' The intron of the same gene removes everything the two sequences share:
  #' mutation rate, linked selection, demography. What remains, if Wright's
  #' model is right, is pi(S) / pi(S_ref) for the bin's mean S. Both sides are
  #' per site over ALL sites (monomorphic included).
  #'
  #' @param wright_genes data.table: Gene, S (S_Wright_signed), pi4 (pi_2allele), n4.
  #' @param intron_two the corrected intron_2_allele.csv.
  #' @param ref_bin label of the reference bin. [0, 0.25) rather than S < 0,
  #'   which holds the genes whose composition sits below the mutational
  #'   equilibrium and is the noisiest bin.
  set.seed(seed)
  wpi <- function(S) vapply(S, function(s) {
    2 * V * U / ((U + V) * (U + V + 1)) *
      exp(log_kummer_1F1(V + 1, U + V + 2, s) - log_kummer_1F1(V, U + V, s))
  }, numeric(1))
  it <- data.table::as.data.table(intron_two)[, .(Gene, piI = pi_2allele_C, nI = n_sites)]
  d <- merge(data.table::as.data.table(wright_genes), it, by = "Gene")
  d <- d[is.finite(S) & is.finite(pi4) & is.finite(piI) & n4 >= min_4fold & nI >= min_intron]
  d[, Sbin := cut(S, S_breaks, right = FALSE)]
  ratio <- function(x) (sum(x$pi4 * x$n4) / sum(x$n4)) / (sum(x$piI * x$nI) / sum(x$nI))
  bins <- split(d, by = "Sbin", sorted = TRUE)
  bins <- bins[vapply(bins, nrow, integer(1)) > 0]
  stopifnot(ref_bin %in% names(bins))
  ref  <- bins[[ref_bin]]
  bs_ref <- replicate(n_boot, ratio(ref[sample.int(nrow(ref), nrow(ref), TRUE)]))
  out <- data.table::rbindlist(lapply(names(bins), function(b) {
    x  <- bins[[b]]
    bs <- replicate(n_boot, ratio(x[sample.int(nrow(x), nrow(x), TRUE)]))
    infl <- bs / bs_ref
    data.table::data.table(
      Sbin = b, n_genes = nrow(x), mean_S = mean(x$S),
      ratio = ratio(x), ratio_lo = stats::quantile(bs, 0.025), ratio_hi = stats::quantile(bs, 0.975),
      observed_inflation = ratio(x) / ratio(ref),
      inflation_lo = stats::quantile(infl, 0.025), inflation_hi = stats::quantile(infl, 0.975),
      predicted_inflation = wpi(mean(x$S)) / wpi(mean(ref$S)))
  }))
  attr(out, "reference_bin") <- ref_bin
  out[]
}

plot_exon_intron_pi_ratio <- function(tab, title = NULL) {
  #' Observed inflation of 4-fold over intronic diversity by S_Wright bin
  #' (points, gene-bootstrap intervals) against Wright's prediction (line).
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  ggplot(d, aes(x = mean_S)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_line(aes(y = predicted_inflation, linetype = "Wright prediction"),
              colour = "grey20", linewidth = 0.7) +
    geom_errorbar(aes(ymin = inflation_lo, ymax = inflation_hi), width = 0.03,
                  colour = "#2166AC", linewidth = 0.5) +
    geom_point(aes(y = observed_inflation, size = n_genes), colour = "#2166AC") +
    scale_size_continuous(name = "Genes", range = c(1.5, 4)) +
    scale_linetype_manual(values = c("Wright prediction" = "solid"), name = NULL) +
    labs(x = expression("Mean " * S[Wright] * " of the bin"),
         y = expression(pi[4 * "-fold"] / pi[intron] * ", relative to the reference bin"),
         title = title) +
    theme_custom() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

# ---------------------------------------------------------------------------
# Shared uptick: total intronic and 4-fold diversity by expression category
# ---------------------------------------------------------------------------

pi_intron_vs_4fold_by_category <- function(pi_feature_file, gene_meta,
                                           expression_var = "Mean_Log10_Exp",
                                           bin_width = 0.2, bulk_max = 1.2,
                                           min_genes = 40L) {
  #' Site-weighted total pi of introns and of 4-fold sites per expression
  #' category, each relative to its own drift-dominated bulk. Whatever part of
  #' the 4-fold uptick the introns share is not specific to codons.
  f <- data.table::fread(pi_feature_file)
  f[, Gene := sub("\\.v[0-9]+\\.[0-9]+$", "", sub("^MgIM767\\.", "", Gene))]
  meta <- data.table::as.data.table(gene_meta)[, .(Gene, Expr = get(expression_var))]
  meta[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  one <- function(sel, lab) {
    x <- merge(f[sel, .(pi_sum = sum(Pi_sum), sites = sum(Sites)), by = Gene],
               meta, by = "Gene")[sites > 0]
    bulk <- x[Exp_cat <= bulk_max, sum(pi_sum) / sum(sites)]
    x[, .(pi = sum(pi_sum) / sum(sites), n_genes = .N), by = Exp_cat][
      , `:=`(rel_to_bulk = pi / bulk, Sequence = lab)][]
  }
  out <- rbind(one(f$Feature_Type == "intron" & f$Degeneracy == "all", "Introns"),
               one(f$Feature_Type == "exon" & f$Degeneracy == "4-fold", "4-fold sites"))
  out[n_genes >= min_genes][order(Sequence, Exp_cat)]
}
