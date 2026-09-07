## ============================================================================
## recombination_expression_control.R
##
## Referee 1's gBGC alternative, addressed with the recombination map.
##
## gBGC acts through recombination, so if it were producing the expression-
## dependent codon signal this paper reports, the opportunity for gBGC would
## have to covary with expression. It does not: the per-gene crossover rate of
## Lovell et al. (2025) is essentially flat across the expression range.
##
## Two cautions are built into how this is reported.
##
##   1. "Flat" means NEGLIGIBLE, not absent. With >22,000 genes a Spearman rho
##      of 0.05 is significant at p ~ 1e-16 while explaining well under 1% of
##      the variance. The effect size is the claim; the p-value is not.
##   2. The companion CDC panel must use the SAME specification as the CDC
##      figure in the main text. The raw CDC-expression relationship is
##      NEGATIVE (rho = -0.18) because CDC is dominated by CDS length
##      (rho = -0.75) and longer genes are more highly expressed here
##      (rho = +0.32). CDC rises with expression only once length and breadth
##      are controlled. Plotting raw CDC beside raw recombination would
##      contradict the paper's own figure.
## ============================================================================

summarise_recombination_by_expression <- function(rec_table, gene_meta,
                                                  expression_var = "Mean_Log10_Exp",
                                                  bin_width = 0.2,
                                                  min_genes = 100L) {
  #' Mean per-gene crossover rate in each expression category
  #'
  #' @param rec_table Lovell Appendix S6, columns geneID and recombRate (cM/Mb).
  #' @param gene_meta per-gene table with Gene_name and the expression variable.
  #' @return data.table Exp_cat, n_genes, rate_mean, ci_low, ci_high

  suppressPackageStartupMessages(require(data.table))
  r <- data.table::as.data.table(rec_table)[, .(Gene_name = geneID, recombRate)]
  m <- data.table::as.data.table(gene_meta)[, c("Gene_name", expression_var), with = FALSE]
  data.table::setnames(m, expression_var, "Expr")
  d <- merge(r, m, by = "Gene_name")[is.finite(recombRate) & is.finite(Expr)]

  d[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  out <- d[, .(n_genes   = .N,
               rate_mean = mean(recombRate),
               se        = stats::sd(recombRate) / sqrt(.N)), by = Exp_cat]
  out <- out[n_genes >= min_genes]
  out[, `:=`(ci_low = rate_mean - 1.96 * se, ci_high = rate_mean + 1.96 * se)]
  data.table::setorder(out, Exp_cat)
  out[]
}


cdc_partial_by_expression <- function(gene_meta,
                                      expression_var = "Mean_Log10_Exp",
                                      bin_width = 0.2, min_genes = 100L) {
  #' CDC across the expression range, with CDS length and breadth held fixed
  #'
  #' The specification follows the CDC figure in the main text,
  #' CDC ~ te(expression, Exp_breadth) + s(CDS_length_nt), with predictions at
  #' the median breadth and median CDS length, so the curve is the partial
  #' effect of expression rather than the raw trend. Reporting the raw trend
  #' instead would show CDC DECLINING, because CDC is dominated by gene length;
  #' see the header.
  #'
  #' `expression_var` defaults to Mean_Log10_Exp, not the Max used by the main
  #' text figure, so that this panel shares an axis with its companion and with
  #' every other expression-category figure in the paper. The partial effect has
  #' the same sign and a similar magnitude either way.
  #'
  #' @return data.table Exp_cat, fit, ci_low, ci_high (link scale = response here)

  suppressPackageStartupMessages({require(data.table); require(mgcv)})
  d <- data.table::as.data.table(gene_meta)
  need <- c("CDC", expression_var, "Exp_breadth", "CDS_length_nt")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("gene_meta lacks: ", paste(miss, collapse = ", "))
  d <- d[stats::complete.cases(d[, ..need])]
  data.table::setnames(d, expression_var, "Expr")

  g <- mgcv::gam(CDC ~ te(Expr, Exp_breadth) + s(CDS_length_nt),
                 data = d, method = "REML")
  d[, Exp_cat := round(round(Expr / bin_width) * bin_width, 1)]
  cats <- d[, .N, by = Exp_cat][N >= min_genes][order(Exp_cat)]

  nd <- data.table::data.table(
    Expr          = cats$Exp_cat,
    Exp_breadth   = stats::median(d$Exp_breadth),
    CDS_length_nt = stats::median(d$CDS_length_nt))
  p <- stats::predict(g, nd, se.fit = TRUE)
  data.table::data.table(Exp_cat = cats$Exp_cat, n_genes = cats$N,
                         fit = as.numeric(p$fit),
                         ci_low  = as.numeric(p$fit - 1.96 * p$se.fit),
                         ci_high = as.numeric(p$fit + 1.96 * p$se.fit))[]
}


plot_recombination_control <- function(rec_table, cdc_table) {
  #' Two panels: the driver of gBGC, and the codon-bias signal it would have to explain
  suppressPackageStartupMessages({require(ggplot2); require(patchwork); require(data.table)})

  brk <- seq(0, max(rec_table$Exp_cat), 0.4)
  pA <- ggplot(rec_table, aes(x = Exp_cat, y = rate_mean)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#4D4D4D", alpha = 0.18) +
    geom_line(linewidth = 0.6, colour = "#4D4D4D") +
    geom_point(size = 2.2, colour = "#4D4D4D") +
    scale_x_continuous(breaks = brk) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "Crossover rate (cM/Mb)") +
    theme_custom()

  pB <- ggplot(cdc_table, aes(x = Exp_cat, y = fit)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#B2182B", alpha = 0.18) +
    geom_line(linewidth = 0.6, colour = "#B2182B") +
    geom_point(size = 2.2, colour = "#B2182B") +
    scale_x_continuous(breaks = brk) +
    labs(x = expression("Expression level category (log"[10] * ")"),
         y = "Codon deviation coefficient\n(CDS length and breadth held constant)") +
    theme_custom()

  (pA / pB) + patchwork::plot_annotation(tag_levels = "A")
}
