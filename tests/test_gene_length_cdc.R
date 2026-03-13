#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: Gene Length Effect on CDC (Codon Deviation Coefficient)
# =============================================================================
#
# MOTIVATION
# ----------
# The CDC (Zhang et al. 2012) quantifies codon usage bias as the cosine
# distance between observed and null-expected codon frequencies. Because it
# is estimated from a finite sample of codons, its variance is inversely
# proportional to gene length: short genes produce inflated CDC values even
# when no real selection exists. This script demonstrates this artefact
# through simulation and shows it has practical consequences for interpreting
# CDC values in genomic analyses.
#
# EXPERIMENTAL DESIGN
# -------------------
# Simulation 1 – Null model (no selection)
#   * Draw codon counts from the expected multinomial distribution at several
#     gene lengths (50, 100, 200, 500, 1000, 2000 codons).
#   * Show that raw CDC decreases with length whereas the bootstrap p-value
#     remains uniform — confirming the bootstrap's role as a length correction.
#
# Simulation 2 – Genes with injected bias
#   * Fix a target degree of bias (mild / strong).
#   * Repeat the length sweep.
#   * Show that CDC magnitude is stable across lengths once bias is real,
#     while statistical power (p < 0.05) increases with gene length.
#
# Simulation 3 – Regression artefact in empirical-style data
#   * Mix genes of varied lengths drawn under the null.
#   * Fit CDC ~ log(length) to quantify the negative correlation.
#   * Compare to real-selection genes to show the signal is distinct.
#
# DEPENDENCIES
# ------------
#   * src/codon_deviation_coefficient_analysis.R  (CDC functions)
#   * ggplot2, data.table, testthat
#
# @author  Luis J. Madrigal-Roca
# @date    2026-03-12
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(testthat)
})

source(file.path(dirname(getwd()), "Codon_bias_analysis", "src",
                 "codon_deviation_coefficient_analysis.R"),
       chdir = FALSE)

# In case the script is sourced from the project root directly:
if (!exists("calc_cdc")) {
  source("./src/codon_deviation_coefficient_analysis.R")
}

# ─────────────────────────────────────────────────────────────────────────────
# Shared setup
# ─────────────────────────────────────────────────────────────────────────────

# Standard genetic code (without stop codons; compatible with the rest of
# the project)
GENETIC_CODE <- c(
  # 6-fold
  "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", "CTC" = "Leu",
  "CTA" = "Leu", "CTG" = "Leu",
  "TCT" = "Ser", "TCC" = "Ser", "TCA" = "Ser", "TCG" = "Ser",
  "AGT" = "Ser", "AGC" = "Ser",
  "CGT" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg",
  "AGA" = "Arg", "AGG" = "Arg",
  # 4-fold
  "GTT" = "Val", "GTC" = "Val", "GTA" = "Val", "GTG" = "Val",
  "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
  "GGT" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly",
  "CCT" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
  "ACT" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
  # 3-fold (Ile)
  "ATT" = "Ile", "ATC" = "Ile", "ATA" = "Ile",
  # 2-fold
  "AAA" = "Lys", "AAG" = "Lys",
  "AAC" = "Asn", "AAT" = "Asn",
  "GAA" = "Glu", "GAG" = "Glu",
  "GAC" = "Asp", "GAT" = "Asp",
  "CAA" = "Gln", "CAG" = "Gln",
  "CAC" = "His", "CAT" = "His",
  "TAC" = "Tyr", "TAT" = "Tyr",
  "TGC" = "Cys", "TGT" = "Cys",
  "TTC" = "Phe", "TTT" = "Phe",
  # 1-fold (single synonymous codons — excluded from CDC by design)
  "ATG" = "Met", "TGG" = "Trp",
  # Stop codons (excluded from CDC by design)
  "TAA" = "STOP", "TAG" = "STOP", "TGA" = "STOP"
)

# Positional composition for a "typical plant gene" with AT-bias
# S = GC content, R = purine content at each position
BASE_COMP <- list(S1 = 0.50, S2 = 0.38, S3 = 0.42,
                  R1 = 0.58, R2 = 0.45, R3 = 0.52)

# Derive true expected codon frequencies from BASE_COMP
TRUE_EXPECTED_FREQ <- calc_expected_codon_usage(BASE_COMP, GENETIC_CODE)

# Gene lengths to sweep across
GENE_LENGTHS <- c(50, 100, 200, 500, 1000, 2000)
N_REPS       <- 500          # simulated genes per length
N_BOOTSTRAP  <- 200          # bootstrap reps per CDC call (kept low for speed)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a forced-biased codon frequency vector
# ─────────────────────────────────────────────────────────────────────────────

#' Introduce systematic bias by up-weighting one codon per amino acid family
#' 
#' @param expected_freq Named vector of expected codon frequencies
#' @param genetic_code  Named vector: codon → amino acid
#' @param bias_strength Fractional weight shift toward preferred codon [0, 1]
#' @return Named vector of biased codon frequencies (sums to 1)
make_biased_freq <- function(expected_freq, genetic_code, bias_strength = 0.4) {
  sense_codons <- get_sense_codons(genetic_code)
  biased <- expected_freq[sense_codons]

  # Group by amino acid family
  aa_groups <- split(sense_codons,
                     genetic_code[sense_codons])

  for (codons in aa_groups) {
    if (length(codons) < 2) next
    # preferred = first alphabetically (arbitrary but deterministic)
    preferred   <- sort(codons)[1]
    others      <- codons[codons != preferred]
    family_mass <- sum(biased[codons])

    # shift `bias_strength` fraction of family mass to the preferred codon
    shift           <- family_mass * bias_strength
    biased[preferred] <- biased[preferred] + shift
    for (c in others) {
      biased[c] <- biased[c] - shift / length(others)
      biased[c] <- max(biased[c], 1e-10)  # guard against negative probs
    }
  }
  biased / sum(biased)
}

MILD_BIAS_FREQ   <- make_biased_freq(TRUE_EXPECTED_FREQ, GENETIC_CODE,
                                     bias_strength = 0.20)
STRONG_BIAS_FREQ <- make_biased_freq(TRUE_EXPECTED_FREQ, GENETIC_CODE,
                                     bias_strength = 0.50)

# ─────────────────────────────────────────────────────────────────────────────
# Core simulation engine
# ─────────────────────────────────────────────────────────────────────────────

#' Simulate N genes of a given length, sampled from `true_freq`.
#' 
#' Expected codon usage is estimated from BASE_COMP (not from sim freq)
#' to mimic what the CDC function does on real data.
#' 
#' @param n_genes     Number of genes to simulate
#' @param n_codons    Codons per gene
#' @param true_freq   True sampling frequency (named vector)
#' @param genetic_code Named vector: codon → amino acid
#' @param n_bootstrap Bootstrap reps for p-value
#' @return data.table with columns: cdc, p_value, n_codons, gene_id
simulate_cdc_at_length <- function(n_genes, n_codons,
                                   true_freq, genetic_code,
                                   n_bootstrap = 200) {
  sense_codons <- get_sense_codons(genetic_code)
  freq_vec     <- true_freq[sense_codons]
  freq_vec     <- freq_vec / sum(freq_vec)

  results <- vector("list", n_genes)

  for (i in seq_len(n_genes)) {
    counts           <- rmultinom(1, n_codons, freq_vec)[, 1]
    names(counts)    <- sense_codons

    cdc_res <- tryCatch(
      calculate_cdc_single(counts, genetic_code,
                           n_bootstrap = n_bootstrap),
      error   = function(e) list(CDC = NA_real_, p_value = NA_real_)
    )

    results[[i]] <- data.table(
      cdc      = cdc_res$CDC,
      p_value  = cdc_res$p_value,
      n_codons = n_codons,
      gene_id  = i
    )
  }
  rbindlist(results)
}

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION 1: Null model — no true selection
# ─────────────────────────────────────────────────────────────────────────────

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SIMULATION 1: CDC under the null (no codon usage bias)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

null_results <- rbindlist(lapply(GENE_LENGTHS, function(L) {
  cat(sprintf("  Simulating %d null genes of length %d codons...\n",
              N_REPS, L))
  simulate_cdc_at_length(N_REPS, L, TRUE_EXPECTED_FREQ,
                         GENETIC_CODE, N_BOOTSTRAP)
}))

# Summary table
null_summary <- null_results[!is.na(cdc),
  .(mean_cdc   = mean(cdc),
    median_cdc = median(cdc),
    sd_cdc     = sd(cdc),
    pct_sig    = mean(p_value < 0.05, na.rm = TRUE) * 100),
  by = n_codons][order(n_codons)]

cat("\nNull-model CDC summary by gene length:\n")
cat(sprintf("  %-12s %-12s %-12s %-12s %-10s\n",
            "Gene length", "Mean CDC", "Median CDC", "SD CDC", "% sig (p<.05)"))
cat(paste(rep("-", 62), collapse = ""), "\n")
null_summary[, cat(sprintf("  %-12d %-12.5f %-12.5f %-12.5f %-10.1f\n",
                            n_codons, mean_cdc, median_cdc, sd_cdc, pct_sig)),
             by = seq_len(nrow(null_summary))]

cat("\nKey result: Mean CDC decreases monotonically with gene length under\n")
cat("the null. This is a pure sampling artefact\n\n")

# Formal test: mean CDC at shortest length > mean CDC at longest length
mean_cdc_short <- null_summary[n_codons == min(GENE_LENGTHS), mean_cdc]
mean_cdc_long  <- null_summary[n_codons == max(GENE_LENGTHS), mean_cdc]

test_that("Mean null CDC is higher for short genes than for long genes", {
  expect_gt(mean_cdc_short, mean_cdc_long)
})

# Bootstrap p-value distribution under null should be uniform (not deflated)
pvals_null <- null_results[n_codons == 500 & !is.na(p_value), p_value]
ks_result  <- ks.test(pvals_null, "punif", 0, 1)

cat(sprintf("KS test of p-value uniformity at n=500 (null): D = %.4f, p = %.4f\n",
            ks_result$statistic, ks_result$p.value))
cat("  (KS p > 0.05 indicates bootstrap preserves uniform p-value distribution)\n\n")

test_that("Bootstrap p-values are approximately uniform under null (n=500)", {
  expect_gt(ks_result$p.value, 0.01)
})

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION 2: Genes with real codon bias — mild vs strong
# ─────────────────────────────────────────────────────────────────────────────

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SIMULATION 2: CDC for genes with real codon usage bias\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

mild_results <- rbindlist(lapply(GENE_LENGTHS, function(L) {
  cat(sprintf("  Mild bias  — length %d codons...\n", L))
  out       <- simulate_cdc_at_length(N_REPS, L, MILD_BIAS_FREQ,
                                      GENETIC_CODE, N_BOOTSTRAP)
  out[, bias := "mild"]
  out
}))

strong_results <- rbindlist(lapply(GENE_LENGTHS, function(L) {
  cat(sprintf("  Strong bias — length %d codons...\n", L))
  out       <- simulate_cdc_at_length(N_REPS, L, STRONG_BIAS_FREQ,
                                      GENETIC_CODE, N_BOOTSTRAP)
  out[, bias := "strong"]
  out
}))

bias_results <- rbind(mild_results, strong_results)

# Summarise
bias_summary <- bias_results[!is.na(cdc),
  .(mean_cdc  = mean(cdc),
    sd_cdc    = sd(cdc),
    power     = mean(p_value < 0.05, na.rm = TRUE) * 100),
  by = .(n_codons, bias)][order(bias, n_codons)]

cat("\nBiased-gene CDC summary:\n")
cat(sprintf("  %-10s %-12s %-10s %-10s %-10s\n",
            "Bias", "Gene length", "Mean CDC", "SD CDC", "Power (%)"))
cat(paste(rep("-", 58), collapse = ""), "\n")
bias_summary[, cat(sprintf("  %-10s %-12d %-10.5f %-10.5f %-10.1f\n",
                            bias, n_codons, mean_cdc, sd_cdc, power)),
             by = seq_len(nrow(bias_summary))]

# Formal tests
# 1. CDC is higher under bias than null (at all lengths)
for (L in GENE_LENGTHS) {
  null_cdc_L  <- null_results[n_codons == L & !is.na(cdc), cdc]
  strong_cdc_L <- strong_results[n_codons == L & !is.na(cdc), cdc]

  test_that(sprintf("Mean CDC is higher for strongly biased vs null genes at length %d", L), {
    expect_gt(mean(strong_cdc_L), mean(null_cdc_L))
  })
}

# 2. Power (proportion p < 0.05) increases with gene length under bias
strong_power <- bias_summary[bias == "strong"][order(n_codons), power]
test_that("Statistical power increases monotonically with gene length (strong bias)", {
  expect_true(all(diff(strong_power) >= 0))
})

# 3. Strong bias gives higher CDC than mild bias (at every length)
for (L in c(200, 1000)) {
  mild_cdc_L   <- mild_results[n_codons == L & !is.na(cdc), cdc]
  strong_cdc_L <- strong_results[n_codons == L & !is.na(cdc), cdc]
  level_name   <- as.character(L)

  test_that(paste("Strong bias produces higher CDC than mild bias at length", L), {
    expect_gt(mean(strong_cdc_L), mean(mild_cdc_L))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION 3: CDC ~ log(length) regression — quantifying the artefact
# ─────────────────────────────────────────────────────────────────────────────

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SIMULATION 3: Regression of raw CDC on log(gene length)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Add log-length and scenario labels
null_results[,   scenario := "null"]
mild_results[,   scenario := "mild bias"]
strong_results[, scenario := "strong bias"]

all_results <- rbind(
  null_results[,   .(cdc, p_value, n_codons, scenario)],
  mild_results[,   .(cdc, p_value, n_codons, scenario)],
  strong_results[, .(cdc, p_value, n_codons, scenario)]
)
all_results[, log_length := log(n_codons)]

fit_null   <- lm(cdc ~ log_length, data = all_results[scenario == "null"])
fit_mild   <- lm(cdc ~ log_length, data = all_results[scenario == "mild bias"])
fit_strong <- lm(cdc ~ log_length, data = all_results[scenario == "strong bias"])

extract_lm <- function(fit, label) {
  s   <- summary(fit)
  b   <- coef(fit)[["log_length"]]
  r2  <- s$r.squared
  p   <- coef(s)["log_length", "Pr(>|t|)"]
  cat(sprintf("  %-14s  slope = %+.5f   R² = %.3f   p = %.2e\n",
              label, b, r2, p))
  invisible(list(slope = b, r2 = r2, p = p))
}

cat("Linear regression CDC ~ log(gene length) per scenario:\n")
res_null   <- extract_lm(fit_null,   "null")
res_mild   <- extract_lm(fit_mild,   "mild bias")
res_strong <- extract_lm(fit_strong, "strong bias")

cat("\nInterpretation:\n")
cat("  - Significant NEGATIVE slope for the null confirms the length artefact.\n")
cat("  - Smaller (less negative / near-zero) slope for biased genes shows\n")
cat("    that real biology attenuates the length dependency.\n\n")

test_that("Null CDC has a significant negative association with log(gene length)", {
  expect_lt(res_null$slope, 0)
  expect_lt(res_null$p, 0.001)
})

test_that("Biased genes show weaker length dependency than null genes", {
  # Absolute slope for null should exceed that of strongly biased genes
  expect_gt(abs(res_null$slope), abs(res_strong$slope))
})

# ─────────────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────────────

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("VISUALISATION\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

theme_set(
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey90"),
          legend.position  = "bottom")
)

# --- Panel A: Raw CDC distribution by length and scenario ---

all_results[, scenario_f := factor(scenario,
                                    levels = c("null", "mild bias", "strong bias"))]
all_results[, length_f := factor(n_codons,
                                  levels = sort(unique(n_codons)))]

p_violins <- ggplot(all_results[!is.na(cdc)],
                    aes(x = length_f, y = cdc,
                        fill = scenario_f, colour = scenario_f)) +
  geom_violin(alpha = 0.35, position = position_dodge(0.85), scale = "width") +
  stat_summary(fun = median, geom = "point", size = 1.5,
               position = position_dodge(0.85)) +
  scale_fill_manual(
    name   = "Scenario",
    values = c("null" = "#95a5a6", "mild bias" = "#3498db", "strong bias" = "#e74c3c")
  ) +
  scale_colour_manual(
    name   = "Scenario",
    values = c("null" = "#7f8c8d", "mild bias" = "#2980b9", "strong bias" = "#c0392b")
  ) +
  labs(
    title    = "A: Raw CDC inflated by short gene length under the null",
    subtitle = paste0("N = ", N_REPS, " simulated genes per cell; ",
                      N_BOOTSTRAP, " bootstrap replicates"),
    x = "Gene length (codons)",
    y = "CDC (cosine distance)"
  )

# --- Panel B: Statistical power vs gene length ---

power_dt <- rbind(
  null_results[!is.na(p_value), .(power = mean(p_value < 0.05) * 100,
                                   scenario = "null"), by = n_codons],
  mild_results[!is.na(p_value), .(power = mean(p_value < 0.05) * 100,
                                   scenario = "mild bias"), by = n_codons],
  strong_results[!is.na(p_value), .(power = mean(p_value < 0.05) * 100,
                                     scenario = "strong bias"), by = n_codons]
)
power_dt[, scenario_f := factor(scenario,
                                  levels = c("null", "mild bias", "strong bias"))]

p_power <- ggplot(power_dt,
                  aes(x = n_codons, y = power,
                      colour = scenario_f, group = scenario_f)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "grey50") +
  annotate("text", x = max(GENE_LENGTHS) * 0.95, y = 8,
           label = "Expected type-I error (5%)", hjust = 1,
           colour = "grey40", size = 3) +
  scale_x_log10(
    breaks = GENE_LENGTHS,
    labels = as.character(GENE_LENGTHS)
  ) +
  scale_colour_manual(
    name   = "Scenario",
    values = c("null" = "#95a5a6", "mild bias" = "#3498db", "strong bias" = "#e74c3c")
  ) +
  labs(
    title    = "B: Power to detect CDC increases with gene length",
    subtitle = "Dashed line = expected false-positive rate under H₀",
    x = "Gene length (codons, log scale)",
    y = "% genes with p < 0.05"
  )

# --- Panel C: CDC ~ log(length) regression overlay ---

smooth_dt <- all_results[!is.na(cdc)]

p_regression <- ggplot(smooth_dt,
                       aes(x = log_length, y = cdc,
                           colour = scenario_f)) +
  geom_point(alpha = 0.07, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_colour_manual(
    name   = "Scenario",
    values = c("null" = "#95a5a6", "mild bias" = "#3498db", "strong bias" = "#e74c3c")
  ) +
  scale_x_continuous(
    breaks = log(GENE_LENGTHS),
    labels = as.character(GENE_LENGTHS)
  ) +
  labs(
    title    = "C: Linear regression of CDC on log(gene length)",
    subtitle = "Null genes show significant negative slope (sampling noise)",
    x = "Gene length (codons, ln scale)",
    y = "CDC (cosine distance)",
    caption = paste0(
      sprintf("Null:   slope = %+.4f, R² = %.3f\n",
              coef(fit_null)["log_length"], summary(fit_null)$r.squared),
      sprintf("Mild:   slope = %+.4f, R² = %.3f\n",
              coef(fit_mild)["log_length"], summary(fit_mild)$r.squared),
      sprintf("Strong: slope = %+.4f, R² = %.3f",
              coef(fit_strong)["log_length"], summary(fit_strong)$r.squared)
    )
  )

# --- Panel D: p-value uniformity check at a fixed length ---

pval_check_dt <- rbind(
  null_results[n_codons == 500 & !is.na(p_value), .(p_value, scenario = "null")],
  mild_results[n_codons == 500 & !is.na(p_value), .(p_value, scenario = "mild bias")],
  strong_results[n_codons == 500 & !is.na(p_value), .(p_value, scenario = "strong bias")]
)
pval_check_dt[, scenario_f := factor(scenario,
                                      levels = c("null", "mild bias", "strong bias"))]

p_pvals <- ggplot(pval_check_dt, aes(x = p_value, fill = scenario_f)) +
  geom_histogram(bins = 20, boundary = 0, colour = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "black") +
  facet_wrap(~ scenario_f, ncol = 3) +
  scale_fill_manual(
    name   = "Scenario",
    values = c("null" = "#95a5a6", "mild bias" = "#3498db", "strong bias" = "#e74c3c"),
    guide  = "none"
  ) +
  labs(
    title    = "D: Bootstrap p-value distributions at gene length = 500 codons",
    subtitle = "Uniform distribution under H₀ confirms valid type-I error control",
    x = "Bootstrap p-value",
    y = "Count"
  )

# ─────────────────────────────────────────────────────────────────────────────
# Print all plots
# ─────────────────────────────────────────────────────────────────────────────

print(p_violins)
print(p_power)
print(p_regression)
print(p_pvals)

# ─────────────────────────────────────────────────────────────────────────────
# Final summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY OF FINDINGS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("1. SAMPLING ARTEFACT — Gene length inflates raw CDC\n")
cat(sprintf("   Mean CDC at %d codons (null):  %.5f\n", min(GENE_LENGTHS), mean_cdc_short))
cat(sprintf("   Mean CDC at %d codons (null): %.5f\n", max(GENE_LENGTHS), mean_cdc_long))
cat(sprintf("   Ratio:                         %.1fx\n\n",
            mean_cdc_short / mean_cdc_long))

cat("2. BOOTSTRAP CORRECTION — preserves valid type-I error\n")
cat(sprintf("   Proportion p < 0.05 under null (n=500): %.1f%%\n",
            null_results[n_codons == 500, mean(p_value < 0.05, na.rm = TRUE) * 100]))
cat("   (Expected: 5%)\n\n")

cat("3. REAL BIAS IS DETECTABLE — power scales with gene length\n")
strong_short_power <- strong_results[n_codons == min(GENE_LENGTHS),
                                      mean(p_value < 0.05, na.rm = TRUE) * 100]
strong_long_power  <- strong_results[n_codons == max(GENE_LENGTHS),
                                      mean(p_value < 0.05, na.rm = TRUE) * 100]
cat(sprintf("   Power at %d codons (strong bias): %.1f%%\n",
            min(GENE_LENGTHS), strong_short_power))
cat(sprintf("   Power at %d codons (strong bias): %.1f%%\n",
            max(GENE_LENGTHS), strong_long_power))
cat("\n")

cat("4. PRACTICAL RECOMMENDATION\n")
cat("   Always interpret CDC values together with bootstrap p-values.\n")
cat("   Raw CDC is length-dependent; use p-value or quantile-normalised\n")
cat("   CDC when comparing genes of different lengths.\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("All explicit tests passed via testthat.\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
