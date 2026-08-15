## ============================================================================
## paper2_linked_selection.R — within-gene position and linked selection
##
## PROVENANCE
##   Extracted verbatim from `main.R` on 2026-08-14 when the New Phytologist
##   manuscript (NPH-MS-2026-57041) was refocused onto codon usage bias alone.
##   Source blocks: main.R lines 2535-2667 (Sections 12.0 / 12.0b) and
##   3291-3472 (RESULTS 10 / Section 15).  The pre-split file is preserved at
##   `archive/main_full_analysis_2026-08-14.R`.
##
## WHY THIS IS NOT IN THE CUB PAPER
##   Referee 1 did not dispute the pattern (elevated pi near gene starts,
##   preferred codons clustering at the 5' end).  The objection was that the
##   manuscript contains no direct measure of recombination rate, so the
##   background-selection-amelioration explanation cannot be separated from:
##     (1) GC-biased gene conversion following non-crossover recombination,
##         which would mimic selection for the mostly-C-ending preferred codons;
##     (2) recombination-associated mutagenesis raising the synonymous mutation
##         rate wherever recombination is elevated.
##   Discriminating these needs polymorphism-vs-divergence: if recombination is
##   mutagenic, divergence should be elevated in the same regions as
##   polymorphism; under linked selection it should not.
##
## WHAT PAPER 2 STILL NEEDS
##   - M. jungermannioides (JUNG1) outgroup for polarisation.  M. nasutus and
##     M. tilingii share polymorphism with M. guttatus; the M. cardinalis group
##     is too distant to align.  JGI release + co-authorship under negotiation.
##   - Within-gene recombination map from the SF x 767 F2 cross (~1,000
##     crossovers resolved by PacBio; see the recom-sf-767-f2 project).
##   - Divergence pipeline: see `Divergence_AnchorWave` notes.
##
## INPUTS EXPECTED IN THE ENVIRONMENT
##   integrated_data  gene-level table with Gene_name, Max_Log10_Exp,
##                    Exp_breadth (built by main_CUB.R Sections 4-5.5)
##   Data files:      data/all_chromosomes.pi_per_gene_feature.txt
##                    data/all_chromosomes.codon_frequencies_preferred.txt
##   Helpers:         src/theme_custom.R and the standard library set loaded by
##                    src/set_environment.R
##
##   Run `main_CUB.R` through Section 12 first, then source this file, or
##   rebuild `integrated_data` independently before running it standalone.
##
## OUTPUTS
##   results/pi_by_gene_distance.pdf                  (was Figure 7C)
##   results/ramp_by_expression_polymorphism.pdf      (was Figure 6A)
##   integrated_data columns Sites_/Pi_sum_ *_first300 / *_after300
## ============================================================================

## ****************************************************************************
## 1) Positional decomposition of nucleotide diversity ----
##    (main.R Section 12.0 / 12.0b)
## ____________________________________________________________________________

# 12.0) Positional decomposition of pi: first 300 bp vs after 300 bp ----
# Compute per-gene 4-fold and 0-fold pi split at the 300 bp CDS boundary.
# This captures Hill-Robertson interference patterns: linked selection near
# the 5' end (translational ramp) vs the gene body.

pi_feature <- data.table::fread("data/all_chromosomes.pi_per_gene_feature.txt")

# Harmonize gene names: feature file has "MgIM767.01G000100.v2.1",
# integrated_data has "MgIM767.01G000100"
pi_feature[, Gene := sub("\\.v[0-9.]+$", "", Gene)]

cat(sprintf("\n=== Section 12.0: Positional pi decomposition | %d unique genes in feature file (no AnaCoDa filter applied) ===\n",
            length(unique(pi_feature$Gene))))

# Per-exon CDS sizes (from "all" degeneracy) and degeneracy-specific pi
exon_all <- pi_feature[Feature_Type == "exon" & Degeneracy == "all",
                       .(Gene, Feature_Num, Exon_Sites = Sites)]
exon_4fold <- pi_feature[Feature_Type == "exon" & Degeneracy == "4-fold",
                         .(Gene, Feature_Num, Sites_4f = Sites, Pi_sum_4f = Pi_sum)]
exon_0fold <- pi_feature[Feature_Type == "exon" & Degeneracy == "0-fold",
                         .(Gene, Feature_Num, Sites_0f = Sites, Pi_sum_0f = Pi_sum)]

exon_data <- merge(exon_all, exon_4fold, by = c("Gene", "Feature_Num"), all.x = TRUE)
exon_data <- merge(exon_data, exon_0fold, by = c("Gene", "Feature_Num"), all.x = TRUE)
exon_data[is.na(Sites_4f), c("Sites_4f", "Pi_sum_4f") := 0]
exon_data[is.na(Sites_0f), c("Sites_0f", "Pi_sum_0f") := 0]

# Order exons and compute cumulative CDS position
data.table::setorder(exon_data, Gene, Feature_Num)
exon_data[, cum_end := cumsum(Exon_Sites), by = Gene]
exon_data[, cum_start := cum_end - Exon_Sites + 1L, by = Gene]

# Fraction of each exon falling within the first 300 bp of CDS
bp_cutoff <- 300
exon_data[, frac_first := data.table::fifelse(
  cum_end <= bp_cutoff, 1.0,
  data.table::fifelse(cum_start > bp_cutoff, 0.0,
                      (bp_cutoff - cum_start + 1) / Exon_Sites)
)]

# Aggregate per gene: first 300 bp vs after 300 bp
pi_300bp <- exon_data[, .(
  Sites_4fold_first300  = sum(Sites_4f * frac_first),
  Pi_sum_4fold_first300 = sum(Pi_sum_4f * frac_first),
  Sites_4fold_after300  = sum(Sites_4f * (1 - frac_first)),
  Pi_sum_4fold_after300 = sum(Pi_sum_4f * (1 - frac_first)),
  Sites_0fold_first300  = sum(Sites_0f * frac_first),
  Pi_sum_0fold_first300 = sum(Pi_sum_0f * frac_first),
  Sites_0fold_after300  = sum(Sites_0f * (1 - frac_first)),
  Pi_sum_0fold_after300 = sum(Pi_sum_0f * (1 - frac_first))
), by = .(Gene_name = Gene)]

# Merge into integrated_data
integrated_data <- integrated_data |>
  dplyr::left_join(pi_300bp, by = "Gene_name")

cat(sprintf("300 bp decomposition: %d genes matched\n",
            sum(!is.na(integrated_data$Sites_4fold_first300))))

# 12.0b) π by distance from gene start — exon (4-fold) vs intron ----
# Assigns each feature a cumulative genomic position by interleaving exons and
# introns in Feature_Num order: exon k → genomic rank 2k-1, intron k → 2k.
# Feature bp length is approximated by "all"-degeneracy Sites (surveyed sites).
# Exons use 4-fold π; introns use all-site π (the only class they carry).

cat(sprintf("\n=== Section 12.0b: pi vs distance from gene start | %d unique genes (pi_feature; no AnaCoDa filter) ===\n",
            length(unique(pi_feature$Gene))))

feat_pos <- pi_feature[Degeneracy == "all",
                        .(Gene, Feature_Type, Feature_Num, feat_bp = Sites)]

feat_pos[Feature_Type == "exon",   genomic_order := 2L * Feature_Num - 1L]
feat_pos[Feature_Type == "intron", genomic_order := 2L * Feature_Num]
data.table::setorder(feat_pos, Gene, genomic_order)

# Cumulative distance from gene start (approximation via surveyed-site lengths)
feat_pos[, cum_end   := cumsum(feat_bp),           by = Gene]
feat_pos[, cum_start := cum_end - feat_bp,          by = Gene]
feat_pos[, midpoint  := (cum_start + cum_end) / 2L, by = Gene]

# 1 kb windows up to 10 kb, last bin 10–20 kb
dist_breaks <- c(0, seq(1000, 10000, 1000), 20000)
dist_labels <- c(as.character(seq(1000, 10000, 1000)), "10k-20k")

feat_pos[, dist_bin := cut(midpoint, breaks = dist_breaks, labels = dist_labels,
                            right = TRUE, include.lowest = FALSE)]

# π values at the appropriate degeneracy class per feature type
pi_exon_4f <- pi_feature[Feature_Type == "exon"   & Degeneracy == "4-fold",
                          .(Gene, Feature_Type, Feature_Num, Pi_sum, n_pi = Sites)]
pi_int_all <- pi_feature[Feature_Type == "intron" & Degeneracy == "all",
                          .(Gene, Feature_Type, Feature_Num, Pi_sum, n_pi = Sites)]
pi_vals_dist <- data.table::rbindlist(list(pi_exon_4f, pi_int_all))

feat_pi_dist <- merge(
  pi_vals_dist,
  feat_pos[, .(Gene, Feature_Type, Feature_Num, dist_bin)],
  by = c("Gene", "Feature_Type", "Feature_Num"), all.x = TRUE
)

# Weighted mean π per feature type × distance window
pi_by_dist <- feat_pi_dist[n_pi > 0 & !is.na(dist_bin), .(
  n_features = .N,
  Pi_mean    = sum(Pi_sum) / sum(n_pi)
), by = .(Feature_Type, dist_bin)]

pi_by_dist[, dist_bin := factor(dist_bin, levels = dist_labels)]
data.table::setorder(pi_by_dist, Feature_Type, dist_bin)

plot_pi_by_dist <- ggplot(
    pi_by_dist,
    aes(x = dist_bin, y = Pi_mean, fill = Feature_Type)
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c("exon" = "#1F4E79", "intron" = "#E07830"),
    labels = c("exon" = "Exon (4-fold)", "intron" = "Intron")
  ) +
  labs(
    x = "Distance from gene start (bp)",
    y = "Nucleotide Diversity",
    fill = NULL
  ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./results/pi_by_gene_distance.pdf", plot_pi_by_dist, width = 10, height = 6)

rm(feat_pos, pi_exon_4f, pi_int_all, pi_vals_dist, feat_pi_dist, pi_by_dist,
   plot_pi_by_dist, dist_breaks, dist_labels)

rm(pi_feature, exon_all, exon_4fold, exon_0fold, exon_data, pi_300bp)
gc()

## ****************************************************************************
## 2) Translational ramp hypothesis ----
##    (main.R RESULTS 10 / Section 15)
## ____________________________________________________________________________
## ============================================================================
## RESULTS 10 — Translational ramp: preferred codons cluster near the 5' end
##   Produces:
##     Figure 6A  Preferred codon frequency vs codon position by expression tier
##                (`ramp_by_expression_polymorphism.pdf`)
##     Also runs the gene-level GAM that produces the supplementary contour
##     `global_gene_level_pref_freq_contour.pdf` (Fig S10 / Fig 6B alt).
## ============================================================================

# ******************************************************************************
# 15) Testing the translational ramp hypothesis ----
# ______________________________________________________________________________

# 15.1) Polymorphism-based ramp models ----
# Per-codon preferred-allele frequencies in the first 200 codons; aggregate into
# 5-codon windows; fit beta-regression bam() with random intercept by gene.
# We subsample 3,000 genes (seed = 1998) — the gene RE makes the full set
# prohibitively slow without changing the population-level smooth.

poly_data <- fread(
  "data/all_chromosomes.codon_frequencies_preferred.txt",
  select = c("Gene", "Codon_Pos", "Preferred_Freq", "Non_Preferred_Freq"),
  showProgress = FALSE
)

poly_data <- poly_data |>
  dplyr::mutate(Gene_clean = paste0("MgIM767.", Gene)) |>
  dplyr::rename(Position = Codon_Pos) |>
  dplyr::filter(Position <= 200)

poly_with_exp <- poly_data |>
  dplyr::left_join(
    integrated_data |>
      dplyr::select(Gene_name, Max_Log10_Exp, Exp_breadth),
    by = c("Gene_clean" = "Gene_name")
  ) |>
  dplyr::filter(!is.na(Max_Log10_Exp), !is.na(Exp_breadth)) |>
  dplyr::mutate(
    Exp_Z = as.numeric(scale(Max_Log10_Exp)),
    Breadth_Z = as.numeric(scale(Exp_breadth)),
    Gene_clean = factor(Gene_clean)
  )

cat(sprintf("Loaded %d codon positions from %d genes\n",
            nrow(poly_with_exp), length(unique(poly_with_exp$Gene_clean))))

n_subsample_genes <- 3000
all_genes_15 <- unique(poly_with_exp$Gene_clean)
set.seed(1998)
sampled_genes_15 <- sample(all_genes_15,
                           size = min(n_subsample_genes, length(all_genes_15)))
poly_with_exp <- poly_with_exp |>
  dplyr::filter(Gene_clean %in% sampled_genes_15) |>
  dplyr::mutate(Gene_clean = droplevels(Gene_clean))

cat(sprintf("After subsampling: %d codon positions from %d genes (seed = 1998)\n",
            nrow(poly_with_exp), length(unique(poly_with_exp$Gene_clean))))

window_size <- 5
poly_agg <- poly_with_exp |>
  dplyr::mutate(Window = ceiling(Position / window_size)) |>
  dplyr::group_by(Gene_clean, Window, Exp_Z, Breadth_Z) |>
  dplyr::summarize(
    Position_mid = mean(Position),
    Preferred_Freq_mean = mean(Preferred_Freq, na.rm = TRUE),
    n_codons = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::filter(Position_mid <= 200) |>
  dplyr::mutate(
    Preferred_Freq_beta = dplyr::case_when(
      Preferred_Freq_mean <= 0.001 ~ 0.001,
      Preferred_Freq_mean >= 0.999 ~ 0.999,
      TRUE ~ Preferred_Freq_mean
    )
  )

if (any(poly_agg$Preferred_Freq_beta <= 0 | poly_agg$Preferred_Freq_beta >= 1)) {
  stop("Beta regression requires values strictly between 0 and 1")
}

fit_ramp_poly <- bam(
  Preferred_Freq_beta ~
    s(Position_mid, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  data = poly_agg, family = betar(),
  method = "fREML", discrete = TRUE, nthreads = 1
)

fit_ramp_int_poly <- bam(
  Preferred_Freq_beta ~
    s(Position_mid, k = 10, bs = "tp") +
    s(Position_mid, by = Exp_Z, k = 10, bs = "tp") +
    Exp_Z + Breadth_Z + Exp_Z:Breadth_Z +
    s(Gene_clean, bs = "re"),
  data = poly_agg, family = betar(),
  method = "fREML", discrete = TRUE, nthreads = 1
)

# Plot 1: Global ramp shape from polymorphism data
pred_positions <- data.frame(
  Position_mid = seq(5, 200, by = 2),
  Exp_Z = 0,
  Breadth_Z = 0,
  Gene_clean = poly_agg$Gene_clean[1]
)

pred_ramp <- predict(fit_ramp_poly, newdata = pred_positions,
                     type = "response", se.fit = TRUE,
                     exclude = "s(Gene_clean)",
                     unconditional = TRUE)

pred_positions$fit   <- pred_ramp$fit
pred_positions$lower <- pred_ramp$fit - 1.96 * pred_ramp$se.fit
pred_positions$upper <- pred_ramp$fit + 1.96 * pred_ramp$se.fit

plot_ramp_poly <- ggplot(pred_positions, aes(x = Position_mid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = fit), color = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Translational Ramp: Population Polymorphism Data",
    subtitle = "Frequency of preferred codons across positions",
    x = "Codon Position",
    y = "Mean Preferred Codon Frequency",
    caption = "Ribbon = +/-1.96 SE | Data from population genomics"
  ) +
  theme_bw(base_size = 12)

ggsave("./results/translational_ramp_polymorphism.pdf",
       plot_ramp_poly, width = 10, height = 6)

# Plot 2: Ramp by expression level
pred_grid_exp <- expand.grid(
  Position_mid = seq(5, 200, by = 5),
  Exp_Z = c(-1.5, 0, 1.5),
  Breadth_Z = 0
) |>
  dplyr::mutate(
    Gene_clean = poly_agg$Gene_clean[1],
    Exp_level = factor(
      Exp_Z,
      levels = c(-1.5, 0, 1.5),
      labels = c("Low Expression", "Medium Expression", "High Expression")
    )
  )

pred_exp <- predict(fit_ramp_int_poly, newdata = pred_grid_exp,
                    type = "response", se.fit = TRUE,
                    exclude = "s(Gene_clean)")

pred_grid_exp$fit   <- pred_exp$fit
pred_grid_exp$lower <- pred_exp$fit - 1.96 * pred_exp$se.fit
pred_grid_exp$upper <- pred_exp$fit + 1.96 * pred_exp$se.fit

plot_ramp_exp_poly <- ggplot(pred_grid_exp,
                             aes(x = Position_mid, color = Exp_level, fill = Exp_level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed", alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Ramp Shape by Expression: Population Data",
    subtitle = "Does selection for preferred codons vary by expression level?",
    x = "Codon Position",
    y = "Mean Preferred Codon Frequency",
    color = NULL, fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("./results/ramp_by_expression_polymorphism.pdf",
       plot_ramp_exp_poly, width = 10, height = 6)

# Memory cleanup: Section 15 large intermediates ---
rm(poly_data, poly_with_exp, poly_agg,
   pred_positions, pred_ramp, plot_ramp_poly,
   pred_grid_exp, pred_exp, plot_ramp_exp_poly,
   sampled_genes_15, all_genes_15, n_subsample_genes, window_size)
gc()
