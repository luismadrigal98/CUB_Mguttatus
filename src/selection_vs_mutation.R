# ---------------------------------------------------------------------------
# Selection vs mutation at 4-fold sites (adversarial review, 2026-11-18/19)
# ---------------------------------------------------------------------------
# The C-vs-rest two-allele model treats G like A and T. The review showed that
# G is favoured too, that part of the C-content gradient is mutational, and that
# the divergence turnover test cannot exclude a mutation-bias mixture. These
# functions separate the pieces with statistics that do not depend on the
# preferred-codon call:
#   polarise_fourfold()          ancestral base = M. jungermannioides = M. glabratus
#   sfs_direction_slopes()       selection from the unfolded spectrum: logistic
#                                slope of P(derived allele is the strong one | DAF);
#                                mutation moves only the intercept (the slope is
#                                ~ S = 4Ns under PRF), shape changes cancel
#   decompose_composition()      log-odds of C and of G against W (A+T) at 4-fold
#                                sites, split into mutational log-odds (rare
#                                variants, DAF <= 5%) and selection (the rest);
#                                weak-mutation equilibrium: ln(X/W) = m_X + S_X
#   turnover_by_context()        Fig 7D split by CpG context, A/G/T and A<->T
#                                yardsticks, exact Wright predictions
#   pi_d_by_site_class()         pi/d for C, C/W, G/W, A/T, A/G/T and introns
#                                (fixed differences), the controls for pi/d
# Grouping: S_Wright bins or expression categories (0.4 wide, top pooled at 2.4).

.acgt <- c("A", "C", "G", "T")

.groups <- function(meta, by) {
  m <- data.table::as.data.table(meta)
  if (by == "S") {
    m[, grp := cut(S, c(-Inf, 0.25, 0.5, 0.75, 1, 1.5, Inf), right = FALSE)]
  } else {
    m[, grp := cut(Mean_Log10_Exp, c(-Inf, 0.2, 0.6, 1.0, 1.4, 1.8, 2.2, Inf),
                   labels = c("0.0", "0.4", "0.8", "1.2", "1.6", "2.0", "2.4+"))]
  }
  m[!is.na(grp), .(Gene, grp)]
}

polarise_fourfold <- function(pop_sites, div_sites) {
  #' Per 4-fold site: line counts, ancestral base (both outgroups agree) and, for
  #' biallelic sites carrying the ancestral base, the derived base and its frequency.
  x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T, n_alleles, next_base)],
             data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref, J = JUNG1_base, K = KONZA_base)],
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r]
  x[, anc := data.table::fifelse(J %chin% .acgt & J == K, J, NA_character_)]
  x[, ac := data.table::fcase(anc == "A", A, anc == "C", C, anc == "G", G, anc == "T", T)]
  x[, der := NA_character_]
  b2 <- x[, !is.na(anc) & n_alleles == 2L & ac > 0L & ac < n]
  x[b2, der := data.table::fcase(anc != "A" & A > 0L, "A", anc != "C" & C > 0L, "C",
                                 anc != "G" & G > 0L, "G", default = "T")]
  x[, daf := data.table::fifelse(is.na(der), NA_real_, (n - ac) / n)]
  x[, c("r", "J", "K") := NULL]
  x[]
}

# ---- selection from the unfolded spectrum -----------------------------------

.slope_pairs <- list(`W->C` = list(fav = c("A>C", "T>C"), unf = c("C>A", "C>T")),
                     `W->G` = list(fav = c("A>G", "T>G"), unf = c("G>A", "G>T")),
                     `G->C` = list(fav = "G>C", unf = "C>G"),
                     `T->A (control)` = list(fav = "T>A", unf = "A>T"))

sfs_direction_slopes <- function(pol, meta, by = c("S", "expression"), n_boot = 200L, seed = 1L) {
  #' Logistic slope of P(derived allele is the favoured one | DAF), per direction
  #' pair and group, with 95% intervals from resampling genes.
  by <- match.arg(by); set.seed(seed)
  s <- merge(pol[!is.na(der), .(Gene, dir = paste0(anc, ">", der), daf)], .groups(meta, by), by = "Gene")
  fit <- function(d) {
    if (nrow(d) < 200 || length(unique(d$y)) < 2) return(NA_real_)
    stats::glm.fit(cbind(1, d$daf), d$y, family = stats::binomial())$coefficients[2]
  }
  out <- list()
  for (k in names(.slope_pairs)) {
    pr <- .slope_pairs[[k]]
    d <- s[dir %chin% c(pr$fav, pr$unf)]
    d[, y := as.integer(dir %chin% pr$fav)]
    for (gname in levels(droplevels(d$grp))) {
      dg <- d[grp == gname]
      gi <- split(seq_len(nrow(dg)), dg$Gene)
      bt <- replicate(n_boot, fit(dg[unlist(gi[sample.int(length(gi), length(gi), TRUE)], use.names = FALSE)]))
      out[[length(out) + 1]] <- data.table::data.table(
        grouping = by, group = gname, pair = k, sites = nrow(dg), slope = fit(dg),
        lo = stats::quantile(bt, 0.025, na.rm = TRUE), hi = stats::quantile(bt, 0.975, na.rm = TRUE))
    }
  }
  r <- data.table::rbindlist(out)
  r[, group := factor(group, levels = levels(.groups(meta, by)$grp))]
  r[order(pair, group)]
}

# ---- composition split into mutation and selection -----------------------------

decompose_composition <- function(pol, meta, intron_comp, by = c("S", "expression"),
                                  rare_max = 0.05, n_boot = 500L, seed = 1L) {
  #' Per group: composition log-odds L_X = ln(X/W) (X = C or G; W = A+T; mean line
  #' frequencies), the intronic baseline ln(X0/W0), the mutational log-odds
  #' m_X = ln(rate W->X / rate X->W) from rare polarised variants (DAF <= rare_max,
  #' per ancestral site), and the selection component S_X = L_X - m_X. Also the
  #' C-over-G versions and the coding-strand asymmetries A->T/T->A and C->T/G->A.
  by <- match.arg(by); set.seed(seed)
  ic <- data.table::fread(intron_comp)
  c0 <- ic[, sum(q_pref_C * n_sites) / sum(n_sites)]; g0 <- ic[, sum(q_pref_G * n_sites) / sum(n_sites)]
  w0 <- 1 - c0 - g0
  x <- merge(pol, .groups(meta, by), by = "Gene")
  x[, rare := !is.na(der) & daf <= rare_max]
  g <- x[, .(grp = grp[1], fC = sum(C / n), fG = sum(G / n), fW = sum((A + T) / n),
             aA = sum(anc == "A", na.rm = TRUE), aC = sum(anc == "C", na.rm = TRUE),
             aG = sum(anc == "G", na.rm = TRUE), aT = sum(anc == "T", na.rm = TRUE),
             WC = sum(rare & anc %chin% c("A", "T") & der == "C"), CW = sum(rare & anc == "C" & der %chin% c("A", "T")),
             WG = sum(rare & anc %chin% c("A", "T") & der == "G"), GW = sum(rare & anc == "G" & der %chin% c("A", "T")),
             GC = sum(rare & anc == "G" & der == "C"), CG = sum(rare & anc == "C" & der == "G"),
             AT = sum(rare & anc == "A" & der == "T"), TA = sum(rare & anc == "T" & der == "A"),
             CT = sum(rare & anc == "C" & der == "T"), GA = sum(rare & anc == "G" & der == "A")), by = Gene]
  st <- function(s) {
    S <- function(v) sum(s[[v]])
    L_C <- log(S("fC") / S("fW")); L_G <- log(S("fG") / S("fW"))
    m_C <- log((S("WC") / (S("aA") + S("aT"))) / (S("CW") / S("aC")))
    m_G <- log((S("WG") / (S("aA") + S("aT"))) / (S("GW") / S("aG")))
    m_CG <- log((S("GC") / S("aG")) / (S("CG") / S("aC")))
    c(L_C = L_C, L_G = L_G, base_C = L_C - log(c0 / w0), base_G = L_G - log(g0 / w0),
      m_C = m_C, m_G = m_G, sel_C = L_C - m_C, sel_G = L_G - m_G,
      L_CG = L_C - L_G, m_CG = m_CG, sel_CG = (L_C - L_G) - m_CG,
      AT_over_TA = (S("AT") / S("aA")) / (S("TA") / S("aT")),
      CT_over_GA = (S("CT") / S("aC")) / (S("GA") / S("aG")))
  }
  lv <- levels(droplevels(g$grp))
  out <- list()
  for (k in lv) {
    s <- g[grp == k]
    e <- st(s)
    b <- replicate(n_boot, st(s[sample.int(nrow(s), nrow(s), TRUE)]))
    out[[k]] <- data.table::data.table(grouping = by, group = k, genes = nrow(s), stat = names(e), est = e,
                                       lo = apply(b, 1, stats::quantile, 0.025, na.rm = TRUE),
                                       hi = apply(b, 1, stats::quantile, 0.975, na.rm = TRUE))
  }
  r <- data.table::rbindlist(out)
  r[, group := factor(group, levels = lv)]
  attr(r, "intron_baseline") <- c(C = c0, G = g0, W = w0)
  r[]
}

# ---- Fig 7D by CpG context -----------------------------------------------------

turnover_by_context <- function(pop_sites, div_sites, meta, U, V, taxon, t = NULL,
                                S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                                min_sites = 20L, n_boot = 500L, seed = 1L) {
  #' The turnover test for all 4-fold sites, for those whose next base is G (C there
  #' forms a CpG) and for the rest, against the A/G/T and A<->T yardsticks. Wright
  #' (exact, fixed differences) is read at the S implied by each context's own
  #' measured C content, scaled to the same S = 0 reference as the data.
  set.seed(seed)
  grid <- seq(-6, 9, by = 0.002); qg <- .wright_Q_kummer(grid, U, V)
  invS <- function(q) stats::approx(qg, grid, xout = q, ties = mean)$y
  x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T, next_base)],
             data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref, o = get(paste0(taxon, "_base")))],
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r & o %chin% .acgt]
  x[, `:=`(oc = data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T), noC = C == 0L,
           W = C == 0L & G == 0L)]
  m <- data.table::as.data.table(meta)[, .(Gene, S)]
  res <- list()
  for (ctx in c("all", "CpG (next base G)", "non-CpG")) {
    y <- if (ctx == "all") x else if (startsWith(ctx, "CpG")) x[next_base == "G"] else x[next_base != "G"]
    g <- y[, .(n4 = .N, QIM = sum(C / n), Qo = sum(o == "C"),
               dC = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
               nN = sum(noC & o != "C"), dN = sum(noC & o != "C" & oc == 0L),
               fA = sum((A / n)[noC & o != "C"]), fG = sum((G / n)[noC & o != "C"]), fT = sum((T / n)[noC & o != "C"]),
               oA = sum(noC & o == "A"), oG = sum(noC & o == "G"), oT = sum(noC & o == "T"),
               nW = sum(W & o %chin% c("A", "T")), dW = sum(W & o %chin% c("A", "T") & oc == 0L),
               wA = sum((A / n)[W & o %chin% c("A", "T")]), wT = sum((T / n)[W & o %chin% c("A", "T")]),
               owA = sum(W & o == "A"), owT = sum(W & o == "T")), by = Gene]
    g <- merge(g, m, by = "Gene")[is.finite(S) & n4 >= min_sites %/% 2]
    g[, Sbin := cut(S, S_breaks, right = FALSE)]
    bins <- levels(droplevels(g$Sbin)); idx <- lapply(bins, function(b) which(g$Sbin == b))
    mS <- vapply(idx, function(i) mean(g$S[i]), 0); w0 <- -mS[1] / (mS[2] - mS[1])
    xS <- vapply(idx, function(i) invS(sum(g$QIM[i]) / sum(g$n4[i])), 0)
    for (yard in c("A/G/T", "A<->T")) {
      st <- function(i) {
        s <- function(v) sum(g[[v]][i])
        QI <- s("QIM") / s("n4"); Qo <- s("Qo") / s("n4")
        tC <- (s("dC") / s("n4")) / (QI * (1 - Qo) + Qo * (1 - QI))
        nr <- if (yard == "A/G/T") {
          fr <- c(s("fA"), s("fG"), s("fT")) / s("nN"); fo <- c(s("oA"), s("oG"), s("oT")) / s("nN")
          (s("dN") / s("nN")) / (1 - sum(fr * fo))
        } else {
          fr <- c(s("wA"), s("wT")) / s("nW"); fo <- c(s("owA"), s("owT")) / s("nW")
          (s("dW") / s("nW")) / (1 - sum(fr * fo))
        }
        tC / nr
      }
      rel <- function(v) v / ((1 - w0) * v[1] + w0 * v[2])
      est <- rel(vapply(idx, st, 0))
      bt <- replicate(n_boot, rel(vapply(idx, function(i) st(i[sample.int(length(i), length(i), TRUE)]), 0)))
      ex <- if (is.null(t)) NA_real_ else rel(vapply(xS, function(s) wright_fixed_per_opportunity_exact(s, t, U, V), 0))
      res[[length(res) + 1]] <- data.table::data.table(
        Outgroup = taxon, context = ctx, yardstick = yard, Sbin = bins, n_genes = lengths(idx),
        S = mS, S_measured = xS, turnover = est, lo = apply(bt, 1, stats::quantile, 0.025),
        hi = apply(bt, 1, stats::quantile, 0.975), Wright_exact = ex)
    }
  }
  data.table::rbindlist(res)
}

# ---- pi/d by site class, the controls -------------------------------------------

pi_d_by_site_class <- function(pop_sites, div_sites, intron_comp, meta, taxon,
                               S_breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, Inf),
                               min_intron = 200L, n_boot = 500L, seed = 1L) {
  #' pi/d (fixed differences) by S bin, relative to S = 0, for the C system, the
  #' allele classes C/W, G/W and A/T (sites carrying only those bases), A/G/T, and
  #' the introns of the same genes (genes with >= min_intron aligned intronic
  #' sites). Under the C-only model only C-involving classes should rise; under
  #' strong-vs-weak selection C/W and G/W rise and A/T does not; a gene-level
  #' process (Ne, depth) moves introns with them.
  set.seed(seed)
  h <- function(k, n) 2 * k * (n - k) / (n * (n - 1))
  x <- merge(data.table::as.data.table(pop_sites)[, .(Gene, pos3, ref3, n, A, C, G, T)],
             data.table::as.data.table(div_sites)[, .(Gene, Pos, r = Ref, o = get(paste0(taxon, "_base")))],
             by.x = c("Gene", "pos3"), by.y = c("Gene", "Pos"))[ref3 == r & o %chin% .acgt]
  x[, oc := data.table::fcase(o == "A", A, o == "C", C, o == "G", G, o == "T", T)]
  x[, `:=`(cw = G == 0L & o != "G", gw = C == 0L & o != "C", at = C == 0L & G == 0L & o %chin% c("A", "T"))]
  g <- x[, .(n4 = .N, piC = sum(h(C, n)), dC = sum((o == "C" & C == 0L) | (o != "C" & C == n)),
             nCW = sum(cw), piCW = sum(h(C, n)[cw]), dCW = sum(cw & ((o == "C" & C == 0L) | (o != "C" & C == n))),
             nGW = sum(gw), piGW = sum(h(G, n)[gw]), dGW = sum(gw & ((o == "G" & G == 0L) | (o != "G" & G == n))),
             nAT = sum(at), piAT = sum(h(A, n)[at]), dAT = sum(at & oc == 0L),
             piAGT = sum((n / (n - 1) * (1 - ((A / n)^2 + (G / n)^2 + (T / n)^2)))[gw]),
             dAGT = sum(gw & oc == 0L)), by = Gene]
  im <- data.table::fread(intron_comp)
  if (!all(c("n_out", "fixed_sub") %in% names(im)))
    stop("[pi_d_by_site_class] ", intron_comp, " lacks fixed-difference columns")
  im <- im[, .(Gene, nI = n_out, nIp = n_sites, piI = pi_total * n_sites, dI = fixed_sub)]
  g <- merge(merge(g, im, by = "Gene"), data.table::as.data.table(meta)[, .(Gene, S)], by = "Gene")
  g <- g[is.finite(S) & n4 >= 20 & nAT >= 5 & nI >= min_intron]
  g[, Sbin := cut(S, S_breaks, right = FALSE)]
  st <- function(s) {
    r <- function(p, np, d, nd) (sum(s[[p]]) / sum(s[[np]])) / (sum(s[[d]]) / sum(s[[nd]]))
    c(`C system` = r("piC", "n4", "dC", "n4"), `C/W` = r("piCW", "nCW", "dCW", "nCW"),
      `G/W` = r("piGW", "nGW", "dGW", "nGW"), `A/T` = r("piAT", "nAT", "dAT", "nAT"),
      `A/G/T` = r("piAGT", "nGW", "dAGT", "nGW"), Introns = r("piI", "nIp", "dI", "nI"))
  }
  bins <- levels(droplevels(g$Sbin)); idx <- lapply(bins, function(b) which(g$Sbin == b))
  mS <- vapply(idx, function(i) mean(g$S[i]), 0); w0 <- -mS[1] / (mS[2] - mS[1])
  rel <- function(m) m / ((1 - w0) * m[, 1] + w0 * m[, 2])
  est <- rel(sapply(idx, function(i) st(g[i])))
  bt <- replicate(n_boot, rel(sapply(idx, function(i) st(g[i[sample.int(length(i), length(i), TRUE)]]))))
  data.table::rbindlist(lapply(rownames(est), function(k) data.table::data.table(
    Outgroup = taxon, class = k, Sbin = bins, n_genes = lengths(idx), S = mS, pi_over_d = est[k, ],
    lo = apply(bt[k, , ], 1, stats::quantile, 0.025, na.rm = TRUE),
    hi = apply(bt[k, , ], 1, stats::quantile, 0.975, na.rm = TRUE))))
}

# ---- plots -----------------------------------------------------------------

plot_sfs_slopes <- function(tab) {
  #' Selection from the unfolded spectrum by group; one line per direction pair.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)
  d[, x := as.numeric(group) + (as.numeric(factor(pair)) - 2.5) * 0.08]
  ggplot(d, aes(x, slope, colour = pair)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.06, linewidth = 0.45) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq_along(levels(d$group)), labels = levels(d$group)) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    labs(x = if (d$grouping[1] == "S") expression(S[Wright] * " bin") else
           expression("Expression level category (log"[10] * ")"),
         y = "Selection from the unfolded spectrum (slope)") +
    theme_custom() + theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 0.5))
}

plot_decomposition <- function(tab) {
  #' Composition log-odds against W split into mutation (rare variants) and
  #' selection, for C, G and C over G, by group.
  suppressPackageStartupMessages(require(ggplot2))
  d <- data.table::as.data.table(tab)[stat %chin% c("L_C", "m_C", "sel_C", "L_G", "m_G", "sel_G", "L_CG", "m_CG", "sel_CG")]
  d[, base := data.table::fcase(stat %chin% c("L_C", "m_C", "sel_C"), "C vs A/T",
                                stat %chin% c("L_G", "m_G", "sel_G"), "G vs A/T", default = "C vs G")]
  d[, part := data.table::fcase(startsWith(stat, "L_"), "Composition (log-odds)",
                                startsWith(stat, "m_"), "Mutation (rare variants)", default = "Selection (difference)")]
  d[, base := factor(base, levels = c("C vs A/T", "G vs A/T", "C vs G"))]
  ggplot(d, aes(as.numeric(group), est, colour = part)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.08, linewidth = 0.45) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.8) +
    facet_wrap(~ base, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(levels(d$group)), labels = levels(d$group)) +
    scale_colour_brewer(palette = "Set1", name = NULL) +
    labs(x = if (d$grouping[1] == "S") expression(S[Wright] * " bin") else
           expression("Expression level category (log"[10] * ")"), y = "Log-odds") +
    theme_custom() + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
}
