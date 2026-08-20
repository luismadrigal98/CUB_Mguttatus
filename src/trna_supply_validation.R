## ============================================================================
## trna_supply_validation.R
##
## Independent validation of the candidate optimal codons, from tRNA gene
## content alone.
##
## ---------------------------------------------------------------------------
## Why this is the primary validity check
## ---------------------------------------------------------------------------
## The preferred-codon call in src/detect_preferred_codons.R is estimated from
## codon counts against expression.  ROC-SEMPPR agrees with it (18/19 families),
## but that agreement is NOT independent: the AnaCoDa run those coefficients
## come from (results_dM_fixed_with_phi_final) was itself given the empirical
## expression data, so both estimators read the same expression-codon
## association.
##
## tRNA gene copy number is orthogonal to both.  It is genomic content: no
## expression data, no codon-usage data.  If the codons called optimal are the
## ones the translational machinery is actually built to read quickly, that is
## corroboration in the proper sense.
##
## ---------------------------------------------------------------------------
## The wobble problem, and why a naive mapping gets the wrong answer
## ---------------------------------------------------------------------------
## A tRNA does not read only its Watson-Crick codon.  Mapping each tRNA gene to
## revcomp(anticodon) and stopping there produces a badly misleading table: in
## four-fold families the C-ending codons appear to have ZERO cognate tRNA genes
## (e.g. Ala GCC:0 against GCT:20), which would look like evidence against the
## call.  It is an artifact.  Eukaryotes largely lack G34 tRNAs for four-fold
## families; the A34 gene is post-transcriptionally deaminated to inosine (I34),
## and I34 reads NNT, NNC and NNA.  A Watson-Crick-only tabulation dumps that
## tRNA's entire copy number onto the T-ending codon.
##
## Decoding rules applied here (standard eukaryotic wobble):
##
##   anticodon 34  reads codons        note
##   -----------   -----------------   ------------------------------------
##   A34 (-> I34)  NNT, NNC, NNA       inosine; assumed for every ANN gene
##   G34           NNC, NNT            G:U wobble
##   T34/U34       NNA, NNG            often modified (mnm5s2U etc.)
##   C34           NNG                 Watson-Crick only
##
## Counting which tRNAs *can* read a codon is still not enough, because wobble
## pairings are slower and less accurate than Watson-Crick ones.  Weighting all
## readers equally makes the two-fold NNY families come out exactly tied (their
## single G34 tRNA reads NNC and NNT alike), which destroys the comparison
## precisely where it should be cleanest.  So supply is computed as the standard
## tAI weight (dos Reis, Savva & Wernisch 2004, NAR 32:5036):
##
##   W(codon) = sum over decoding anticodons of (1 - s_pairing) * gene_copies
##
##   codon   Watson-Crick reader      wobble reader (penalty s)
##   -----   ---------------------    -------------------------
##   NNT     A34 (I:U),     s = 0     G34 (G:U),  s = 0.41
##   NNC     G34 (G:C),     s = 0     A34 (I:C),  s = 0.28
##   NNA     T34 (U:A),     s = 0     A34 (I:A),  s = 0.9999
##   NNG     C34 (C:G),     s = 0     T34 (U:G),  s = 0.68
##
## ---------------------------------------------------------------------------
## Read the output as SUPPORTED / NOT SUPPORTED, never as agree / disagree
## ---------------------------------------------------------------------------
## Whether a family can say anything at all is decided by its anticodon
## inventory, not by its codon counts.  For the C-vs-T contrast:
##
##   G34 absent (A34 only)  ->  W(NNC)/W(NNT) = 1 - s(I:C) = 0.720 EXACTLY,
##                              whatever the copy numbers.  The C-ending codon
##                              can never win.  NOT evidence against the call.
##   A34 absent (G34 only)  ->  ratio = 1/(1 - s(G:U)) = 1.695 exactly.  Here
##                              the evidence is the inventory itself: the genome
##                              encodes a dedicated Watson-Crick reader for the
##                              C-ending codon, and did not take the inosine
##                              route it uses elsewhere.  SUPPORTED (inventory).
##   both present           ->  ratio varies with the A34:G34 copy ratio, so the
##                              comparison is genuinely quantitative.
##
## In M. guttatus this splits 13 supported / 8 not supported / 0 contradicted.
## Ala, Arg_4, Ile, Thr and Val have zero G34 genes and sit at exactly 0.720;
## Leu_4, Pro and Ser_4 carry a single token G34 gene against 12-22 A34 genes
## (0.745-0.766).  Those eight are the inosine-dominated set: tRNA gene content
## does not support their C-ending call, and it is also the regime where the
## generic s(I:C) weight is least trustworthy, so neither reading is safe.
## Gly is the informative four-fold case — 35 G34 genes and no A34 — and it
## supports GGC.  The two-fold NNR families (Gln, Glu, Lys, Arg_2, Leu_2) are
## genuinely copy-number driven via C34 vs T34.
##
## Note this is where the paper must be careful: the "mostly C-ending" headline
## rests largely on four-fold families, and those are exactly the eight with no
## tRNA support.  Say so rather than reporting a bare 13/21.
##
## This is a known property of eukaryotic tRNA gene sets, not a quirk of this
## genome: the tRNA(GNN) species of the four-fold families are largely extinct
## in eukaryotes, and NNC codons there are read by I34.  Note also that generic
## tAI s-values are fitted to S. cerevisiae, and the literature holds I:C to be
## an EFFICIENT pairing, so s(I:C) = 0.28 probably understates it.  Species-
## specific weights (stAIcalc / gtAI) are the standard remedy and would be the
## right upgrade if the four-fold families need to be tested properly.
## ============================================================================

trna_supply_by_codon <- function(trna_file,
                                 genetic_code,
                                 wobble = TRUE,
                                 min_score = NA_real_) {
  #' Per-codon tRNA gene supply under eukaryotic wobble decoding
  #'
  #' @param trna_file tRNAscan-style table with `Anticodon` and (optionally)
  #'   `Score` columns — e.g.
  #'   `data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt`.
  #' @param genetic_code named character vector, codon -> family label.
  #' @param wobble if FALSE, use Watson-Crick pairing only (kept so the naive
  #'   tabulation can be reproduced and shown to be misleading).
  #' @param min_score optional tRNAscan score cutoff.
  #' @return data.table with Codon, Family, Third_base, `N_cognate_WC`
  #'   (Watson-Crick gene copies), `N_readers` (number of tRNA genes that can
  #'   read the codon at all) and `Supply` (the tAI weight: gene copies summed
  #'   with the dos Reis wobble penalties applied).

  suppressPackageStartupMessages({
    require(data.table)
    require(Biostrings)
  })

  tr <- data.table::fread(trna_file)
  if (!"Anticodon" %in% names(tr)) {
    stop("tRNA file must contain an `Anticodon` column; got: ",
         paste(names(tr), collapse = ", "))
  }
  tr[, Anticodon := toupper(gsub("U", "T", Anticodon))]
  tr <- tr[nchar(Anticodon) == 3L & !grepl("[^ACGT]", Anticodon)]
  if (!is.na(min_score) && "Score" %in% names(tr)) tr <- tr[Score >= min_score]

  # Watson-Crick codon of each anticodon
  tr[, Codon_WC := as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(Anticodon))
  )]

  # Which codons can each anticodon read, and with what tAI weight (1 - s)?
  # s-values from dos Reis, Savva & Wernisch (2004).
  codons_read <- function(anticodon, wc_codon) {
    base34 <- substr(anticodon, 1, 1)          # 5' base of the anticodon
    stem   <- substr(wc_codon, 1, 2)
    if (!wobble) return(data.frame(Codon = wc_codon, Weight = 1))
    switch(base34,
      # A34 -> I34: I:U (WC-like), I:C, I:A
      "A" = data.frame(Codon  = paste0(stem, c("T", "C", "A")),
                       Weight = c(1, 1 - 0.28, 1 - 0.9999)),
      # G34: G:C (WC), G:U wobble
      "G" = data.frame(Codon  = paste0(stem, c("C", "T")),
                       Weight = c(1, 1 - 0.41)),
      # U34: U:A (WC), U:G wobble
      "T" = data.frame(Codon  = paste0(stem, c("A", "G")),
                       Weight = c(1, 1 - 0.68)),
      # C34: C:G only
      "C" = data.frame(Codon = paste0(stem, "G"), Weight = 1),
      data.frame(Codon = wc_codon, Weight = 1)
    )
  }

  reader_rows <- data.table::rbindlist(lapply(seq_len(nrow(tr)), function(i) {
    cr <- codons_read(tr$Anticodon[i], tr$Codon_WC[i])
    data.table::data.table(
      Anticodon = tr$Anticodon[i],
      Codon     = cr$Codon,
      Weight    = cr$Weight
    )
  }))

  supply <- reader_rows[, .(N_readers = .N, Supply = sum(Weight)), by = Codon]
  wc     <- tr[, .(N_cognate_WC = .N), by = .(Codon = Codon_WC)]

  code <- genetic_code[genetic_code != "STOP"]
  out <- data.table::data.table(Codon = names(code), Family = unname(code))
  out <- merge(out, supply, by = "Codon", all.x = TRUE)
  out <- merge(out, wc,     by = "Codon", all.x = TRUE)
  out[is.na(N_readers),    N_readers := 0L]
  out[is.na(Supply),       Supply := 0]
  out[is.na(N_cognate_WC), N_cognate_WC := 0L]
  out[, Third_base := substr(Codon, 3, 3)]

  # Anticodon inventory per family, carried along because it — not the codon
  # counts — determines whether a family's C-vs-T contrast carries any
  # information (see header).
  tr[, Family := unname(genetic_code[Codon_WC])]
  tr[, Base34 := substr(Anticodon, 1, 1)]
  inv <- data.table::dcast(
    tr[!is.na(Family) & Family != "STOP", .N, by = .(Family, Base34)],
    Family ~ Base34, value.var = "N", fill = 0L
  )
  for (b in c("A", "C", "G", "T")) if (!b %in% names(inv)) inv[[b]] <- 0L
  data.table::setnames(inv, c("A", "C", "G", "T"),
                       c("n_A34", "n_C34", "n_G34", "n_T34"), skip_absent = TRUE)
  data.table::setattr(out, "anticodon_inventory",
                      inv[, .(Family, n_A34, n_C34, n_G34, n_T34)])
  out[]
}


validate_preferred_by_trna <- function(preferred, trna_supply, genetic_code) {
  #' Concordance between the preferred-codon call and tRNA supply
  #'
  #' @param preferred the `preferred` element of `detect_preferred_codons()`.
  #' @param trna_supply output of `trna_supply_by_codon()`.
  #' @param genetic_code named character vector, codon -> family label.
  #' @return list with `by_family` (per family: the called codon, whether it has
  #'   the highest tAI supply, and the margin over the runner-up),
  #'   `clean_subset` (the two-fold families, the clean test) and `summary`
  #'   counts split by family size — four-fold families are NOT a clean test,
  #'   see the header.

  suppressPackageStartupMessages(require(data.table))

  pref <- data.table::as.data.table(preferred)[, .(Family, Preferred_Codon)]
  sup  <- data.table::as.data.table(trna_supply)
  sup  <- sup[Family %in% pref$Family]

  by_family <- sup[order(Family, -Supply), {
    top <- .SD[1L]
    runner <- if (.N >= 2L) .SD[2L] else NULL
    list(
      Top_supply_codon = top$Codon,
      Top_supply       = top$Supply,
      Runner_supply    = if (is.null(runner)) NA_real_ else runner$Supply,
      Margin           = if (is.null(runner)) NA_real_ else top$Supply - runner$Supply
    )
  }, by = Family]

  by_family <- merge(pref, by_family, by = "Family")
  by_family <- merge(
    by_family,
    sup[, .(Family, Preferred_Codon = Codon,
            Preferred_supply = Supply, Preferred_WC = N_cognate_WC)],
    by = c("Family", "Preferred_Codon"), all.x = TRUE
  )
  by_family[, Tops_supply := Preferred_Codon == Top_supply_codon]

  fam_codons <- split(names(genetic_code), unname(genetic_code))
  by_family[, Family_size := vapply(Family, function(f) length(fam_codons[[f]]), 1L)]

  # Anticodon inventory per family — this is what decides whether a family can
  # carry evidence at all.
  inv <- attr(trna_supply, "anticodon_inventory")
  if (is.null(inv)) {
    stop("trna_supply must carry an `anticodon_inventory` attribute; ",
         "regenerate it with trna_supply_by_codon().")
  }
  by_family <- merge(by_family, inv, by = "Family", all.x = TRUE)
  for (cc in c("n_A34", "n_G34", "n_T34", "n_C34")) {
    by_family[is.na(get(cc)), (cc) := 0L]
  }

  # Classify the evidence.  The third base of the CALLED codon decides which
  # contrast matters.
  by_family[, Called_base := substr(Preferred_Codon, 3, 3)]
  # Exactly pinned: with no G34 gene at all, W(NNC)/W(NNT) collapses to
  # 1 - s(I:C) = 0.72 whatever the copy numbers, so no data enters.
  by_family[, Ratio_pinned := Called_base == "C" & n_G34 == 0L]

  by_family[, Evidence := data.table::fcase(
    # C-ending call in an inosine-dominated family. Includes both the exactly
    # pinned case (no G34) and families with a token single G34 gene against
    # 12-22 A34 genes, where the outcome is likewise not in doubt. These are
    # precisely the families where generic tAI s-values are least trustworthy:
    # s(I:C) = 0.28 was fitted to S. cerevisiae, and I:C is held in the
    # literature to be an efficient pairing. Not evidence against the call.
    Called_base == "C" & n_A34 >= n_G34,
      "not supported (inosine-dominated)",
    # Dedicated Watson-Crick reader, no inosine route: the inventory itself is
    # the evidence — the genome built a G34 reader instead of using I34.
    Called_base == "C" & n_G34 > 0L & n_A34 == 0L,
      "supported (inventory)",
    # Both routes materially represented, or a G-ending call decided by C34 vs
    # T34: genuinely quantitative.
    default = data.table::fifelse(Tops_supply, "supported (copy number)",
                                  "contradicted")
  )]

  by_family[, Evidence_basis := data.table::fcase(
    Evidence == "not supported (inosine-dominated)" & Ratio_pinned,
      sprintf("no G34 gene (A34=%d); W(NNC)/W(NNT) pinned at 1-s(I:C)=0.72 regardless of copy number", n_A34),
    Evidence == "not supported (inosine-dominated)",
      sprintf("inosine route dominates (A34=%d vs G34=%d); outcome not in doubt and s(I:C) is the least reliable tAI weight", n_A34, n_G34),
    Evidence == "supported (inventory)",
      sprintf("genome encodes a dedicated Watson-Crick G34 reader (%d copies) and no A34", n_G34),
    default =
      sprintf("copy-number driven (A34=%d G34=%d T34=%d C34=%d)", n_A34, n_G34, n_T34, n_C34)
  )]

  informative <- by_family[Evidence != "not supported (inosine-dominated)"]

  list(
    by_family     = by_family[],
    informative   = informative[],
    summary = list(
      n_families        = nrow(by_family),
      n_supported       = sum(grepl("^supported", by_family$Evidence)),
      n_contradicted    = sum(by_family$Evidence == "contradicted"),
      n_not_supported   = sum(by_family$Evidence == "not supported (inosine-dominated)"),
      n_supported_copy  = sum(by_family$Evidence == "supported (copy number)"),
      n_supported_inv   = sum(by_family$Evidence == "supported (inventory)"),
      n_ratio_pinned    = sum(by_family$Ratio_pinned),
      inosine_dominated_families =
        by_family$Family[by_family$Evidence == "not supported (inosine-dominated)"]
    )
  )
}


plot_trna_supply_vs_preference <- function(trna_supply, preferred,
                                           output_file = NULL) {
  #' Per-family tRNA supply, with the called optimal codon marked
  #'
  #' Reports the wobble-aware supply so the figure cannot be read as the naive
  #' Watson-Crick tabulation, which is misleading in four-fold families.

  suppressPackageStartupMessages({ require(ggplot2); require(data.table) })

  d <- data.table::as.data.table(trna_supply)
  pref <- data.table::as.data.table(preferred)[, .(Family, Codon = Preferred_Codon)]
  pref[, Is_preferred := TRUE]
  d <- merge(d, pref, by = c("Family", "Codon"), all.x = TRUE)
  d <- d[!is.na(Family)]
  d[is.na(Is_preferred), Is_preferred := FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(Codon, Supply),
                                       y = Supply, fill = Is_preferred)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ Family, scales = "free_y", ncol = 5) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1B7837", "FALSE" = "grey70"),
      labels = c("TRUE" = "Called optimal", "FALSE" = "Other synonym"),
      name = NULL) +
    ggplot2::labs(
      title = "tRNA gene supply per codon (tAI-weighted) vs the called optimal codon",
      subtitle = "tAI weight: gene copies summed with dos Reis (2004) wobble penalties. Two-fold families are the clean test; four-fold C-ending codons are inosine-decoded and penalised by construction.",
      x = NULL, y = "tRNA gene supply") +
    theme_custom() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 14, height = 11, dpi = 300)
  }
  p
}
