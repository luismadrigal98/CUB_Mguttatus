#' GC-biased Gene Conversion (gBGC) Diagnostic Analysis
#'
#' Separates the effect of gBGC from translational selection on codon usage
#' by comparing polymorphism patterns at W↔S vs S↔S synonymous sites.
#'
#' Key insight: gBGC only affects W↔S mutations (A/T ↔ G/C), pushing G/C
#' alleles to higher frequency. S↔S mutations (G ↔ C) are immune to gBGC.
#' By comparing selection signals between these classes, we can:
#'   - Quantify gBGC contribution to observed SFS patterns
#'   - Isolate genuine CUB selection from gBGC contamination
#'
#' @references
#' Galtier, N. et al. (2009) Trends Genet 25:1-5.
#' Ratnakumar, A. et al. (2010) Genome Res 20:1001-1009.
#' Rousselle, M. et al. (2019) Mol Biol Evol 36:1092-1103.
#'
#' @author Generated for codon bias analysis in Mimulus guttatus
# ______________________________________________________________________________

#' Classify polymorphic synonymous sites by mutation strength class
#'
#' Parses the codon_frequencies_preferred.txt file and classifies each
#' polymorphic 3rd-position site as W↔S, S↔S, or W↔W based on the
#' segregating nucleotides at codon position 3.
#'
#' @param gene_list Character vector of gene names to include
#' @param codon_freq_file Path to codon_frequencies_preferred.txt
#' @param target_n Sample size for SFS projection
#' @return List with:
#'   \describe{
#'     \item{WS}{data.table with n, k (count of S allele), count for W↔S sites}
#'     \item{SS}{data.table with n, k (count of ROC-preferred allele), count for S↔S sites}
#'     \item{summary}{data.frame with site counts by class}
#'   }
classify_sites_by_mutation_class <- function(
    gene_list,
    codon_freq_file = "./data/all_chromosomes.codon_frequencies_preferred.txt",
    target_n = 90
) {
  require(data.table)
  require(stringr)
  
  cat("\n=== Classifying Sites by Mutation Strength Class ===\n")
  cat(sprintf("Genes: %d | Target n: %d\n", length(gene_list), target_n))
  
  # Read data
  raw <- fread(codon_freq_file, header = TRUE, sep = "\t") |>
    _[, Gene := paste0("MgIM767.", Gene)] |>
    _[Gene %in% gene_list]
  
  cat(sprintf("Sites in gene set: %s\n", format(nrow(raw), big.mark = ",")))
  
  # Classify S (strong: G,C) and W (weak: A,T) bases
  is_strong <- function(base) base %in% c("G", "C")
  is_weak   <- function(base) base %in% c("A", "T")
  
  # ---- Parse each site ----
  # Pre-allocate for speed
  n_rows <- nrow(raw)
  
  site_class   <- character(n_rows)
  site_n       <- integer(n_rows)
  site_k       <- integer(n_rows)  # k depends on class (see below)
  site_valid   <- logical(n_rows)
  site_pref_base <- character(n_rows)
  
  cat("Parsing variants and classifying sites...\n")
  
  for (i in seq_len(n_rows)) {
    variants_str <- raw$Codon_Variants[i]
    pref_codon   <- raw$Preferred_Codon[i]
    pref_base    <- substr(pref_codon, 3, 3)
    
    site_pref_base[i] <- pref_base
    
    # Parse "AAC:150;AAT:37;NNN:3" into codons and counts
    entries <- strsplit(variants_str, ";")[[1]]
    codons  <- sub(":.*", "", entries)
    counts  <- as.integer(sub(".*:", "", entries))
    
    # Filter out NNN and invalid codons
    valid_mask <- codons != "NNN" & nchar(codons) == 3
    codons <- codons[valid_mask]
    counts <- counts[valid_mask]
    
    if (length(codons) < 1) {
      site_valid[i] <- FALSE
      next
    }
    
    # Extract 3rd-position bases
    bases_3 <- substr(codons, 3, 3)
    unique_bases <- unique(bases_3)
    
    total_n <- sum(counts)
    
    # Must be polymorphic (>=2 distinct bases at position 3)
    if (length(unique_bases) < 2 || total_n == 0) {
      site_valid[i] <- FALSE
      next
    }
    
    site_valid[i] <- TRUE
    site_n[i] <- total_n
    
    # Classify by strength of segregating alleles
    n_strong <- sum(is_strong(unique_bases))
    n_weak   <- sum(is_weak(unique_bases))
    
    if (n_strong > 0 && n_weak > 0) {
      # W↔S site: gBGC affected
      site_class[i] <- "WS"
      # k = total count of S (G/C) alleles
      s_mask <- is_strong(bases_3)
      site_k[i] <- sum(counts[s_mask])
      
    } else if (n_strong >= 2 && n_weak == 0) {
      # S↔S site (G↔C): gBGC immune
      site_class[i] <- "SS"
      # k = count of the ROC-preferred allele
      pref_mask <- (codons == pref_codon)
      if (any(pref_mask)) {
        site_k[i] <- sum(counts[pref_mask])
      } else {
        # Preferred codon not segregating; use C-ending as default orientation
        c_mask <- (bases_3 == "C")
        site_k[i] <- sum(counts[c_mask])
      }
      
    } else if (n_weak >= 2 && n_strong == 0) {
      # W↔W site (A↔T): gBGC immune
      site_class[i] <- "WW"
      # k = count of preferred allele
      pref_mask <- (codons == pref_codon)
      if (any(pref_mask)) {
        site_k[i] <- sum(counts[pref_mask])
      } else {
        site_k[i] <- sum(counts[bases_3 == "A"])
      }
    }
  }
  
  # Build filtered data.table
  dt <- data.table(
    class = site_class[site_valid],
    n     = site_n[site_valid],
    k     = site_k[site_valid],
    pref_base = site_pref_base[site_valid]
  )
  
  # Summary
  class_counts <- dt[, .N, by = class]
  cat("\nSite classification:\n")
  for (r in seq_len(nrow(class_counts))) {
    cat(sprintf("  %s: %s sites\n", 
                class_counts$class[r],
                format(class_counts$N[r], big.mark = ",")))
  }
  
  # Build SFS summary tables (grouped (n,k) → count)
  WS_sfs <- dt[class == "WS", .(.N), by = .(n, k)]
  setnames(WS_sfs, "N", "count")
  
  SS_sfs <- dt[class == "SS", .(.N), by = .(n, k)]
  setnames(SS_sfs, "N", "count")
  
  WW_sfs <- dt[class == "WW", .(.N), by = .(n, k)]
  if (nrow(WW_sfs) > 0) setnames(WW_sfs, "N", "count")
  
  # Compute mean preferred frequency for quick diagnostic
  dt[, freq := k / n]
  mean_freqs <- dt[, .(mean_pref_freq = mean(freq, na.rm = TRUE),
                        median_pref_freq = median(freq, na.rm = TRUE),
                        n_sites = .N), by = class]
  
  cat("\nMean preferred-allele frequency by class:\n")
  print(mean_freqs)
  
  return(list(
    WS = WS_sfs,
    SS = SS_sfs,
    WW = WW_sfs,
    summary = class_counts,
    mean_freqs = mean_freqs,
    site_data = dt
  ))
}


#' Run gBGC diagnostic across expression quantiles
#'
#' For each expression quantile, classifies 3rd-position sites by mutation
#' class and compares the preferred-allele frequency spectrum between
#' W↔S (gBGC-affected) and S↔S (gBGC-immune) sites.
#'
#' @param integrated_data Data frame with Gene_name and Geom_Mean_CPM
#' @param codon_freq_file Path to codon_frequencies_preferred.txt
#' @param target_n Projection sample size
#' @param n_quantiles Number of expression quantiles (default 20)
#' @param neutral_params List with alpha_G, beta_G, alpha_C, beta_C from introns
#' @return List with comparison data frame and diagnostic plots
run_gbgc_diagnostic <- function(
    integrated_data,
    codon_freq_file = "./data/all_chromosomes.codon_frequencies_preferred.txt",
    target_n = 90,
    n_quantiles = 20,
    neutral_params = NULL
) {
  require(data.table)
  require(ggplot2)
  require(future.apply)
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("gBGC DIAGNOSTIC: W↔S vs S↔S COMPARISON\n")
  cat(strrep("=", 80), "\n\n")
  
  cat("Logic: If gBGC inflates the SFS signal:\n")
  cat("  → W↔S sites (gBGC-affected) will have elevated preferred-allele freq\n")
  cat("  → S↔S sites (gBGC-immune) will show ONLY genuine selection signal\n")
  cat("  → The W↔S excess over S↔S quantifies gBGC contribution\n\n")
  
  # ---- Step 1: Read ALL variant data once (expensive I/O) ----
  cat("Step 1: Loading full variant data...\n")
  raw_all <- fread(codon_freq_file, header = TRUE, sep = "\t")
  raw_all[, Gene := paste0("MgIM767.", Gene)]
  cat(sprintf("  Total sites loaded: %s\n", format(nrow(raw_all), big.mark = ",")))
  
  # ---- Step 2: Pre-parse all sites ----
  cat("Step 2: Pre-parsing variant data...\n")
  
  parsed <- parse_sites_for_gbgc(raw_all)
  
  cat(sprintf("  Polymorphic sites: %s\n", format(nrow(parsed), big.mark = ",")))
  cat(sprintf("  W↔S: %s | S↔S: %s | W↔W: %s\n",
              format(sum(parsed$class == "WS"), big.mark = ","),
              format(sum(parsed$class == "SS"), big.mark = ","),
              format(sum(parsed$class == "WW"), big.mark = ",")))
  
  # ---- Step 3: Compute per-quantile statistics ----
  cat("\nStep 3: Computing statistics by expression quantile...\n")
  
  cutoffs <- quantile(integrated_data$Geom_Mean_CPM, 
                      seq(0, 1, length.out = n_quantiles + 1), na.rm = TRUE)
  
  quantile_results <- lapply(seq_len(n_quantiles), function(idx) {
    low  <- ifelse(idx == 1, cutoffs[idx] - 0.001, cutoffs[idx])
    high <- cutoffs[idx + 1]
    
    gene_set <- integrated_data |>
      dplyr::filter(Geom_Mean_CPM > low & Geom_Mean_CPM <= high) |>
      dplyr::pull(Gene_name)
    
    sub_parsed <- parsed[Gene %in% gene_set]
    
    if (nrow(sub_parsed) == 0) {
      return(data.frame(
        quantile = idx, mean_exp = NA,
        mean_freq_WS = NA, mean_freq_SS = NA,
        n_WS = 0, n_SS = 0,
        delta_freq = NA
      ))
    }
    
    ws_data <- sub_parsed[class == "WS"]
    ss_data <- sub_parsed[class == "SS"]
    
    data.frame(
      quantile   = idx,
      mean_exp   = mean(integrated_data$Max_Log10_Exp[
        integrated_data$Gene_name %in% gene_set], na.rm = TRUE),
      mean_freq_WS = if (nrow(ws_data) > 0) mean(ws_data$freq, na.rm = TRUE) else NA,
      mean_freq_SS = if (nrow(ss_data) > 0) mean(ss_data$freq, na.rm = TRUE) else NA,
      median_freq_WS = if (nrow(ws_data) > 0) median(ws_data$freq, na.rm = TRUE) else NA,
      median_freq_SS = if (nrow(ss_data) > 0) median(ss_data$freq, na.rm = TRUE) else NA,
      n_WS = nrow(ws_data),
      n_SS = nrow(ss_data)
    )
  })
  
  comparison_df <- do.call(rbind, quantile_results)
  comparison_df$delta_freq <- comparison_df$mean_freq_WS - comparison_df$mean_freq_SS
  
  cat("\n=== Per-Quantile Comparison ===\n")
  print(comparison_df[, c("quantile", "mean_exp", "mean_freq_WS", 
                           "mean_freq_SS", "delta_freq", "n_WS", "n_SS")])
  
  # ---- Step 4: Estimate gamma for W↔S and S↔S separately ----
  gamma_results <- NULL
  
  if (!is.null(neutral_params)) {
    cat("\nStep 4: Estimating gamma by mutation class...\n")
    
    # For W↔S: use existing neutral params (alpha_C, beta_C averaged with alpha_G, beta_G)
    # For S↔S: approximate neutral params as symmetric-ish G↔C
    # Best approximation: alpha_SS ≈ mean(alpha_G, alpha_C), beta_SS ≈ mean(beta_G, beta_C)
    alpha_WS <- (neutral_params$alpha_C + neutral_params$alpha_G) / 2
    beta_WS  <- (neutral_params$beta_C + neutral_params$beta_G) / 2
    
    # For S↔S, G→C and C→G should be roughly symmetric transversions
    # Use the same averaged params as a first approximation
    alpha_SS <- alpha_WS
    beta_SS  <- beta_WS
    
    cat(sprintf("  W↔S neutral params: alpha=%.6f, beta=%.6f\n", alpha_WS, beta_WS))
    cat(sprintf("  S↔S neutral params: alpha=%.6f, beta=%.6f (averaged approx)\n", 
                alpha_SS, beta_SS))
    
    gamma_by_quantile <- lapply(seq_len(n_quantiles), function(idx) {
      low  <- ifelse(idx == 1, cutoffs[idx] - 0.001, cutoffs[idx])
      high <- cutoffs[idx + 1]
      
      gene_set <- integrated_data |>
        dplyr::filter(Geom_Mean_CPM > low & Geom_Mean_CPM <= high) |>
        dplyr::pull(Gene_name)
      
      sub_parsed <- parsed[Gene %in% gene_set]
      
      ws_sub <- sub_parsed[class == "WS"]
      ss_sub <- sub_parsed[class == "SS"]
      
      # Build (n,k,count) tables for SFS projection
      ws_sfs_tab <- ws_sub[, .(.N), by = .(n, k)]
      setnames(ws_sfs_tab, "N", "count")
      
      ss_sfs_tab <- ss_sub[, .(.N), by = .(n, k)]
      setnames(ss_sfs_tab, "N", "count")
      
      gamma_WS <- tryCatch({
        obs_ws <- project_sfs(ws_sfs_tab, target_n)
        estimate_gamma(obs_ws, 
                       list(alpha_C = alpha_WS, beta_C = beta_WS,
                            alpha_G = alpha_WS, beta_G = beta_WS), 
                       "C", target_n)
      }, error = function(e) list(gamma = NA, p_value = NA))
      
      gamma_SS <- tryCatch({
        obs_ss <- project_sfs(ss_sfs_tab, target_n)
        estimate_gamma(obs_ss,
                       list(alpha_C = alpha_SS, beta_C = beta_SS,
                            alpha_G = alpha_SS, beta_G = beta_SS),
                       "C", target_n)
      }, error = function(e) list(gamma = NA, p_value = NA))
      
      data.frame(
        quantile = idx,
        gamma_WS = gamma_WS$gamma,
        gamma_WS_p = gamma_WS$p_value,
        gamma_SS = gamma_SS$gamma,
        gamma_SS_p = gamma_SS$p_value
      )
    })
    
    gamma_results <- do.call(rbind, gamma_by_quantile)
    gamma_results$delta_gamma <- gamma_results$gamma_WS - gamma_results$gamma_SS
    
    comparison_df <- merge(comparison_df, gamma_results, by = "quantile")
    
    cat("\n=== Gamma Estimates by Mutation Class ===\n")
    print(gamma_results[, c("quantile", "gamma_WS", "gamma_SS", "delta_gamma")])
    
    cat(sprintf("\nMean gamma_WS: %.3f | Mean gamma_SS: %.3f | Mean delta (gBGC): %.3f\n",
                mean(gamma_results$gamma_WS, na.rm = TRUE),
                mean(gamma_results$gamma_SS, na.rm = TRUE),
                mean(gamma_results$delta_gamma, na.rm = TRUE)))
  }
  
  # ---- Step 5: Diagnostic plots ----
  cat("\nStep 5: Creating diagnostic plots...\n")
  
  # Plot A: Mean preferred frequency by class × expression
  plot_long <- comparison_df |>
    tidyr::pivot_longer(
      cols = c(mean_freq_WS, mean_freq_SS),
      names_to = "Class",
      values_to = "Mean_Pref_Freq"
    ) |>
    dplyr::mutate(
      Class = dplyr::case_when(
        Class == "mean_freq_WS" ~ "W↔S (gBGC affected)",
        Class == "mean_freq_SS" ~ "S↔S (gBGC immune)"
      )
    )
  
  p_freq_comparison <- ggplot(plot_long, 
                               aes(x = mean_exp, y = Mean_Pref_Freq, 
                                   color = Class, shape = Class)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("W↔S (gBGC affected)" = "#E41A1C",
                                   "S↔S (gBGC immune)" = "#377EB8")) +
    labs(
      title = "gBGC Diagnostic: Preferred Allele Frequency by Mutation Class",
      subtitle = paste(
        "W↔S sites are subject to gBGC; S↔S (G↔C) sites are immune.",
        "\nThe gap between curves estimates the gBGC contribution.",
        "\nA rising S↔S trend indicates genuine CUB selection."
      ),
      x = "Mean Expression (Log10 CPM)",
      y = "Mean Preferred-Allele Frequency",
      color = "Mutation Class",
      shape = "Mutation Class"
    ) +
    theme_custom() +
    theme(legend.position = "bottom")
  
  # Plot B: Delta (W↔S - S↔S) across expression
  p_delta <- ggplot(comparison_df, aes(x = mean_exp, y = delta_freq)) +
    geom_point(size = 3, color = "#E41A1C") +
    geom_smooth(method = "loess", se = TRUE, color = "#E41A1C", fill = "#E41A1C",
                alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "gBGC Excess: W↔S minus S↔S Preferred Frequency",
      subtitle = paste(
        "Positive values = gBGC inflating preferred-allele frequency.",
        "\nConstant across expression = gBGC acts independently of selection.",
        "\nIncreasing = gBGC correlated with expression (e.g., via recombination)."
      ),
      x = "Mean Expression (Log10 CPM)",
      y = "Δ Preferred Frequency (W↔S − S↔S)"
    ) +
    theme_custom()
  
  # Plot C: Gamma comparison (if computed)
  p_gamma <- NULL
  if (!is.null(gamma_results)) {
    gamma_long <- comparison_df |>
      tidyr::pivot_longer(
        cols = c(gamma_WS, gamma_SS),
        names_to = "Class",
        values_to = "Gamma"
      ) |>
      dplyr::mutate(
        Class = dplyr::case_when(
          Class == "gamma_WS" ~ "W↔S (gBGC affected)",
          Class == "gamma_SS" ~ "S↔S (gBGC immune)"
        )
      )
    
    p_gamma <- ggplot(gamma_long,
                       aes(x = mean_exp, y = Gamma, color = Class, shape = Class)) +
      geom_point(size = 3) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("W↔S (gBGC affected)" = "#E41A1C",
                                     "S↔S (gBGC immune)" = "#377EB8")) +
      labs(
        title = "Selection Coefficient by Mutation Class",
        subtitle = paste(
          "gamma(W↔S) − gamma(S↔S) estimates the gBGC contribution.",
          "\ngamma(S↔S) represents genuine selection (CUB) free from gBGC."
        ),
        x = "Mean Expression (Log10 CPM)",
        y = "Estimated γ (4Nes)",
        color = "Mutation Class",
        shape = "Mutation Class"
      ) +
      theme_custom() +
      theme(legend.position = "bottom")
  }
  
  # ---- Step 6: Extreme group SFS comparison ----
  cat("\nStep 6: Building SFS for extreme expression groups...\n")
  
  bottom_genes <- integrated_data |>
    dplyr::filter(Geom_Mean_CPM <= cutoffs[2]) |>
    dplyr::pull(Gene_name)
  
  top_genes <- integrated_data |>
    dplyr::filter(Geom_Mean_CPM > cutoffs[n_quantiles]) |>
    dplyr::pull(Gene_name)
  
  # Build projected SFS for top and bottom, by class
  sfs_extremes <- list()
  for (grp_name in c("Bottom", "Top")) {
    genes <- if (grp_name == "Bottom") bottom_genes else top_genes
    sub <- parsed[Gene %in% genes]
    
    for (cls in c("WS", "SS")) {
      cls_sub <- sub[class == cls]
      tab <- cls_sub[, .(.N), by = .(n, k)]
      if (nrow(tab) > 0) {
        setnames(tab, "N", "count")
        proj <- project_sfs(tab, target_n)
        sfs_extremes[[paste0(grp_name, "_", cls)]] <- proj
      }
    }
  }
  
  # Build comparison plot: folded SFS for extremes × class
  freq_bins <- 0:target_n
  sfs_plot_data <- do.call(rbind, lapply(names(sfs_extremes), function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    data.frame(
      Freq_bin = freq_bins,
      Count = sfs_extremes[[nm]],
      Expression = parts[1],
      Mutation_Class = ifelse(parts[2] == "WS", "W↔S (gBGC affected)", 
                               "S↔S (gBGC immune)")
    )
  }))
  
  # Remove monomorphic bins (0 and target_n) for clearer polymorphism visualization
  sfs_plot_poly <- sfs_plot_data |>
    dplyr::filter(Freq_bin > 0, Freq_bin < target_n, Count > 0)
  
  p_sfs_class <- ggplot(sfs_plot_poly, 
                          aes(x = Freq_bin, y = Count + 1, 
                              fill = interaction(Expression, Mutation_Class))) +
    geom_col(position = "dodge", width = 0.8) +
    scale_y_log10() +
    facet_wrap(~ Mutation_Class, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Bottom.W↔S (gBGC affected)" = "#377EB8",
        "Top.W↔S (gBGC affected)" = "#E41A1C",
        "Bottom.S↔S (gBGC immune)" = "#7FCDBB",
        "Top.S↔S (gBGC immune)" = "#FC9272"
      ),
      labels = c("Bottom 5%", "Top 5%", "Bottom 5%", "Top 5%"),
      name = "Expression Group"
    ) +
    labs(
      title = "SFS by Mutation Class: Bottom 5% vs Top 5% Expression",
      subtitle = "Right-shift in S↔S = genuine selection | Right-shift only in W↔S = gBGC",
      x = "Preferred Allele Count",
      y = "Number of Sites (log scale)"
    ) +
    theme_custom() +
    theme(legend.position = "bottom")
  
  return(list(
    comparison = comparison_df,
    gamma_results = gamma_results,
    plots = list(
      freq_comparison = p_freq_comparison,
      delta = p_delta,
      gamma = p_gamma,
      sfs_by_class = p_sfs_class
    ),
    sfs_extremes = sfs_extremes
  ))
}


#' Parse all sites from codon frequency file for gBGC classification
#'
#' Internal helper that pre-parses the full variant file into a data.table
#' with columns: Gene, class (WS/SS/WW), n, k, freq.
#'
#' @param raw_dt data.table from fread of codon_frequencies_preferred.txt
#' @return data.table with parsed site data
parse_sites_for_gbgc <- function(raw_dt) {
  require(data.table)
  
  n_rows <- nrow(raw_dt)
  
  # Pre-allocate vectors
  out_gene  <- character(n_rows)
  out_class <- character(n_rows)
  out_n     <- integer(n_rows)
  out_k     <- integer(n_rows)
  out_valid <- logical(n_rows)
  
  is_S <- function(b) b == "G" | b == "C"
  is_W <- function(b) b == "A" | b == "T"
  
  for (i in seq_len(n_rows)) {
    variants_str <- raw_dt$Codon_Variants[i]
    pref_codon   <- raw_dt$Preferred_Codon[i]
    
    # Parse "AAC:150;AAT:37" format
    entries <- strsplit(variants_str, ";", fixed = TRUE)[[1]]
    codons  <- sub(":.*", "", entries)
    counts  <- as.integer(sub(".*:", "", entries))
    
    # Remove NNN / invalid
    ok <- codons != "NNN" & nchar(codons) == 3
    codons <- codons[ok]
    counts <- counts[ok]
    
    if (length(codons) < 2) {
      out_valid[i] <- FALSE
      next
    }
    
    bases_3 <- substr(codons, 3, 3)
    uniq_b  <- unique(bases_3)
    
    if (length(uniq_b) < 2) {
      out_valid[i] <- FALSE
      next
    }
    
    total <- sum(counts)
    if (total == 0) { out_valid[i] <- FALSE; next }
    
    out_valid[i] <- TRUE
    out_gene[i]  <- raw_dt$Gene[i]
    out_n[i]     <- total
    
    has_S <- any(is_S(uniq_b))
    has_W <- any(is_W(uniq_b))
    
    if (has_S && has_W) {
      # W↔S
      out_class[i] <- "WS"
      out_k[i] <- sum(counts[is_S(bases_3)])
      
    } else if (has_S && !has_W) {
      # S↔S (G↔C)
      out_class[i] <- "SS"
      # Orient by ROC-preferred codon
      pref_mask <- (codons == pref_codon)
      if (any(pref_mask)) {
        out_k[i] <- sum(counts[pref_mask])
      } else {
        # Preferred codon not segregating; orient by C allele
        out_k[i] <- sum(counts[bases_3 == "C"])
      }
      
    } else {
      # W↔W (A↔T)
      out_class[i] <- "WW"
      pref_mask <- (codons == pref_codon)
      if (any(pref_mask)) {
        out_k[i] <- sum(counts[pref_mask])
      } else {
        out_k[i] <- sum(counts[bases_3 == "A"])
      }
    }
  }
  
  dt <- data.table(
    Gene  = out_gene[out_valid],
    class = out_class[out_valid],
    n     = out_n[out_valid],
    k     = out_k[out_valid]
  )
  dt[, freq := k / n]
  
  return(dt)
}
