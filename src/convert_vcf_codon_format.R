#' Convert Codon VCF Format to Gamma Estimation Input
#' 
#' Transforms codon variant data from your VCF processing pipeline
#' into the format required by estimate_gamma_by_gene_with_neutral_params()
#' 
#' @author Luis Javier Madrigal-Roca
#' _____________________________________________________________________________

prepare_vcf_for_gamma_estimation <- function(vcf_codon_dt, genetic_code_df) {
  #' Convert VCF codon data to format needed for gamma estimation
  #' 
  #' @param vcf_codon_dt data.table with columns:
  #'   - Gene
  #'   - Codon_Pos
  #'   - AA
  #'   - Preferred_Codon
  #'   - Codon_Variants (format: "AAC:0;AAT:187")
  #' @param genetic_code_df data.table with columns: Codon, AA
  #'   Mapping of codons to amino acids for filtering synonymous variants
  #' 
  #' @return data.table with columns:
  #'   - Gene
  #'   - AA
  #'   - Preferred_Codon
  #'   - k (count of preferred alleles - SYNONYMOUS ONLY)
  #'   - n (total sample size - SYNONYMOUS ONLY)
  #'   - p (frequency of preferred)
  #' 
  #' @details
  #' The function parses Codon_Variants, filters out non-synonymous changes,
  #' then calculates k and n using only synonymous variants.
  #' This is critical because non-synonymous variants would bias the sample size.
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(vcf_codon_dt)) setDT(vcf_codon_dt)
  if (!is.data.table(genetic_code_df)) setDT(genetic_code_df)
  
  cat("Preparing VCF codon data for gamma estimation...\n")
  cat("Filtering for SYNONYMOUS variants only...\n\n")
  
  # Create codon-to-AA lookup
  codon_to_aa <- setNames(genetic_code_df$AA, genetic_code_df$Codon)
  
  # Expand VCF to long format (one row per variant)
  dt_expanded <- vcf_codon_dt[, {
    # Parse variant string: "AAC:0;AAT:187"
    variant_parts <- unlist(strsplit(Codon_Variants, ";", fixed = TRUE))
    
    # Extract codon and count pairs
    codon_list <- character()
    count_list <- integer()
    
    for (v in variant_parts) {
      parts <- strsplit(v, ":", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        codon_list <- c(codon_list, parts[1])
        count_list <- c(count_list, as.integer(parts[2]))
      }
    }
    
    list(
      Variant_Codon = codon_list,
      Count = count_list
    )
  }, by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
  
  # Add amino acid for each variant codon
  dt_expanded[, Variant_AA := codon_to_aa[Variant_Codon]]
  
  # Filter: Keep only synonymous variants (Variant_AA == Site_AA)
  dt_syn <- dt_expanded[Variant_AA == AA]
  
  cat(sprintf("Total variant observations: %d\n", nrow(dt_expanded)))
  cat(sprintf("Synonymous variants: %d (%.1f%%)\n", 
              nrow(dt_syn), 
              100 * nrow(dt_syn) / nrow(dt_expanded)))
  cat(sprintf("Non-synonymous filtered out: %d\n\n", 
              nrow(dt_expanded) - nrow(dt_syn)))
  
  # Aggregate: Calculate k and n from synonymous variants only
  # CRITICAL: Group by Codon_Pos to preserve independent sites
  # DO NOT collapse across positions - that destroys variance!
  result <- dt_syn[, {
    
    # Skip Met and Trp (single codon)
    if (AA[1] %in% c("M", "W")) {
      list(k = NA_integer_, n = NA_integer_, p = NA_real_)
    } else {
      
      pref_codon <- Preferred_Codon[1]
      
      # Count of preferred codon (synonymous only)
      k_count <- sum(Count[Variant_Codon == pref_codon])
      
      # Total sample size (synonymous only)
      total_n <- sum(Count)
      
      # Calculate frequency
      freq_pref <- if (total_n > 0) k_count / total_n else 0
      
      list(
        k = k_count,
        n = total_n,
        p = freq_pref
      )
    }
  }, by = .(Gene, Codon_Pos, AA, Preferred_Codon)]  # <-- FIXED: Added Codon_Pos]
  
  # Remove sites with no data or monomorphic
  result <- result[!is.na(k) & n > 0]
  result <- result[k > 0 & k < n]  # Remove fixed sites (need polymorphism)
  
  # Quality control
  cat("=== Quality Control ===\n")
  cat(sprintf("Total polymorphic sites (one row = one codon position): %d\n", nrow(result)))
  cat(sprintf("Sites removed (Met/Trp): %d\n", 
              sum(vcf_codon_dt$AA %in% c("M", "W"))))
  cat(sprintf("Sites removed (monomorphic or fixed): %d\n\n", 
              nrow(vcf_codon_dt) - nrow(result) - sum(vcf_codon_dt$AA %in% c("M", "W"))))
  
  # CRITICAL VALIDATION: Check that we have MULTIPLE sites per Gene×AA
  sites_per_gene_aa <- result[, .N, by = .(Gene, AA)]
  
  cat("=== CRITICAL VALIDATION ===\n")
  cat(sprintf("Mean sites per Gene×AA: %.1f\n", mean(sites_per_gene_aa$N)))
  cat(sprintf("Median sites per Gene×AA: %.0f\n", median(sites_per_gene_aa$N)))
  cat(sprintf("Gene×AA with single site: %d (%.1f%% - BAD if high!)\n",
              sum(sites_per_gene_aa$N == 1),
              100 * mean(sites_per_gene_aa$N == 1)))
  cat(sprintf("Gene×AA with ≥5 sites: %d (%.1f%% - GOOD!)\n\n",
              sum(sites_per_gene_aa$N >= 5),
              100 * mean(sites_per_gene_aa$N >= 5)))
  
  if (mean(sites_per_gene_aa$N) < 2) {
    warning("⚠️  CRITICAL: Most Gene×AA have only 1 site! Check Codon_Pos grouping!")
  }
  
  # Check sample sizes (should be ~187 for inbred lines, not diploid count)
  cat("Sample size distribution (n = homozygous genotypes per site):\n")
  cat(sprintf("  Mean: %.1f\n", mean(result$n)))
  cat(sprintf("  Median: %.0f\n", median(result$n)))
  cat(sprintf("  Range: %d - %d\n\n", min(result$n), max(result$n)))
  
  if (mean(result$n) > 500) {
    warning("⚠️  CRITICAL: Sample sizes are too large! You may be summing across sites!")
  }
  
  if (mean(result$n) > 250) {
    warning("⚠️  WARNING: Sample sizes suggest diploid counting. For inbred lines, expect n≈187.")
  }
  
  # Check amino acid coverage
  cat("Sites per amino acid:\n")
  aa_counts <- result[, .N, by = AA][order(-N)]
  print(aa_counts)
  cat("\n")
  
  return(result)
}


add_mutation_rates_to_vcf <- function(vcf_prepared_dt, aa_mut_rates_dt) {
  #' Add mutation rate columns to prepared VCF data
  #' 
  #' @param vcf_prepared_dt Output from prepare_vcf_for_gamma_estimation()
  #' @param aa_mut_rates_dt data.table with columns: AA, u, v
  #' @return Merged data.table
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(vcf_prepared_dt)) setDT(vcf_prepared_dt)
  if (!is.data.table(aa_mut_rates_dt)) setDT(aa_mut_rates_dt)
  
  setkey(vcf_prepared_dt, AA)
  setkey(aa_mut_rates_dt, AA)
  
  result <- aa_mut_rates_dt[vcf_prepared_dt, nomatch = 0]
  
  cat(sprintf("Merged %d sites with mutation rates\n", nrow(result)))
  
  return(result)
}


validate_vcf_format <- function(vcf_codon_dt) {
  #' Check if VCF codon data has expected format
  #' 
  #' @param vcf_codon_dt Input data.table
  #' @return TRUE if valid, stops with error message if not
  #' ___________________________________________________________________________
  
  required_cols <- c("Gene", "Codon_Pos", "AA", "Preferred_Codon", 
                     "Codon_Variants")
  
  missing_cols <- setdiff(required_cols, names(vcf_codon_dt))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", 
                 paste(missing_cols, collapse = ", ")))
  }
  
  # Check Codon_Variants format
  sample_variant <- vcf_codon_dt$Codon_Variants[1]
  
  if (!grepl(":", sample_variant)) {
    stop("Codon_Variants column should have format 'CODON:COUNT;CODON:COUNT'")
  }
  
  cat("✓ VCF format validation passed\n")
  return(TRUE)
}


# Example usage workflow:
# 
# # 1. Load genetic code
# genetic_code_df <- data.table(
#   Codon = names(Biostrings::GENETIC_CODE),
#   AA = as.character(Biostrings::GENETIC_CODE)
# )
# 
# # 2. Validate format
# validate_vcf_format(vcf_codon)
# 
# # 3. Convert to gamma estimation format (WITH synonymous filtering)
# vcf_prepared <- prepare_vcf_for_gamma_estimation(vcf_codon, genetic_code_df)
# 
# # 4. Estimate gamma
# gamma_results <- estimate_gamma_by_gene_with_neutral_params(
#   codon_vcf_data = vcf_prepared,
#   neutral_params = neutral_params,
#   preferred_codons_df = preferred_codons
# )
