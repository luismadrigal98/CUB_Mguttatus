fit_bi_multinom_family_model <- function(family_name, genetic_code, 
                                         usage_dt, meta_dt, preferred_codons_df)
{
  #' Fits binomial or multinomial GAM for codon usage within an amino acid family
  #'
  #' - For 2-codon families: binomial model (preferred vs non-preferred)
  #' - For 3+ codon families: multinomial model (all codons, preferred as baseline)
  #' 
  #' Models codon choice as a function of:
  #' - High_exp_log2 (selection effect)
  #' - s(CDS_length_nt) (mutational bias - length)
  #' - s(GC3s) (mutational bias - GC content)
  #'
  #' @param family_name The name of the family (e.g., "Ala", "Leu_4", "Lys")
  #' @param genetic_code Named vector: names = codons, values = amino acids
  #' @param usage_dt The codon_usage data.table or data.frame
  #' @param meta_dt The integrated_data data.table or data.frame
  #' @param preferred_codons_df Data frame with Codon and Amino_Acid columns
  #' @return A data.frame with model coefficients
  
  suppressPackageStartupMessages({
    require(mgcv)
    require(dplyr)
  })
  
  # Get all synonymous codons for this family
  codons_in_family <- names(genetic_code[genetic_code == family_name])
  
  # Validate family has multiple codons
  if (length(codons_in_family) < 2) {
    return(data.frame(
      Family = family_name,
      Preferred_Codon = NA,
      Intercept = 0,
      Slope_Expression = 0,
      SE = NA,
      p_value = 1.0,
      Converged = FALSE
    ))
  }
  
  # Identify the preferred codon for this family
  # preferred_codons_df has columns: Codon, Amino_Acid, relative_adaptiveness
  family_preferred <- preferred_codons_df %>%
    dplyr::filter(Amino_Acid == family_name)
  
  if (nrow(family_preferred) == 0) {
    warning(sprintf("No preferred codon found for family %s", family_name))
    return(data.frame(
      Family = family_name,
      Preferred_Codon = NA,
      Slope_Expression = 0,
      SE = NA,
      p_value = 1.0,
      p_CDS_length = NA,
      p_GC3s = NA,
      Converged = FALSE,
      N_genes = 0
    ))
  }
  
  # Use the Codon column (already the corrected preferred codon)
  preferred_codon <- family_preferred$Codon[1]
  non_preferred_codons <- setdiff(codons_in_family, preferred_codon)
  K <- length(non_preferred_codons)  # Number of non-baseline codons
  
  # Prepare the data for this family
  family_data <- as.data.frame(meta_dt) %>%
    dplyr::select(Gene_name, High_exp_log2, CDS_length_nt, GC3s) %>%
    dplyr::left_join(
      as.data.frame(usage_dt) %>% dplyr::select(dplyr::all_of(c("Gene_name", codons_in_family))),
      by = "Gene_name"
    ) %>%
    dplyr::mutate(
      Total_Codons = rowSums(dplyr::across(dplyr::all_of(codons_in_family)), na.rm = TRUE)
    )
  
  # Clean data: remove genes with zero counts or missing values
  family_data_clean <- family_data %>%
    dplyr::filter(Total_Codons > 0,
                  !is.na(High_exp_log2),
                  !is.na(CDS_length_nt),
                  !is.na(GC3s))
  
  # Check if we have enough data
  if (nrow(family_data_clean) < 50) {
    return(data.frame(
      Family = family_name,
      Preferred_Codon = preferred_codon,
      N_codons = length(codons_in_family),
      Slope_Expression = 0,
      SE = NA,
      p_value = 1.0,
      p_CDS_length = NA,
      p_GC3s = NA,
      Converged = FALSE,
      N_genes = 0
    ))
  }
  
  # Fit model: binomial for 2-codon families, multinomial for 3+
  tryCatch({
    if (K == 1) {
      # Binomial: 2-codon family (preferred vs non-preferred)
      family_data_clean$Preferred_Count <- family_data_clean[[preferred_codon]]
      family_data_clean$NonPreferred_Count <- family_data_clean[[non_preferred_codons[1]]]
      
      model <- mgcv::gam(
        cbind(Preferred_Count, NonPreferred_Count) ~ High_exp_log2 + s(CDS_length_nt, k = 5) + s(GC3s, k = 5),
        family = binomial(),
        data = family_data_clean,
        na.action = na.exclude
      )
      
      # Extract coefficients
      coef_summary <- summary(model)$p.table
      slope_exp <- coef_summary["High_exp_log2", "Estimate"]
      se_exp <- coef_summary["High_exp_log2", "Std. Error"]
      p_exp <- coef_summary["High_exp_log2", "Pr(>|z|)"]
      
    } else {
      # Multinomial: 3+ codon family (preferred as baseline)
      # Build response matrix: cbind(preferred, codon1, codon2, ...)
      ordered_codons <- c(preferred_codon, non_preferred_codons)
      Y <- as.matrix(family_data_clean[, ordered_codons])
      
      suppressPackageStartupMessages(require(VGAM))
      
      model <- VGAM::vgam(
        Y ~ High_exp_log2 + s(CDS_length_nt, df = 4) + s(GC3s, df = 4),
        family = VGAM::multinomial(refLevel = 1),  # First codon (preferred) is baseline
        data = family_data_clean,
        na.action = na.exclude
      )
      
      # Extract High_exp_log2 coefficients for all non-baseline codons
      # VGAM uses coef() and vcov() directly, not coef(summary())
      coef_est <- coef(model)
      se_all <- sqrt(diag(vcov(model)))
      z_values <- coef_est / se_all
      p_values <- 2 * pnorm(-abs(z_values))
      
      # Find High_exp_log2 terms
      slope_indices <- grep("High_exp_log2", names(coef_est))
      slope_exp <- mean(coef_est[slope_indices])
      se_exp <- mean(se_all[slope_indices])
      p_exp <- min(p_values[slope_indices])  # Most significant
      
      # For VGAM, smoothers are in the model but not easily extracted
      # Set to NA for multinomial models
      p_CDS_length <- NA
      p_GC3s <- NA
    }
    
    # Get smoother statistics (mutation effects) - only for mgcv models
    if (K == 1) {
      smooth_summary <- summary(model)$s.table
      p_CDS_length <- smooth_summary["s(CDS_length_nt)", "p-value"]
      p_GC3s <- smooth_summary["s(GC3s)", "p-value"]
    }
    
    return(data.frame(
      Family = family_name,
      Preferred_Codon = preferred_codon,
      N_codons = length(codons_in_family),
      Slope_Expression = slope_exp,
      SE = se_exp,
      p_value = p_exp,
      p_CDS_length = p_CDS_length,
      p_GC3s = p_GC3s,
      Converged = TRUE,
      N_genes = nrow(family_data_clean)
    ))
    
  }, error = function(e) {
    warning(sprintf("Error fitting model for family %s: %s", family_name, e$message))
    return(data.frame(
      Family = family_name,
      Preferred_Codon = preferred_codon,
      N_codons = length(codons_in_family),
      Slope_Expression = 0,
      SE = NA,
      p_value = 1.0,
      p_CDS_length = NA,
      p_GC3s = NA,
      Converged = FALSE,
      N_genes = 0
    ))
  })
}
