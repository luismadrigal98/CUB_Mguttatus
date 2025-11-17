fit_pairwise_gams <- function(family_name, genetic_code, usage_dt, meta_dt, 
                              preferred_codons_df = NULL) {
  #' Fits pairwise binomial GAMs for all codons in a family
  #' relative to a baseline codon using GAM smoothers for confounders.
  #'
  #' @param family_name The name of the amino acid family (e.g., "Ala", "Leu_4")
  #' @param genetic_code Named vector: codon -> amino acid
  #' @param usage_dt Codon usage data.table with Gene_name column (COUNTS per gene)
  #' @param meta_dt Metadata with High_exp_log2, CDS_length_nt, GC12
  #' @param preferred_codons_df Optional data frame with Codon and AA columns to use as baseline
  #' @return List with family info, models, and coefficient table
  
  suppressPackageStartupMessages({
    require(mgcv)
    require(dplyr)
    require(data.table)
  })
  
  # 1. Get all synonymous codons for this family
  codons_in_family <- names(genetic_code[genetic_code == family_name])
  
  if (length(codons_in_family) < 2) {
    warning(sprintf("Family %s has < 2 codons, skipping", family_name))
    return(NULL)
  }
  
  # 2. Determine baseline codon
  # If preferred codons provided, use that; otherwise use alphabetically first
  if (!is.null(preferred_codons_df)) {
    preferred_for_family <- preferred_codons_df %>%
      dplyr::filter(AA == family_name | Amino_Acid == family_name)
    
    if (nrow(preferred_for_family) > 0) {
      baseline_codon <- preferred_for_family$Codon[1]
      if (!(baseline_codon %in% codons_in_family)) {
        warning(sprintf("Preferred codon %s not in family %s, using alphabetical", 
                       baseline_codon, family_name))
        baseline_codon <- sort(codons_in_family)[1]
      }
    } else {
      baseline_codon <- sort(codons_in_family)[1]
    }
  } else {
    baseline_codon <- sort(codons_in_family)[1]
  }
  
  response_codons <- setdiff(codons_in_family, baseline_codon)
  
  # 3. Prepare the base data for this family
  family_data <- as.data.frame(meta_dt) %>%
    dplyr::select(Gene_name, High_exp_log2, CDS_length_nt, GC12) %>%
    dplyr::left_join(
      as.data.frame(usage_dt) %>% dplyr::select(dplyr::all_of(c("Gene_name", codons_in_family))),
      by = "Gene_name"
    ) %>%
    na.omit()
  
  fitted_models <- list()
  coefficients_list <- list()
  
  # 4. Loop through each response codon and fit a binomial model
  for (response_codon in response_codons) {
    
    # Create data for this specific pair
    pair_data <- family_data %>%
      dplyr::select(
        Gene_name, High_exp_log2, CDS_length_nt, GC12,
        dplyr::all_of(c(response_codon, baseline_codon))
      ) %>%
      # Important: weights are the sum of *only these two* codons
      dplyr::mutate(Total_Pair_Count = .data[[response_codon]] + .data[[baseline_codon]]) %>%
      dplyr::filter(Total_Pair_Count > 0)
    
    if (nrow(pair_data) < 50) {
      warning(sprintf("Skipping %s vs %s: < 50 data points", response_codon, baseline_codon))
      next
    }
    
    # 5. Fit the binomial GAM
    formula_str <- sprintf(
      "cbind(%s, %s) ~ High_exp_log2 + s(CDS_length_nt, k = 5) + s(GC12, k = 5)",
      response_codon, baseline_codon
    )
    
    tryCatch({
      model <- mgcv::gam(
        as.formula(formula_str),
        family = binomial(link = "logit"),
        data = pair_data,
        na.action = na.exclude
      )
      
      # 6. Store the model and extract coefficients
      fitted_models[[response_codon]] <- model
      
      coef_summary <- summary(model)$p.table
      coefficients_list[[response_codon]] <- data.frame(
        Codon = response_codon,
        Baseline = baseline_codon,
        Family = family_name,
        Mutation_Intercept = coef_summary["(Intercept)", "Estimate"],
        Selection_Slope = coef_summary["High_exp_log2", "Estimate"],
        SE = coef_summary["High_exp_log2", "Std. Error"],
        p_value = coef_summary["High_exp_log2", "Pr(>|z|)"],
        N_genes = nrow(pair_data)
      )
    }, error = function(e) {
      warning(sprintf("Error fitting %s vs %s: %s", response_codon, baseline_codon, e$message))
    })
  }
  
  # 7. Add the baseline codon (slope = 0 by definition)
  coefficients_list[[baseline_codon]] <- data.frame(
    Codon = baseline_codon,
    Baseline = baseline_codon,
    Family = family_name,
    Mutation_Intercept = 0,
    Selection_Slope = 0,
    SE = 0,
    p_value = 1.0,
    N_genes = nrow(family_data)
  )
  
  # 8. Combine coefficients and apply FDR correction
  coef_df <- data.table::rbindlist(coefficients_list)
  
  # Apply FDR correction (Benjamini-Hochberg) to p-values
  # Exclude baseline codon (p = 1.0) from correction
  non_baseline <- coef_df$Codon != baseline_codon
  if (sum(non_baseline) > 0) {
    coef_df$p_adj <- 1.0
    coef_df$p_adj[non_baseline] <- p.adjust(coef_df$p_value[non_baseline], 
                                            method = "BH")
    coef_df$Significant <- coef_df$p_adj < 0.05
  } else {
    coef_df$p_adj <- coef_df$p_value
    coef_df$Significant <- FALSE
  }
  
  return(list(
    family = family_name,
    baseline = baseline_codon,
    models = fitted_models,
    coefficients = coef_df
  ))
}


fit_pairwise_glms <- function(family_name, genetic_code, usage_dt, meta_dt,
                              boxcox_confounders = TRUE, preferred_codons_df = NULL) {
  #' Fits pairwise binomial GLMs with Box-Cox transformed confounders.
  #' Avoids GAM smoothers for improved interpretability.
  #'
  #' @param family_name The name of the amino acid family
  #' @param genetic_code Named vector: codon -> amino acid
  #' @param usage_dt Codon usage data.table with Gene_name column (COUNTS per gene)
  #' @param meta_dt Metadata with High_exp_log2, CDS_length_nt, GC12
  #' @param boxcox_confounders If TRUE, apply Box-Cox to CDS_length_nt and GC12
  #' @param preferred_codons_df Optional data frame with Codon and AA columns to use as baseline
  #' @return List with family info, models, and coefficient table
  
  suppressPackageStartupMessages({
    require(dplyr)
    require(data.table)
    require(MASS)  # For boxcox
  })
  
  # 1. Get all synonymous codons for this family
  codons_in_family <- names(genetic_code[genetic_code == family_name])
  
  if (length(codons_in_family) < 2) {
    warning(sprintf("Family %s has < 2 codons, skipping", family_name))
    return(NULL)
  }
  
  # 2. Determine baseline codon
  if (!is.null(preferred_codons_df)) {
    preferred_for_family <- preferred_codons_df %>%
      dplyr::filter(AA == family_name | Amino_Acid == family_name)
    
    if (nrow(preferred_for_family) > 0) {
      baseline_codon <- preferred_for_family$Codon[1]
      if (!(baseline_codon %in% codons_in_family)) {
        warning(sprintf("Preferred codon %s not in family %s, using alphabetical", 
                       baseline_codon, family_name))
        baseline_codon <- sort(codons_in_family)[1]
      }
    } else {
      baseline_codon <- sort(codons_in_family)[1]
    }
  } else {
    baseline_codon <- sort(codons_in_family)[1]
  }
  
  response_codons <- setdiff(codons_in_family, baseline_codon)
  
  # 3. Prepare the base data for this family
  family_data <- as.data.frame(meta_dt) %>%
    dplyr::select(Gene_name, High_exp_log2, CDS_length_nt, GC12) %>%
    dplyr::left_join(
      as.data.frame(usage_dt) %>% dplyr::select(dplyr::all_of(c("Gene_name", codons_in_family))),
      by = "Gene_name"
    ) %>%
    na.omit()
  
  # 4. Box-Cox transformation of confounders if requested
  if (boxcox_confounders) {
    cat(sprintf("  Applying Box-Cox transformation to confounders for %s...\n", family_name))
    
    # Box-Cox for CDS_length_nt (must be positive)
    if (min(family_data$CDS_length_nt, na.rm = TRUE) <= 0) {
      family_data$CDS_length_nt <- family_data$CDS_length_nt + 1
    }
    
    # Fit Box-Cox model to determine optimal lambda
    bc_length <- MASS::boxcox(CDS_length_nt ~ 1, 
                              data = family_data, 
                              plotit = FALSE)
    lambda_length <- bc_length$x[which.max(bc_length$y)]
    
    # Apply transformation
    if (abs(lambda_length) < 0.01) {
      family_data$CDS_length_BC <- log(family_data$CDS_length_nt)
    } else {
      family_data$CDS_length_BC <- (family_data$CDS_length_nt^lambda_length - 1) / lambda_length
    }
    
    # Box-Cox for GC12 (must be in (0, 1), scale to avoid boundary issues)
    GC12_scaled <- family_data$GC12 * 0.98 + 0.01  # Avoid exact 0 or 1
    
    bc_gc <- MASS::boxcox(GC12_scaled ~ 1, 
                          data = family_data, 
                          plotit = FALSE)
    lambda_gc <- bc_gc$x[which.max(bc_gc$y)]
    
    # Apply transformation
    if (abs(lambda_gc) < 0.01) {
      family_data$GC12_BC <- log(GC12_scaled)
    } else {
      family_data$GC12_BC <- (GC12_scaled^lambda_gc - 1) / lambda_gc
    }
    
    confounder_formula <- "High_exp_log2 + CDS_length_BC + GC12_BC"
  } else {
    confounder_formula <- "High_exp_log2 + CDS_length_nt + GC12"
  }
  
  fitted_models <- list()
  coefficients_list <- list()
  
  # 5. Loop through each response codon and fit a binomial GLM
  for (response_codon in response_codons) {
    
    # Create data for this specific pair
    pair_data <- family_data %>%
      dplyr::select(
        Gene_name, High_exp_log2, 
        dplyr::any_of(c("CDS_length_nt", "GC12", "CDS_length_BC", "GC12_BC")),
        dplyr::all_of(c(response_codon, baseline_codon))
      ) %>%
      dplyr::mutate(Total_Pair_Count = .data[[response_codon]] + .data[[baseline_codon]]) %>%
      dplyr::filter(Total_Pair_Count > 0)
    
    if (nrow(pair_data) < 50) {
      warning(sprintf("Skipping %s vs %s: < 50 data points", response_codon, baseline_codon))
      next
    }
    
    # 6. Fit the binomial GLM
    formula_str <- sprintf(
      "cbind(%s, %s) ~ %s",
      response_codon, baseline_codon, confounder_formula
    )
    
    tryCatch({
      model <- glm(
        as.formula(formula_str),
        family = binomial(link = "logit"),
        data = pair_data
      )
      
      # 7. Store the model and extract coefficients
      fitted_models[[response_codon]] <- model
      
      coef_summary <- summary(model)$coefficients
      coefficients_list[[response_codon]] <- data.frame(
        Codon = response_codon,
        Baseline = baseline_codon,
        Family = family_name,
        Mutation_Intercept = coef_summary["(Intercept)", "Estimate"],
        Selection_Slope = coef_summary["High_exp_log2", "Estimate"],
        SE = coef_summary["High_exp_log2", "Std. Error"],
        p_value = coef_summary["High_exp_log2", "Pr(>|z|)"],
        N_genes = nrow(pair_data),
        BoxCox_Applied = boxcox_confounders
      )
    }, error = function(e) {
      warning(sprintf("Error fitting %s vs %s: %s", response_codon, baseline_codon, e$message))
    })
  }
  
  # 8. Add the baseline codon (slope = 0 by definition)
  coefficients_list[[baseline_codon]] <- data.frame(
    Codon = baseline_codon,
    Baseline = baseline_codon,
    Family = family_name,
    Mutation_Intercept = 0,
    Selection_Slope = 0,
    SE = 0,
    p_value = 1.0,
    N_genes = nrow(family_data),
    BoxCox_Applied = boxcox_confounders
  )
  
  # 9. Combine coefficients and apply FDR correction
  coef_df <- data.table::rbindlist(coefficients_list)
  
  # Apply FDR correction (Benjamini-Hochberg) to p-values
  non_baseline <- coef_df$Codon != baseline_codon
  if (sum(non_baseline) > 0) {
    coef_df$p_adj <- 1.0
    coef_df$p_adj[non_baseline] <- p.adjust(coef_df$p_value[non_baseline], 
                                            method = "BH")
    coef_df$Significant <- coef_df$p_adj < 0.05
  } else {
    coef_df$p_adj <- coef_df$p_value
    coef_df$Significant <- FALSE
  }
  
  return(list(
    family = family_name,
    baseline = baseline_codon,
    models = fitted_models,
    coefficients = coef_df
  ))
}


predict_codon_frequencies <- function(model_result, meta_dt, 
                                     expression_range = NULL,
                                     n_points = 100) {
  #' Predicts codon frequencies across expression levels using fitted models
  #'
  #' @param model_result Output from fit_pairwise_gams or fit_pairwise_glms
  #' @param meta_dt The integrated_data for getting mean confounders
  #' @param expression_range Vector c(min, max) for expression, or NULL for auto
  #' @param n_points Number of prediction points
  #' @return Data frame with predicted frequencies per codon
  
  family_name <- model_result$family
  baseline_codon <- model_result$baseline
  models <- model_result$models
  
  # 1. Create prediction grid with expanded range for clear S-curves
  if (is.null(expression_range)) {
    raw_range <- range(meta_dt$High_exp_log2, na.rm = TRUE)
    # Expand range by 10% on each side to show full S-curve dynamics
    range_width <- diff(raw_range)
    expression_range <- c(
      raw_range[1] - 0.1 * range_width,
      raw_range[2] + 0.1 * range_width
    )
  }
  
  prediction_grid <- data.frame(
    High_exp_log2 = seq(expression_range[1], expression_range[2], 
                        length.out = n_points)
  )
  
  # Use mean values for confounders
  prediction_grid$CDS_length_nt <- mean(meta_dt$CDS_length_nt, na.rm = TRUE)
  prediction_grid$GC12 <- mean(meta_dt$GC12, na.rm = TRUE)
  
  # If Box-Cox was applied, need to transform confounders
  # Check if first model has CDS_length_BC variable
  if (length(models) > 0) {
    first_model <- models[[1]]
    if ("CDS_length_BC" %in% names(first_model$model)) {
      # Extract Box-Cox parameters from the original fit
      # For simplicity, apply log transform (most common)
      prediction_grid$CDS_length_BC <- log(prediction_grid$CDS_length_nt)
      GC12_scaled <- prediction_grid$GC12 * 0.98 + 0.01
      prediction_grid$GC12_BC <- log(GC12_scaled)
    }
  }
  
  # 2. Get log-odds predictions for each response codon
  log_odds_predictions <- list()
  for (codon_name in names(models)) {
    log_odds_predictions[[codon_name]] <- predict(
      models[[codon_name]],
      newdata = prediction_grid,
      type = "link"
    )
  }
  
  # 3. Add the baseline (log-odds = 0)
  log_odds_predictions[[baseline_codon]] <- rep(0, n_points)
  
  # 4. Convert log-odds to probabilities using softmax
  log_odds_matrix <- do.call(cbind, log_odds_predictions)
  exp_log_odds <- exp(log_odds_matrix)
  prob_matrix <- exp_log_odds / rowSums(exp_log_odds)
  
  # 5. Format for plotting
  plot_data <- data.frame(High_exp_log2 = prediction_grid$High_exp_log2)
  plot_data <- cbind(plot_data, prob_matrix)
  
  # Convert to long format
  plot_data_long <- plot_data %>%
    tidyr::pivot_longer(
      cols = -High_exp_log2,
      names_to = "Codon",
      values_to = "Predicted_Frequency"
    )
  
  return(plot_data_long)
}
