#!/usr/bin/env Rscript
# model_selection_diversity.R
# 
# Systematic model exploration for testing:
# H1: Does Pi at 4-fold sites decrease with expression?
# H2: Does Tajima's D decrease with expression?
#
# Approach:
# 1. Define candidate models with varying complexity
# 2. Compare models using AIC/BIC
# 3. Detrend response from confounders (excluding expression)
# 4. Test clean expression-diversity relationship
#
# Author: Luis J. Madrigal-Roca
# Date: December 2025

#' Explore model space for diversity ~ expression relationship
#' 
#' @param data Data frame with integrated gene data
#' @param response Response variable ("Pi_mean_4fold" or "TajimaD_4fold")
#' @param expression_var Expression variable ("High_exp_log2" or "phi_estimate")
#' @param confounders Character vector of confounder names
#' @param include_interactions Include expression:confounder interactions?
#' @param include_quadratic Include quadratic terms?
#' @param scale_predictors Scale continuous predictors to mean=0, sd=1?
#' @param output_dir Directory to save results
#' 
#' @return List with model comparison, best model, and detrended analysis
explore_diversity_expression_models <- function(
    data,
    response = "Pi_mean_4fold",
    expression_var = "High_exp_log2",
    confounders = c("CDS_length_nt", "GC3s"),
    include_interactions = TRUE,
    include_quadratic = TRUE,
    scale_predictors = FALSE,
    output_dir = "./results/model_selection"
) {
  
  require(mgcv)
  require(dplyr)
  require(ggplot2)
  
  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  
  # Check that all required columns exist
  vars_needed <- c(response, expression_var, confounders)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing columns in data: ", paste(missing_vars, collapse = ", "),
         "\nAvailable columns: ", paste(head(names(data), 20), collapse = ", "), "...")
  }
  
  # Check for infinite/NA values in expression variable (log(0) protection)
  expr_values <- data[[expression_var]]
  n_inf <- sum(is.infinite(expr_values), na.rm = TRUE)
  n_na <- sum(is.na(expr_values))
  if (n_inf > 0) {
    warning(sprintf("Found %d infinite values in %s (possibly from log(0)). These will be excluded.",
                    n_inf, expression_var))
  }
  if (n_na > 0) {
    cat(sprintf("Note: %d NA values in %s will be excluded.\n", n_na, expression_var))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n", strrep("=", 70), "\n")
  cat("MODEL SELECTION FOR DIVERSITY-EXPRESSION RELATIONSHIP\n")
  cat(strrep("=", 70), "\n\n")
  cat("Response:", response, "\n")
  cat("Expression variable:", expression_var, "\n")
  cat("Confounders:", paste(confounders, collapse = ", "), "\n\n")
  
  # Prepare clean data (remove NAs and infinite values for consistent model comparison)
  vars_needed <- c(response, expression_var, confounders)
  clean_data <- data %>%
    dplyr::select(all_of(vars_needed)) %>%
    dplyr::filter(if_all(everything(), ~ !is.infinite(.))) %>%
    na.omit()
  
  n_obs <- nrow(clean_data)
  cat("Sample size (complete cases):", n_obs, "\n")
  
  # Optional: Scale predictors for coefficient comparability
  scaling_info <- NULL
  if (scale_predictors) {
    cat("Scaling predictors to mean=0, sd=1...\n")
    scaling_info <- list()
    
    # Scale expression variable
    scaling_info[[expression_var]] <- list(
      mean = mean(clean_data[[expression_var]], na.rm = TRUE),
      sd = sd(clean_data[[expression_var]], na.rm = TRUE)
    )
    clean_data[[expression_var]] <- scale(clean_data[[expression_var]])[,1]
    
    # Scale confounders
    for (conf in confounders) {
      scaling_info[[conf]] <- list(
        mean = mean(clean_data[[conf]], na.rm = TRUE),
        sd = sd(clean_data[[conf]], na.rm = TRUE)
      )
      clean_data[[conf]] <- scale(clean_data[[conf]])[,1]
    }
    
    cat("  Scaling applied. Original means/SDs stored in output.\n")
  }
  
  cat("\n")
  
  # =========================================================================
  # STEP 1: Define candidate models
  # =========================================================================
  cat("STEP 1: Building candidate models...\n\n")
  
  # Base formula components
  conf_smooth <- paste0("s(", confounders, ")", collapse = " + ")
  conf_linear <- paste(confounders, collapse = " + ")
  
  # Build model formulas
  model_formulas <- list()
  model_descriptions <- list()
  
  # Model 1: Null model (confounders only, smooth)
  model_formulas[["M1_null_smooth"]] <- as.formula(
    paste(response, "~", conf_smooth)
  )
  model_descriptions[["M1_null_smooth"]] <- "Null: confounders only (smooth)"
  
  # Model 2: Null model (confounders only, linear)
  model_formulas[["M2_null_linear"]] <- as.formula(
    paste(response, "~", conf_linear)
  )
  model_descriptions[["M2_null_linear"]] <- "Null: confounders only (linear)"
  
  # Model 3: Expression linear + confounders smooth
  model_formulas[["M3_expr_linear"]] <- as.formula(
    paste(response, "~", expression_var, "+", conf_smooth)
  )
  model_descriptions[["M3_expr_linear"]] <- "Expression linear + confounders smooth"
  
  # Model 4: Expression smooth + confounders smooth
  model_formulas[["M4_expr_smooth"]] <- as.formula(
    paste(response, "~ s(", expression_var, ") +", conf_smooth)
  )
  model_descriptions[["M4_expr_smooth"]] <- "Expression smooth + confounders smooth"
  
  # Model 5: All linear
  model_formulas[["M5_all_linear"]] <- as.formula(
    paste(response, "~", expression_var, "+", conf_linear)
  )
  model_descriptions[["M5_all_linear"]] <- "All predictors linear"
  
  if (include_quadratic) {
    # Model 6: Expression quadratic + confounders smooth
    model_formulas[["M6_expr_quad"]] <- as.formula(
      paste(response, "~ poly(", expression_var, ", 2) +", conf_smooth)
    )
    model_descriptions[["M6_expr_quad"]] <- "Expression quadratic + confounders smooth"
    
    # Model 7: All quadratic
    quad_terms <- paste0("poly(", c(expression_var, confounders), ", 2)", collapse = " + ")
    model_formulas[["M7_all_quad"]] <- as.formula(
      paste(response, "~", quad_terms)
    )
    model_descriptions[["M7_all_quad"]] <- "All predictors quadratic"
  }
  
  if (include_interactions && length(confounders) > 0) {
    # Model 8: Expression + confounders + first-order interactions
    interaction_terms <- paste0(expression_var, ":", confounders, collapse = " + ")
    model_formulas[["M8_interactions"]] <- as.formula(
      paste(response, "~", expression_var, "+", conf_linear, "+", interaction_terms)
    )
    model_descriptions[["M8_interactions"]] <- "Linear with expression:confounder interactions"
    
    # Model 9: Expression smooth + tensor interactions with confounders
    if (length(confounders) >= 1) {
      tensor_terms <- paste0("te(", expression_var, ", ", confounders, ")", collapse = " + ")
      model_formulas[["M9_tensor"]] <- as.formula(
        paste(response, "~", tensor_terms)
      )
      model_descriptions[["M9_tensor"]] <- "Tensor product smooth interactions"
    }
    
    # Model 10: Expression smooth + main effects smooth + interactions
    model_formulas[["M10_smooth_interact"]] <- as.formula(
      paste(response, "~ s(", expression_var, ") +", conf_smooth, "+", 
            paste0("ti(", expression_var, ", ", confounders, ")", collapse = " + "))
    )
    model_descriptions[["M10_smooth_interact"]] <- "Smooth main effects + smooth interactions"
  }
  
  cat("Built", length(model_formulas), "candidate models\n\n")
  
  # =========================================================================
  # STEP 2: Fit and compare models
  # =========================================================================
  cat("STEP 2: Fitting models and comparing...\n\n")
  
  model_results <- data.frame(
    Model = character(),
    Description = character(),
    AIC = numeric(),
    BIC = numeric(),
    GCV = numeric(),
    Deviance_explained = numeric(),
    Adj_R2 = numeric(),
    edf = numeric(),
    stringsAsFactors = FALSE
  )
  
  fitted_models <- list()
  
  for (model_name in names(model_formulas)) {
    
    formula <- model_formulas[[model_name]]
    description <- model_descriptions[[model_name]]
    
    tryCatch({
      # Fit GAM (works for both smooth and linear terms)
      fit <- gam(formula, data = clean_data, method = "REML")
      fitted_models[[model_name]] <- fit
      
      # Extract metrics
      summ <- summary(fit)
      
      model_results <- rbind(model_results, data.frame(
        Model = model_name,
        Description = description,
        AIC = AIC(fit),
        BIC = BIC(fit),
        GCV = fit$gcv.ubre,
        Deviance_explained = summ$dev.expl * 100,
        Adj_R2 = summ$r.sq,
        edf = sum(fit$edf),
        stringsAsFactors = FALSE
      ))
      
      cat(sprintf("  ✓ %s: AIC=%.1f, R²=%.3f\n", 
                  model_name, AIC(fit), summ$r.sq))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s: Failed - %s\n", model_name, e$message))
    })
  }
  
  # Sort by AIC
  model_results <- model_results %>%
    arrange(AIC) %>%
    mutate(
      Delta_AIC = AIC - min(AIC),
      AIC_weight = exp(-0.5 * Delta_AIC) / sum(exp(-0.5 * Delta_AIC))
    )
  
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("MODEL COMPARISON (sorted by AIC)\n")
  cat(strrep("-", 70), "\n\n")
  
  print(model_results %>% 
          dplyr::select(Model, Description, AIC, Delta_AIC, AIC_weight, Adj_R2) %>%
          mutate(across(where(is.numeric), ~round(., 3))))
  
  # Best model
  best_model_name <- model_results$Model[1]
  best_model <- fitted_models[[best_model_name]]
  
  cat("\n→ Best model:", best_model_name, "\n")
  cat("  ", model_results$Description[1], "\n\n")
  
  # =========================================================================
  # STEP 3: Examine expression effect in best model
  # =========================================================================
  cat("STEP 3: Expression effect in best model\n")
  cat(strrep("-", 50), "\n\n")
  
  print(summary(best_model))
  
  # =========================================================================
  # STEP 4: Detrend from confounders and test pure expression effect
  # =========================================================================
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("STEP 4: DETRENDING APPROACH\n")
  cat(strrep("=", 70), "\n\n")
  cat("Removing confounder effects to isolate expression-diversity relationship\n\n")
  
  # Fit confounder-only model (null model with smooth terms)
  null_formula <- as.formula(paste(response, "~", conf_smooth))
  null_model <- gam(null_formula, data = clean_data, method = "REML")
  
  # Get residuals = response with confounders removed
  clean_data$response_detrended <- residuals(null_model)
  
  # Also detrend expression from confounders (for clean visualization)
  expr_null_formula <- as.formula(paste(expression_var, "~", conf_smooth))
  expr_null_model <- gam(expr_null_formula, data = clean_data, method = "REML")
  clean_data$expression_detrended <- residuals(expr_null_model)
  
  cat("Detrended", response, "from:", paste(confounders, collapse = ", "), "\n")
  cat("Detrended", expression_var, "from:", paste(confounders, collapse = ", "), "\n\n")
  
  # Test detrended relationship
  detrend_formula <- as.formula("response_detrended ~ expression_detrended")
  detrend_model <- lm(detrend_formula, data = clean_data)
  detrend_summary <- summary(detrend_model)
  
  cat("DETRENDED RELATIONSHIP TEST:\n")
  cat(strrep("-", 50), "\n")
  print(detrend_summary)
  
  # Extract key statistics
  coef_estimate <- coef(detrend_model)[2]
  coef_se <- detrend_summary$coefficients[2, 2]
  t_value <- detrend_summary$coefficients[2, 3]
  p_value <- detrend_summary$coefficients[2, 4]
  r_squared <- detrend_summary$r.squared
  
  cat("\n")
  cat("SUMMARY OF DETRENDED ANALYSIS:\n")
  cat(strrep("-", 50), "\n")
  cat(sprintf("Coefficient (expression): %.6f ± %.6f\n", coef_estimate, coef_se))
  cat(sprintf("t-value: %.3f\n", t_value))
  cat(sprintf("p-value: %.2e\n", p_value))
  cat(sprintf("R²: %.4f\n", r_squared))
  cat("\n")
  
  if (p_value < 0.05) {
    if (coef_estimate < 0) {
      cat("✓ SIGNIFICANT NEGATIVE relationship\n")
      cat("  Higher expression → LOWER diversity (consistent with selection)\n")
    } else {
      cat("⚠ SIGNIFICANT POSITIVE relationship\n")
      cat("  Higher expression → HIGHER diversity (opposite of expectation)\n")
    }
  } else {
    cat("✗ NO SIGNIFICANT relationship after detrending\n")
    cat("  Expression does not predict diversity beyond confounders\n")
  }
  
  # =========================================================================
  # STEP 5: Create visualizations
  # =========================================================================
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("STEP 5: CREATING VISUALIZATIONS\n")
  cat(strrep("=", 70), "\n\n")
  
  # Plot 1: Model comparison
  p_comparison <- ggplot(model_results, aes(x = reorder(Model, -AIC), y = AIC)) +
    geom_bar(stat = "identity", aes(fill = AIC_weight), color = "black") +
    geom_hline(yintercept = min(model_results$AIC) + 2, 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_viridis_c(name = "AIC Weight", option = "plasma") +
    coord_flip() +
    labs(
      title = paste("Model Comparison for", response, "~", expression_var),
      subtitle = "Red line = Delta AIC = 2 threshold",
      x = "Model",
      y = "AIC (lower = better)"
    ) +
    theme_custom() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, paste0("model_comparison_", response, ".pdf")),
         p_comparison, width = 10, height = 6)
  
  # Plot 2: Detrended relationship
  p_detrended <- ggplot(clean_data, aes(x = expression_detrended, y = response_detrended)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("Detrended Relationship:", response, "vs", expression_var),
      subtitle = sprintf("Confounders removed: %s\nCoef = %.2e, p = %.2e, R2 = %.4f",
                        paste(confounders, collapse = ", "),
                        coef_estimate, p_value, r_squared),
      x = paste(expression_var, "(residuals, confounders removed)"),
      y = paste(response, "(residuals, confounders removed)")
    ) +
    theme_custom() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, paste0("detrended_relationship_", response, ".pdf")),
         p_detrended, width = 8, height = 6)
  
  # Plot 3: Raw relationship with GC3s coloring (confound illustration)
  if ("GC3s" %in% names(clean_data)) {
    p_confound <- ggplot(clean_data, aes(x = .data[[expression_var]], y = .data[[response]])) +
      geom_point(aes(color = GC3s), alpha = 0.4, size = 1) +
      geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
      scale_color_viridis_c(name = "GC3s", option = "viridis") +
      labs(
        title = paste("Raw Relationship:", response, "vs", expression_var),
        subtitle = "Color shows GC3s confounding structure",
        x = expression_var,
        y = response
      ) +
      theme_custom() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(output_dir, paste0("raw_with_confound_", response, ".pdf")),
           p_confound, width = 8, height = 6)
  }
  
  # Plot 4: Partial effect from best model
  if (grepl("s\\(", as.character(model_formulas[[best_model_name]])[3])) {
    # If best model has smooth terms, plot partial effects
    pdf(file.path(output_dir, paste0("partial_effects_", response, ".pdf")),
        width = 10, height = 8)
    plot(best_model, pages = 1, residuals = FALSE, rug = TRUE,
         main = paste("Partial Effects -", best_model_name))
    dev.off()
  }
  
  cat("✓ Plots saved to:", output_dir, "\n\n")
  
  # =========================================================================
  # Save results
  # =========================================================================
  
  # Save model comparison table
  write.csv(model_results, 
            file.path(output_dir, paste0("model_comparison_", response, ".csv")),
            row.names = FALSE)
  
  # Save detrended data
  write.csv(clean_data,
            file.path(output_dir, paste0("detrended_data_", response, ".csv")),
            row.names = FALSE)
  
  # Return results
  results <- list(
    model_comparison = model_results,
    fitted_models = fitted_models,
    best_model_name = best_model_name,
    best_model = best_model,
    null_model = null_model,
    detrend_model = detrend_model,
    detrended_data = clean_data,
    detrend_summary = list(
      coefficient = coef_estimate,
      se = coef_se,
      t_value = t_value,
      p_value = p_value,
      r_squared = r_squared,
      direction = ifelse(coef_estimate < 0, "negative", "positive"),
      significant = p_value < 0.05
    ),
    settings = list(
      response = response,
      expression_var = expression_var,
      confounders = confounders,
      n_observations = n_obs
    )
  )
  
  return(results)
}


#' Run sensitivity analysis comparing empirical expression vs phi estimates
#' 
#' @param data Data frame with integrated gene data
#' @param response Response variable
#' @param confounders Confounder variables
#' @param output_dir Directory to save results
#' 
#' @return List with comparison results
run_expression_sensitivity <- function(
    data,
    response = "Pi_mean_4fold",
    phi_column = "Log_Phi",
    empirical_column = NULL,  # Will auto-detect if NULL
    confounders = c("CDS_length_nt", "GC3s"),
    output_dir = "./results/model_selection/sensitivity"
) {
  
  require(dplyr)
  require(ggplot2)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n", strrep("=", 70), "\n")
  cat("SENSITIVITY ANALYSIS: Empirical Expression vs Phi Estimates\n")

  cat(strrep("=", 70), "\n\n")
  
  # Auto-detect empirical expression column if not specified
  if (is.null(empirical_column)) {
    possible_cols <- c("High_exp_log10", "High_exp_log2", "High_exp")
    for (col in possible_cols) {
      if (col %in% names(data)) {
        empirical_column <- col
        cat("Auto-detected empirical expression column:", col, "\n")
        break
      }
    }
    if (is.null(empirical_column)) {
      stop("Could not auto-detect empirical expression column. ",
           "Please specify empirical_column parameter.\n",
           "Available columns: ", paste(names(data)[1:20], collapse = ", "), "...")
    }
  }
  
  # Validate empirical column exists
  if (!empirical_column %in% names(data)) {
    stop("Empirical column '", empirical_column, "' not found in data.\n",
         "Available columns: ", paste(names(data)[1:20], collapse = ", "), "...")
  }
  
  # Check if phi column exists
  if (!phi_column %in% names(data)) {
    cat("Warning: phi column '", phi_column, "' not found in data.\n")
    cat("Available columns:", paste(names(data)[1:20], collapse = ", "), "...\n")
    cat("Proceeding with empirical expression only.\n\n")
    
    # Run only empirical analysis
    results_empirical <- explore_diversity_expression_models(
      data = data,
      response = response,
      expression_var = empirical_column,
      confounders = confounders,
      output_dir = file.path(output_dir, "empirical_expression")
    )
    
    return(list(
      empirical = results_empirical,
      phi = NULL,
      comparison = NULL,
      robust = results_empirical$detrend_summary$significant && 
               results_empirical$detrend_summary$direction == "Negative"
    ))
  }
  
  # Run analysis with empirical expression
  cat("\n--- Analysis 1: Using Empirical Expression ---\n")
  results_empirical <- explore_diversity_expression_models(
    data = data,
    response = response,
    expression_var = empirical_column,
    confounders = confounders,
    output_dir = file.path(output_dir, "empirical_expression")
  )
  
  # Run analysis with phi estimates
  cat("\n--- Analysis 2: Using Phi Estimates (AnaCoDa) ---\n")
  results_phi <- explore_diversity_expression_models(
    data = data,
    response = response,
    expression_var = phi_column,
    confounders = confounders,
    output_dir = file.path(output_dir, "phi_estimates")
  )
  
  # Compare results
  cat("\n", strrep("=", 70), "\n")
  cat("SENSITIVITY COMPARISON\n")
  cat(strrep("=", 70), "\n\n")
  
  comparison <- data.frame(
    Measure = c("Coefficient", "t-value", "p-value", "R²", "Direction", "Significant"),
    Empirical_Expression = c(
      sprintf("%.2e", results_empirical$detrend_summary$coefficient),
      sprintf("%.3f", results_empirical$detrend_summary$t_value),
      sprintf("%.2e", results_empirical$detrend_summary$p_value),
      sprintf("%.4f", results_empirical$detrend_summary$r_squared),
      results_empirical$detrend_summary$direction,
      as.character(results_empirical$detrend_summary$significant)
    ),
    Phi_Estimates = c(
      sprintf("%.2e", results_phi$detrend_summary$coefficient),
      sprintf("%.3f", results_phi$detrend_summary$t_value),
      sprintf("%.2e", results_phi$detrend_summary$p_value),
      sprintf("%.4f", results_phi$detrend_summary$r_squared),
      results_phi$detrend_summary$direction,
      as.character(results_phi$detrend_summary$significant)
    )
  )
  
  print(comparison)
  
  # Interpretation
  cat("\n")
  cat("INTERPRETATION:\n")
  cat(strrep("-", 50), "\n")
  
  both_sig <- results_empirical$detrend_summary$significant && 
              results_phi$detrend_summary$significant
  same_dir <- results_empirical$detrend_summary$direction == 
              results_phi$detrend_summary$direction
  
  if (both_sig && same_dir) {
    cat("✓ ROBUST: Both measures show significant effect in same direction\n")
    cat("  Conclusion is robust to expression estimation method\n")
  } else if (both_sig && !same_dir) {
    cat("⚠ CONFLICTING: Both significant but opposite directions\n")
    cat("  Results depend on expression estimation method - investigate further\n")
  } else if (!both_sig) {
    cat("⚠ SENSITIVE: Significance depends on expression measure\n")
    if (results_empirical$detrend_summary$significant) {
      cat("  Empirical expression shows effect, phi estimates do not\n")
    } else if (results_phi$detrend_summary$significant) {
      cat("  Phi estimates show effect, empirical expression does not\n")
    } else {
      cat("  Neither measure shows significant effect\n")
    }
  }
  
  # Save comparison
  write.csv(comparison,
            file.path(output_dir, paste0("sensitivity_comparison_", response, ".csv")),
            row.names = FALSE)
  
  # Create comparison plot
  plot_data <- data.frame(
    Method = rep(c("Empirical Expression", "Phi Estimates"), each = 2),
    Measure = rep(c("Coefficient", "R²"), 2),
    Value = c(
      results_empirical$detrend_summary$coefficient,
      results_empirical$detrend_summary$r_squared,
      results_phi$detrend_summary$coefficient,
      results_phi$detrend_summary$r_squared
    )
  )
  
  return(list(
    empirical = results_empirical,
    phi = results_phi,
    comparison = comparison,
    robust = both_sig && same_dir
  ))
}


#' Run complete diversity-expression analysis pipeline
#' 
#' @param data Integrated data frame
#' @param output_base Base directory for outputs
#' 
#' @return List with all analysis results
run_complete_diversity_analysis <- function(
    data,
    phi_column = "Log_Phi",
    empirical_column = NULL,  # Will auto-detect if NULL
    output_base = "./results/diversity_expression_analysis"
) {
  
  cat("\n")
  cat(strrep("#", 80), "\n")
  cat("COMPLETE DIVERSITY-EXPRESSION ANALYSIS PIPELINE\n")
  cat(strrep("#", 80), "\n\n")
  
  results <- list()
  
  # Analysis 1: Pi at 4-fold sites
  cat("\n", strrep("*", 70), "\n")
  cat("ANALYSIS 1: Nucleotide Diversity (Pi) at 4-fold sites\n")
  cat(strrep("*", 70), "\n")
  
  if ("Pi_mean_4fold" %in% names(data)) {
    results$pi_4fold <- run_expression_sensitivity(
      data = data,
      response = "Pi_mean_4fold",
      phi_column = phi_column,
      empirical_column = empirical_column,
      confounders = c("CDS_length_nt", "GC3s"),
      output_dir = file.path(output_base, "pi_4fold")
    )
  } else {
    cat("Warning: Pi_mean_4fold not found in data\n")
  }

  
  # Analysis 2: Tajima's D at 4-fold sites
  cat("\n", strrep("*", 70), "\n")
  cat("ANALYSIS 2: Tajima's D at 4-fold sites\n")
  cat(strrep("*", 70), "\n")
  
  if ("TajimaD_4fold" %in% names(data)) {
    results$tajima_4fold <- run_expression_sensitivity(
      data = data,
      response = "TajimaD_4fold",
      phi_column = phi_column,
      empirical_column = empirical_column,
      confounders = c("CDS_length_nt", "GC3s"),
      output_dir = file.path(output_base, "tajima_4fold")
    )
  } else {
    cat("Warning: TajimaD_4fold not found in data\n")
  }
  
  # Analysis 3: Pi at all sites (for comparison)
  cat("\n", strrep("*", 70), "\n")
  cat("ANALYSIS 3: Nucleotide Diversity (Pi) at all sites\n")
  cat(strrep("*", 70), "\n")
  
  if ("Pi_mean_all" %in% names(data)) {
    results$pi_all <- run_expression_sensitivity(
      data = data,
      response = "Pi_mean_all",
      phi_column = phi_column,
      empirical_column = empirical_column,
      confounders = c("CDS_length_nt", "GC3s"),
      output_dir = file.path(output_base, "pi_all")
    )
  } else {
    cat("Warning: Pi_mean_all not found in data\n")
  }
  
  # Summary
  cat("\n")
  cat(strrep("#", 80), "\n")
  cat("ANALYSIS COMPLETE - SUMMARY\n")
  cat(strrep("#", 80), "\n\n")
  
  cat("Results saved to:", output_base, "\n\n")
  
  cat("Key findings:\n")
  if (!is.null(results$pi_4fold)) {
    cat("  Pi (4-fold): Robust =", results$pi_4fold$robust, "\n")
  }
  if (!is.null(results$tajima_4fold)) {
    cat("  Tajima's D (4-fold): Robust =", results$tajima_4fold$robust, "\n")
  }
  if (!is.null(results$pi_all)) {
    cat("  Pi (all sites): Robust =", results$pi_all$robust, "\n")
  }
  
  return(results)
}

# =============================================================================
# WRAPPER FUNCTIONS (aliases for main.R compatibility)
# =============================================================================

#' Run diversity analysis pipeline
#' 
#' Wrapper function that provides a simpler interface matching main.R expectations
#' 
#' @param data Data frame with integrated gene data
#' @param response Response variable (e.g., "Pi_mean_4fold")
#' @param expression_var Expression variable (e.g., "High_exp_log10")
#' @param predictors Character vector of predictor/confounder names
#' @param include_quadratic Include quadratic terms?
#' @param include_interactions Include interaction terms?
#' @param scale_predictors Scale predictors to mean=0, sd=1 for comparability?
#' @param output_dir Directory to save results
#' 
#' @return List with model_comparison, best_model, and analysis results
run_diversity_analysis_pipeline <- function(
    data,
    response,
    expression_var,
    predictors = c("GC3s", "CDS_length_nt"),
    include_quadratic = TRUE,
    include_interactions = TRUE,
    scale_predictors = FALSE,
    output_dir = "./results/diversity_modeling"
) {
  
  # Call the main exploration function
  results <- explore_diversity_expression_models(
    data = data,
    response = response,
    expression_var = expression_var,
    confounders = predictors,
    include_interactions = include_interactions,
    include_quadratic = include_quadratic,
    scale_predictors = scale_predictors,
    output_dir = output_dir
  )
  
  return(results)
}


#' Explore model robustness via detrending
#' 
#' Residualize response from confounders, then test expression effect
#' 
#' @param data Data frame with integrated gene data
#' @param response Response variable
#' @param expression_var Expression variable
#' @param confounders Character vector of confounder names
#' 
#' @return List with detrending model, expression model, and robustness assessment
explore_model_robustness <- function(
    data,
    response,
    expression_var,
    confounders = c("GC3s", "CDS_length_nt")
) {
  
  require(mgcv)
  require(dplyr)
  
  cat("\n--- Robustness Check: Detrending Approach ---\n")
  cat("Response:", response, "\n")
  cat("Expression:", expression_var, "\n")
  cat("Confounders:", paste(confounders, collapse = ", "), "\n\n")
  
  # Prepare clean data
  vars_needed <- c(response, expression_var, confounders)
  clean_data <- data %>%
    dplyr::select(all_of(vars_needed)) %>%
    na.omit()
  
  n_obs <- nrow(clean_data)
  cat("Sample size:", n_obs, "\n")
  
  # Step 1: Fit GAM with confounders only (excluding expression)
  conf_formula <- as.formula(
    paste(response, "~", paste0("s(", confounders, ")", collapse = " + "))
  )
  
  confounder_model <- gam(conf_formula, data = clean_data)
  
  # Step 2: Extract residuals (response detrended from confounders)
  clean_data$residuals_detrended <- residuals(confounder_model)
  
  # Step 3: Regress residuals on expression
  expr_formula <- as.formula(paste("residuals_detrended ~", expression_var))
  expression_model <- lm(expr_formula, data = clean_data)
  
  expr_summary <- summary(expression_model)
  
  # Extract results
  expr_coef <- coef(expression_model)[expression_var]
  expr_se <- expr_summary$coefficients[expression_var, "Std. Error"]
  expr_t <- expr_summary$coefficients[expression_var, "t value"]
  expr_p <- expr_summary$coefficients[expression_var, "Pr(>|t|)"]
  expr_r2 <- expr_summary$r.squared
  
  cat("\n--- Detrended Analysis Results ---\n")
  cat("Expression coefficient:", round(expr_coef, 6), "\n")
  cat("Standard error:", round(expr_se, 6), "\n")
  cat("t-value:", round(expr_t, 3), "\n")
  cat("p-value:", format.pval(expr_p, digits = 4), "\n")
  cat("R² (expression alone):", round(expr_r2, 4), "\n")
  
  # Interpretation
  if (expr_p < 0.05 && expr_coef < 0) {
    interpretation <- "SUPPORTED: Diversity decreases with expression (after controlling for confounders)"
    robust <- TRUE
  } else if (expr_p < 0.05 && expr_coef > 0) {
    interpretation <- "OPPOSITE: Diversity increases with expression"
    robust <- FALSE
  } else {
    interpretation <- "NOT SIGNIFICANT: No clear relationship after controlling for confounders"
    robust <- FALSE
  }
  
  cat("\nInterpretation:", interpretation, "\n")
  
  return(list(
    confounder_model = confounder_model,
    expression_model = expression_model,
    detrended_data = clean_data,
    expression_effect = list(
      coefficient = expr_coef,
      se = expr_se,
      t_value = expr_t,
      p_value = expr_p,
      r_squared = expr_r2
    ),
    interpretation = interpretation,
    robust = robust
  ))
}
