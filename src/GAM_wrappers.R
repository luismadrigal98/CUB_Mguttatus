fit_codon_gam_suite <- function(data, 
                                model_list = alist(
                                  Null        = ~ 1,
                                  Additive    = ~ s(Max_Log10_Exp) + s(Exp_breadth) + s(CDS_length_nt),
                                  Interaction = ~ te(Max_Log10_Exp, Exp_breadth) + s(CDS_length_nt),
                                  Complex     = ~ te(Max_Log10_Exp, Exp_breadth, CDS_length_nt)
                                ),
                                response_var = "CDC", 
                                family = betar(link = "logit")) {
  require(mgcv)
  
  # Ensure directory exists
  if(!dir.exists("./results")) dir.create("./results")
  
  message(sprintf("\n--- Fitting GAM suite for response: %s ---", response_var))
  
  # Map over the alist and fit models
  fitted_models <- lapply(names(model_list), function(m_name) {
    # Combine response and predictor formula safely
    f_str <- paste(response_var, paste(deparse(model_list[[m_name]]), collapse = ""))
    f <- as.formula(f_str)
    
    message(sprintf("  Fitting %s...", m_name))
    
    # Using tryCatch to handle convergence issues in complex models
    tryCatch({
      gam(f, data = data, family = family, method = "REML", select = TRUE)
    }, error = function(e) {
      message(sprintf("  !! Model %s failed: %s", m_name, e$message))
      return(NULL)
    })
  })
  
  names(fitted_models) <- names(model_list)
  return(fitted_models)
}

get_model_selection_stats <- function(model_list) {
  stats <- lapply(names(model_list), function(n) {
    m <- model_list[[n]]
    data.frame(
      Model = n,
      AIC = AIC(m),
      Deviance_Expl = summary(m)$dev.expl * 100, # Percentage
      R_sq = summary(m)$r.sq,
      Converged = m$converged
    )
  })
  return(do.call(rbind, stats) |> dplyr::arrange(AIC))
}

run_slope_analysis <- function(model, data, breadths = c(1, 15, 29)) {
  require(marginaleffects)
  
  # 1. Calculate slopes at specific breadths
  # Holding length at the mean automatically via datagrid
  slopes <- avg_slopes(
    model,
    variables = "Max_Log10_Exp",
    by = "Exp_breadth",
    newdata = datagrid(Exp_breadth = breadths, 
                       CDS_length_nt = mean(data$CDS_length_nt, na.rm = TRUE))
  )
  
  # 2. Formal Hypothesis: Does the slope increase with breadth?
  # Compares the last breadth in the vector to the first
  hyp_text <- sprintf("b%d - b1 = 0", length(breadths))
  test_result <- hypotheses(slopes, hypothesis = hyp_text)
  
  return(list(slopes = slopes, test = test_result))
}

visualize_gam_results <- function(model, data, output_prefix = "CDC") {
  require(marginaleffects)
  require(ggplot2)
  
  # Prediction Plot
  p1 <- plot_predictions(model, 
                         condition = c("Max_Log10_Exp", "Exp_breadth"),
                         newdata = datagrid(CDS_length_nt = mean(data$CDS_length_nt))) +
    geom_rug(data = data, aes(x = Max_Log10_Exp), sides = "b", alpha = 0.05) +
    theme_custom() + 
    labs(title = paste("Predicted", output_prefix), subtitle = "Effect of Exp x Breadth")
  
  ggsave(sprintf("./results/GAM_%s_Predictions.pdf", output_prefix), p1, width = 10, height = 6)
  
  return(p1)
}

analyze_evolutionary_constraint <- function(model, data, predictor = "Exp_breadth", output_prefix = "CDC") {
  require(ggplot2)
  require(mgcv)
  
  message("--- Analyzing Evolutionary Constraint (Residual Variance) ---")
  
  # 1. Extract absolute residuals (deviation from predicted optimum)
  data$bias_deviation <- abs(residuals(model, type = "response"))
  
  # 2. Model the deviation (Gamma is great for strictly positive noise)
  # Using a simple smooth to see if 'noise' decreases as breadth increases
  f_noise <- as.formula(paste("bias_deviation ~ s(", predictor, ", bs = 'cs')"))
  m_noise <- gam(f_noise, data = data, family = Gamma(link = "log"))
  
  # 3. Plotting
  p <- ggplot(data, aes_string(x = predictor, y = "bias_deviation")) +
    geom_point(alpha = 0.1, color = "gray60") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
                color = "firebrick", fill = "firebrick") +
    theme_custom() +
    labs(title = "Evolutionary Constraint Analysis",
         subtitle = paste("Deviation from predicted", output_prefix, "across", predictor),
         y = "|Residuals| (Unexplained Variation)",
         x = predictor)
  
  ggsave(sprintf("./results/Constraint_%s_%s.pdf", output_prefix, predictor), 
         p, width = 8, height = 6)
  
  return(list(model = m_noise, plot = p))
}

run_posteriori_gam_analysis <- function(model, data, response_name = "Pi", 
                                        focal_pred = "Max_Log10_Exp", 
                                        interact_pred = "Exp_breadth",
                                        third_pred = NULL,
                                        prefix = "Exp") {
  require(marginaleffects)
  require(ggplot2)
  require(dplyr)
  
  message(sprintf("--- Post-hoc Analysis: %s ---", response_name))
  
  # 1. Define slices for interaction variables
  interact_levels <- quantile(data[[interact_pred]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  
  # 2. Build datagrid list dynamically
  grid_params <- list(model = model, grid_type = "mean_or_mode")
  grid_params[[interact_pred]] <- interact_levels
  
  # Only slice third_pred if it actually exists in the model
  if (!is.null(third_pred)) {
    length_slices <- quantile(data[[third_pred]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
    grid_params[[third_pred]] <- length_slices
    conditions <- c(focal_pred, interact_pred, third_pred)
  } else {
    conditions <- c(focal_pred, interact_pred)
  }
  
  # 3. INTERACTION PLOT
  p_interaction <- plot_predictions(
    model, 
    condition = conditions,
    newdata = do.call(datagrid, grid_params)
  ) +
    theme_custom() +
    scale_color_viridis_d(option = "plasma", name = interact_pred) +
    scale_fill_viridis_d(option = "plasma", name = interact_pred) +
    labs(title = paste("Interaction Analysis:", response_name),
         x = focal_pred, y = paste("Predicted", response_name))
  
  if (!is.null(third_pred)) {
    p_interaction <- p_interaction + 
      facet_wrap(as.formula(paste0("~", third_pred)), labeller = label_both)
  }
  
  # 4. PARTIAL EFFECT PLOT (Link Scale)
  p_partial <- plot_slopes(
    model, 
    variables = focal_pred, 
    condition = focal_pred,
    type = "link" 
  ) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_custom() +
    labs(title = paste("Partial Effect of", focal_pred),
         subtitle = "Controlled for other model predictors",
         y = "Partial Effect (Log Scale)")
  
  # 5. SLOPE ANALYSIS & NUMERIC HYPOTHESIS
  slope_grid_params <- list(model = model, grid_type = "mean_or_mode")
  slope_grid_params[[interact_pred]] <- interact_levels
  
  slopes_base <- avg_slopes(
    model, 
    variables = focal_pred, 
    by = interact_pred,
    newdata = do.call(datagrid, slope_grid_params)
  )
  
  # Use a numeric vector: (1 * Highest Level) - (1 * Lowest Level)
  # This avoids all "pairwise" or string-based checkmate errors
  hyp_matrix <- c(-1, 0, 1) 
  slope_hyp <- hypotheses(slopes_base, hypothesis = hyp_matrix)
  
  message("--- Formal Test: High vs Low Levels of Interaction Predicate ---")
  print(slope_hyp)
  
  # 6. CONSTRAINT ANALYSIS
  data$resid_abs <- abs(residuals(model, type = "response"))
  p_constraint <- ggplot(data, aes(x = .data[[interact_pred]], y = resid_abs)) +
    geom_point(alpha = 0.05, color = "gray60") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "firebrick") +
    theme_custom() +
    labs(title = "Evolutionary Constraint", 
         y = "|Residuals| (Variation)", x = interact_pred)
  
  # Output handling
  if(!dir.exists("./results")) dir.create("./results")
  ggsave(sprintf("./results/%s_%s_Interaction_Plot.pdf", prefix, response_name), p_interaction, width = 10, height = 6)
  ggsave(sprintf("./results/%s_%s_Partial_Effect.pdf", prefix, response_name), p_partial, width = 8, height = 6)
  
  return(list(interaction_plot = p_interaction, 
              partial_plot = p_partial,
              constraint_plot = p_constraint, 
              slopes = slopes_base, 
              hyp_test = slope_hyp))
}

#' Check Confounding and Multicollinearity (L_ROC vs GC)
#'
#' @param data The integrated data frame.
#' @param focal_pred The selection metric (e.g., "L_ROC").
#' @param comp_pred The composition metric (e.g., "GC" or "GC3").
#' @return A list containing correlation stats and a diagnostic plot.
check_confounding_vif <- function(data, focal_pred = "L_ROC", comp_pred = "GC") {
  require(ggplot2)
  require(dplyr)
  
  message(sprintf("--- Analyzing Confounding: %s vs %s ---", focal_pred, comp_pred))
  
  # 1. Spearman Correlation (Non-linear relationship check)
  cor_test <- cor.test(data[[focal_pred]], data[[comp_pred]], method = "spearman",
                       exact = F)
  rho <- cor_test$estimate
  p_val <- cor_test$p.value
  
  # 2. Variance Inflation Factor (VIF)
  # We fit a simple linear model to see how well one predicts the other
  # VIF = 1 / (1 - R^2)
  fit_vif <- lm(as.formula(paste(comp_pred, "~", focal_pred)), data = data)
  r_sq <- summary(fit_vif)$r.squared
  vif_val <- 1 / (1 - r_sq)
  
  # 3. Visualization
  p <- ggplot(data, aes_string(x = focal_pred, y = comp_pred)) +
    geom_bin2d(bins = 100) + # Heatmap style to handle high density
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
    scale_fill_viridis_c(option = "mako") +
    theme_minimal() +
    labs(title = "Confounding Diagnostic",
         subtitle = sprintf("Spearman rho = %.3f | VIF = %.2f", rho, vif_val),
         x = paste("Selection Intensity (", focal_pred, ")"),
         y = paste("Composition (", comp_pred, ")"))
  
  # 4. Recommendation Logic
  recommendation <- if(vif_val > 5) {
    "HIGH COLLINEARITY: Including both may lead to circularity/variance masking."
  } else if(rho > 0.7) {
    "STRONG CORRELATION: Variables are likely confounded; use residuals or pick one."
  } else {
    "MODERATE/LOW: Potentially safe to include GC as a control."
  }
  
  message(sprintf("Recommendation: %s", recommendation))
  
  return(list(rho = rho, p_value = p_val, vif = vif_val, plot = p, 
              recommendation = recommendation))
}

plot_selection_surface <- function(model, data, response_name = "Preferred_Freq",
                                   x_var = "Max_Log10_Exp", 
                                   y_var = "Exp_breadth", 
                                   facet_var = "Total_Codons") {
  require(marginaleffects)
  require(ggplot2)
  require(viridis)
  
  message(sprintf("--- Generating Selection Surface: %s ---", response_name))
  
  # 1. Define specific slices for Total_Codons (the "facets")
  # We use the 10th, 50th, and 90th percentiles to represent Small, Medium, and Large genes
  facet_levels <- quantile(data[[facet_var]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  
  # 2. Generate a dense grid for the surface (e.g., 50x50 per facet)
  grid_args <- list(model = model)
  grid_args[[x_var]] <- seq(min(data[[x_var]]), max(data[[x_var]]), length.out = 50)
  grid_args[[y_var]] <- seq(min(data[[y_var]]), max(data[[y_var]]), length.out = 50)
  grid_args[[facet_var]] <- facet_levels
  
  surface_grid <- do.call(datagrid, grid_args)
  
  # 3. Predict across the surface
  surface_preds <- predictions(model, newdata = surface_grid, type = "response")
  
  # 4. Plot the Heatmap Surface
  p_surface <- ggplot(surface_preds, aes(x = .data[[x_var]], y = .data[[y_var]], fill = estimate)) +
    geom_tile() +
    geom_contour(aes(z = estimate), color = "white", alpha = 0.2) +
    facet_wrap(as.formula(paste0("~", facet_var)), 
               labeller = label_both, ncol = 3) +
    scale_fill_viridis_c(option = "magma", name = "Predicted\nFreq") +
    theme_custom() +
    labs(title = paste("Selection Surface:", response_name),
         subtitle = "Visualizing the 3-way interaction: Intensity x Breadth x Size",
         x = "Expression Intensity (Log10)",
         y = "Expression Breadth (Tissues)") +
    theme(legend.position = "right",
          panel.spacing = unit(1, "lines"))
  
  # Save result
  if(!dir.exists("./results")) dir.create("./results")
  ggsave(sprintf("./results/Surface_%s.pdf", response_name), p_surface, width = 14, height = 5)
  
  return(p_surface)
}