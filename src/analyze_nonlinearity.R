#' Analyze Non-Linearity of Predictor
#' 
#' Evaluates linear vs non-linear behavior by comparing a GLM-style GAM 
#' to a smooth GAM using the same family.
#' 
#' @param resp Character string of the response variable.
#' @param pred Character string of the predictor variable.
#' @param data Data frame containing variables.
#' @param family The mgcv family object (e.g., betar(), Gamma(link="log"), nb()).
#' @return A tibble with AIC comparison and model recommendation.

analyze_nonlinearity <- function(resp, pred, data, family = betar(link = "logit")) {
  require(mgcv)
  require(tibble)
  
  tryCatch({
    # 1. Construct Formulas
    # Linear model (GLM equivalent in GAM framework)
    form_linear <- as.formula(paste0(resp, " ~ ", pred))
    # Non-linear model (Smooth spline)
    form_smooth <- as.formula(paste0(resp, " ~ s(", pred, ")"))
    
    # 2. Fit Models using the SAME engine (mgcv::gam)
    # This ensures AIC is calculated using the same likelihood definitions
    model_linear <- gam(form_linear, data = data, family = family, method = "REML")
    model_smooth <- gam(form_smooth, data = data, family = family, method = "REML")
    
    # 3. Compare AIC
    aic_lin <- AIC(model_linear)
    aic_sm  <- AIC(model_smooth)
    delta_aic <- aic_lin - aic_sm # Positive value means Smooth is better
    
    # 4. Determine recommendation
    # Threshold of 2 is standard for "significant" improvement in parsimony
    is_nonlinear <- delta_aic > 2.0
    
    return(tibble(
      Response       = resp,
      Predictor      = pred,
      Family         = family$family,
      Linear_AIC     = round(aic_lin, 2),
      Smooth_AIC     = round(aic_sm, 2),
      Delta_AIC      = round(delta_aic, 2),
      Recommendation = ifelse(is_nonlinear, "GAM (Non-Linear)", "GLM (Linear)")
    ))
    
  }, error = function(e) {
    return(tibble(
      Response       = resp,
      Predictor      = pred,
      Family         = NA_character_,
      Linear_AIC     = NA_real_,
      Smooth_AIC     = NA_real_,
      Delta_AIC      = NA_real_,
      Recommendation = "ERROR",
      Note           = as.character(e$message)
    ))
  })
}

#' Batch Analyze Non-Linearity for Multiple Predictors
#' 
#' Iterates through a vector of predictors and compares Linear vs Smooth fits.
#'
#' @param resp Character string of the response variable.
#' @param predictors Character vector of predictor names.
#' @param data Data frame containing variables.
#' @param family The mgcv family object.
#' @return A consolidated tibble with results for all predictors.

analyze_nonlinearity_suite <- function(resp, predictors, data, family = betar()) {
  require(dplyr)
  require(purrr)
  
  # Map the analysis function over the vector of predictors
  results <- map_df(predictors, function(p) {
    
    message(sprintf("Testing non-linearity: %s ~ %s", resp, p))
    
    tryCatch({
      # 1. Formulas
      f_lin <- as.formula(paste0(resp, " ~ ", p))
      f_sm  <- as.formula(paste0(resp, " ~ s(", p, ")"))
      
      # 2. Fit with mgcv
      m_lin <- mgcv::gam(f_lin, data = data, family = family, method = "REML")
      m_sm  <- mgcv::gam(f_sm,  data = data, family = family, method = "REML")
      
      # 3. Extract Stats
      aic_lin <- AIC(m_lin)
      aic_sm  <- AIC(m_sm)
      delta   <- aic_lin - aic_sm
      
      # EDF tells us how 'wiggly' the spline is (1 = linear)
      edf <- summary(m_sm)$s.table[1, "edf"]
      
      tibble::tibble(
        Response = resp,
        Predictor = p,
        Delta_AIC = delta,
        EDF = round(edf, 2),
        Recommendation = ifelse(delta > 2 & edf > 1.1, "GAM (Non-Linear)", "GLM (Linear)"),
        Status = "Success"
      )
      
    }, error = function(e) {
      tibble::tibble(Predictor = p, Status = "Error", Note = e$message)
    })
  })
  
  return(results)
}