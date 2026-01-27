analyze_nonlinearity <- function(pred, data) {
  
  tryCatch({
    # Formulas
    form_glm <- as.formula(paste0("CDC_beta ~ ", pred))
    form_gam <- as.formula(paste0("CDC_beta ~ s(", pred, ")"))
    
    # A. Fit Models
    # GLM (Linear Beta) - Note: link.phi=NULL uses constant precision
    model_glm <- betareg(form_glm, data = data, link = "logit")
    
    # GAM (Smooth Beta)
    model_gam <- gam(form_gam, data = data, family = betar(link = "logit"), 
                     method = "REML")
    
    # B. Compare AIC (The Gold Standard for Beta)
    aic_glm <- AIC(model_glm)
    aic_gam <- AIC(model_gam)
    
    delta_aic <- aic_glm - aic_gam # Positive = GLM is worse
    
    # C. Logic
    # If GAM is more than 2 AIC units better, we prefer it.
    is_nonlinear <- delta_aic > 2.0
    
    return(tibble::tibble(
      Predictor = pred,
      GLM_AIC = round(aic_glm, 2),
      GAM_AIC = round(aic_gam, 2),
      Delta_AIC = round(delta_aic, 2),
      Recommendation = ifelse(is_nonlinear, "GAM (Non-Linear)", "GLM (Linear)"),
      Note = "AIC Selection"
    ))
    
  }, error = function(e) {
    return(tibble::tibble(
      Predictor = pred,
      GLM_AIC = NA_real_,
      GAM_AIC = NA_real_,
      Delta_AIC = NA_real_,
      Recommendation = "ERROR",
      Note = as.character(e$message)
    ))
  })
}