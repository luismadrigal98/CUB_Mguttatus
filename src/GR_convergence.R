get_mcmc_trace <- function(dir, pattern = "R_objects/parameter.Rda") 
{
  #' Auxiliary function to get the traces for every parameter of interest.
  #' 
  #' @param dir Directory that stores individual runs of the mcmc
  #' @param pattern Pattern to identify the parameter object inside each of the
  #' run directories.
  #' 
  #' @return Traces in independent lists. Each trace is represented in a matrix,
  #' where each column correspond to an iteration and each row a codon.
  
  # Load parameters
  param <- admisc::listRDA((paste(dir, pattern, sep="/")))
  
  # Get and format individual traces
  sel <- do.call('cbind', param$selectionTrace[[1]])
  mut <- do.call('cbind', param$mutationTrace[[1]])
  
  trace <- list(selection = sel,
                mutation = mut)
  
  return(trace)
}

GR_convergence <- function(dirs_run, parameter = 'all', burn_in = 1000) {
  #' Function to calculate the Gelman and Rubin's convergence diagnostic metric
  #' 
  #' @param dirs_run Vector of directory paths.
  #' @param parameter "all", "selection", or "mutation".
  #' @param burn_in Number of iterations to discard from start of chain.
  
  assertthat::assert_that(parameter %in% c("all", "selection", "mutation"),
                          msg = "Undefined parameter used")
  
  # 1. Extract chains for the requested parameter
  raw_chains <- lapply(dirs_run, function(x) {
    get_mcmc_trace(x)
  })
  
  # 2. Identify which trace types we have (selection, mutation, or both)
  trace_types <- names(raw_chains[[1]])
  if(parameter != 'all')
  {
    trace_types <- trace_types[trace_types == parameter]
  }
  
  # 3. Calculate GR Stat for each trace type
  GR_stats_list <- lapply(trace_types, function(type) {
    
    # Extract the specific type (e.g., selection) from all runs
    # And convert to coda::mcmc objects immediately
    mcmc_chain_list <- lapply(raw_chains, function(run) {
      
      data <- run[[type]]
      
      # Remove burn-in if rows are sufficient
      if(nrow(data) > burn_in) {
        data <- data[(burn_in + 1):nrow(data), ]
      }
      
      # Convert to MCMC object for coda
      coda::as.mcmc(data)
    })
    
    # Create the mcmc.list required by gelman.diag
    combined_chains <- coda::mcmc.list(mcmc_chain_list)
    
    # Calculate Gelman-Rubin
    # multivariate=FALSE gives us the stat for each parameter individually
    diag <- coda::gelman.diag(combined_chains, multivariate = FALSE)
    
    return(diag$psrf)
  })
  
  names(GR_stats_list) <- trace_types
  return(GR_stats_list)
}