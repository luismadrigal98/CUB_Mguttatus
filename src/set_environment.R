set_environment <- function(required_pckgs,
                            automatic_download = FALSE,
                            personal_seed = as.numeric(Sys.time()),
                            parallel_backend = FALSE,
                            n_cores = NULL,
                            src_dir = './src')
{
  #' This fucntion will set up the working environment for performing all the
  #' analysis. It will load all the required packages, set the seed for
  #' reproducibility and set the parallel backend if required.
  #' 
  #' @param required_pckgs A character vector of the required packages.
  #' 
  #' @param automatic_download A logical value to set the automatic download of
  #' the required packages. Default set to FALSE.
  #' 
  #' @param personal_seed An integer to set the seed for reproducibility. Default
  #' set to system time.
  #' 
  #' @param parallel_backend A logical value to set the parallel backend. Default
  #' set to FALSE.
  #' 
  #' @return invisible
  #' ___________________________________________________________________________
  
  ## Loading the required libraries ----
  
  message("Loading the required libraries")
  
  tryCatch(
    {
      for(pckg in required_pckgs)
      {
        
        if(!require(pckg, character.only = TRUE))
        {
          if (automatic_download == TRUE)
          {
            install.packages(pckg)
            library(pckg, character.only = TRUE)
          }
          else
          {
            message(paste0("The automatic download of the required package ", 
                           pckg,  " is disabled"))
            message("Install it manually and run the script again")
            stop()
          }
        }
        else
        {
          library(pckg, character.only = TRUE)
        }
      }
      
      message("The required libraries have been loaded")
    }, error = function(e) {
      message("Some packages cannot be installed through CRAN, trying with remotes")
      print(e)
    })
  
  ## Setting the seed ----
  set.seed(personal_seed)
  
  ## Setting the parallel backend ----
  if(parallel_backend == TRUE)
  {
    require(future)
    require(future.apply)
    require(doParallel)
    require(doFuture)
    
    ## 1.2) Setting the parallel backend
    options(doFuture.rng.onMisuse = "ignore")
    registerDoFuture()
    plan(multisession, workers = ifelse(!is.null(n_cores), n_cores, parallelly::availableCores() - 1))
  }
  
  # Source required functions
  message("Sourcing the required functions")
  source_files <- list.files(path = src_dir, pattern = "*.R", full.names = TRUE)
  source_files <- setdiff(source_files, './src/set_environment.R')
  
  sapply(source_files, function(x) source(x))
  
  return(invisible())
}