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
  #' @param n_cores An integer to set the number of cores to be used for the
  #' parallel backend. Default set to NULL, which will use all the available
  #' cores minus one.
  #' 
  #' @param src_dir A character string with the path to the source directory
  #' containing the required functions. Default set to './src'.
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
  source_files <- list.files(path = src_dir, pattern = "\\.R$", full.names = TRUE)
  source_files <- setdiff(source_files, './src/set_environment.R')
  
  is_executable_script <- function(file_path) {
    file_text <- tryCatch(
      paste(readLines(file_path, warn = FALSE, n = 80), collapse = "\n"),
      error = function(e) ""
    )
    grepl("commandArgs\\s*\\(", file_text, perl = TRUE) ||
      grepl("Usage: Rscript", file_text, fixed = TRUE)
  }
  
  source_files <- source_files[!vapply(source_files, is_executable_script, logical(1))]
  
  sapply(source_files, function(x) source(x))
  
  return(invisible())
}