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
  # Exclude set_environment.R using basename so the match works with any path form
  source_files <- source_files[basename(source_files) != "set_environment.R"]

  is_executable_script <- function(file_path) {
    file_text <- tryCatch(
      paste(readLines(file_path, warn = FALSE, n = 80), collapse = "\n"),
      error = function(e) ""
    )
    # Any CLI marker triggers exclusion
    grepl("^#!.*Rscript", file_text, perl = TRUE) ||
      grepl("commandArgs\\s*\\(", file_text, perl = TRUE) ||
      grepl("Usage:\\s*Rscript", file_text, perl = TRUE) ||
      grepl("\\bquit\\s*\\(", file_text, perl = TRUE) ||
      grepl("\\bq\\s*\\(\\s*['\"]", file_text, perl = TRUE)  # q("no"), q('yes'), etc.
  }

  exec_flags <- vapply(source_files, is_executable_script, logical(1))
  if (any(exec_flags)) {
    message(sprintf("Skipping %d executable script(s): %s",
                    sum(exec_flags),
                    paste(basename(source_files[exec_flags]), collapse = ", ")))
  }
  source_files <- source_files[!exec_flags]

  # Defensive override: mask base::quit / base::q in globalenv for the duration
  # of sourcing, so a stray quit() in a sourced file warns instead of tearing
  # down the R session (which triggers the "Save workspace image?" prompt).
  had_quit <- exists("quit", envir = globalenv(), inherits = FALSE)
  had_q    <- exists("q",    envir = globalenv(), inherits = FALSE)
  old_quit <- if (had_quit) get("quit", envir = globalenv()) else NULL
  old_q    <- if (had_q)    get("q",    envir = globalenv()) else NULL
  on.exit({
    if (had_quit) assign("quit", old_quit, envir = globalenv()) else suppressWarnings(rm("quit", envir = globalenv()))
    if (had_q)    assign("q",    old_q,    envir = globalenv()) else suppressWarnings(rm("q",    envir = globalenv()))
  }, add = TRUE)

  safe_quit <- function(...) {
    warning("Suppressed quit()/q() call from a sourced file", call. = FALSE)
    invisible(NULL)
  }
  assign("quit", safe_quit, envir = globalenv())
  assign("q",    safe_quit, envir = globalenv())

  for (f in source_files) {
    tryCatch(
      source(f, local = FALSE),
      error = function(e) {
        warning(sprintf("Error sourcing %s: %s", f, conditionMessage(e)),
                call. = FALSE)
      }
    )
  }

  return(invisible())
}