extract_traces <- function(objects, codon_idx) {
  df_list <- list()
  
  for (run_name in names(objects)) {
    # Access logic: Run -> selectionTrace -> Mixture 1 -> Codon Index
    trace_data <- objects[[run_name]]$selectionTrace[[1]][[codon_idx]]
    
    # Create temporary DF
    temp_df <- data.frame(
      Iteration = 1:length(trace_data),
      Value = trace_data,
      Run = run_name
    )
    df_list[[run_name]] <- temp_df
  }
  
  return(do.call(rbind, df_list))
}