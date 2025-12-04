# =============================================================================
# ROC Model Trajectory Plotting Functions
# =============================================================================
# Functions to visualize codon frequency trajectories across expression levels
# using CSP parameters (dM, dEta) from AnaCoDa ROC model
# =============================================================================

library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Load CSP parameters from AnaCoDa MCMC output
# -----------------------------------------------------------------------------
# AnaCoDa outputs separate files for mutation (dM) and selection (dEta)
# Format: columns are codons, single row with parameter estimates
# -----------------------------------------------------------------------------
load_csp_parameters <- function(mutation_file, selection_file) {
  
  dM_raw <- read.csv(mutation_file, stringsAsFactors = FALSE)
  dEta_raw <- read.csv(selection_file, stringsAsFactors = FALSE)
  
  # Get codon names from column headers (remove any prefix)
  codon_cols_dM <- colnames(dM_raw)
  codon_cols_dEta <- colnames(dEta_raw)
  
  # Extract values
  dM_values <- as.numeric(dM_raw[1, ])
  dEta_values <- as.numeric(dEta_raw[1, ])
  
  # Create data frame
  csp <- data.frame(
    Codon = codon_cols_dM,
    dM = dM_values,
    dEta = dEta_values,
    stringsAsFactors = FALSE
  )
  
  # Standard genetic code for AA mapping
  codon_table <- c(
    "TTT" = "F", "TTC" = "F",
    "TTA" = "L", "TTG" = "L", "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I",
    "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "AGT" = "S", "AGC" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y",
    "CAT" = "H", "CAC" = "H",
    "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N",
    "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D",
    "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C",
    "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  # Add AA column
  csp$AA <- codon_table[csp$Codon]
  
  # AnaCoDa convention: "Z" for Ser AGC/AGT codons
  csp$AA_anacoda <- ifelse(csp$Codon %in% c("AGC", "AGT"), "Z", csp$AA)
  
  # Handle reference codons (NA becomes 0)
  csp$dM[is.na(csp$dM)] <- 0
  csp$dEta[is.na(csp$dEta)] <- 0
  
  return(csp)
}

# -----------------------------------------------------------------------------
# Map amino acids to AnaCoDa convention in codon frequency data frame
# -----------------------------------------------------------------------------
map_aa_to_anacoda <- function(codon_freq_df) {
  
  codon_freq_df$AA_anacoda <- ifelse(
    codon_freq_df$Codon %in% c("AGC", "AGT"), 
    "Z", 
    codon_freq_df$AA
  )
  
  return(codon_freq_df)
}

# -----------------------------------------------------------------------------
# Calculate observed codon frequencies per gene from CDS sequences
# -----------------------------------------------------------------------------
calculate_observed_codon_frequencies <- function(cds_seqs) {
  
  # Standard genetic code
  codon_table <- c(
    "TTT" = "F", "TTC" = "F",
    "TTA" = "L", "TTG" = "L", "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I",
    "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "AGT" = "S", "AGC" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y",
    "CAT" = "H", "CAC" = "H",
    "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N",
    "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D",
    "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C",
    "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G",
    "TAA" = "Stop", "TAG" = "Stop", "TGA" = "Stop"
  )
  
  result_list <- lapply(seq_along(cds_seqs), function(i) {
    gene_name <- names(cds_seqs)[i]
    seq_str <- as.character(cds_seqs[[i]])
    
    # Split into codons
    seq_len <- nchar(seq_str)
    n_codons <- seq_len %/% 3
    
    if (n_codons == 0) return(NULL)
    
    codons <- substring(seq_str, 
                        seq(1, n_codons * 3, by = 3), 
                        seq(3, n_codons * 3, by = 3))
    
    # Count codons
    codon_counts <- table(codons)
    
    # Create data frame
    df <- data.frame(
      Gene = gene_name,
      Codon = names(codon_counts),
      Count = as.numeric(codon_counts),
      stringsAsFactors = FALSE
    )
    
    # Add AA
    df$AA <- codon_table[df$Codon]
    
    # Remove stop codons and any unknown
    df <- df[!is.na(df$AA) & df$AA != "Stop", ]
    
    # Calculate frequency within AA family
    df <- df %>%
      group_by(Gene, AA) %>%
      mutate(
        AA_total = sum(Count),
        Observed_freq = Count / AA_total
      ) %>%
      ungroup()
    
    return(df)
  })
  
  result <- do.call(rbind, result_list)
  return(result)
}

# -----------------------------------------------------------------------------
# Predict codon probabilities for a batch of genes
# -----------------------------------------------------------------------------
# ROC multinomial model: P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
# where Z normalizes within each amino acid family
# -----------------------------------------------------------------------------
predict_codon_probs_batch <- function(gene_expression_df, csp, phi_col = "Exp_log10") {
  
  aa_families <- unique(csp$AA_anacoda)
  
  result_list <- lapply(1:nrow(gene_expression_df), function(i) {
    gene <- gene_expression_df$Gene[i]
    phi <- gene_expression_df[[phi_col]][i]
    
    # Predict for each AA family
    aa_results <- lapply(aa_families, function(aa) {
      csp_aa <- csp[csp$AA_anacoda == aa, ]
      
      if (nrow(csp_aa) <= 1) return(NULL)  # Skip non-degenerate
      
      # Calculate unnormalized log probabilities
      log_unnorm <- -csp_aa$dM - csp_aa$dEta * phi
      
      # Softmax normalization
      max_log <- max(log_unnorm)
      unnorm <- exp(log_unnorm - max_log)
      probs <- unnorm / sum(unnorm)
      
      data.frame(
        Gene = gene,
        AA = aa,
        Codon = csp_aa$Codon,
        predicted_prob = probs,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, aa_results)
  })
  
  result <- do.call(rbind, result_list)
  return(result)
}
