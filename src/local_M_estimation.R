#! /usr/bin/env Rscript
#' Local estimation of M across genomic windows using introns as proxy for unconstrained
#' sequences
#' 
#' @author Luis Javier Madrigal-Roca
#' 
#' @date 2025/11/20
#' _____________________________________________________________________________ 

# ******************************************************************************
# STEP 1: Extract and Filter Introns ----
# ______________________________________________________________________________

get_intron_sequences <- function(fasta_file, ann_file,
                                 trim_bp = 30,
                                 min_width = 86,
                                 organism = NA)
{
  #' @title Extract Trimmed Deep Intron Sequences
  #'
  #' @description
  #' Loads an annotation file (GFF3) and a genome FASTA file, filters for primary 
  #' transcripts on main chromosomes, trims the splice site regions of introns,
  #' and extracts the resulting deep intronic sequences.
  #'
  #' @param fasta_file Character string. Path to the genome reference FASTA file.
  #' @param ann_file Character string. Path to the gene annotation file (GFF3).
  #' @param trim_bp Numeric. Number of base pairs to trim from both the 5' and 3' ends of each intron
  #'   (to remove splice site signals). Default is 30.
  #' @param min_width Numeric. Minimum required intron width *after* trimming (i.e., total required length 
  #'   is min_width + 2*trim_bp). Introns shorter than this are discarded. Default is 86.
  #' @param organism Character string. The organism name for the TxDb object metadata 
  #'   (e.g., "Mimulus guttatus"). Default is NA.
  #'
  #' @details
  #' \strong{Filtering Assumptions:}
  #' \itemize{
  #'   \item \strong{Primary Transcript:} The transcript ID ends with ".1" (e.g., "GeneID.1").
  #'   \item \strong{Main Chromosomes:} Chromosome names follow the pattern "^Chr_[0-9]{1,2}$".
  #' }
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \strong{intron_seqs:} A \code{DNAStringSet} of the extracted, trimmed deep intronic sequences.
  #'   \item \strong{trimmed_introns:} A \code{GRanges} object with the coordinates of the extracted sequences.
  #' }
  #' @import GenomicFeatures txdbmaker Rsamtools GenomicRanges Biostrings S4Vectors
  #' @export
  #' ___________________________________________________________________________
  
  message("Step 1: Loading Annotation and Extracting Introns...")
  
  # --- 1. Load Annotation and Create TxDb ---
  txdb <- txdbmaker::makeTxDbFromGFF(file = ann_file, 
                                     format = "gff3",
                                     organism = organism)
  
  introns_list <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
  
  # --- 2. Filter for Primary Transcripts (.1) ---
  transcript_ids <- names(introns_list)
  # Look for IDs ending in ".1"
  is_primary_transcript <- sapply(strsplit(transcript_ids, "\\."), 
                                  function(x) tail(x, 1) == "1")
  primary_introns_list <- introns_list[is_primary_transcript]
  
  message(paste("Original number of transcripts:", length(introns_list)))
  message(paste("Number of primary transcripts (.1):", length(primary_introns_list)))
  
  # --- 3. Filter for Main Chromosomes (e.g., Chr_01) ---
  
  # Define main chromosome pattern based on common annotation formats
  main_chr_pattern <- "^Chr_[0-9]{1,2}$"
  
  # Extract the single chromosome name for each transcript
  seq_names <- sapply(primary_introns_list, function(gr) {
    return(as.character(GenomeInfoDb::seqlevelsInUse(gr)))
  })
  
  # Filter the GRangesList
  is_main_chromosome <- grepl(main_chr_pattern, seq_names)
  final_introns_list <- primary_introns_list[is_main_chromosome]
  
  message(paste("Final list contains", length(final_introns_list), 
                "intron sets for primary transcripts on main chromosomes."))
  
  # --- 4. Intron Width Filtering and Trimming ---
  
  all_introns <- unlist(final_introns_list)
  
  # Filtering: Check if the intron is long enough to survive trimming
  required_min_width <- min_width + (2 * trim_bp)
  clean_introns <- all_introns[width(all_introns) > required_min_width]
  
  # Trimming: Remove splice sites from both ends
  trimmed_introns <- GenomicRanges::narrow(clean_introns, 
                                           start = trim_bp + 1, 
                                           end = width(clean_introns) - trim_bp)
  
  # --- 5. Prepare Genome FASTA for Sequence Extraction (FaFile) ---
  
  # Use FaFile for robust coordinate-based sequence extraction (required by getSeq)
  if (!file.exists(fasta_file)) stop("FASTA file not found.")
  
  dna <- readDNAStringSet(filepath = fasta_file)
  
  original_names <- names(dna)
  
  # Simplify the names
  names(dna) <- sub("^(\\S+)\\s.*", "\\1", original_names)
  
  main_chroms_to_keep <- names(dna)[grepl(main_chr_pattern, 
                                                 names(dna))]
  
  dna <- dna[main_chroms_to_keep]
  
  # --- 6. Final Synchronization and Sequence Extraction ---
  
  # Crucial: Synchronize the GRanges seqlevels to the cleaned FASTA names
  trimmed_introns <- GenomeInfoDb::keepSeqlevels(
    trimmed_introns, 
    main_chroms_to_keep, 
    pruning.mode = "coarse"
  )
  
  # --- CUSTOM DIRECT SEQUENCE EXTRACTION (Your requested method) ---
  
  # 1. Get components for extraction
  r <- GenomicRanges::ranges(trimmed_introns)
  seq_names <- as.character(GenomicRanges::seqnames(trimmed_introns))
  
  intron_seqs_raw <- Biostrings::DNAStringSet(
    sapply(1:length(r), function(i) {
      # Extract the sequence for the current chromosome (seq_names[i])
      chrom_seq <- dna[seq_names[i]]
      
      # Subseq the single-chromosome DNAStringSet
      subseq_obj <- Biostrings::subseq(chrom_seq, start(r[i]), end(r[i]))

      return(as.character(subseq_obj))
    })
  )
  
  # 3. Handle Reverse Complementation for Negative Strand
  neg_strand_idx <- GenomicRanges::strand(trimmed_introns) == "-"
  
  intron_seqs_final <- intron_seqs_raw
  intron_seqs_final[neg_strand_idx] <- Biostrings::reverseComplement(intron_seqs_final[neg_strand_idx])
  
  names(intron_seqs_final) <- names(trimmed_introns)
  
  # Get the Seqinfo object for the next function call
  genome_seqinfo <- GenomeInfoDb::seqinfo(dna) 
  
  return(list(
    # The full genome DNAStringSet (optional, but requested earlier)
    dna = dna, 
    intron_seqs = intron_seqs_final, 
    trimmed_introns = trimmed_introns,
    genome_seqinfo = genome_seqinfo
  ))
}

# ******************************************************************************
# STEP 2: Measure Base Composition per Window ----
# ______________________________________________________________________________

# Create a mapping function to aggregate counts

calculate_window_metrics <- function(window_idx, 
                                     all_windows, 
                                     all_seqs, 
                                     hit_list,
                                     return_Ns = TRUE) 
{
  #' @title Calculate Base Composition for a Single Genomic Window
  #' @description Internal function to extract intron sequences within a window 
  #' and calculate A, C, G, T frequencies.
  #' @param window_idx Numeric. Index of the current window.
  #' @param all_windows GRanges object defining all genomic windows.
  #' @param all_seqs DNAStringSet containing all trimmed intron sequences.
  #' @param hit_list Hits object from findOverlaps mapping windows to introns.
  #' @param return_Ns Whether to return also the count of Ns if working with a 
  #' masked fasta file
  #' @return A named vector of window metadata, total base pairs, and base frequencies, 
  #' or NULL if empty.
  #' ___________________________________________________________________________
  
  # Identify which introns belong to this specific window index
  intron_indices <- S4Vectors::subjectHits(hit_list)[S4Vectors::queryHits(hit_list) == window_idx]
  
  if (length(intron_indices) == 0) {
    return(NULL) # Skip empty windows
  }
  
  # Extract sequences for this window
  local_seqs <- all_seqs[intron_indices]
  
  # Count bases (A, C, G, T)
  # baseOnly = TRUE returns 5 columns: A, C, G, T, and "other" (which includes N)
  counts <- Biostrings::alphabetFrequency(local_seqs, baseOnly = TRUE, as.prob = FALSE)
  
  # Sum counts across all introns in this window
  total_counts <- colSums(counts)
  
  # Calculate Total Valid Base Pairs (Strictly A + C + G + T)
  # This effectively ignores 'N's from the hard masking
  total_bp <- sum(total_counts[1:4])
  
  # CRITICAL CHECK: If the window was mostly 'N's, total_bp might be tiny or zero.
  # We return NULL if there are no valid bases to avoid division by zero or garbage data.
  if(total_bp == 0) return(NULL)
  
  # Calculate Frequencies (Pi vector) based on VALID bases only
  freqs <- total_counts[1:4] / total_bp
  names(freqs) <- c("pi_A", "pi_C", "pi_G", "pi_T")
  
  # Get window metadata
  current_window <- all_windows[window_idx]
  window_data <- c(
    seqnames = as.character(GenomicRanges::seqnames(current_window)),
    start = GenomicRanges::start(current_window),
    end = GenomicRanges::end(current_window),
    window_idx = window_idx
  )
  
  # Optionally return N count
  if (return_Ns) N_count <- as.numeric(total_counts["other"]) else N_count <- NA
  
  # Return data including the valid base pair count
  return(c(window_data, total_bp = total_bp, freqs, N_count = N_count))
}

get_base_composition_per_windows <- function(genome_seqinfo, 
                                             trimmed_introns,
                                             intron_seqs,
                                             window_size = 100000)
{
  #' @title Aggregate Intron Base Composition per Window
  #'
  #' @description Tiles the genome into windows and calculates the A, C, G, T 
  #' frequencies of all trimmed deep introns falling within each window.
  #'
  #' @param genome_seqinfo GenomeInfoDb::Seqinfo object with chromosome names and lengths.
  #' @param trimmed_introns GenomicRanges::GRanges object of trimmed introns.
  #' @param intron_seqs Biostrings::DNAStringSet of trimmed intron sequences.
  #' @param window_size Numeric. The size of the genomic windows (default: 100000 bp).
  #' @return A data frame with window coordinates and base frequencies (pi_A, pi_C, pi_G, pi_T).
  #' @export
  #' ___________________________________________________________________________
  
  message("Step 1: Calculating Base Frequencies per Window...")
  
  # Define genomic windows
  # FIX: Use qualified function call and correct input name
  seq_lengths <- GenomeInfoDb::seqlengths(genome_seqinfo)
  windows <- GenomicRanges::tileGenome(seq_lengths, 
                                       tilewidth = window_size, 
                                       cut.last.tile.in.chrom = TRUE)
  
  # Map introns to windows
  # We use 'findOverlaps' to see which introns fall into which 100kb window
  overlaps <- GenomicRanges::findOverlaps(windows, trimmed_introns)
  
  # Apply calculation across windows (Use lapply for list, then bind)
  results_list <- lapply(1:length(windows), calculate_window_metrics, 
                         all_windows = windows, 
                         all_seqs = intron_seqs, 
                         hit_list = overlaps)
  
  # Remove NULL (empty windows) and Convert to Data Frame
  df_results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  df_results <- as.data.frame(df_results)
  
  # Ensure numeric columns are cast correctly (do.call(rbind) converts to matrix, 
  # making them character/factor)
  numeric_cols <- c("start", "end", "window_idx", "pi_A", "pi_C", "pi_G", "pi_T",
                    "total_bp", "N_count")
  for(col in numeric_cols) {
    df_results[[col]] <- as.numeric(df_results[[col]])
  }
  
  return(df_results)
}

refine_windows_for_genes <- function(nuc_composition, min_bp = 1000) 
{
  #' @title Filter Low-Coverage Windows and Perform Gap-Filling
  #'
  #' @description Filters out genomic windows with insufficient coverage (total_bp < min_bp),
  #' then collapses the remaining, high-confidence windows into a new GRanges object
  #' that maintains the original genomic coverage structure for mapping to genes.
  #'
  #' @param nuc_composition Data frame containing window composition and total_bp (output of get_base_composition_per_windows).
  #' @param min_bp Numeric. The minimum required total base pairs of intron sequence within a window to be retained. Default is 1000.
  #' @return A GRanges object representing the high-confidence windows, with original composition data attached as metadata.
  #' @export
  #' ___________________________________________________________________________
  
  message(paste("Original number of windows:", nrow(nuc_composition)))
  
  # --- 1. Filter Windows by Coverage ---
  
  # Filter out windows with total_bp less than the threshold
  high_coverage_df <- nuc_composition[nuc_composition$total_bp >= min_bp, ]
  
  message(paste("Windows retained (total_bp >= ", min_bp, "):", nrow(high_coverage_df)))
  
  if (nrow(high_coverage_df) == 0) {
    stop("No windows met the minimum coverage threshold.")
  }
  
  # --- 2. Convert Filtered Data Frame to GRanges ---
  # We must operate on coordinates using Bioconductor objects.
  
  # Ensure seqnames are treated as factors/Rle-compatible strings
  high_coverage_df$seqnames <- as.character(high_coverage_df$seqnames)
  
  # Create the GRanges object using the start and end coordinates
  high_coverage_gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = high_coverage_df,
    keep.extra.columns = TRUE, # Keep all pi_X, total_bp, etc. as metadata
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end"
  )
  
  # --- 3. Perform Gap-Filling (Merging Adjacent Windows) ---
  
  # The 'reduce' function merges adjacent or overlapping ranges.
  # Since the original windows were tiled (adjacent), reducing them 
  # will merge contiguous high-coverage windows into longer segments, 
  # effectively filling the gaps left by the discarded low-coverage windows.
  
  # NOTE: The metadata (pi_A, total_bp) is LOST here, as 'reduce' only 
  # combines coordinates. This result is used only for the *final* mapping step.
  # The actual composition data for a gene will be derived by mapping the gene 
  # back to the *original* high_coverage_gr object.
  
  # Resulting object is for determining gene assignments (locus)
  gap_filled_gr <- GenomicRanges::reduce(high_coverage_gr)
  
  message(paste("Final number of merged/reduced windows:", length(gap_filled_gr)))
  
  # Return the full set of high-coverage original windows AND the gap-filled object
  return(list(
    high_coverage_windows = high_coverage_gr,
    gap_filled_locus = gap_filled_gr
  ))
}

# ******************************************************************************
# STEP 3: Solve for Mutation Rate Matrix (Q) ----
# ______________________________________________________________________________

message("Step 3: Solving for Mutation Matrix Q...")

# Function to generate the Q matrix from Equilibrium Frequencies
# This assumes the HKY85 model structure, which is robust for this data type.
# It solves Q such that pi * Q = 0.
# Q_ij = pi_j (if transversion)
# Q_ij = kappa * pi_j (if transition)

solve_Q_matrix <- function(pi_A, pi_C, pi_G, pi_T, kappa = 2) {
  
  # Define the matrix rows/cols: A, C, G, T
  # Transitions: A<->G, C<->T
  # Transversions: All others
  
  # Initialize 4x4 matrix
  Q <- matrix(0, nrow = 4, ncol = 4, 
              dimnames = list(c("A","C","G","T"), c("A","C","G","T")))
  
  freqs <- c(A=pi_A, C=pi_C, G=pi_G, T=pi_T)
  
  # Fill off-diagonals
  bases <- c("A", "C", "G", "T")
  for (i in bases) {
    for (j in bases) {
      if (i == j) next
      
      is_transition <- (i=="A" & j=="G") | (i=="G" & j=="A") | 
        (i=="C" & j=="T") | (i=="T" & j=="C")
      
      rate <- freqs[j] # Basic F81 assumption: rate depends on target freq
      
      if (is_transition) {
        rate <- rate * kappa # HKY85 adjustment
      }
      
      Q[i, j] <- rate
    }
  }
  
  # Fill diagonals (Rows must sum to 0)
  diag(Q) <- -rowSums(Q)
  
  # Scale matrix so that -sum(pi_i * Q_ii) = 1 (Standard normalization)
  # This ensures unit time represents 1 substitution per site
  scaling_factor <- -sum(freqs * diag(Q))
  Q_normalized <- Q / scaling_factor
  
  return(Q_normalized)
}

apply_q_matrix_to_windows <- function(nuc_composition_df, kappa = 2) {
  
  message(paste("Calculating Q matrices for", nrow(nuc_composition_df), 
                "windows..."))
  
  # Use dplyr to iterate over rows and apply the function
  df_with_q_matrix <- nuc_composition_df %>%
    # 1. Start row-wise operation
    rowwise() %>%
    
    # 2. Mutate to create a new list column named 'Q_matrix'
    mutate(
      Q_matrix = list(solve_Q_matrix(
        pi_A = pi_A, 
        pi_C = pi_C, 
        pi_G = pi_G, 
        pi_T = pi_T, 
        kappa = kappa
      ))
    ) %>%
    
    # 3. Stop row-wise operation (returns standard data frame structure)
    ungroup()
  
  return(df_with_q_matrix)
}

plot_genomic_rate_variation <- function(df_plot) {
  
  #' @title Plot Mutational Pressure and GC Content (with Legend)
  
  require(ggplot2)
  require(scales)
  
  # Define nice labels for the facet rows
  facet_labels <- c(
    "GC_content" = "Intron GC Content",
    "Q_AG_rate" = "A -> G Transition Rate"
  )
  
  p <- ggplot(df_plot, aes(x = midpoint, y = Rate_Value, color = Variable)) +
    
    # 1. Raw Data: Use faint lines
    geom_line(alpha = 0.4, linewidth = 0.2) +
    
    # 2. Trend: Add a smoothing line
    geom_smooth(se = FALSE, method = "loess", span = 0.3, 
                color = "black", linewidth = 0.8) +
    
    # 3. Faceting
    facet_grid(Variable ~ seqnames, 
               scales = "free", 
               space = "free_x", 
               labeller = labeller(Variable = facet_labels)) +
    
    # 4. Scaling
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), 
                       breaks = pretty_breaks(n = 3)) +
    
    # 5. Colors & Legend Keys
    # FIX: Added 'labels' and 'name' to make the legend readable
    scale_color_manual(
      values = c("GC_content" = "#1f78b4", "Q_AG_rate" = "#e31a1c"),
      labels = c("GC_content" = "Intron GC Content", "Q_AG_rate" = "A -> G Transition Rate"),
      name = "Metric"
    ) +
    
    labs(
      title = "Genomic Landscape of Mutational Pressure",
      subtitle = "Raw window values (faint lines) and smoothed trends (black lines)",
      x = "Genomic Position (Mb)",
      y = "Rate / Frequency"
    ) +
    
    theme_bw() +
    theme(
      # FIX: Changed from "none" to "bottom" (or "right")
      legend.position = "bottom", 
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.spacing.x = unit(0.1, "lines"), 
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}
# ******************************************************************************
# STEP 4: Cluster genomic windows by mutational pressure ----
# ______________________________________________________________________________

make_clusters <- function(data, G = 1:5, 
                          check_convergence = T,
                          retrieve_full_model = T)
{
  #' Functions to cluster genomic windows as a funciton of mutational spectrum.
  #' The clustering exercise is based on Gaussian Mixed Models (GMM)
  #' 
  #' @param data Data frame with windows id in the first column and features for
  #' clustering (e. g. PC as summary variables)
  #' @param G Number of groups or clusters.
  #' @param check_convergence Whether to check if the GMM converged if the G selected
  #' is the ceiling defined by the user.
  #' @param retrieve_full_model Whether to retrive the full model or not
  #' 
  #' @return Cluster assignations.
  #' ___________________________________________________________________________
  
  model <- mclust::Mclust(data, G = G)
  G_empirical <- model$G
  
  if (check_convergence) 
  {
    # G_selected must be lower than G_empirical (this is my convergence criterium)
    while(max(G) == G_empirical)
    {
      # Create new cealing
      G = G + 10
      
      # Fir the GMM
      model <- mclust::Mclust(data_ts, G = G)
      G_empirical <- model$G
    }
  }
  
  if (!retrieve_full_model)
  {
    cluster_res <- model$classification
  }
  else
  {
    cluster_res <- model
  }
  
  return(cluster_res)
}

# ******************************************************************************
# STEP 5: Generate AnaCoDa dM File ----
# ______________________________________________________________________________

generate_anacoda_dM <- function(pi_A, pi_C, pi_G, pi_T, output_file,
                                splitS = TRUE) 
{
  message("Generating AnaCoDa dM (Mutation Bias) file...")
  
  # 1. Define Standard Genetic Code
  # Exclude Stop codons (*), Methionine (M), and Tryptophan (W)
  # M and W have only 1 codon, so dM is undefined/irrelevant for them.
  genetic_code <- Biostrings::GENETIC_CODE
  valid_codons <- names(genetic_code)[!genetic_code %in% c("*", "M", "W")]
  
  if(splitS)
  {
    # Assign Z to the 2-fold codon family of S (AGY)
    
    codon_idx <- names(genetic_code) %in% c("AGT", "AGC")
    genetic_code[codon_idx] <- "Z"
  }
  
  # 2. Create Data Frame for calculations
  dM_df <- data.frame(
    AA = genetic_code[valid_codons],
    Codon = valid_codons,
    stringsAsFactors = FALSE
  )
  
  # 3. Calculate Expected Frequency based on Nucleotide Pi
  # (Assuming independence of positions, which is standard for intron-derived priors)
  
  # Helper to get prob of a triplet string
  get_codon_prob <- function(codon) {
    bases <- strsplit(codon, "")[[1]]
    probs <- c(A=pi_A, C=pi_C, G=pi_G, T=pi_T)
    return(prod(probs[bases]))
  }
  
  dM_df$Expected_Freq <- sapply(dM_df$Codon, get_codon_prob)
  
  # 4. Calculate dM relative to reference codon
  # AnaCoDa typically uses the last codon alphabetically as reference, 
  # or simply requires log(freq / ref_freq).
  
  dM_df$dM <- NA
  
  for (aa in unique(dM_df$AA)) {
    # Subset for this Amino Acid
    aa_indices <- which(dM_df$AA == aa)
    sub_data <- dM_df[aa_indices, ]
    
    # Identify Reference Codon (Last alphabetically is standard convention)
    # e.g., for Alanine: GCA, GCC, GCG, GCT -> Reference is GCT
    ref_codon_row <- sub_data[order(sub_data$Codon, decreasing = FALSE), ][nrow(sub_data), ]
    ref_freq <- ref_codon_row$Expected_Freq
    
    # Calculate Log Ratio
    # dM = log( frequency_ref / frequency_current )
    dM_df$dM[aa_indices] <- log(ref_freq / dM_df$Expected_Freq[aa_indices])
  }
  
  # 5. Write to CSV
  # AnaCoDa expects columns: AA, Codon, dM
  write.csv(dM_df[, c("AA", "Codon", "dM")], file = output_file, 
            row.names = FALSE, quote = FALSE)
  
  message(paste("dM file written to:", output_file))
  return(dM_df)
}