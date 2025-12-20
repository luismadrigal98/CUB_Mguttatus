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

get_inergenic_sequences <- function(fasta_file, ann_file,
                                 trim_bp = 1000,
                                 width = 10000,
                                 organism = NA)
{
  #' @title Extract Trimmed Intergenic Sequences
  #'
  #' @description
  #' Loads an annotation file (GFF3) and a genome FASTA file, filters for primary 
  #' transcripts on main chromosomes, get width / 2 bp from each side of the 
  #' gene after removing the trim_bp proximal to it (main regulatory core).
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
  
  message("Step 1: Loading Annotation...")
  
  # --- 1. Load Annotation and Create TxDb ---
  txdb <- txdbmaker::makeTxDbFromGFF(file = ann_file, 
                                     format = "gff3",
                                     organism = organism)
  
  all_genes <- GenomicFeatures::genes(txdb)
  
  # --- 2. Filter for Main Chromosomes ---
  main_chr_pattern <- "^Chr_[0-9]{1,2}$"
  genes_on_main <- all_genes[grepl(main_chr_pattern, seqnames(all_genes))]
  
  message(paste("Original genes:", length(all_genes)))
  message(paste("Genes on main chromosomes:", length(genes_on_main)))
  
  # --- 3. Define Upstream and Downstream Regions ---
  message("Step 2: Calculating Flanking Coordinates...")
  
  # A. Upstream (Promoter-ish side)
  # flank(start=TRUE) gets the 5' side.
  upstream_raw <- GenomicRanges::flank(genes_on_main, width = width + trim_bp, start = TRUE)
  # Remove the 'trim_bp' part closest to the gene.
  upstream_ranges <- GenomicRanges::resize(upstream_raw, width = width, fix = "start")
  
  # B. Downstream (Terminator-ish side)
  # flank(start=FALSE) gets the 3' side.
  downstream_raw <- GenomicRanges::flank(genes_on_main, width = width + trim_bp, start = FALSE)
  # Remove the 'trim_bp' part closest to the gene.
  downstream_ranges <- GenomicRanges::resize(downstream_raw, width = width, fix = "end")
  
  # --- 4. Load Genome and Harmonize Seqlevels ---
  message("Step 3: Loading Genome Sequence...")
  
  if (!file.exists(fasta_file)) stop("FASTA file not found.")
  dna <- Biostrings::readDNAStringSet(filepath = fasta_file)
  
  # Simplify FASTA names
  names(dna) <- sub("^(\\S+)\\s.*", "\\1", names(dna))
  
  # Harmonize
  fasta_chroms <- names(dna)[grepl(main_chr_pattern, names(dna))]
  upstream_ranges <- GenomeInfoDb::keepSeqlevels(upstream_ranges, fasta_chroms, pruning.mode = "coarse")
  downstream_ranges <- GenomeInfoDb::keepSeqlevels(downstream_ranges, fasta_chroms, pruning.mode = "coarse")
  
  # Ensure ranges fit within chromosome limits
  seqlengths(upstream_ranges) <- width(dna)[match(seqlevels(upstream_ranges), names(dna))]
  seqlengths(downstream_ranges) <- width(dna)[match(seqlevels(downstream_ranges), names(dna))]
  
  upstream_ranges <- GenomicRanges::trim(upstream_ranges)
  downstream_ranges <- GenomicRanges::trim(downstream_ranges)
  
  # Remove truncated ranges
  upstream_ranges <- upstream_ranges[width(upstream_ranges) == width]
  downstream_ranges <- downstream_ranges[width(downstream_ranges) == width]
  
  # --- 4b. Handle Overlapping Regions ---
  # When genes are close together, their flanking regions can overlap.
  # We use `reduce()` to merge overlapping ranges, ensuring each genomic
  # position is counted only ONCE in downstream nucleotide composition analysis.
  
  message("Step 3b: Resolving overlapping regions...")
  
  n_upstream_before <- length(upstream_ranges)
  n_downstream_before <- length(downstream_ranges)
  
  # Combine upstream and downstream, then reduce overlaps
  # This creates non-overlapping intervals that cover all intergenic space
  all_intergenic <- c(upstream_ranges, downstream_ranges)
  all_intergenic_reduced <- GenomicRanges::reduce(all_intergenic)
  
  # Also remove any regions that overlap with genes themselves
  # (in case flanking regions from adjacent genes intrude)
  all_intergenic_reduced <- GenomicRanges::setdiff(all_intergenic_reduced, genes_on_main)
  
  message(sprintf("  Upstream regions: %d -> reduced", n_upstream_before))
  message(sprintf("  Downstream regions: %d -> reduced", n_downstream_before))
  message(sprintf("  Combined non-overlapping intergenic regions: %d", length(all_intergenic_reduced)))
  message(sprintf("  Total bp coverage: %s", format(sum(width(all_intergenic_reduced)), big.mark = ",")))
  
  # --- 5. Extract Sequences (Manual Method) ---
  message("Step 4: Extracting Sequences...")
  
  # Helper function to extract and reverse complement based on strand
  extract_and_orient <- function(granges_obj, genome_dna) {
    
    # 1. Map chromosome names to integer indices for faster extraction
    chr_names <- as.character(GenomicRanges::seqnames(granges_obj))
    starts <- GenomicRanges::start(granges_obj)
    ends <- GenomicRanges::end(granges_obj)
    
    # 2. Vectorized extraction using lapply over the unique chromosomes
    # (This is faster than looping over every single gene)
    unique_chrs <- unique(chr_names)
    
    # Create a list to store results
    extracted_list <- character(length(granges_obj))
    
    # Extract sequences per chromosome to minimize subsets
    for (chr in unique_chrs) {
      # Identify indices for this chromosome
      idx <- which(chr_names == chr)
      
      # Extract from the specific chromosome DNAString
      # Views is efficient for multiple ranges on one sequence
      if (length(idx) > 0) {
        chr_seq <- genome_dna[[chr]]
        v <- Biostrings::Views(chr_seq, start = starts[idx], end = ends[idx])
        extracted_list[idx] <- as.character(v)
      }
    }
    
    # Convert to DNAStringSet
    dss <- Biostrings::DNAStringSet(extracted_list)
    names(dss) <- names(granges_obj)
    
    # # 3. Handle Strand (Reverse Complement Negative Strand)
    # # For reduced intergenic regions, strand is "*" (unstranded)
    # # Only apply RC if explicitly negative strand
    # neg_strand <- as.character(GenomicRanges::strand(granges_obj)) == "-"
    # if (any(neg_strand)) {
    #   dss[neg_strand] <- Biostrings::reverseComplement(dss[neg_strand])
    # } >> Not necessary for intergeinc regions
    
    return(dss)
  }
  
  # Extract sequences from the REDUCED (non-overlapping) intergenic regions
  # This ensures each genomic position is counted only once
  intergenic_seqs <- extract_and_orient(all_intergenic_reduced, dna)
  
  message("Extraction complete.")
  message(sprintf("  Extracted %d non-overlapping intergenic sequences", length(intergenic_seqs)))
  
  # Also keep the original upstream/downstream for backward compatibility if needed
  upstream_seqs <- extract_and_orient(upstream_ranges, dna)
  downstream_seqs <- extract_and_orient(downstream_ranges, dna)
  
  return(list(
    # PRIMARY: Use these for unbiased nucleotide composition
    intergenic_seqs = intergenic_seqs,
    intergenic_ranges = all_intergenic_reduced,
    # SECONDARY: Original (potentially overlapping) for other analyses
    upstream_seqs = upstream_seqs,
    downstream_seqs = downstream_seqs,
    upstream_ranges = upstream_ranges,
    downstream_ranges = downstream_ranges,
    genome_seqinfo = GenomeInfoDb::seqinfo(dna)
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
  #' @description Internal function to extract sequences within a window 
  #' and calculate A, C, G, T frequencies.
  #' @param window_idx Numeric. Index of the current window.
  #' @param all_windows GRanges object defining all genomic windows.
  #' @param all_seqs DNAStringSet containing the sequences (upstream or downstream).
  #' @param hit_list Hits object from findOverlaps mapping windows to features.
  #' @param return_Ns Whether to return also the count of Ns.
  #' @return A named vector of window metadata, total base pairs, and base frequencies.
  #' ___________________________________________________________________________
  
  # Identify which features (upstream/downstream regions) belong to this specific window index
  feature_indices <- S4Vectors::subjectHits(hit_list)[S4Vectors::queryHits(hit_list) == window_idx]
  
  if (length(feature_indices) == 0) {
    return(NULL) # Skip empty windows
  }
  
  # Extract sequences for this window
  local_seqs <- all_seqs[feature_indices]
  
  # Count bases (A, C, G, T)
  # baseOnly = TRUE returns 5 columns: A, C, G, T, and "other" (which includes N)
  counts <- Biostrings::alphabetFrequency(local_seqs, baseOnly = TRUE, as.prob = FALSE)
  
  # Sum counts across all features in this window
  total_counts <- colSums(counts)
  
  # Calculate Total Valid Base Pairs (Strictly A + C + G + T)
  total_bp <- sum(total_counts[1:4])
  
  # CRITICAL CHECK: If the window was mostly 'N's, total_bp might be tiny or zero.
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
  
  # Return data
  return(c(window_data, total_bp = total_bp, freqs, N_count = N_count))
}

get_base_composition_per_windows <- function(input_data, 
                                             window_size = 100000)
{
  #' @title Aggregate Base Composition per Window (Generic)
  #'
  #' @description Tiles the genome into windows and calculates A, C, G, T frequencies
  #' for any valid genomic regions present in the input object (introns, upstream, or downstream).
  #'
  #' @param input_data A list output from either `get_intron_sequences()` or `get_intergenic_sequences()`.
  #'   Must contain 'genome_seqinfo' and paired ranges/sequences (e.g., 'trimmed_introns'/'intron_seqs' 
  #'   OR 'upstream_ranges'/'upstream_seqs', etc.).
  #' @param window_size Numeric. The size of the genomic windows (default: 100000 bp).
  #' @return A combined data frame with window coordinates, base frequencies, and a 'region_type' column.
  #' @export
  #' ___________________________________________________________________________
  
  # --- 1. Validation and Setup ---
  
  if (!"genome_seqinfo" %in% names(input_data)) {
    stop("input_data must contain a 'genome_seqinfo' object to define chromosome lengths.")
  }
  
  genome_seqinfo <- input_data$genome_seqinfo
  
  message("Step 1: Tiling Genome...")
  
  # Define genomic windows
  seq_lengths <- GenomeInfoDb::seqlengths(genome_seqinfo)
  windows <- GenomicRanges::tileGenome(seq_lengths, 
                                       tilewidth = window_size, 
                                       cut.last.tile.in.chrom = TRUE)
  
  # --- 2. Define Processing Helper ---
  
  process_region <- function(ranges_obj, seqs_obj, region_label) {
    
    message(paste("  Processing", region_label, "regions..."))
    
    # Map features to windows
    overlaps <- GenomicRanges::findOverlaps(windows, ranges_obj)
    
    # Apply calculation across windows
    results_list <- lapply(1:length(windows), calculate_window_metrics, 
                           all_windows = windows, 
                           all_seqs = seqs_obj, 
                           hit_list = overlaps)
    
    # Filter empty results
    non_null_results <- results_list[!sapply(results_list, is.null)]
    
    if (length(non_null_results) == 0) {
      warning(paste("No overlapping windows found for", region_label))
      return(NULL)
    }
    
    # Combine and Format
    df_results <- do.call(rbind, non_null_results)
    df_results <- as.data.frame(df_results)
    
    # Numeric conversion for relevant columns
    numeric_cols <- c("start", "end", "window_idx", "pi_A", "pi_C", "pi_G", "pi_T",
                      "total_bp", "N_count")
    for(col in numeric_cols) {
      if(col %in% colnames(df_results)) {
        df_results[[col]] <- as.numeric(df_results[[col]])
      }
    }
    
    # Add identifier column
    df_results$region_type <- region_label
    
    return(df_results)
  }
  
  # --- 3. Dynamic Execution based on Input Content ---
  
  results_list <- list()
  
  # Case A: Introns (from get_intron_sequences)
  if (!is.null(input_data$trimmed_introns) && !is.null(input_data$intron_seqs)) {
    results_list[["intron"]] <- process_region(input_data$trimmed_introns, 
                                               input_data$intron_seqs, 
                                               "intron")
  }
  
  # Case B: PREFERRED - Reduced/merged intergenic regions (non-overlapping)
  # This should be used for unbiased nucleotide composition estimation
  if (!is.null(input_data$intergenic_ranges) && !is.null(input_data$intergenic_seqs)) {
    results_list[["intergenic"]] <- process_region(input_data$intergenic_ranges, 
                                                   input_data$intergenic_seqs, 
                                                   "intergenic")
    message("  Using REDUCED (non-overlapping) intergenic regions for unbiased composition")
  } else {
    # Fallback: Use separate upstream/downstream (may have overlaps - backward compatibility)
    # Case B-alt: Upstream Intergenic (from get_intergenic_sequences)
    if (!is.null(input_data$upstream_ranges) && !is.null(input_data$upstream_seqs)) {
      results_list[["upstream"]] <- process_region(input_data$upstream_ranges, 
                                                   input_data$upstream_seqs, 
                                                   "upstream")
      warning("Using separate upstream regions - may contain overlaps with downstream!")
    }
    
    # Case C-alt: Downstream Intergenic (from get_intergenic_sequences)
    if (!is.null(input_data$downstream_ranges) && !is.null(input_data$downstream_seqs)) {
      results_list[["downstream"]] <- process_region(input_data$downstream_ranges, 
                                                     input_data$downstream_seqs, 
                                                     "downstream")
      warning("Using separate downstream regions - may contain overlaps with upstream!")
    }
  }
  
  # --- 4. Final Combination ---
  
  if (length(results_list) == 0) {
    stop("Input data does not contain recognizable sequence data (introns, intergenic, upstream, or downstream).")
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL # Clean up row names
  
  return(final_df)
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
    dplyr::mutate(
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

# ******************************************************************************
# STEP 6: Complete dM Estimation Pipeline (Wrapper) ----
# ______________________________________________________________________________

estimate_dM_from_neutral_regions <- function(fasta_file, 
                                              ann_file,
                                              output_dir = "./data",
                                              output_prefix = "Mguttatus",
                                              source = c("introns", "intergenic", "both"),
                                              window_size = 100000,
                                              min_bp = 1000,
                                              max_N_freq = 0.25,
                                              organism = "Mimulus guttatus",
                                              # Intron-specific parameters
                                              intron_trim_bp = 30,
                                              intron_min_width = 86,
                                              # Intergenic-specific parameters
                                              intergenic_trim_bp = 1000,
                                              intergenic_width = 10000,
                                              return_intermediates = FALSE) 
{
  #' @title Complete dM Estimation Pipeline from Neutral Genomic Regions
  #'

  #' @description

  #' Wrapper function that extracts neutral sequences (introns and/or intergenic regions),
  #' calculates nucleotide composition per genomic window, filters for quality,
  #' and generates AnaCoDa-compatible dM (mutation bias) files.
  #'
  #' This consolidates the workflow that was previously duplicated in main.R for

  #' introns vs intergenic regions.
  #'
  #' @param fasta_file Character. Path to the genome reference FASTA file (preferably hard-masked).
  #' @param ann_file Character. Path to the gene annotation GFF3 file.
  #' @param output_dir Character. Directory to save output dM files. Default: "./data"
  #' @param output_prefix Character. Prefix for output filenames. Default: "Mguttatus"
  #' @param source Character. Which neutral regions to use: "introns", "intergenic", or "both".
  #' @param window_size Numeric. Size of genomic windows for composition analysis. Default: 100000
  #' @param min_bp Numeric. Minimum base pairs per window to retain. Default: 1000
  #' @param max_N_freq Numeric. Maximum frequency of N's allowed per window. Default: 0.25
  #' @param organism Character. Organism name for TxDb metadata. Default: "Mimulus guttatus"
  #' @param intron_trim_bp Numeric. BP to trim from splice sites in introns. Default: 30
  #' @param intron_min_width Numeric. Minimum intron width after trimming. Default: 86
  #' @param intergenic_trim_bp Numeric. BP to trim near genes in intergenic regions. Default: 1000
  #' @param intergenic_width Numeric. Width of intergenic flanking regions to extract. Default: 10000
  #' @param return_intermediates Logical. Return intermediate data objects? Default: FALSE
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \code{dM_introns}: dM data frame from introns (if source includes "introns")
  #'   \item \code{dM_intergenic}: dM data frame from intergenic (if source includes "intergenic")
  #'   \item \code{global_stats_introns}: Global nucleotide frequencies from introns
  #'   \item \code{global_stats_intergenic}: Global nucleotide frequencies from intergenic
  #'   \item \code{output_files}: Paths to generated dM CSV files
  #'   \item \code{intermediates}: (Optional) Intermediate data objects if requested
  #' }
  #'
  #' @examples
  #' \dontrun{
  #' # Generate dM from both introns and intergenic regions
  #' dM_results <- estimate_dM_from_neutral_regions(
  #'   fasta_file = "./data/genome.hardmasked.fa",
  #'   ann_file = "./data/genes.gff3",
  #'   source = "both",
  #'   output_prefix = "Mguttatus"
  #' )
  #' 
  #' # Access the dM files
  #' print(dM_results$output_files)
  #' }
  #'
  #' @export
  #' ___________________________________________________________________________
  
  source <- match.arg(source)
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list(
    dM_introns = NULL,
    dM_intergenic = NULL,
    global_stats_introns = NULL,
    global_stats_intergenic = NULL,
    output_files = character(),
    intermediates = list()
  )
  
  # ---------------------------------------------------------------------------
  # Helper: Process a single source type
  # ---------------------------------------------------------------------------
  process_source <- function(source_type) {
    
    cat(sprintf("\n%s\n", paste(rep("=", 70), collapse = "")))
    cat(sprintf("Processing: %s\n", toupper(source_type)))
    cat(sprintf("%s\n\n", paste(rep("=", 70), collapse = "")))
    
    # 1. Extract sequences based on source type
    if (source_type == "introns") {
      cat("Step 1: Extracting intronic sequences...\n")
      seq_data <- get_intron_sequences(
        fasta_file = fasta_file,
        ann_file = ann_file,
        trim_bp = intron_trim_bp,
        min_width = intron_min_width,
        organism = organism
      )
    } else {
      cat("Step 1: Extracting intergenic sequences...\n")
      seq_data <- get_inergenic_sequences(
        fasta_file = fasta_file,
        ann_file = ann_file,
        trim_bp = intergenic_trim_bp,
        width = intergenic_width,
        organism = organism
      )
    }
    
    # 2. Calculate nucleotide composition per window
    cat("\nStep 2: Calculating nucleotide composition per window...\n")
    nuc_composition <- get_base_composition_per_windows(
      input_data = seq_data,
      window_size = window_size
    )
    
    cat(sprintf("  Total windows: %d\n", nrow(nuc_composition)))
    
    # 3. Filter windows
    cat("\nStep 3: Filtering windows...\n")
    
    # Filter by minimum base pairs
    nuc_filtered <- nuc_composition |>
      dplyr::filter(total_bp >= min_bp)
    
    cat(sprintf("  After min_bp filter (>= %d bp): %d windows\n", 
                min_bp, nrow(nuc_filtered)))
    
    # Calculate and filter by N frequency
    nuc_filtered <- nuc_filtered |>
      dplyr::mutate(
        N_freq = N_count / (total_bp + N_count)
      ) |>
      dplyr::filter(N_freq < max_N_freq)
    
    cat(sprintf("  After N frequency filter (< %.2f): %d windows\n", 
                max_N_freq, nrow(nuc_filtered)))
    
    if (nrow(nuc_filtered) == 0) {
      warning(sprintf("No windows passed filtering for %s!", source_type))
      return(list(dM = NULL, global_stats = NULL, output_file = NULL, 
                  intermediates = list(seq_data = seq_data, 
                                       nuc_composition = nuc_composition)))
    }
    
    # 4. Calculate global weighted nucleotide frequencies
    cat("\nStep 4: Calculating global nucleotide frequencies (weighted by bp)...\n")
    
    global_stats <- nuc_filtered |>
      dplyr::summarize(
        total_genome_bp = sum(total_bp),
        n_windows = dplyr::n(),
        avg_pi_A = sum(pi_A * total_bp) / sum(total_bp),
        avg_pi_C = sum(pi_C * total_bp) / sum(total_bp),
        avg_pi_G = sum(pi_G * total_bp) / sum(total_bp),
        avg_pi_T = sum(pi_T * total_bp) / sum(total_bp)
      )
    
    cat(sprintf("  Total BP analyzed: %s\n", format(global_stats$total_genome_bp, big.mark = ",")))
    cat(sprintf("  Windows used: %d\n", global_stats$n_windows))
    cat(sprintf("  π(A) = %.4f, π(C) = %.4f, π(G) = %.4f, π(T) = %.4f\n",
                global_stats$avg_pi_A, global_stats$avg_pi_C, 
                global_stats$avg_pi_G, global_stats$avg_pi_T))
    cat(sprintf("  GC content = %.2f%%\n", 
                100 * (global_stats$avg_pi_G + global_stats$avg_pi_C)))
    
    # 5. Generate dM file
    cat("\nStep 5: Generating AnaCoDa dM file...\n")
    
    output_file <- file.path(output_dir, 
                              sprintf("%s_%s_derived_dM.csv", output_prefix, source_type))
    
    dM_data <- generate_anacoda_dM(
      pi_A = global_stats$avg_pi_A,
      pi_C = global_stats$avg_pi_C,
      pi_G = global_stats$avg_pi_G,
      pi_T = global_stats$avg_pi_T,
      output_file = output_file
    )
    
    cat(sprintf("\n✓ %s processing complete!\n", toupper(source_type)))
    cat(sprintf("  Output: %s\n", output_file))
    
    return(list(
      dM = dM_data,
      global_stats = global_stats,
      output_file = output_file,
      intermediates = list(
        seq_data = seq_data,
        nuc_composition = nuc_composition,
        nuc_filtered = nuc_filtered
      )
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Execute based on source selection
  # ---------------------------------------------------------------------------
  
  if (source %in% c("introns", "both")) {
    intron_results <- process_source("introns")
    results$dM_introns <- intron_results$dM
    results$global_stats_introns <- intron_results$global_stats
    results$output_files <- c(results$output_files, introns = intron_results$output_file)
    if (return_intermediates) {
      results$intermediates$introns <- intron_results$intermediates
    }
  }
  
  if (source %in% c("intergenic", "both")) {
    intergenic_results <- process_source("intergenic")
    results$dM_intergenic <- intergenic_results$dM
    results$global_stats_intergenic <- intergenic_results$global_stats
    results$output_files <- c(results$output_files, intergenic = intergenic_results$output_file)
    if (return_intermediates) {
      results$intermediates$intergenic <- intergenic_results$intermediates
    }
  }
  
  # ---------------------------------------------------------------------------
  # Summary comparison (if both sources were processed)
  # ---------------------------------------------------------------------------
  
  if (source == "both" && 
      !is.null(results$global_stats_introns) && 
      !is.null(results$global_stats_intergenic)) {
    
    cat("\n")
    cat(sprintf("%s\n", paste(rep("=", 70), collapse = "")))
    cat("COMPARISON: Introns vs Intergenic\n")
    cat(sprintf("%s\n", paste(rep("=", 70), collapse = "")))
    
    cat("\nNucleotide Frequencies:\n")
    cat(sprintf("  %-12s  %8s  %8s  %8s\n", "", "Introns", "Intergenic", "Diff"))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %+8.4f\n", "π(A)", 
                results$global_stats_introns$avg_pi_A,
                results$global_stats_intergenic$avg_pi_A,
                results$global_stats_intergenic$avg_pi_A - results$global_stats_introns$avg_pi_A))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %+8.4f\n", "π(C)", 
                results$global_stats_introns$avg_pi_C,
                results$global_stats_intergenic$avg_pi_C,
                results$global_stats_intergenic$avg_pi_C - results$global_stats_introns$avg_pi_C))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %+8.4f\n", "π(G)", 
                results$global_stats_introns$avg_pi_G,
                results$global_stats_intergenic$avg_pi_G,
                results$global_stats_intergenic$avg_pi_G - results$global_stats_introns$avg_pi_G))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %+8.4f\n", "π(T)", 
                results$global_stats_introns$avg_pi_T,
                results$global_stats_intergenic$avg_pi_T,
                results$global_stats_intergenic$avg_pi_T - results$global_stats_introns$avg_pi_T))
    
    gc_introns <- results$global_stats_introns$avg_pi_G + results$global_stats_introns$avg_pi_C
    gc_intergenic <- results$global_stats_intergenic$avg_pi_G + results$global_stats_intergenic$avg_pi_C
    cat(sprintf("\n  %-12s  %7.2f%%  %7.2f%%  %+7.2f%%\n", "GC Content", 
                100 * gc_introns, 100 * gc_intergenic, 100 * (gc_intergenic - gc_introns)))
    
    cat("\n")
  }
  
  cat("\n✓ dM estimation complete!\n")
  cat("Output files:\n")
  for (nm in names(results$output_files)) {
    cat(sprintf("  - %s: %s\n", nm, results$output_files[nm]))
  }
  
  return(results)
}