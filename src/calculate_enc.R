calculate_enc <- function(codon_counts, genetic_code, have_F6 = FALSE)
{
  #' Calculate Effective Number of Codons (ENC)
  #' 
  #' @description ENC quantifies codon bias, ranging from 20 (extreme bias, 
  #' one codon per amino acid) to 61 (no bias, all codons used equally).
  #' Based on Wright (1990) method.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param have_F6 Logical indicating if F6 values are precomputed (default FALSE)ve_F6 Wheter families of 6 codons are being considered. Default
  #' is FALSE, because the recommended approach is to split them in the F2 and F4
  #' subfamilies.
  #' 
  #' @return Data frame with Gene_name and ENC values
  #' ___________________________________________________________________________
  
  # --- 1. Prepare Metadata ---
  # Create a data.table of the genetic code and add degeneracy class
  gc_dt <- data.table(Codon = names(genetic_code), AA = genetic_code)
  
  if (!have_F6) {
    # Split 6-codon families into 4-fold and 2-fold degenerate families
    # This is the recommended approach (Wright 1990)
    
    # Define the split families for 6-codon amino acids
    # Arg: CGN (4-fold) + AGA/AGG (2-fold)
    # Leu: CTN (4-fold) + TTA/TTG (2-fold)
    # Ser: TCN (4-fold) + AGT/AGC (2-fold)
    
    # Create family assignments
    gc_dt[, Family := AA]  # Default: family = amino acid
    
    # Split Arg into Arg_4fold and Arg_2fold
    gc_dt[Codon %in% c("CGT", "CGC", "CGA", "CGG"), Family := "Arg_4fold"]
    gc_dt[Codon %in% c("AGA", "AGG"), Family := "Arg_2fold"]
    
    # Split Leu into Leu_4fold and Leu_2fold
    gc_dt[Codon %in% c("CTT", "CTC", "CTA", "CTG"), Family := "Leu_4fold"]
    gc_dt[Codon %in% c("TTA", "TTG"), Family := "Leu_2fold"]
    
    # Split Ser into Ser_4fold and Ser_2fold
    gc_dt[Codon %in% c("TCT", "TCC", "TCA", "TCG"), Family := "Ser_4fold"]
    gc_dt[Codon %in% c("AGT", "AGC"), Family := "Ser_2fold"]
    
    # Get degeneracy for each family (excluding STOP, Met, Trp)
    family_groups <- split(gc_dt[!(AA %in% c("STOP", "Met", "Trp"))]$Codon, 
                          gc_dt[!(AA %in% c("STOP", "Met", "Trp"))]$Family)
    
    family_degeneracy <- sapply(family_groups, length)
    family_to_deg_dt <- data.table(Family = names(family_degeneracy), 
                                    Degeneracy = family_degeneracy)
    
    # Join degeneracy info back to the main genetic code table
    gc_dt <- gc_dt[family_to_deg_dt, on = "Family"]
    
  } else {
    # Keep 6-codon families together (old approach)
    # Get degeneracy for AAs (excluding STOP, Met, Trp)
    aa_groups <- split(gc_dt[!(AA %in% c("STOP", "Met", "Trp"))]$Codon, 
                       gc_dt[!(AA %in% c("STOP", "Met", "Trp"))]$AA)
    
    aa_degeneracy <- sapply(aa_groups, length)
    aa_to_deg_dt <- data.table(AA = names(aa_degeneracy), Degeneracy = aa_degeneracy)
    
    # Join degeneracy info back to the main genetic code table
    gc_dt <- gc_dt[aa_to_deg_dt, on = "AA"]
  }
  
  # --- 2. Melt Codon Counts to Long Format ---
  # We only need to melt codons that are part of the calculation
  codon_cols <- gc_dt$Codon
  melted_counts <- melt(codon_counts, 
                        id.vars = "Gene_name", 
                        measure.vars = intersect(codon_cols, names(codon_counts)), 
                        variable.name = "Codon", 
                        value.name = "Count")
  
  # --- 3. Join Metadata with Counts ---
  melted_counts <- gc_dt[melted_counts, on = "Codon"]
  
  # --- 4. Calculate F_aa (Homozygosity) for all Gene/AA (or Family) pairs ---
  
  if (!have_F6) {
    # When splitting families, group by Family instead of AA
    # a. Get total count for each Family in each gene
    melted_counts[, Total := sum(Count), by = .(Gene_name, Family)]
    
    # b. Calculate sum of squared proportions (p_squared_sum)
    #    Only for Families that are present (Total > 0)
    aa_f_values <- melted_counts[Total > 0, 
                                 .(p_squared_sum = sum((Count / Total)^2)), 
                                 by = .(Gene_name, Family, Degeneracy, Total)]
    
    # c. Calculate F_family for each Family
    #    Only where Total > 1 (to avoid 0/0 division)
    aa_f_values[Total > 1, F_aa := (Total * p_squared_sum - 1) / (Total - 1)]
    
  } else {
    # When keeping 6-codon families together, group by AA
    # a. Get total count for each AA in each gene
    melted_counts[, Total := sum(Count), by = .(Gene_name, AA)]
    
    # b. Calculate sum of squared proportions (p_squared_sum)
    #    Only for AAs that are present (Total > 0)
    aa_f_values <- melted_counts[Total > 0, 
                                 .(p_squared_sum = sum((Count / Total)^2)), 
                                 by = .(Gene_name, AA, Degeneracy, Total)]
    
    # c. Calculate F_aa for each AA
    #    Only where Total > 1 (to avoid 0/0 division)
    aa_f_values[Total > 1, F_aa := (Total * p_squared_sum - 1) / (Total - 1)]
  }
  
  # d. Filter out NA/Inf values (where Total was 1)
  aa_f_values <- aa_f_values[!is.na(F_aa)]
  
  # --- 5. Calculate Average F (F_bar_k) for each Degeneracy Class ---
  f_bar_k <- aa_f_values[, .(F_bar = mean(F_aa)), by = .(Gene_name, Degeneracy)]
  
  # --- 6. Dcast to Wide Format ---
  # Create a table with one row per gene and columns: Gene_name, F2, F3, F4, F6
  f_bar_wide <- dcast(f_bar_k, Gene_name ~ paste0("F", Degeneracy), value.var = "F_bar")
  
  # --- 7. Join with All Genes ---
  # This ensures genes with no valid F values are included
  all_genes <- data.table(Gene_name = codon_counts$Gene_name)
  f_bar_wide <- f_bar_wide[all_genes, on = "Gene_name"]
  
  # --- 8. Calculate Final ENC (Corrected) ---
  # This section fixes the bug that caused ENC > 61.
  # We cap the contribution of each class at its theoretical maximum ("no bias" value).
  
  # ENC Formula (Wright 1990):
  # When splitting 6-codon families into 4-fold and 2-fold:
  # - 2-fold families: Asn, Asp, Cys, Gln, Glu, His, Lys, Phe, Tyr + Arg(2), Leu(2), Ser(2) = 12 families
  # - 3-fold family: Ile = 1 family
  # - 4-fold families: Ala, Gly, Pro, Thr, Val + Arg(4), Leu(4), Ser(4) = 8 families
  # - 1-fold: Met, Trp = 2 amino acids
  #
  # Max contributions (no bias):
  # - Met + Trp: 2
  # - 2-fold (12 families): 12 × 2 = 24
  # - 3-fold (1 family): 1 × 3 = 3
  # - 4-fold (8 families): 8 × 4 = 32
  # Total: 2 + 24 + 3 + 32 = 61 ✓
  #
  # When using 6-codon families (not split):
  # - 2-fold: 9 families → 18
  # - 3-fold: 1 family → 3
  # - 4-fold: 5 families → 20
  # - 6-fold: 3 families → 18
  # Total: 2 + 18 + 3 + 20 + 18 = 61 ✓
  
  # Start with 2 for Met and Trp
  f_bar_wide[, ENC := 2.0] 
  
  if(have_F6) {
    # Using 6-codon families (old approach)
    # --- Class 2: 9 families (Max contribution = 18) ---
    f_bar_wide[, ENC_2 := 9.0]  # Default (extreme bias)
    f_bar_wide[!is.na(F2) & F2 > 0, ENC_2 := 9.0 / F2]
    f_bar_wide[, ENC_2 := pmin(ENC_2, 18.0)]  # Cap at 18
    
    # --- Class 3: 1 family (Max contribution = 3) ---
    f_bar_wide[, ENC_3 := 1.0]  # Default
    f_bar_wide[!is.na(F3) & F3 > 0, ENC_3 := 1.0 / F3]
    f_bar_wide[, ENC_3 := pmin(ENC_3, 3.0)]  # Cap at 3
    
    # --- Class 4: 5 families (Max contribution = 20) ---
    f_bar_wide[, ENC_4 := 5.0]  # Default
    f_bar_wide[!is.na(F4) & F4 > 0, ENC_4 := 5.0 / F4]
    f_bar_wide[, ENC_4 := pmin(ENC_4, 20.0)]  # Cap at 20
    
    # --- Class 6: 3 families (Max contribution = 18) ---
    f_bar_wide[, ENC_6 := 3.0]  # Default
    f_bar_wide[!is.na(F6) & F6 > 0, ENC_6 := 3.0 / F6]
    f_bar_wide[, ENC_6 := pmin(ENC_6, 18.0)]  # Cap at 18
    
    # Sum all contributions
    f_bar_wide[, ENC := ENC + ENC_2 + ENC_3 + ENC_4 + ENC_6]
    
  } else {
    # Splitting 6-codon families into 4-fold and 2-fold (recommended approach)
    # --- Class 2: 12 families (Max contribution = 24) ---
    # 9 original 2-fold + 3 new 2-fold from splitting (Arg, Leu, Ser)
    f_bar_wide[, ENC_2 := 12.0]  # Default (extreme bias)
    f_bar_wide[!is.na(F2) & F2 > 0, ENC_2 := 12.0 / F2]
    f_bar_wide[, ENC_2 := pmin(ENC_2, 24.0)]  # Cap at 24
    
    # --- Class 3: 1 family (Max contribution = 3) ---
    f_bar_wide[, ENC_3 := 1.0]  # Default
    f_bar_wide[!is.na(F3) & F3 > 0, ENC_3 := 1.0 / F3]
    f_bar_wide[, ENC_3 := pmin(ENC_3, 3.0)]  # Cap at 3
    
    # --- Class 4: 8 families (Max contribution = 32) ---
    # 5 original 4-fold + 3 new 4-fold from splitting (Arg, Leu, Ser)
    f_bar_wide[, ENC_4 := 8.0]  # Default
    f_bar_wide[!is.na(F4) & F4 > 0, ENC_4 := 8.0 / F4]
    f_bar_wide[, ENC_4 := pmin(ENC_4, 32.0)]  # Cap at 32
    
    # Sum all contributions
    f_bar_wide[, ENC := ENC + ENC_2 + ENC_3 + ENC_4]
  }   
  
  # --- 9. Return Final Result ---
  # Select only the columns we need
  result <- f_bar_wide[, .(Gene_name, ENC)]
  
  return(result)
}