merge_2_and_4_to_6_fold <- function(preference_df, AA_family_col)
{
  #' This function will condense together the 4 and 2 fold families from a 6-fold
  #' aminoacid family back to the one preferred codon per amino acid. It will take
  #' the aminoacid with the greater adaptiveness.
  #' 
  #' FOR COMPATIBILITY WITH OTHER PLANT STUDIES
  #' 
  #' Args:
  #' preference_df: Data frame with columns for Amino Acid, Codon_RNA, relative_adaptiveness
  #' AA_family_col: Column name indicating the root amino acid family (e.g., "Leu" for "Leu_2" and "Leu_4")
  #' 
  #' ___________________________________________________________________________
  
  condensed_preferences <- preference_df |>
    dplyr::group_by(!!sym(AA_family_col)) |>
    dplyr::arrange(!!sym(AA_family_col), 
                   eta) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(Amino_Acid = !!sym(AA_family_col), 
                  Codon_RNA, eta)
  
  return(condensed_preferences)
}