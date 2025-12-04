"""
Script to extract high-confidence ribosomal and elongation factor genes from annotation files for CAI calculation.
The script identifies genes based on:
  1. Keywords in their annotations (ribosomal proteins, elongation factors)
  2. Empirical expression data (must be in top X% expressed genes)

This dual-filtering approach ensures the reference set contains genes that are both
functionally appropriate AND empirically highly expressed.

@author: Luis Javier Madrigal-Roca

@date: 2024-12-4

"""

import pandas as pd
import numpy as np
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

# Expression percentile threshold (genes must be in top X% by expression)
EXPRESSION_PERCENTILE = 5  # Top 5% - can adjust to 5, 15, 20, etc.

# Whether to use max expression across tissues or require high in all
USE_MAX_EXPRESSION = True  # True = high in at least one tissue

# Relaxed exclusion mode - set to True to include "putative" and "family" genes
RELAXED_MODE = True

# ============================================================================
# FUNCTIONS
# ============================================================================

def load_expression_data(expression_file):
    """
    Load expression data and identify highly expressed genes.
    Returns a set of gene IDs that are in the top EXPRESSION_PERCENTILE.
    """
    try:
        expr_df = pd.read_csv(expression_file)
        print(f"Loaded expression data: {len(expr_df)} genes")
        
        # Calculate max or mean expression across tissues
        expr_cols = [col for col in expr_df.columns if col.startswith('Exp_')]
        
        if USE_MAX_EXPRESSION:
            expr_df['max_expr'] = expr_df[expr_cols].max(axis=1)
            expr_metric = 'max_expr'
        else:
            expr_df['mean_expr'] = expr_df[expr_cols].mean(axis=1)
            expr_metric = 'mean_expr'
        
        # Calculate percentile threshold
        threshold = np.percentile(expr_df[expr_metric], 100 - EXPRESSION_PERCENTILE)
        print(f"Expression threshold (top {EXPRESSION_PERCENTILE}%): {threshold:.2f} CPM")
        
        # Get highly expressed genes
        high_expr_genes = set(expr_df[expr_df[expr_metric] >= threshold]['GeneID'])
        print(f"Genes in top {EXPRESSION_PERCENTILE}%: {len(high_expr_genes)}")
        
        return high_expr_genes, expr_df
        
    except FileNotFoundError:
        print(f"Error: Expression file {expression_file} not found.")
        return set(), None

def extract_reference_genes(file_path, gene_type, high_expr_genes=None):
    """
    Parses the annotation file to extract gene IDs based on high-expression keywords.
    Optionally filters by empirical expression data.
    """
    
    # Keywords for ribosomal proteins and elongation factors
    high_expression_keywords = [
        # Ribosomal Keywords (structural proteins)
        "ribosomal protein l", 
        "ribosomal protein s",
        "ribosomal protein rps",
        "ribosomal protein rpl",
        "60s ribosomal",
        "40s ribosomal",
        "50s ribosomal",
        "30s ribosomal",
        "cytosolic ribosomal",
        
        # Elongation Keywords (core translation machinery)
        "elongation factor ef",
        "elongation factor tu",
        "elongation factor 1",
        "elongation factor g",
        "elongation factor ts",
        "translation elongation factor",
        "eef1a", 
        "ef-tu",
        "ef-g",
        "ef-ts",
    ]

    # Exclusion keywords - always exclude organellar genes
    exclusion_keywords_strict = [
        "mitochondrial", 
        "chloroplast", 
        "biogenesis", 
        "processing", 
        "regulator", 
        "recycling",
    ]
    
    # Additional exclusions for strict mode
    exclusion_keywords_relaxed = [
        "like",       # "like" proteins may be low expression
        "family",     # "family protein" is often generic
        "putative",   # Less confident annotations
    ]
    
    # Choose exclusion list based on mode
    if RELAXED_MODE:
        exclusion_keywords = exclusion_keywords_strict
    else:
        exclusion_keywords = exclusion_keywords_strict + exclusion_keywords_relaxed
    
    unique_genes = set()
    annotation_matched = set()  # Genes matching annotation criteria

    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip() or "PAC:" not in line:
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                gene_id = parts[1]
                base_gene_id = gene_id.split('.')[0] + "." + gene_id.split('.')[1] 

                full_text = line.lower()
                
                # Check annotation keywords
                if any(k in full_text for k in high_expression_keywords):
                    if not any(ex.lower() in full_text for ex in exclusion_keywords):
                        annotation_matched.add(base_gene_id)
                        
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return set()

    print(f"[{gene_type}] Annotation matches: {len(annotation_matched)}")
    
    # Filter by expression if provided
    if high_expr_genes is not None:
        unique_genes = annotation_matched & high_expr_genes
        print(f"[{gene_type}] After expression filter (top {EXPRESSION_PERCENTILE}%): {len(unique_genes)}")
    else:
        unique_genes = annotation_matched
        
    return unique_genes, annotation_matched

# ============================================================================
# MAIN EXECUTION
# ============================================================================

print("=" * 60)
print("CAI Reference Gene Selection")
print(f"Mode: {'RELAXED' if RELAXED_MODE else 'STRICT'}")
print(f"Expression filter: Top {EXPRESSION_PERCENTILE}%")
print("=" * 60)

# 1. Load expression data
high_expr_genes, expr_df = load_expression_data('data/observed_expression_multitissue.csv')

print("-" * 60)

# 2. Process Ribosomal Genes
ribo_genes, ribo_all = extract_reference_genes(
    'data/candidate_ribosome_associated_genes.txt', 
    'Ribosomal',
    high_expr_genes
)

# 3. Process Elongation Genes
elong_genes, elong_all = extract_reference_genes(
    'data/candidate_elongation_associated_genes.txt', 
    'Elongation',
    high_expr_genes
)

print("-" * 60)

# 4. Combine results
final_reference_set = sorted(list(ribo_genes | elong_genes))
all_annotation_matches = ribo_all | elong_all

# 5. Save filtered set
output_filename = "./data/CAI_Reference_Set_Mguttatus.txt"
with open(output_filename, "w") as f:
    for gene in final_reference_set:
        f.write(gene + "\n")

# 6. Also save annotation-only matches for comparison
annotation_only_file = "./data/CAI_Reference_Set_annotation_only.txt"
with open(annotation_only_file, "w") as f:
    for gene in sorted(all_annotation_matches):
        f.write(gene + "\n")

print(f"Total unique genes (annotation + expression): {len(final_reference_set)}")
print(f"Total annotation-only matches: {len(all_annotation_matches)}")
print("-" * 60)
print(f"Final reference set saved to: {output_filename}")
print(f"Annotation-only set saved to: {annotation_only_file}")

# 7. Show expression stats for selected genes
if expr_df is not None and len(final_reference_set) > 0:
    print("-" * 60)
    print("Expression summary for selected genes:")
    selected_expr = expr_df[expr_df['GeneID'].isin(final_reference_set)]
    if len(selected_expr) > 0:
        for col in [c for c in expr_df.columns if c.startswith('Exp_')]:
            print(f"  {col}: median={selected_expr[col].median():.1f}, range=[{selected_expr[col].min():.1f}, {selected_expr[col].max():.1f}]")
        print(f"\nSelected genes:")
        for gene in final_reference_set:
            gene_data = expr_df[expr_df['GeneID'] == gene]
            if len(gene_data) > 0:
                expr_vals = gene_data[[c for c in expr_df.columns if c.startswith('Exp_')]].values[0]
                print(f"  {gene}: {', '.join([f'{v:.1f}' for v in expr_vals])}")