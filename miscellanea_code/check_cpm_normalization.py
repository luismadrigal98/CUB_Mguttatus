#!/usr/bin/env python3
"""
Diagnostic script to verify CPM normalization and summarize gene-level distributions.

This script checks:
1. Per-sample CPM sums (should be ~1e6 for each plant)
2. Gene-level summary statistics (mean, median, max CPM)
3. Distribution skewness and proportion of genes with median == 0

Usage:
    python check_cpm_normalization.py <cpm_file.csv>
    
Example:
    python check_cpm_normalization.py results/leaf_expression_cpm.csv
"""

import pandas as pd
import numpy as np
import sys
import os


def check_cpm_normalization(filename):
    """
    Check CPM normalization and report gene-level summary statistics.
    """
    if not os.path.exists(filename):
        print(f"ERROR: File not found: {filename}")
        sys.exit(1)
    
    print(f"Reading CPM data from: {filename}")
    print("=" * 80)
    
    # Read CPM file
    df = pd.read_csv(filename, sep='\t')
    
    print(f"\nDataset dimensions: {df.shape[0]} genes × {df.shape[1]-1} samples")
    print("-" * 80)
    
    # ========================================================================
    # 1. CHECK PER-SAMPLE CPM SUMS
    # ========================================================================
    print("\n[1] PER-SAMPLE CPM NORMALIZATION CHECK")
    print("-" * 80)
    
    # Sum across genes for each sample (excluding gene_id column)
    sample_sums = df.iloc[:, 1:].sum(axis=0)
    
    print(f"Per-sample CPM sums:")
    print(f"  Min:  {sample_sums.min():.10f}")
    print(f"  Mean: {sample_sums.mean():.10f}")
    print(f"  Max:  {sample_sums.max():.10f}")
    
    # Check deviations from 1e6
    deviations = (sample_sums - 1e6).abs()
    max_dev = deviations.max()
    n_bad = (deviations > 1e-6).sum()
    
    print(f"\nDeviation from 1,000,000:")
    print(f"  Max absolute deviation: {max_dev:.2e}")
    print(f"  Samples with deviation > 1e-6: {n_bad}")
    
    if max_dev < 1e-6:
        print("  ✓ CPM normalization is CORRECT (all samples sum to 1e6)")
    else:
        print(f"  ⚠ WARNING: {n_bad} samples have deviations > 1e-6")
    
    print("\nSample sum quantiles:")
    quantiles = sample_sums.quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
    for q, val in quantiles.items():
        print(f"  {q*100:5.1f}%: {val:.6f}")
    
    # ========================================================================
    # 2. GENE-LEVEL SUMMARY STATISTICS
    # ========================================================================
    print("\n[2] GENE-LEVEL EXPRESSION DISTRIBUTIONS")
    print("-" * 80)
    
    # Calculate gene-level statistics across samples
    gene_means = df.iloc[:, 1:].mean(axis=1)
    gene_medians = df.iloc[:, 1:].median(axis=1)
    gene_maxs = df.iloc[:, 1:].max(axis=1)
    
    # Summary function
    def print_quantiles(series, name):
        print(f"\n{name} quantiles:")
        q = series.quantile([0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0])
        for perc, val in q.items():
            print(f"  {perc*100:5.1f}%: {val:>12.2f}")
    
    print_quantiles(gene_means, "Mean CPM per gene")
    print_quantiles(gene_medians, "Median CPM per gene")
    print_quantiles(gene_maxs, "Max CPM per gene")
    
    # ========================================================================
    # 3. DISTRIBUTION SKEWNESS ANALYSIS
    # ========================================================================
    print("\n[3] DISTRIBUTION SKEWNESS")
    print("-" * 80)
    
    # Proportion with median == 0
    prop_median_zero = (gene_medians == 0).mean()
    print(f"Proportion of genes with median CPM == 0: {prop_median_zero:.3f} ({prop_median_zero*100:.1f}%)")
    print(f"  → {int(prop_median_zero * len(gene_medians))} / {len(gene_medians)} genes")
    
    # Genes expressed in few vs many samples
    n_samples = df.shape[1] - 1
    samples_detected = (df.iloc[:, 1:] > 0).sum(axis=1)
    
    print(f"\nGenes detected (CPM > 0) in:")
    print(f"  All {n_samples} samples: {(samples_detected == n_samples).sum()}")
    print(f"  ≥ 95% of samples: {(samples_detected >= 0.95*n_samples).sum()}")
    print(f"  ≥ 50% of samples: {(samples_detected >= 0.5*n_samples).sum()}")
    print(f"  < 10% of samples: {(samples_detected < 0.1*n_samples).sum()}")
    print(f"  0 samples (not detected): {(samples_detected == 0).sum()}")
    
    # ========================================================================
    # 4. TOP GENES BY DIFFERENT METRICS
    # ========================================================================
    print("\n[4] TOP 10 GENES BY EACH METRIC")
    print("-" * 80)
    
    # Top by max CPM
    print("\nTop 10 by MAX CPM (peak expression):")
    top_max_idx = gene_maxs.sort_values(ascending=False).head(10).index
    top_max = pd.DataFrame({
        'gene_id': df.loc[top_max_idx, 'gene_id'].values,
        'mean': gene_means[top_max_idx].values,
        'median': gene_medians[top_max_idx].values,
        'max': gene_maxs[top_max_idx].values
    })
    print(top_max.to_string(index=False))
    
    # Top by mean CPM
    print("\nTop 10 by MEAN CPM (average expression):")
    top_mean_idx = gene_means.sort_values(ascending=False).head(10).index
    top_mean = pd.DataFrame({
        'gene_id': df.loc[top_mean_idx, 'gene_id'].values,
        'mean': gene_means[top_mean_idx].values,
        'median': gene_medians[top_mean_idx].values,
        'max': gene_maxs[top_mean_idx].values
    })
    print(top_mean.to_string(index=False))
    
    # Examples of highly skewed genes (high max, low median)
    print("\nExamples of HIGHLY SKEWED genes (high max, median ≈ 0):")
    skew_mask = (gene_maxs > gene_maxs.quantile(0.9)) & (gene_medians < 10)
    if skew_mask.sum() > 0:
        skewed_idx = gene_maxs[skew_mask].sort_values(ascending=False).head(10).index
        skewed = pd.DataFrame({
            'gene_id': df.loc[skewed_idx, 'gene_id'].values,
            'mean': gene_means[skewed_idx].values,
            'median': gene_medians[skewed_idx].values,
            'max': gene_maxs[skewed_idx].values
        })
        print(skewed.to_string(index=False))
    else:
        print("  No genes meet criteria")
    
    print("=" * 80)
    print("Analysis complete!")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_cpm_normalization.py <cpm_file.csv>")
        print("\nExample:")
        print("  python check_cpm_normalization.py results/leaf_expression_cpm.csv")
        sys.exit(1)
    
    filename = sys.argv[1]
    check_cpm_normalization(filename)
