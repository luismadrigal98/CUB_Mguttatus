#!/usr/bin/env python3
"""
Format leaf expression data to match bud expression file format.

Takes the leaf expression summary file and extracts a single metric (max_CPM, mean_CPM, or median_CPM)
to create a file compatible with bud_gene_expression_cpm_remapped.txt format (without Remapped_Gene field).

Output format: Gene | Expression

Usage:
    python format_leaf_expression.py <summary_file> <metric> <output_file>
    
Arguments:
    summary_file: Path to leaf_expression_summary.csv
    metric: Which metric to use (max_cpm, mean_cpm, or median_cpm)
    output_file: Path for output file
    
Example:
    python format_leaf_expression.py results/leaf_expression_summary.csv max_cpm data/leaf_gene_expression_max_cpm.txt
"""

import pandas as pd
import sys
import os


def format_expression_data(summary_file, metric, output_file):
    """
    Extract expression metric and format as Gene | Expression.
    """
    # Validate inputs
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        sys.exit(1)
    
    valid_metrics = ['max_cpm', 'mean_cpm', 'median_cpm']
    metric_lower = metric.lower()
    if metric_lower not in valid_metrics:
        print(f"ERROR: Invalid metric '{metric}'. Must be one of: {', '.join(valid_metrics)}")
        sys.exit(1)
    
    # Map user input to column name
    metric_col_map = {
        'max_cpm': 'max_CPM',
        'mean_cpm': 'mean_CPM',
        'median_cpm': 'median_CPM'
    }
    metric_col = metric_col_map[metric_lower]
    
    print(f"Reading summary file: {summary_file}")
    print(f"Extracting metric: {metric_col}")
    print("-" * 60)
    
    # Read summary file
    df = pd.read_csv(summary_file, sep='\t')
    
    # Check that required columns exist
    if 'gene_id' not in df.columns:
        print(f"ERROR: 'gene_id' column not found in {summary_file}")
        sys.exit(1)
    if metric_col not in df.columns:
        print(f"ERROR: '{metric_col}' column not found in {summary_file}")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # Extract gene and expression columns
    output_df = df[['gene_id', metric_col]].copy()
    output_df.columns = ['Gene', 'Expression']
    
    # Sort by expression (descending) for easier inspection
    output_df = output_df.sort_values('Expression', ascending=False)
    
    # Report statistics
    print(f"\nDataset statistics:")
    print(f"  Total genes: {len(output_df)}")
    print(f"  Genes with Expression > 0: {(output_df['Expression'] > 0).sum()}")
    print(f"  Genes with Expression == 0: {(output_df['Expression'] == 0).sum()}")
    print(f"  Min expression: {output_df['Expression'].min():.6f}")
    print(f"  Max expression: {output_df['Expression'].max():.2f}")
    print(f"  Mean expression: {output_df['Expression'].mean():.2f}")
    print(f"  Median expression: {output_df['Expression'].median():.2f}")
    
    print(f"\nTop 10 genes by {metric_col}:")
    print(output_df.head(10).to_string(index=False))
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\nCreated directory: {output_dir}")
    
    # Write output file (tab-delimited)
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ Output written to: {output_file}")
    print(f"  Format: Gene | Expression")
    print(f"  Rows: {len(output_df)}")
    print(f"  Sorted: by Expression (descending)")
    
    return output_df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python format_leaf_expression.py <summary_file> <metric> <output_file>")
        print("\nArguments:")
        print("  summary_file: Path to leaf_expression_summary.csv")
        print("  metric: Which metric to use (max_cpm, mean_cpm, or median_cpm)")
        print("  output_file: Path for output file")
        print("\nExample:")
        print("  python format_leaf_expression.py results/leaf_expression_summary.csv max_cpm data/leaf_gene_expression_max_cpm.txt")
        sys.exit(1)
    
    summary_file = sys.argv[1]
    metric = sys.argv[2]
    output_file = sys.argv[3]
    
    format_expression_data(summary_file, metric, output_file)