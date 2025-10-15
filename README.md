# Codon Usage Bias Analysis in Mimulus guttatus

A comprehensive R-based pipeline for analyzing codon usage bias (CUB) in *Mimulus guttatus* transcriptomes.

## Overview

This repository provides tools to analyze codon usage patterns, evaluate selection/mutation/drift balance, and create visualizations for understanding evolutionary forces shaping codon preferences.

## Features

### Core Analyses
- **Codon quantification**: Count codon occurrences across all genes
- **RSCU (Relative Synonymous Codon Usage)**: Measure relative usage of synonymous codons
- **ENC (Effective Number of Codons)**: Quantify overall codon bias (20-61 scale)
- **GC content metrics**: Calculate GC, GC1, GC2, GC3, GC12, and GC3s

### Evolutionary Analysis
- **Neutrality Plot**: Distinguish mutation pressure from selection (GC12 vs GC3)
- **ENC Plot**: Identify genes under selection for codon bias (ENC vs GC3s)
- **PR2 Bias Plot**: Analyze purine/pyrimidine bias at 3rd codon position

### Visualizations
- **Codon usage heatmaps**: Genome-wide RSCU patterns
- **Codon usage barplots**: RSCU by amino acid
- **Sequence logos**: Nucleotide preferences at each codon position per amino acid

## Requirements

```r
required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel',
                        'doFuture', 'ggplot2')
```

## Quick Start

```r
# Run the complete analysis pipeline
source("main.R")
```

The pipeline will:
1. Load transcript data from FASTA file
2. Count codons across all genes
3. Calculate all CUB metrics
4. Generate all plots and visualizations
5. Save results to `./results/` directory

## Input Data

The pipeline expects a FASTA file with coding sequences (CDS):
- `./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa`

Each sequence should:
- Start with canonical ATG start codon
- Have length divisible by 3 (valid reading frame)
- Contain only coding sequence (introns removed)

## Output Files

All results are saved to `./results/`:

### Data Files
- `cub_analysis_complete.csv`: Complete dataset with all metrics
- `enc_values.csv`: ENC values per gene
- `gc_content.csv`: GC content metrics per gene
- `summary_statistics.csv`: Genome-wide summary statistics

### Plots
- `codon_usage_heatmap.pdf`: Heatmap of RSCU values
- `codon_usage_barplot.pdf`: Barplot of RSCU by amino acid
- `neutrality_plot.pdf`: GC12 vs GC3 plot
- `enc_plot.pdf`: ENC vs GC3s plot
- `pr2_plot.pdf`: Parity Rule 2 bias plot
- `codon_logos/*.pdf`: Individual logos for each amino acid

## Function Reference

### Data Processing
- `codon_quant()`: Count codons across transcriptome
- `calculate_rscu()`: Calculate relative synonymous codon usage
- `calculate_enc()`: Calculate effective number of codons
- `calculate_gc_content()`: Calculate GC content metrics

### Visualization
- `visualize_codon_usage()`: Create heatmaps and barplots
- `create_codon_logo()`: Create sequence logo for amino acid

### Evolutionary Analysis
- `neutrality_plot()`: Create neutrality plot (mutation vs selection)
- `enc_plot()`: Create ENC plot (identify genes under selection)
- `pr2_bias_plot()`: Create PR2 bias plot (purine/pyrimidine bias)

### Comprehensive Analysis
- `cub_summary()`: Run complete CUB analysis pipeline
- `create_aa_specific_logos()`: Generate logos for all amino acids

## Interpretation Guide

### RSCU Values
- **RSCU = 1**: Codon used at expected frequency
- **RSCU > 1**: Codon used more than expected (positive bias)
- **RSCU < 1**: Codon used less than expected (negative bias)

### ENC Values
- **ENC = 20**: Extreme bias (one codon per amino acid)
- **ENC = 61**: No bias (all codons used equally)
- Lower values indicate stronger codon bias

### Neutrality Plot
- **Slope ≈ 1**: Mutation pressure dominates
- **Slope ≈ 0**: Selection dominates
- Strong correlation suggests mutation bias

### ENC Plot
- **Genes on curve**: Mutation-drift equilibrium
- **Genes below curve**: Under selection for codon bias
- **Genes above curve**: Other factors (e.g., amino acid composition)

## Authors

Luis Javier Madrigal-Roca & John K. Kelly

## License

[Add license information]
