# Codon Usage Bias Analysis in Mimulus guttatus

A comprehensive R-based pipeline for analyzing codon usage bias (CUB) in *Mimulus guttatus* transcriptomes, integrating population genomics, expression data, and evolutionary modeling.

## Overview

This repository provides advanced tools to analyze codon usage patterns, evaluate selection/mutation/drift balance, integrate polymorphism data, and model translational selection. The pipeline supports both classical metrics (RSCU, ENC, CAI) and modern Bayesian approaches (AnaCoDa framework).

## Features

### Core Codon Usage Metrics
- **Codon quantification**: Count codon occurrences across all genes
- **RSCU (Relative Synonymous Codon Usage)**: Measure relative usage of synonymous codons
- **ENC (Effective Number of Codons)**: Quantify overall codon bias (20-61 scale)
  - Wright (1990) method with family splitting
  - Sun et al. (2012) pseudocount correction
- **CDC (Codon Deviation Coefficient)**: Quantify deviation from expected usage patterns
- **CAI (Codon Adaptation Index)**: Measure adaptation to optimal codons
- **RF (Relative Frequency)**: Alternative bias metric
- **GC content metrics**: Calculate GC, GC1, GC2, GC3, GC12, and GC3s

### Expression-Based Analysis
- **Integration with RNA-seq data**: Multi-tissue expression profiles (leaf, bud)
- **Expression-stratified analysis**: Compare high vs. low expression genes
- **CAI vs. Expression correlation**: Validate translational selection
- **Detrended analysis**: Control for confounding factors (gene length, GC content)

### Population Genomics Integration
- **Polymorphism-based selection inference**: Estimate selection coefficients from diversity patterns
- **Intron-based mutation rate estimation**: Calculate dM from neutral regions
- **Site frequency spectrum (SFS) analysis**: Analyze allele frequency distributions
- **π (nucleotide diversity) vs. selection**: Validate "hump effect" predictions
- **Preferred vs. non-preferred codon diversity**: Compare polymorphism levels

### Bayesian Modeling (AnaCoDa Framework)
- **ROC (Ribosome Overhead Cost) model**: Infer selection on codon usage
- **Multi-tissue expression integration**: Model tissue-specific codon preferences
- **Mutation bias estimation**: Fix or estimate dM from neutral regions
- **MCMC convergence diagnostics**: Gelman-Rubin statistics across chains
- **Codon frequency trajectories**: Validate model predictions vs. empirical data

### Statistical Testing
- **Multinomial models**: Test codon proportion differences
- **Binomial models**: Pairwise codon comparisons
- **G-test for independence**: Family-level codon usage tests
- **FDR correction**: Multiple testing correction (Benjamini-Hochberg)
- **Effect size calculation**: Cohen's d for biological significance

### Advanced Visualizations
- **Enhanced biplots**: PCA/CA with codon loadings and gene scores
- **ENC plots with CDC highlighting**: Identify genes with significant deviation
- **Neutrality plots**: Mutation vs. selection balance (GC12 vs. GC3)
- **PR2 bias plots**: Purine/pyrimidine patterns at 3rd position
- **Codon frequency heatmaps**: Genome-wide RSCU patterns
- **Amino acid-specific logos**: Position-specific nucleotide preferences
- **3D PCA videos**: Dynamic visualization of codon usage space
- **ROC trajectory plots**: Model validation with empirical data

## Requirements

```r
required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel', 'doFuture',
                        'ggplot2', 'dplyr', 'tidyr', 'mgcv', 'MASS',
                        'dunn.test', 'ggrepel', 'RColorBrewer',
                        'gridExtra', 'cowplot')
```

### System Requirements
- R >= 4.0
- Sufficient RAM for large transcriptomes (>8GB recommended)
- Multi-core processor recommended for parallel processing

## Quick Start

```r
# Run the complete analysis pipeline
source("main.R")
```

The pipeline will:
1. Load transcript data from FASTA file
2. Count codons across all genes
3. Calculate all CUB metrics (RSCU, ENC, CDC, CAI, GC content)
4. Integrate expression data (if available)
5. Perform CDC-based selection analysis
6. Calculate CAI and identify optimal codons
7. Estimate mutation rates from neutral regions
8. Run AnaCoDa Bayesian modeling (optional)
9. Perform polymorphism-based selection inference
10. Generate all plots and visualizations
11. Save results to `./results/` directory

## Input Data

### Required Files
- **CDS sequences**: `./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnlyClean.fa`
  - FASTA file with coding sequences
  - Must have canonical ATG start codon
  - Length divisible by 3 (valid reading frame)
  - Stop codons removed

### Optional Files (for advanced analyses)
- **Expression data**:
  - `./data/bud_gene_expression_cpm_remapped.txt`
  - `./data/leaf_gene_expression_mean_cpm_renamed.txt`
- **VCF files**: For polymorphism analysis (one per chromosome)
- **GFF3 annotation**: `./data/Mguttatusvar_IM767_887_v2.1.gene.gff3`
- **tRNA data**: `./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt`
- **Neutral regions**:
  - Intron sequences for mutation rate estimation
  - Intergenic sequences (alternative)

## Output Files

All results are saved to `./results/`:

### Core Analysis Data
- `cub_analysis_complete.csv`: Complete dataset with all metrics
- `enc_values.csv`: ENC values per gene (Wright and/or Sun methods)
- `cdc_results.csv`: Codon Deviation Coefficient analysis
- `gc_content.csv`: GC content metrics per gene
- `cai_values.csv`: Codon Adaptation Index per gene
- `integrated_data.csv`: Merged CUB metrics with expression data
- `summary_statistics.csv`: Genome-wide summary statistics

### CAI Analysis
- `optimal_codons_relative_adaptiveness.pdf`: W-values for optimal codons
- `CAI_by_expression_group.pdf`: CAI distribution across expression levels
- `CAI_vs_ENC_scatter.pdf`: CAI vs. CDC correlation
- `preferred_codon_usage_comparison.pdf`: Enriched codon usage in highly expressed genes
- `codon_proportion_test_results.csv`: Statistical tests for codon enrichment
- `codon_classification_corrected.csv`: Corrected preferred codon classification

### Selection Analysis Plots
- `ENC_plot_CDC_highlighted.pdf`: ENC vs. GC3s with CDC-significant genes
- `ENC_deviation_by_CDC.pdf`: Distribution of ENC deviations
- `neutrality_plot.pdf`: GC12 vs. GC3 plot (mutation vs. selection)
- `pr2_plot.pdf`: Parity Rule 2 bias plot
- `Detrended_ENC_by_expression_group.pdf`: Length-corrected CDC analysis

### Codon Usage Visualizations
- `codon_usage_heatmap.pdf`: Heatmap of RSCU values
- `codon_usage_barplot.pdf`: Barplot of RSCU by amino acid
- `codon_frequency_top5_vs_rest.pdf`: Absolute frequency comparison
- `codon_logos/*.pdf`: Individual sequence logos for each amino acid

### Polymorphism Analysis
- `diversity_modeling/`: Directory with polymorphism-based selection inference
  - `pi_vs_selection_binned.pdf`: Nucleotide diversity vs. selection coefficient
  - `sfs_comparison_*.pdf`: Site frequency spectra for preferred/non-preferred sites
  - `selection_coefficient_estimates.csv`: Estimated selection coefficients

### AnaCoDa Results (if run)
- `MCMC_results/`: Bayesian modeling outputs
  - Parameter traces, convergence diagnostics
  - Selection coefficients (dEta) per codon
  - Gene expression estimates (phi)
- `ROC_codon_trajectories.pdf`: Model validation plot

## Function Reference

### Core CUB Metrics
- `codon_quant()`: Count codons across transcriptome
- `calculate_rscu()`: Calculate relative synonymous codon usage
- `calculate_enc()`: Calculate effective number of codons (Wright or Sun methods)
- `calculate_cdc_all()`: Calculate codon deviation coefficient
- `calculate_cai()`: Calculate codon adaptation index
- `calculate_rf()`: Calculate relative frequency
- `calculate_gc_content()`: Calculate GC content metrics

### Expression Integration
- `integrate_cdc_analysis()`: Combine CDC with expression data
- `diagnose_cai_vs_proportion()`: Reconcile CAI w-values with statistical tests

### Statistical Testing
- `test_codon_proportions()`: Per-codon proportion tests with FDR correction
- `fit_pairwise_binomial_models()`: Pairwise codon comparisons
- `fit_bi_multinom_family_model()`: Multinomial family-level tests
- `CUB_g_test()`: G-test for codon usage independence
- `cohens_d_calc()`: Effect size calculation

### Polymorphism Analysis
- `derivation_gamma_from_polymorphism()`: Estimate selection from diversity
- `integrate_intronic_polymorphism()`: Process VCF data for introns
- `local_M_estimation()`: Estimate mutation rates from neutral regions
- `null_and_empirical_sfs()`: Compare observed vs. expected SFS

### Visualization
- `visualize_codon_usage()`: Create heatmaps and barplots
- `create_aa_logo()`: Create sequence logo for amino acid
- `enhanced_biplot()`: PCA/CA biplots with loadings
- `plot_scurves()`: S-curve visualization for selection inference
- `create_3d_pca_video()`: Animated 3D PCA visualization

### AnaCoDa Integration
- `GR_convergence()`: Gelman-Rubin MCMC diagnostics
- `roc_model_validation()`: Validate ROC model predictions
- `run_trajectory_analysis()`: Codon frequency vs. expression trajectories

### Utility Functions
- `cub_summary()`: Run complete CUB analysis pipeline
- `check_canonical_start()`: Validate start codons
- `check_cds()`: Validate coding sequences
- `trim_uninformative()`: Remove non-informative codons (Met, Trp, STOP)
- `count_preferred_by_aa()`: Count preferred codon usage per amino acid
- `count_preferred_by_family()`: Count preferred codon usage per family

## Interpretation Guide

### RSCU Values
- **RSCU = 1**: Codon used at expected frequency (no bias)
- **RSCU > 1**: Codon used more than expected (positive bias)
- **RSCU < 1**: Codon used less than expected (negative bias)

### ENC Values
- **ENC = 20**: Extreme bias (one codon per amino acid)
- **ENC = 61**: No bias (all synonymous codons used equally)
- **20-40**: High codon bias (strong selection)
- **40-50**: Moderate codon bias
- **50-61**: Low codon bias (weak selection)

### CDC (Codon Deviation Coefficient)
- Measures deviation from expected codon frequencies
- Statistical significance tested via chi-square or G-test
- p < 0.05 indicates significant deviation from neutral expectations

### CAI (Codon Adaptation Index)
- **CAI = 0-1**: Ranges from 0 (not adapted) to 1 (perfectly adapted)
- **CAI > 0.7**: Highly adapted to optimal codons (typical for highly expressed genes)
- **CAI < 0.5**: Poorly adapted (typical for lowly expressed genes)

### Neutrality Plot (GC12 vs. GC3)
- **Slope ≈ 1**: Mutation pressure dominates
- **Slope ≈ 0**: Selection dominates
- **Intermediate slope (0.3-0.7)**: Both forces contribute

### ENC Plot (ENC vs. GC3s)
- **Genes on Wright curve**: Mutation-drift equilibrium
- **Genes below curve**: Under selection for codon bias
- **Genes above curve**: Unusual amino acid composition or other factors

### Selection Coefficients from Polymorphism
- **Ns > 1**: Strong selection (effective selection)
- **0.1 < Ns < 1**: Weak selection
- **Ns < 0.1**: Nearly neutral
- Negative Ns indicates purifying selection against the focal codon

## Advanced Features

### CDC Analysis with Statistical Testing
The pipeline performs comprehensive statistical testing of codon usage deviations:
- Chi-square goodness-of-fit tests per gene
- Multiple testing correction (FDR)
- Integration with expression data to identify selection signatures

### Expression-Stratified Analysis
Compare codon usage patterns across expression levels:
- Top 5% highly expressed genes
- Middle 90% genes
- Bottom 5% lowly expressed genes
- Statistical tests (Kruskal-Wallis, Dunn's post-hoc)
- Effect sizes (Cohen's d)

### Polymorphism-Based Selection Inference
Estimate selection coefficients using the "hump effect":
- Higher diversity at weakly selected sites (0.1 < Ns < 1)
- Lower diversity at strongly selected sites (Ns > 1)
- Compare preferred vs. non-preferred codons
- Validate with site frequency spectra

### AnaCoDa Bayesian Framework
Run Markov Chain Monte Carlo (MCMC) to infer:
- Selection coefficients (dEta) per codon
- Gene expression levels (phi)
- Mutation bias parameters (dM)
- Multiple chain convergence diagnostics
- Model validation with empirical data

### Multivariate Analysis
- Principal Component Analysis (PCA) of codon usage
- Correspondence Analysis (CA)
- Enhanced biplots with codon loadings
- 3D visualization with animated GIFs/videos

## Pipeline Workflow

### Phase 1: Data Preparation
1. Load and validate CDS sequences
2. Filter for canonical start codons
3. Check reading frames
4. Count codons per gene

### Phase 2: Core Metrics
1. Calculate RSCU values
2. Compute ENC (Wright and/or Sun methods)
3. Calculate CDC with statistical tests
4. Estimate GC content metrics

### Phase 3: Expression Integration (Optional)
1. Load RNA-seq data
2. Normalize expression values (CPM, log2 transformation)
3. Merge with codon usage data
4. Stratify genes by expression level

### Phase 4: CAI Analysis
1. Define reference set (highly expressed genes)
2. Calculate relative adaptiveness (w-values)
3. Compute CAI for all genes
4. Test for enrichment of preferred codons
5. Statistical validation with proportion tests

### Phase 5: Polymorphism Analysis (Optional)
1. Process VCF files
2. Estimate mutation rates from neutral regions
3. Calculate nucleotide diversity (π) per site
4. Classify sites by preference
5. Estimate selection coefficients
6. Validate "hump effect"

### Phase 6: AnaCoDa Modeling (Optional - requires external execution)
1. Prepare input files (sequences, expression data, mutation rates)
2. Run MCMC chains (command-line execution)
3. Check convergence (Gelman-Rubin diagnostics)
4. Extract selection coefficients
5. Validate model predictions

### Phase 7: Visualization and Reporting
1. Generate all plots
2. Create summary statistics
3. Export results to CSV files
4. Produce publication-quality figures

## Testing

Unit tests are provided for core functions:

```bash
# Run tests for calculate_enc function
cd tests
Rscript test_calculate_enc.R
```

Tests cover:
- Extreme bias scenarios (ENC ≈ 20)
- No bias scenarios (ENC ≈ 61)
- Multiple genes
- Both Wright and Sun methods
- Edge cases and error handling

## Project Structure

```
Codon_bias_analysis/
├── main.R                    # Main analysis pipeline
├── README.md                 # This file
├── data/                     # Input data files
│   ├── *cds*.fa             # CDS sequences
│   ├── *expression*.txt     # RNA-seq data
│   ├── *.vcf                # Polymorphism data
│   └── *.gff3               # Genome annotation
├── src/                      # Source R scripts (50+ functions)
│   ├── calculate_*.R        # Core metric calculations
│   ├── *_analysis.R         # Analysis workflows
│   ├── plot_*.R             # Visualization functions
│   └── *.R                  # Utility functions
├── tests/                    # Unit tests
│   └── test_*.R
├── doc/                      # Documentation
│   ├── METHODS.md           # Detailed methodology
│   ├── VALIDATION_FRAMEWORK.md
│   └── POLYMORPHISM_BASED_SELECTION_INFERENCE.md
├── results/                  # Output directory
│   ├── *.csv                # Data outputs
│   ├── *.pdf                # Plots
│   ├── diversity_modeling/  # Polymorphism analysis
│   └── MCMC_results/        # AnaCoDa outputs
└── examples/                 # Example scripts
    └── example_*.R
```

## Citation

If you use this pipeline in your research, please cite:

- Wright, F. (1990). The 'effective number of codons' used in a gene. Gene, 87(1), 23-29.
- Sun, X., et al. (2012). An improved implementation of effective number of codons (Nc). Molecular Biology and Evolution, 30(1), 191-196.
- Sharp, P.M., & Li, W.H. (1987). The codon adaptation index. Nucleic Acids Research, 15(3), 1281-1295.

## Authors

**Luis Javier Madrigal-Roca** & **John K. Kelly**

## License

[Specify license]

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

## Acknowledgments

This pipeline integrates methods from multiple seminal papers in codon usage research and builds upon the AnaCoDa framework for Bayesian modeling of translational selection.
