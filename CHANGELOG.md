# Changelog

## [1.0.0] - 2024-10-15

### Added - Complete CUB Analysis Pipeline

#### Core Analysis Functions
- **calculate_rscu.R**: Relative Synonymous Codon Usage calculation
  - Calculates RSCU values for all codons per gene
  - Groups codons by amino acid
  - Excludes STOP codons from analysis
  - Handles zero counts gracefully

- **calculate_enc.R**: Effective Number of Codons calculation
  - Implements Wright (1990) method
  - Calculates homozygosity F for each degeneracy class
  - Handles 2-fold, 3-fold, 4-fold, and 6-fold degenerate codons
  - Returns ENC values ranging from 20 to 61

- **calculate_gc_content.R**: GC content metrics
  - Calculates 6 different GC metrics: GC, GC1, GC2, GC3, GC12, GC3s
  - GC3s excludes non-degenerate codons (Met and Trp)
  - Position-specific GC content analysis

#### Visualization Functions
- **visualize_codon_usage.R**: Codon usage visualization
  - Heatmap visualization with RSCU values
  - Barplot visualization by amino acid
  - Faceted plots for easy comparison
  - Color-coded by RSCU values
  - Includes sequence logo creation for specific amino acids
  - Custom color scheme for nucleotides (A=green, T=red, G=orange, C=blue)

#### Selection/Mutation/Drift Analysis
- **neutrality_analysis.R**: Evolutionary analysis plots
  - **Neutrality plot**: GC12 vs GC3 with regression line
    - Distinguishes mutation pressure from selection
    - Calculates correlation and slope
    - Includes interpretation guidelines
  - **ENC plot**: ENC vs GC3s with expected curve
    - Identifies genes under selection
    - Wright (1990) expected curve for mutation-drift equilibrium
    - Highlights genes below/above expected values
  - **PR2 bias plot**: Parity Rule 2 analysis
    - Analyzes purine/pyrimidine bias at 3rd position
    - Shows A/(A+T) vs G/(G+C)
    - Reveals strand asymmetry and mutational patterns

#### Comprehensive Analysis
- **cub_summary.R**: Integrated analysis pipeline
  - Runs all CUB analyses in one function
  - Generates all visualizations
  - Calculates summary statistics
  - Saves all results to organized output directory
  - Includes function to create logos for all amino acids

#### Integration
- **main.R**: Complete workflow script
  - Loads and sources all required functions
  - Reads transcript data
  - Performs codon quantification
  - Runs comprehensive CUB analysis
  - Generates all outputs
  - Includes commented examples for custom analyses

### Fixed
- **check_canonical_start.R**: Fixed variable name bug
  - Changed `trans` to `transcript_set` (line 16)
  - Now correctly references the input parameter

### Documentation
- **README.md**: Complete user guide
  - Overview of all features
  - Function reference
  - Interpretation guidelines
  - Output file descriptions
  - Requirements and installation

- **METHODS.md**: Theoretical background
  - Detailed methodology for each analysis
  - Mathematical formulas
  - Biological interpretation
  - References to key papers
  - Integration guidelines

- **QUICKSTART.md**: Quick reference guide
  - Installation instructions
  - Basic usage examples
  - Common tasks
  - Troubleshooting tips
  - Performance optimization

- **example_usage.R**: Practical examples
  - Complete automated analysis
  - Step-by-step manual analysis
  - Individual metric calculations
  - Custom visualizations
  - Amino acid-specific logos
  - Result exploration examples

- **test_functions.R**: Testing suite
  - Unit tests for all functions
  - Synthetic data generation
  - Validation of outputs
  - Creates test visualizations
  - Verifies all plots work correctly

### Dependencies
Updated to use only essential packages:
- data.table (efficient data manipulation)
- Biostrings (sequence handling)
- assertthat (input validation)
- stringi (string operations)
- foreach (iteration)
- doParallel (parallel processing)
- doFuture (advanced parallel backend)
- ggplot2 (visualization)

Removed:
- coRdon (not used in implementation)

### Performance
- Parallel processing support for large datasets
- Configurable number of cores
- Efficient data structures (data.table)
- Handles datasets with >50,000 genes

### Output Files
All results organized in `./results/` directory:
- CSV files with all metrics
- PDF plots for publication
- Summary statistics
- Amino acid-specific logos in subdirectory

### Quality
- All ggplot2 deprecation warnings fixed (size -> linewidth)
- Functions properly sourced in correct order
- Comprehensive error handling
- Informative messages during execution
- Code review feedback addressed

## References
Implementation follows these key publications:
1. Wright, F. (1990). Gene, 87(1), 23-29.
2. Sharp, P.M., & Li, W.H. (1987). NAR, 15(3), 1281-1295.
3. Sueoka, N. (1988). PNAS, 85(8), 2653-2657.
4. Sueoka, N. (1995). JME, 40(3), 318-325.
