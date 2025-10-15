# Quick Start Guide

## Installation

### 1. Install R packages
```r
# Required packages
install.packages(c(
  "data.table",
  "assertthat", 
  "stringi",
  "foreach",
  "doParallel",
  "doFuture",
  "ggplot2"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "coRdon"))
```

### 2. Clone the repository
```bash
git clone https://github.com/luismadrigal98/CUB_Mguttatus.git
cd CUB_Mguttatus
```

## Basic Usage

### Option 1: Complete Analysis (Recommended)
Run the entire pipeline with one command:

```r
source("main.R")
```

This will:
- Load your transcript data
- Count all codons
- Calculate all CUB metrics
- Generate all plots
- Save results to `./results/`

### Option 2: Step-by-Step Analysis
For more control, follow the examples in `example_usage.R`:

```r
source("example_usage.R")
```

### Option 3: Test with Synthetic Data
Verify installation without needing the full dataset:

```r
source("test_functions.R")
```

## Output Files

After running the analysis, check the `./results/` directory:

### Data Files (CSV)
- `cub_analysis_complete.csv` - All metrics for all genes
- `enc_values.csv` - ENC values per gene
- `gc_content.csv` - GC metrics per gene
- `summary_statistics.csv` - Genome-wide statistics

### Plots (PDF)
- `codon_usage_heatmap.pdf` - RSCU heatmap
- `codon_usage_barplot.pdf` - RSCU barplot by amino acid
- `neutrality_plot.pdf` - Mutation vs selection
- `enc_plot.pdf` - Genes under selection
- `pr2_plot.pdf` - Base composition bias
- `codon_logos/*.pdf` - One logo per amino acid

## Common Tasks

### Analyze Your Own Data
1. Place your FASTA file in `./data/`
2. Edit `main.R` line 66 to point to your file:
```r
trans <- Biostrings::readDNAStringSet(
  filepath = "./data/YOUR_FILE.fa", 
  format = 'fasta'
)
```
3. Run: `source("main.R")`

### Calculate Only Specific Metrics

```r
# Load your data first
source("./src/set_environment.R")
set_environment(required_pckgs = c('data.table', 'Biostrings', 'ggplot2'), 
                parallel_backend = FALSE)
source("./src/codon_quant.R")
source("./src/calculate_enc.R")

# Calculate ENC only
enc_values <- calculate_enc(codon_usage, genetic_code_dna_long)
```

### Create Custom Visualizations

```r
# After running analysis
source("./src/visualize_codon_usage.R")

# Custom heatmap
visualize_codon_usage(codon_usage, genetic_code_dna_long, 
                     "my_custom_heatmap.pdf", type = "heatmap")

# Logo for specific amino acid
create_codon_logo(codon_usage, genetic_code_dna_long, "Leu", 
                 "leucine_usage.pdf")
```

### Analyze Subset of Genes

```r
# Filter genes by some criterion
high_expression_genes <- codon_usage[1:1000, ]

# Run analysis on subset
subset_results <- cub_summary(high_expression_genes, genetic_code_dna_long,
                             output_dir = "./results/high_expression")
```

## Troubleshooting

### Error: "Cannot find data file"
- Check that your FASTA file exists in `./data/`
- Use absolute path if relative path doesn't work

### Error: "Package not found"
- Install missing packages (see Installation section)
- Check that Bioconductor packages are installed correctly

### Error: "Cannot allocate memory" (large datasets)
- Reduce number of cores: change `n_cores = 10` to `n_cores = 2` in main.R
- Process genes in batches
- Disable parallel processing: `parallel = FALSE`

### Plots look weird
- Check that ggplot2 is up to date: `update.packages("ggplot2")`
- Adjust plot size in ggsave() calls if needed

### Long computation time
- Enable parallel processing: `parallel = TRUE`
- Increase cores: `n_cores = 10` (or more if available)
- Use filtered dataset (remove low-quality genes first)

## Performance Tips

### For Large Datasets (>50,000 genes)
```r
# Use parallel processing
set_environment(parallel_backend = TRUE, n_cores = 10)

# Run codon quantification in parallel
codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = TRUE)
```

### For Small Datasets (<1,000 genes)
```r
# Sequential processing may be faster
codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = FALSE)
```

## Next Steps

- Read `README.md` for detailed documentation
- Check `METHODS.md` for theoretical background
- Explore `example_usage.R` for advanced examples
- Modify functions in `./src/` to customize analyses

## Getting Help

- Check function documentation with `?function_name`
- Review code comments in `./src/` files
- See `example_usage.R` for working examples
- Consult `METHODS.md` for interpretation guidance

## Citation

If you use this pipeline, please cite:

```
Madrigal-Roca, L.J. & Kelly, J.K. (2024). 
CUB_Mguttatus: Comprehensive Codon Usage Bias Analysis Pipeline.
GitHub repository: https://github.com/luismadrigal98/CUB_Mguttatus
```

And the key methodological papers listed in `METHODS.md`.
