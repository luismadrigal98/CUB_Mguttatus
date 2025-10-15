# CUB Analysis Workflow

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT: FASTA FILE                        │
│           (Coding sequences - transcripts)                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│              QUALITY CONTROL & FILTERING                    │
│  • Check canonical start codon (ATG)                        │
│  • Validate reading frame (length % 3 == 0)                 │
│  • Extract gene names                                       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                 CODON QUANTIFICATION                        │
│  • Split sequences into codons                              │
│  • Count each of 64 codons per gene                         │
│  • Create gene × codon count matrix                         │
└────────────────────┬────────────────────────────────────────┘
                     │
         ┌───────────┴───────────┐
         │                       │
         ▼                       ▼
┌──────────────────┐   ┌──────────────────┐
│   BIAS METRICS   │   │  COMPOSITION     │
│                  │   │                  │
│  • RSCU          │   │  • GC content    │
│  • ENC           │   │  • GC1, GC2, GC3 │
│  • CAI (future)  │   │  • GC12, GC3s    │
└────────┬─────────┘   └────────┬─────────┘
         │                      │
         └──────────┬───────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│            EVOLUTIONARY ANALYSIS                            │
│                                                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌────────────┐ │
│  │ NEUTRALITY PLOT │  │   ENC PLOT      │  │  PR2 PLOT  │ │
│  │  (GC12 vs GC3)  │  │  (ENC vs GC3s)  │  │  (A/T G/C) │ │
│  │                 │  │                 │  │            │ │
│  │ Mutation vs     │  │ Genes under     │  │ Strand     │ │
│  │ Selection       │  │ Selection       │  │ Asymmetry  │ │
│  └─────────────────┘  └─────────────────┘  └────────────┘ │
└─────────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                   VISUALIZATIONS                            │
│                                                             │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐ │
│  │   HEATMAP    │  │   BARPLOT    │  │  SEQUENCE LOGOS  │ │
│  │              │  │              │  │                  │ │
│  │ Genome-wide  │  │ By amino     │  │ Per amino acid   │ │
│  │ RSCU pattern │  │ acid group   │  │ codon preference │ │
│  └──────────────┘  └──────────────┘  └──────────────────┘ │
└─────────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                  OUTPUT FILES                               │
│                                                             │
│  DATA (CSV):                    PLOTS (PDF):                │
│  • cub_analysis_complete.csv    • codon_usage_heatmap.pdf  │
│  • enc_values.csv               • codon_usage_barplot.pdf  │
│  • gc_content.csv               • neutrality_plot.pdf      │
│  • summary_statistics.csv       • enc_plot.pdf             │
│                                 • pr2_plot.pdf             │
│                                 • codon_logos/*.pdf        │
└─────────────────────────────────────────────────────────────┘
```

## Function Call Flow

### Main Workflow (main.R)
```r
main.R
  │
  ├─> set_environment()              [src/set_environment.R]
  │     └─> Load libraries, set parallel backend
  │
  ├─> readDNAStringSet()             [Biostrings package]
  │     └─> Load FASTA file
  │
  ├─> codon_quant()                  [src/codon_quant.R]
  │     ├─> check_canonical_start()  [src/check_canonical_start.R]
  │     ├─> splitInPartsAux()        [src/splitInPartsAux.R]
  │     └─> codons_counter()         [src/codons_counter.R]
  │           └─> gene_name_extractor() [src/gene_name_extractor.R]
  │
  ├─> cub_summary()                  [src/cub_summary.R]
  │     ├─> calculate_rscu()         [src/calculate_rscu.R]
  │     ├─> calculate_enc()          [src/calculate_enc.R]
  │     ├─> calculate_gc_content()   [src/calculate_gc_content.R]
  │     ├─> visualize_codon_usage()  [src/visualize_codon_usage.R]
  │     ├─> neutrality_plot()        [src/neutrality_analysis.R]
  │     ├─> enc_plot()               [src/neutrality_analysis.R]
  │     └─> pr2_bias_plot()          [src/neutrality_analysis.R]
  │
  └─> create_aa_specific_logos()     [src/cub_summary.R]
        └─> create_codon_logo()      [src/visualize_codon_usage.R]
```

## Data Flow

### Input Data Structure
```
DNAStringSet object
├─> Gene1: "ATGTTCGCA..."
├─> Gene2: "ATGAAAGGT..."
└─> Gene3: "ATGCCCTAT..."
```

### Codon Count Matrix
```
              TTT TTC TTA TTG ... GGG
Gene1         3   5   2   1   ... 4
Gene2         1   2   6   3   ... 2
Gene3         4   1   1   2   ... 5
```

### Analysis Results
```
RSCU Matrix          ENC Values        GC Content
(Gene × Codon)       (Gene × 1)        (Gene × 6)
┌──────────┐         ┌──────┐          ┌──────────┐
│ 0.8 1.2  │         │ 45.2 │          │ GC  0.52 │
│ 1.5 0.5  │         │ 38.7 │          │ GC3 0.48 │
│ 0.9 1.1  │         │ 52.1 │          │ ... ...  │
└──────────┘         └──────┘          └──────────┘
```

## Parallel Processing

```
Sequential Mode              Parallel Mode
─────────────────           ─────────────────────────
Gene1 → Process             ┌─> Gene1 → Process
  ↓                         │
Gene2 → Process             ├─> Gene2 → Process
  ↓                    ────>│
Gene3 → Process             ├─> Gene3 → Process
  ↓                         │
...                         └─> ... → Process
                                    ↓
                              Combine Results
```

## Usage Examples

### Quick Start (One Command)
```r
source("main.R")
```

### Custom Analysis
```r
# 1. Load data
trans <- readDNAStringSet("data/your_file.fa")

# 2. Count codons
codon_counts <- codon_quant(trans, names(genetic_code), parallel = TRUE)

# 3. Calculate specific metric
enc <- calculate_enc(codon_counts, genetic_code)

# 4. Create specific plot
neutrality_plot(gc_content, "my_plot.pdf")
```

### Batch Processing
```r
# Process multiple datasets
datasets <- c("species1.fa", "species2.fa", "species3.fa")

for(dataset in datasets) {
  trans <- readDNAStringSet(paste0("data/", dataset))
  codon_counts <- codon_quant(trans, names(genetic_code))
  results <- cub_summary(codon_counts, genetic_code,
                        output_dir = paste0("results/", dataset))
}
```

## Decision Tree for Analysis

```
Start
  │
  ├─ Have full dataset? ──Yes──> Use main.R (full pipeline)
  │                              └─> All analyses + plots
  │
  └─ No/Want specific analysis?
      │
      ├─ Need RSCU? ──Yes──> calculate_rscu()
      │
      ├─ Need ENC? ──Yes──> calculate_enc()
      │
      ├─ Need GC content? ──Yes──> calculate_gc_content()
      │
      ├─ Need selection analysis? ──Yes──> neutrality_plot()
      │                                     enc_plot()
      │                                     pr2_bias_plot()
      │
      └─ Need visualizations? ──Yes──> visualize_codon_usage()
                                       create_codon_logo()
```

## Performance Considerations

| Dataset Size | Recommended Settings        | Expected Time   |
|-------------|----------------------------|-----------------|
| < 1,000     | parallel = FALSE           | < 1 minute      |
| 1,000-10K   | parallel = TRUE, cores = 4 | 1-5 minutes     |
| 10K-50K     | parallel = TRUE, cores = 8 | 5-20 minutes    |
| > 50K       | parallel = TRUE, cores = 10| 20+ minutes     |

## File Organization

```
CUB_Mguttatus/
├── main.R                      # Main workflow script
├── example_usage.R             # Usage examples
├── test_functions.R            # Testing suite
├── README.md                   # User guide
├── METHODS.md                  # Theory & methods
├── QUICKSTART.md               # Quick reference
├── WORKFLOW.md                 # This file
├── CHANGELOG.md                # Version history
│
├── src/                        # Source functions
│   ├── set_environment.R       # Setup
│   ├── codon_quant.R          # Counting
│   ├── calculate_rscu.R       # Metrics
│   ├── calculate_enc.R        # Metrics
│   ├── calculate_gc_content.R # Metrics
│   ├── visualize_codon_usage.R# Plots
│   ├── neutrality_analysis.R  # Evolutionary
│   ├── cub_summary.R          # Integration
│   └── [helper functions]     # Utilities
│
├── data/                       # Input data
│   └── *.fa                   # FASTA files
│
└── results/                    # Output (created)
    ├── *.csv                  # Data files
    ├── *.pdf                  # Plots
    └── codon_logos/           # AA-specific logos
```

## Troubleshooting Flow

```
Problem?
  │
  ├─ Installation issue?
  │   └─> Check QUICKSTART.md → Installation section
  │
  ├─ Data loading issue?
  │   └─> Verify file path, check FASTA format
  │
  ├─ Function error?
  │   └─> Check function documentation, see example_usage.R
  │
  ├─ Plot not generated?
  │   └─> Check ggplot2 version, verify output directory exists
  │
  └─ Performance issue?
      └─> Adjust parallel settings, reduce dataset size
```
