# Codon Usage Bias Analysis Methods

## Overview

This document describes the methodology and theoretical background for each analysis implemented in the CUB pipeline.

## 1. Codon Quantification

### Method
Counts the occurrence of each of the 64 possible codons across all genes in the transcriptome.

### Implementation
- Reads coding sequences from FASTA file
- Filters sequences with canonical start codon (ATG)
- Validates reading frame (length divisible by 3)
- Splits sequences into codons
- Counts frequency of each codon per gene

### Output
Matrix with genes as rows and 64 codon columns, plus gene names.

---

## 2. Relative Synonymous Codon Usage (RSCU)

### Theory
RSCU measures the relative frequency of a codon compared to the expected frequency if all synonymous codons were used equally. It normalizes for amino acid composition and genetic code degeneracy.

### Formula
```
RSCU = (observed frequency of codon) / (expected frequency)
     = (Xij / Σj Xij) × ni

Where:
- Xij = count of codon j for amino acid i
- ni = number of synonymous codons for amino acid i
```

### Interpretation
- **RSCU = 1**: Codon used at expected frequency (no bias)
- **RSCU > 1**: Positive codon bias (overused)
- **RSCU < 1**: Negative codon bias (underused)

### Reference
Sharp, P.M., & Li, W.H. (1987). The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research, 15(3), 1281-1295.

---

## 3. Effective Number of Codons (ENC)

### Theory
ENC quantifies the degree of codon bias in a gene, independent of gene length and amino acid composition. It represents the effective number of codons used if they were used with equal probability.

### Method (Wright 1990)
1. Calculate homozygosity F for each degeneracy class (2-fold, 3-fold, 4-fold, 6-fold)
2. Compute ENC from F values

### Formula
```
ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6

Where Fk = average homozygosity for k-fold degenerate amino acids
Fk = Σ(n × p² - 1) / (n - 1)

n = total count of amino acid
p = proportion of each synonymous codon
```

### Interpretation
- **ENC = 20**: Extreme bias (one codon per amino acid)
- **ENC = 61**: No bias (all synonymous codons used equally)
- **20 < ENC < 40**: High codon bias
- **40 < ENC < 50**: Moderate codon bias
- **50 < ENC < 61**: Low codon bias

### Reference
Wright, F. (1990). The 'effective number of codons' used in a gene. Gene, 87(1), 23-29.

---

## 4. GC Content Metrics

### Metrics Calculated
- **GC**: Overall GC content across all three codon positions
- **GC1**: GC content at first codon position
- **GC2**: GC content at second codon position
- **GC3**: GC content at third codon position
- **GC12**: Combined GC content at first and second positions
- **GC3s**: GC content at synonymous third positions (excludes Met and Trp)

### Biological Significance
- **GC12**: Primarily constrained by amino acid sequence (under selection)
- **GC3**: Mostly synonymous changes (may reflect mutational bias)
- **GC3s**: Best indicator of mutational pressure vs selection

---

## 5. Neutrality Plot (GC12 vs GC3)

### Theory
The neutrality plot helps distinguish between mutation pressure and natural selection in shaping codon usage.

### Method
- Plot GC12 (first two positions) against GC3 (third position)
- Calculate linear regression and correlation
- Compare slope to theoretical expectations

### Interpretation

**Slope ≈ 1 (strong positive correlation)**
- Mutation pressure dominates
- GC content determined primarily by mutational bias
- Little selective constraint on synonymous sites

**Slope ≈ 0 (weak or no correlation)**
- Selection dominates
- GC3 independent of GC12
- Strong selective pressure on synonymous codon usage

**Intermediate slope (0.3-0.7)**
- Both mutation and selection contribute
- Balance between mutational bias and selection

### Reference
Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. Proceedings of the National Academy of Sciences, 85(8), 2653-2657.

---

## 6. ENC Plot (ENC vs GC3s)

### Theory
The ENC plot identifies genes under selection for codon usage bias by comparing observed ENC to expected ENC under mutation-drift equilibrium.

### Expected ENC Curve (Wright 1990)
```
ENC_expected = 2 + GC3s + (29 / (GC3s² + (1-GC3s)²))
```

This curve represents the expected relationship between ENC and GC3s if codon usage were determined solely by mutational bias.

### Interpretation

**Genes on the expected curve**
- Codon usage determined by mutation-drift equilibrium
- No selection for codon bias

**Genes below the curve**
- Lower ENC than expected (stronger bias)
- Under selection for codon usage bias
- May indicate translational selection

**Genes above the curve**
- Higher ENC than expected (weaker bias)
- May reflect unusual amino acid composition
- Other factors affecting codon usage

### Reference
Wright, F. (1990). The 'effective number of codons' used in a gene. Gene, 87(1), 23-29.

---

## 7. PR2 Bias Plot (Parity Rule 2)

### Theory
PR2 analyzes the balance between purines (A, G) and pyrimidines (T, C) at the third codon position, revealing biases in mutation and selection patterns.

### Method
Calculate two ratios for each gene:
```
AU3 = A3 / (A3 + T3)  # A vs T at 3rd position
GC3 = G3 / (G3 + C3)  # G vs C at 3rd position
```

Plot AU3 vs GC3.

### Interpretation

**Center point (0.5, 0.5)**
- No bias (Parity Rule 2 holds)
- Equal usage of complementary bases

**Deviation from center**
- Indicates bias in base composition
- Reveals patterns of mutation pressure

**Diagonal patterns**
- AT bias: movement toward (1, 0) or (0, 1)
- GC bias: movement toward (0, 0) or (1, 1)

**Quadrant interpretation**
- Upper right: G and A preferred (purine rich)
- Lower left: C and T preferred (pyrimidine rich)
- Upper left or lower right: Strand asymmetry

### Reference
Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. Journal of Molecular Evolution, 40(3), 318-325.

---

## 8. Codon Logo Visualization

### Theory
Sequence logos visually represent the information content at each position within codons for a specific amino acid.

### Method
For each amino acid:
1. Collect all synonymous codons and their frequencies
2. Calculate nucleotide frequency at each of the three positions
3. Visualize as stacked bar chart with color-coded nucleotides

### Color Scheme (Standard)
- **Adenine (A)**: Green
- **Thymine (T)**: Red
- **Guanine (G)**: Orange
- **Cytosine (C)**: Blue

### Interpretation
- Height indicates relative usage frequency
- Allows quick visual identification of preferred nucleotides at each position
- Facilitates comparison across amino acids

---

## Integration of Multiple Analyses

### Comprehensive CUB Assessment
Combining multiple metrics provides a complete picture:

1. **RSCU**: Identifies specific preferred codons
2. **ENC**: Quantifies overall bias strength
3. **GC metrics**: Reveals compositional patterns
4. **Neutrality plot**: Distinguishes mutation vs selection
5. **ENC plot**: Identifies genes under selection
6. **PR2 plot**: Reveals strand asymmetry and mutational patterns
7. **Codon logos**: Visualizes position-specific preferences

### Biological Conclusions
By integrating these analyses, one can:
- Determine if codon bias is driven by mutation or selection
- Identify highly expressed genes (often show strong bias)
- Detect recent changes in mutational patterns
- Understand evolutionary forces shaping the genome
- Optimize codons for heterologous expression

---

## Software Implementation Details

### Parallel Processing
- Uses `doFuture` for parallel computation
- Significant speedup for large transcriptomes (>10,000 genes)
- Can be disabled for debugging

### Data Structures
- `data.table` for efficient data manipulation
- `Biostrings` for sequence handling
- `ggplot2` for publication-quality graphics

### Quality Control
- Filters non-canonical start codons
- Validates reading frames
- Handles missing data gracefully
- Provides informative error messages

---

## References

1. Sharp, P.M., & Li, W.H. (1987). The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295.

2. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29.

3. Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *Proceedings of the National Academy of Sciences*, 85(8), 2653-2657.

4. Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325.

5. Plotkin, J.B., & Kudla, G. (2011). Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12(1), 32-42.

6. Novembre, J.A. (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390-1394.
