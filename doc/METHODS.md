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

## 8. Codon Deviation Coefficient (CDC)

### Theory
CDC quantifies the deviation of observed codon frequencies from expected frequencies based on amino acid composition. It provides a complementary measure to ENC for detecting genes under selection for codon usage.

### Method
1. Calculate expected codon frequencies based on amino acid composition
2. Compare observed vs. expected using chi-square or G-test
3. Compute deviation coefficient
4. Test statistical significance (p-value)
5. Apply FDR correction for multiple testing

### Statistical Framework
```
CDC = sum((observed - expected)^2 / expected)
```

Under the null hypothesis of no selection, CDC follows a chi-square distribution with degrees of freedom equal to the number of synonymous codon families minus constraints.

### Interpretation
- **CDC ≈ 0**: Observed matches expected (no selection)
- **CDC > 0**: Deviation from expected (potential selection)
- **p < 0.05**: Statistically significant deviation
- Genes with significant CDC below the ENC curve indicate selection for translational efficiency

### Integration with ENC Plot
CDC analysis enhances the ENC plot by providing statistical rigor:
- Genes significantly below the Wright curve have low ENC and significant CDC
- This combination strongly suggests selection for codon bias

---

## 9. Codon Adaptation Index (CAI)

### Theory
CAI measures the degree to which a gene's codon usage matches the codon usage of a reference set of highly expressed genes. It assumes that highly expressed genes have optimized codon usage for translational efficiency.

### Method (Sharp & Li 1987)
1. Define reference set (e.g., top 5% expressed genes)
2. Calculate relative adaptiveness (w) for each codon:
   ```
   w_i = (frequency of codon i in reference set) / 
         (frequency of most common codon for that amino acid)
   ```
3. For each gene, calculate geometric mean of w values:
   ```
   CAI = exp(sum(log(w_i)) / L)
   ```
   where L is the number of codons in the gene

### Reference Set Selection
Critical for accurate CAI calculation:
- Should include genes under strong selection for translation
- Common choices: ribosomal proteins, elongation factors, highly expressed genes
- In this pipeline: top 5% expressed genes with relevant functional annotations

### Interpretation
- **CAI = 0-1**: Ranges from no adaptation (0) to perfect adaptation (1)
- **CAI > 0.7**: Highly adapted (typical for housekeeping genes)
- **CAI = 0.5-0.7**: Moderately adapted
- **CAI < 0.5**: Poorly adapted (typical for lowly expressed genes)

### Validation
- Strong positive correlation with expression level
- Higher CAI in genes encoding ribosomal/translation machinery
- Lower CAI in tissue-specific or conditionally expressed genes

### Reference
Sharp, P.M., & Li, W.H. (1987). The codon adaptation index-a measure of directional synonymous codon usage bias. Nucleic Acids Research, 15(3), 1281-1295.

---

## 10. Polymorphism-Based Selection Inference

### Theory
The "hump effect" predicts that nucleotide diversity (π) peaks at sites under weak selection (0.1 < Ns < 1) due to the balance between selection and drift. This pattern can be used to estimate selection coefficients.

### Method
1. Process VCF files to extract polymorphism data
2. Classify sites by codon preference (preferred vs. non-preferred)
3. Calculate nucleotide diversity (π) per site
4. Estimate local mutation rate (M) from neutral regions
5. Apply selection inference model to estimate Ns
6. Validate with site frequency spectrum (SFS) analysis

### Theoretical Framework
Under weak selection, the expected diversity is:
```
π_selected ≈ 4Nμ × h(Ns)
```
where h(Ns) is a function that peaks at Ns ≈ 1.

For neutral sites:
```
π_neutral = 4Nμ
```

The ratio π_selected / π_neutral reveals the strength of selection.

### Mutation Rate Estimation
Critical for accurate selection inference:
- Use introns or intergenic regions (putatively neutral)
- Calculate mutation rate (M) from sequence divergence or polymorphism
- Account for local variation in mutation rate
- Options: intron-based (preferred) or intergenic-based

### Site Classification
- **Preferred codons**: Identified via CAI (w = 1) or statistical enrichment
- **Non-preferred codons**: All other synonymous codons
- **Unpreferred codons**: Significantly depleted in highly expressed genes

### Interpretation
- **π_preferred < π_non-preferred**: Selection against non-preferred codons
- **Peak diversity at intermediate Ns**: Validates weak selection model
- **Low diversity at preferred sites**: Strong effective selection (Ns >> 1)
- **SFS skew toward rare alleles**: Indicates purifying selection

### Validation
- Compare observed vs. expected SFS under neutral model
- Test for correlation between selection coefficient and expression
- Verify hump-shaped relationship between π and Ns

### References
- Kimura, M. (1983). The Neutral Theory of Molecular Evolution. Cambridge University Press.
- Bulmer, M. (1991). The selection-mutation-drift theory of synonymous codon usage. Genetics, 129(3), 897-907.

---

## 11. AnaCoDa: Bayesian Modeling of Codon Usage

### Theory
AnaCoDa (Analyzing Codon Data) implements a Bayesian framework to simultaneously infer:
- Selection coefficients per codon (dEta)
- Gene expression levels (phi)
- Mutation bias parameters (dM)

The ROC (Ribosome Overhead Cost) model assumes codon fitness depends on ribosome availability and elongation rate.

### ROC Model
The probability of observing codon i in gene g is:
```
P(codon_i | phi_g) = exp(-dM_i - dEta_i × phi_g) / Z
```

where:
- dM_i: mutation bias for codon i
- dEta_i: selection coefficient for codon i
- phi_g: expression level of gene g
- Z: normalization constant

### Method
1. Prepare input files:
   - CDS sequences (FASTA)
   - Expression data (optional)
   - Mutation rate matrix (dM, from neutral regions)
2. Run MCMC sampling (command-line execution):
   ```bash
   Rscript AnaCoDa_pipeline.R \
     -i sequences.fa \
     -o output_dir \
     --fix_dM --dM mutation_rates.csv \
     -s 10000 -d 4000 -n 10
   ```
3. Check convergence (Gelman-Rubin diagnostic)
4. Extract posterior estimates
5. Validate model predictions

### MCMC Configuration
- **Sampling iterations**: 10,000+ (more for complex models)
- **Burn-in**: 4,000-5,000 iterations
- **Thinning**: Every 10th or 25th sample
- **Chains**: 3-6 independent chains for convergence assessment
- **Adaptation**: First 20-30% of samples

### Convergence Diagnostics
Gelman-Rubin statistic (R̂) for each parameter:
- **R̂ < 1.1**: Good convergence
- **R̂ = 1.1-1.2**: Acceptable convergence
- **R̂ > 1.2**: Poor convergence, run longer or adjust priors

Additional checks:
- Visual inspection of trace plots
- Effective sample size (ESS > 100)
- Autocorrelation < 0.1 at lag 10

### Model Variants
1. **Naive model**: Estimate all parameters (dM, dEta, phi)
2. **Fixed dM**: Use mutation rates from neutral regions (recommended)
3. **Fixed dM + empirical phi**: Use observed expression data
4. **Multi-tissue**: Model tissue-specific selection

### Interpretation
- **dEta < 0**: Codon is preferred (lower cost)
- **dEta > 0**: Codon is disfavored (higher cost)
- **|dEta| > 0.5**: Strong selection
- **|dEta| < 0.1**: Weak or no selection

### Validation
Compare model predictions with empirical patterns:
1. **Codon frequency trajectories**: Plot codon frequency vs. expression
2. **Expected vs. observed**: Model should reproduce observed patterns
3. **Preferred codon identification**: Agreement with CAI-based methods
4. **Expression estimates**: Correlation with RNA-seq data (r > 0.7)

### Advantages over Traditional Methods
- Integrates mutation and selection simultaneously
- Accounts for expression level
- Provides uncertainty estimates (credible intervals)
- Can incorporate prior biological knowledge
- Models multiple tissues/conditions

### Limitations
- Computationally intensive (hours to days)
- Requires careful convergence checking
- Sensitive to prior specifications
- Assumes steady-state evolution

### Software
AnaCoDa R package: https://github.com/clandere/AnaCoDa

---

## 12. Statistical Testing Framework

### Per-Codon Proportion Tests
Test whether codon proportions differ between groups (e.g., high vs. low expression).

**Method**: Two-sample test for proportions
```
H0: p_high = p_low (no difference in codon usage)
H1: p_high ≠ p_low (significant difference)
```

**Test statistic**:
```
z = (p1 - p2) / sqrt(p_pooled × (1 - p_pooled) × (1/n1 + 1/n2))
```

**Multiple testing correction**: Benjamini-Hochberg FDR < 0.05

### Family-Level Multinomial Tests
Test overall codon usage patterns within amino acid families.

**Method**: Multinomial logistic regression or chi-square test
- Tests simultaneous differences across all codons in a family
- More powerful than multiple pairwise tests
- Accounts for compositional constraint (proportions sum to 1)

### Effect Size Calculation
Cohen's d for pairwise comparisons:
```
d = (mean1 - mean2) / pooled_SD
```

**Interpretation**:
- |d| < 0.2: Negligible
- 0.2 ≤ |d| < 0.5: Small
- 0.5 ≤ |d| < 0.8: Medium
- |d| ≥ 0.8: Large

---

## 13. Multivariate Analysis of Codon Usage

### Principal Component Analysis (PCA)
Reduces dimensionality of codon usage data (59 dimensions) to principal components.

**Method**:
1. Create codon usage matrix (genes × codons)
2. Center and scale data
3. Compute covariance matrix
4. Extract eigenvectors (principal components)
5. Project genes onto PC space

**Interpretation**:
- PC1 often correlates with GC content or expression
- PC2 may reflect selection vs. mutation
- Genes cluster by expression level or functional category

### Correspondence Analysis (CA)
Alternative to PCA, optimized for compositional data (proportions).

**Advantages over PCA**:
- Preserves distances between codons and genes simultaneously
- Creates symmetric maps (biplot)
- No negative values in transformed space
- Better suited for count data

### Enhanced Biplots
Visualize both genes (points) and codons (vectors) simultaneously:
- **Gene points**: Position reflects overall codon usage pattern
- **Codon vectors**: Direction indicates codon contribution to PCs
- **Vector length**: Importance of codon in explaining variance
- **Color coding**: By expression, GC3s, or other attributes

**Biological insights**:
- Preferred codons point toward highly expressed genes
- GC-rich codons point toward high-GC genes
- Clustering reveals functional gene groups

---

## Integration of Multiple Analyses

### Comprehensive CUB Assessment
Combining multiple metrics provides a complete picture:

1. **RSCU + ENC**: Identify strength and specifics of bias
2. **CAI + Expression**: Validate translational selection hypothesis
3. **CDC + ENC plot**: Statistical rigor for selection detection
4. **GC metrics + Neutrality plot**: Distinguish mutation vs. selection
5. **Polymorphism + AnaCoDa**: Validate selection coefficients
6. **Multivariate analysis**: Reveal global patterns and gene clusters

### Biological Conclusions Framework

**Evidence for translational selection**:
✓ High CAI in highly expressed genes
✓ Genes below ENC curve with significant CDC
✓ Lower diversity at preferred codon sites
✓ Positive dEta values from AnaCoDa
✓ Codon frequency increases with expression

**Evidence for mutational bias**:
✓ Strong GC12 vs. GC3 correlation (neutrality plot)
✓ Genes on ENC curve
✓ Similar diversity at preferred and non-preferred sites
✓ dM dominates over dEta in AnaCoDa

**Evidence for balanced forces**:
✓ Intermediate slopes in neutrality plot
✓ Some genes below, some on ENC curve
✓ Hump-shaped π vs. Ns relationship
✓ Both dM and dEta contribute significantly

---

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
- Automatically detects available cores

### Data Structures
- `data.table` for efficient data manipulation (100x faster than data.frame)
- `Biostrings` for sequence handling and validation
- `ggplot2` for publication-quality graphics
- Native R matrices for numerical operations

### Quality Control
- Filters non-canonical start codons
- Validates reading frames (length divisible by 3)
- Removes sequences with internal stop codons
- Handles missing data gracefully
- Provides informative error messages
- Extensive unit testing for core functions

### Performance Optimization
- Vectorized operations where possible
- Memory-efficient algorithms for large datasets
- Incremental saving of intermediate results
- Progress bars for long-running operations

---

## References

### Core Methods
1. Sharp, P.M., & Li, W.H. (1987). The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295.

2. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29. https://doi.org/10.1016/0378-1119(90)90491-9

3. Sun, X., Yang, Q., & Xia, X. (2012). An improved implementation of effective number of codons (Nc). *Molecular Biology and Evolution*, 30(1), 191-196. https://doi.org/10.1093/molbev/mss201

4. Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *Proceedings of the National Academy of Sciences*, 85(8), 2653-2657.

5. Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325.

### Population Genetics & Selection
6. Kimura, M. (1983). *The Neutral Theory of Molecular Evolution*. Cambridge University Press.

7. Bulmer, M. (1991). The selection-mutation-drift theory of synonymous codon usage. *Genetics*, 129(3), 897-907.

8. McVean, G.A., & Charlesworth, B. (1999). A population genetic model for the evolution of synonymous codon usage: patterns and predictions. *Genetics Research*, 74(2), 145-158.

### Bayesian Modeling
9. Gilchrist, M.A., et al. (2015). Estimating gene expression and codon-specific translational efficiencies, mutation biases, and selection coefficients from genomic data alone. *Genome Biology and Evolution*, 7(6), 1559-1579.

10. Landerer, C., et al. (2018). AnaCoDa: analyzing codon data with Bayesian mixture models. *Bioinformatics*, 34(14), 2496-2498.

### General Reviews
11. Plotkin, J.B., & Kudla, G. (2011). Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12(1), 32-42.

12. Novembre, J.A. (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390-1394.

13. Hershberg, R., & Petrov, D.A. (2008). Selection on codon bias. *Annual Review of Genetics*, 42, 287-299.

---

## Appendix: Formulas and Notation

### Key Notation
- **N**: Effective population size
- **s**: Selection coefficient
- **μ**: Mutation rate
- **π**: Nucleotide diversity
- **θ**: Population mutation rate (4Nμ)
- **Ns**: Scaled selection coefficient (determines fixation probability)
- **w**: Relative adaptiveness (CAI)
- **φ (phi)**: Gene expression level
- **dM**: Mutation bias parameter
- **dEta**: Selection coefficient (AnaCoDa)

### Wright's ENC Formula
```
ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
```
where F_k is homozygosity for k-fold degenerate families.

With family splitting (recommended):
```
ENC = 2 + 12/F2 + 1/F3 + 8/F4
```

### Expected ENC Curve
```
ENC_expected = 2 + GC3s + 29/(GC3s² + (1-GC3s)²)
```

### CAI Formula
```
CAI = exp((1/L) × Σ ln(w_i))
```
where L is gene length in codons, w_i is relative adaptiveness of codon i.

### Selection-Diversity Relationship
Under weak selection (Ns ≈ 1):
```
π/π_neutral ≈ h(Ns)
```
where h(Ns) is maximized near Ns = 1 (the "hump effect").

---

## Glossary

**Codon Usage Bias (CUB)**: Non-uniform usage of synonymous codons

**Synonymous codons**: Codons encoding the same amino acid

**Degeneracy**: Number of synonymous codons for an amino acid (2, 3, 4, or 6-fold)

**Preferred codon**: Codon enriched in highly expressed genes

**Optimal codon**: Synonym for preferred codon

**Translational selection**: Natural selection favoring codons for translation efficiency/accuracy

**Mutation bias**: Tendency of mutation process to favor certain nucleotides

**GC3s**: GC content at synonymous third codon positions

**ENC**: Effective Number of Codons (20-61 scale)

**RSCU**: Relative Synonymous Codon Usage (ratio to expected frequency)

**CAI**: Codon Adaptation Index (0-1 scale)

**CDC**: Codon Deviation Coefficient (deviation from expected usage)

**ROC model**: Ribosome Overhead Cost model (AnaCoDa framework)

**dM**: Mutation bias parameter (log-scale)

**dEta**: Selection coefficient parameter (log-scale)

**MCMC**: Markov Chain Monte Carlo (Bayesian sampling method)

**FDR**: False Discovery Rate (multiple testing correction)

**SFS**: Site Frequency Spectrum (distribution of allele frequencies)

**π (pi)**: Nucleotide diversity (average pairwise differences)

**Ns**: Scaled selection coefficient (N × s, determines fate of mutations)

---
