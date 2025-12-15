# Polymorphism-Based Selection Inference for Codon Usage Bias

## Overview

This document describes the methodology for estimating selection coefficients (γ = 4Nes) on codon usage using polymorphism data from intronic and coding sequences. The approach leverages neutral intronic sites to empirically estimate mutation rate parameters, eliminating the need for *a priori* assumptions about mutation-drift equilibrium.

---

## Theoretical Framework

### Wright-Fisher Diffusion Model

Under the Wright-Fisher model with mutation and selection, the equilibrium frequency distribution of a biallelic locus follows:

```
φ(x) ∝ x^(α-1) · (1-x)^(β-1) · exp(S·x)
```

Where:
- **x**: Frequency of the preferred allele
- **α = 4N·u**: Scaled mutation rate (unpreferred → preferred)
- **β = 4N·v**: Scaled mutation rate (preferred → unpreferred)
- **S = 4N·s**: Scaled selection coefficient (gamma, γ)
- **N**: Effective population size
- **u, v**: Per-generation mutation rates
- **s**: Selection coefficient per generation

### Site Frequency Spectrum (SFS)

The probability of observing **k** preferred alleles in a sample of **n** individuals is:

```
P(k | n, α, β, S) = (n choose k) · [Beta(k+α, n-k+β) / Beta(α, β)] · [₁F₁(k+α, n+α+β, S) / ₁F₁(α, α+β, S)]
```

Where ₁F₁ is the confluent hypergeometric function.

---

## Key Innovation: Using Intronic Sites

### Why Introns?

1. **Selectively Neutral**: Intronic sites (excluding splice boundaries) evolve under mutation-drift equilibrium (S = 0)
2. **Same Mutation Process**: G↔A and C↔T transitions share the same molecular mechanisms in introns and coding sequences
3. **Independent Estimation**: Allows empirical calibration of α and β without circular reasoning

### Nucleotide-Specific Parameters

Since all preferred codons in your analysis end in **G** or **C**, we estimate:

- **α_G, β_G**: From intronic sites where we track G vs non-G
- **α_C, β_C**: From intronic sites where we track C vs non-C

These parameters capture:
- Population size (N)
- Mutation rates (u, v)
- Demographic history
- Any GC-biased gene conversion

---

## Computational Pipeline

### Phase 1: Extract Intronic Site Frequency Spectra

**Python Script**: `python_scripts/filter_vcf_for_introns.py`

**Purpose**: Parse VCF and GFF3 to identify intronic variants, accounting for strand orientation.

**Key Steps**:

1. **Parse GFF3**:
   - Extract gene boundaries and exon positions
   - Calculate introns as gaps between exons
   - Trim 30bp from exon boundaries (remove splice sites)

2. **Filter VCF**:
   - Keep only variants within trimmed intronic regions
   - Correct for strand orientation (complement on minus strand)
   - Count G vs non-G and C vs non-C separately

3. **Output**:
   - `sfs_introns_G.csv`: Columns (n, k, count) for G-process
   - `sfs_introns_C.csv`: Columns (n, k, count) for C-process

**Usage**:
```bash
# Streaming mode (recommended for large VCFs)
zcat your_variants.vcf.gz | python python_scripts/filter_vcf_for_introns.py \
  --stream \
  --gff data/Mguttatusvar_IM767_887_v2.1.gene.gff3 \
  --workers 8
```

---

### Phase 2: Estimate Neutral Mutation Parameters

**R Function**: `solve_alpha_and_beta_from_introns()`

**Method**: Beta-Binomial Maximum Likelihood

For neutral sites (S = 0), the Wright distribution simplifies to Beta(α, β). We fit this model to the intronic SFS using:

```
L(α, β | data) = ∏ (n choose k) · Beta(k+α, n-k+β) / Beta(α, β)
```

**Implementation**:
```r
neutral_params <- load_and_estimate_neutral_params(
  sfs_G_file = "data/sfs_introns_G.csv",
  sfs_C_file = "data/sfs_introns_C.csv"
)

# Returns:
# - alpha_G: 4N·u for G sites
# - beta_G: 4N·v for G sites
# - alpha_C: 4N·u for C sites
# - beta_C: 4N·v for C sites
```

**Validation**: Calculate expected nucleotide diversity (π) and compare to observed:

```
E[π] = ∫₀¹ 2x(1-x) · φ(x) dx / ∫₀¹ φ(x) dx
```

For S=0, this has an analytical solution using hypergeometric functions.

---

### Phase 3: Infer Selection at Coding Sites

**R Function**: `estimate_gamma_for_AA()`

**Input**:
- Codon variant frequencies at synonymous sites
- Pre-computed α and β from introns
- Terminal nucleotide of preferred codon (G or C)

**Process**:

1. For each gene and amino acid family:
   - Extract counts of preferred vs non-preferred codons
   - Determine if preferred codon ends in G or C
   - Select appropriate (α, β) parameters

2. Optimize gamma by maximizing likelihood:
   ```r
   NLL(γ) = -∑ log[P(k_i | n_i, α, β, γ)]
   ```

3. Statistical significance:
   - Likelihood ratio test: 2·ln(L_selection / L_neutral) ~ χ²(1)
   - Threshold: |γ| > 1.92 ≈ p < 0.05

**Implementation**:
```r
gamma_results <- estimate_gamma_by_gene_with_neutral_params(
  codon_vcf_data = codon_polymorphism_data,
  neutral_params = neutral_params,
  preferred_codons_df = preferred_codons_table
)
```

---

## Biological Interpretation

### Sign of Gamma

| γ Value | Interpretation |
|---------|---------------|
| γ > 0   | **Positive selection**: Preferred codons increase fitness (translational efficiency, accuracy, mRNA stability) |
| γ ≈ 0   | **Neutral evolution**: No selection on codon choice |
| γ < 0   | **Purifying selection**: Preferred codons are actually deleterious (rare but possible if "preferred" definition is wrong) |

### Magnitude of Gamma

| \|γ\| | Effect Size |
|------|------------|
| < 1  | Weak selection (drift dominates) |
| 1-5  | Moderate selection (mutation-selection balance) |
| > 5  | Strong selection (rapid fixation of preferred alleles) |

### Expected Patterns

If codon usage bias is adaptive:
1. **High-expression genes**: γ > 0 (selection for translational efficiency)
2. **Low-expression genes**: γ ≈ 0 (relaxed selection)
3. **Correlation**: γ should correlate positively with CAI, CDC, and expression level

---

## Advantages Over Traditional Methods

### Compared to CAI/ENC

| Method | Strength | Limitation |
|--------|----------|------------|
| **CAI/ENC** | Simple, no population data needed | Assumes selection, can't measure strength |
| **Polymorphism γ** | Directly estimates selection strength | Requires high-quality VCF and annotation |

### Compared to dN/dS

| Method | Scope | Time Scale |
|--------|-------|-----------|
| **dN/dS** | Measures protein evolution | Between-species divergence |
| **γ estimation** | Measures codon preference | Within-population dynamics |

### Compared to McDonald-Kreitman

| Method | Target | Interpretation |
|--------|--------|---------------|
| **MK test** | Fixed vs polymorphic ratio | Detects selection but not magnitude |
| **γ estimation** | Full SFS shape | Quantifies selection coefficient |

---

## Validation Strategy

### 1. Internal Consistency
- Compare α_G vs α_C (should be similar if mutation process is uniform)
- Verify E[π] from model matches observed π in introns

### 2. Cross-Validation with Expression
```r
cor.test(gamma_per_gene, expression_level, method = "spearman")
```
**Expectation**: ρ > 0.3 if selection is expression-driven

### 3. Comparison with CAI
```r
cor.test(gamma_per_gene, CAI, method = "spearman")
```
**Expectation**: ρ > 0.4 (both measure adaptive codon usage)

### 4. Gene Ontology Enrichment
- High-γ genes should be enriched for:
  - Ribosomal proteins
  - Translation machinery
  - Highly abundant enzymes

---

## Assumptions and Limitations

### Assumptions

1. **Wright-Fisher Model Applies**
   - Random mating
   - Constant population size (or smooth changes)
   - Weak selection (|s| << 1)

2. **Introns Are Neutral**
   - No regulatory elements
   - No selection on RNA secondary structure
   - Splice sites properly excluded

3. **Mutation Process Is Consistent**
   - Same u, v in introns and coding sequences
   - No codon-context effects (CpG islands, etc.)

4. **Synonymous Sites Are Biallelic**
   - Preferred vs non-preferred (valid for your data)
   - No complex multi-allelic dynamics

### Limitations

1. **Demographic History**
   - α, β absorb population size changes
   - Recent bottlenecks can mimic selection

2. **Gene Conversion**
   - GC-biased gene conversion inflates α for G/C
   - Can be detected if α_G, α_C >> α_A, α_T

3. **Sample Size**
   - Need sufficient polymorphism (π > 0.01)
   - Low diversity reduces power

4. **Linkage**
   - Closely linked sites violate independence
   - Use gene-level aggregation to minimize

---

## Output Files

### Primary Results

1. **`neutral_mutation_parameters.csv`**
   - α_G, β_G, α_C, β_C estimates
   - Expected π for validation

2. **`gamma_estimates_by_gene_and_aa.csv`**
   - Columns: Gene, AA, Gamma, N_Sites, Significant, Terminal_Nuc
   - One row per gene × amino acid combination

3. **`gamma_validation_vs_cai_cdc.csv`**
   - Merged gamma, CAI, CDC, expression data
   - For cross-validation analyses

### Diagnostic Plots

1. **`gamma_distribution_by_nucleotide.pdf`**
   - Histogram of γ values for G-ending vs C-ending codons

2. **`gamma_vs_expression.pdf`**
   - Scatter plot: Mean γ per gene vs expression level

3. **`gamma_vs_CAI.pdf`**
   - Scatter plot: Mean γ per gene vs CAI

---

## Troubleshooting

### Low Convergence Rate

**Problem**: Many genes return NA for gamma
**Solutions**:
- Increase S_interval range (try `c(-20, 50)`)
- Aggregate across genes (family-level analysis)
- Filter out low-diversity amino acids

### Unrealistic Parameter Estimates

**Problem**: α or β > 1 (implies θ > 1, very high mutation rate)
**Checks**:
- Verify VCF filtering (are invariant sites included?)
- Check for reference bias in genotype calls
- Inspect GFF parsing (are introns correctly identified?)

### Negative Gamma Everywhere

**Problem**: All γ < 0 suggests model misspecification
**Diagnoses**:
- Preferred codon definition may be inverted
- Check strand orientation in VCF processing
- Verify ancestral state inference

---

## Future Enhancements

1. **Multiallelic Extension**
   - Dirichlet-Multinomial framework
   - Estimate γ for each codon separately

2. **Codon Context Effects**
   - Model CpG dinucleotides separately
   - Account for tRNA abundance variation

3. **Hierarchical Bayesian Model**
   - Share information across amino acids
   - Regularize estimates for rare codons

4. **Time-Varying Selection**
   - Integrate with ancient DNA
   - Detect recent adaptation

---

## References

### Theoretical Foundation

- **Wright, S. (1938)**. "The distribution of gene frequencies under irreversible mutation." *PNAS* 24(7): 253-259.

- **Sawyer, S.A. & Hartl, D.L. (1992)**. "Population genetics of polymorphism and divergence." *Genetics* 132: 1161-1176.

### Codon Usage Bias

- **Bulmer, M. (1991)**. "The selection-mutation-drift theory of synonymous codon usage." *Genetics* 129: 897-907.

- **Akashi, H. (1995)**. "Inferring weak selection from patterns of polymorphism and divergence at 'silent' sites in *Drosophila* DNA." *Genetics* 139: 1067-1076.

### Hypergeometric Functions in SFS

- **Steinrücken, M., et al. (2014)**. "A novel spectral method for inferring general diploid selection from time series genetic data." *Ann. Appl. Stat.* 8(4): 2203-2222.

- **Tataru, P., et al. (2017)**. "Statistical inference in the Wright–Fisher model using allele frequency data." *Syst. Biol.* 66(1): e30-e46.

---

## Contact

For questions about this methodology:
- **Luis Javier Madrigal-Roca** (Implementation)
- **John K. Kelly** (Theoretical framework)

*University of Kansas, Department of Ecology and Evolutionary Biology*
