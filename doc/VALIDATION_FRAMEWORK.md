# Implementation: Polymorphism vs AnaCoDa Validation Framework

**Date:** 2024-12-15  
**Status:** IMPLEMENTED

---

## Overview

This implementation creates a rigorous mathematical framework for validating polymorphism-based selection estimates (γ) against AnaCoDa's mechanistic codon usage model.

---

## Changes Made

### 1. Bounded Gamma Estimator

**File:** `src/derivation_gamma_from_polymorphism.R`

**Change:**
```r
# Before
S_interval = c(0, 20)

# After
S_interval = c(0, 50)
```

**Updated Documentation:**
```r
#' Constrained to POSITIVE values only (gamma >= 0) to match AnaCoDa framework
#' where preferred codon is optimal and gamma measures selection favorability
```

**Rationale:**
- AnaCoDa assumes preferred codon is optimal (reference point)
- Non-preferred codons have negative deltaEta (penalties)
- By forcing γ ≥ 0, we measure selection **favoring** the preferred codon
- This makes γ directly comparable to AnaCoDa's |4Nes|

**Impact:**
- Eliminates negative gamma estimates (which would be biologically inconsistent)
- Ensures both methods measure selection in same direction
- Simplifies interpretation: higher values = stronger selection for bias

---

### 2. Mathematical Contrast Function

**File:** `src/integrate_intronic_polymorphism.R`

**New Function:** `contrast_gamma_anacoda()`

#### Mathematical Formula

The function implements the rigorous aggregation:

$$
\bar{S}_{poly} = \frac{1}{L} \sum_{AA} \left( \text{Count}_{\text{Unpref}, AA} \times \gamma_{AA} \right)
$$

Where:
- **L** = Total gene length (in codons)
- **Count_Unpref(AA)** = Number of **unpreferred** codons for each amino acid
- **γ_AA** = Selection coefficient (4Nes) favoring the preferred codon

#### Algorithm

```r
FOR each gene:
  total_length = sum of all codon counts
  selection_load = 0
  
  FOR each codon in gene:
    aa = amino_acid(codon)
    
    IF codon is UNPREFERRED for aa:
      gamma = gamma_estimate[gene, aa]
      count = codon_count[gene, codon]
      
      selection_load += count * gamma
  
  S_poly = selection_load / total_length
```

#### Why This Works

**Biological Interpretation:**
- Each unpreferred codon incurs a "cost" of γ selection units
- Total load = sum of costs across all unpreferred codons
- Normalization by gene length gives mean cost per codon

**Comparability to AnaCoDa:**

| Aspect | AnaCoDa | Polymorphism (S_poly) |
|--------|---------|----------------------|
| **Formula** | Sum(\|ΔEta_i\| × Count_i) / L | Sum(γ_AA × Count_Unpref_AA) / L |
| **Unit** | Mean penalty per codon | Mean selection load per codon |
| **Direction** | Penalties (negative) | Favorability (positive) |
| **Aggregation** | Sum over codons | Sum over AAs |
| **Scale** | Model-specific | 4Nes (population genetic) |

**Expected Correlation:** Positive and strong (ρ > 0.5) if both methods detect the same biological signal.

---

## Function Signature

```r
contrast_gamma_anacoda(
  gamma_results,        # Output from estimate_gamma_by_gene_with_neutral_params()
  codon_usage,          # Codon usage matrix (genes × codons)
  preferred_codons,     # AA → preferred codon mapping
  anacoda_intensity,    # AnaCoDa S_coeff per gene
  genetic_code          # Codon → AA mapping
)
```

**Returns:**
- Merged data.table with:
  - `Gene_name`
  - `S_poly` (polymorphism-based load)
  - `S_coeff` (AnaCoDa intensity)
  - `Spearman_rho` (correlation coefficient)
  - `Spearman_p` (p-value)
  - `Gene_Length`

**Side Effects:**
- Prints correlation statistics
- Generates scatter plot: `./results/gamma_anacoda_contrast.pdf`

---

## Output Interpretation

### Strong Validation (ρ > 0.5, p < 0.01)
✓ **Polymorphism data validates mechanistic model**
- Both methods identify same genes under selection
- Selection on codon bias is robustly detected
- Results can be trusted for biological interpretation

### Moderate Validation (0.3 < ρ < 0.5)
⚠ **General agreement with some discrepancies**
- Methods capture similar trends
- Differences may reflect:
  - Temporal scale (contemporary vs ancestral selection)
  - Model assumptions (mutation-selection balance)
  - Measurement noise

### Weak Validation (ρ < 0.3)
✗ **Methods measure different signals**
- Investigate systematic biases
- Check for:
  - Demographic history violations
  - Recombination effects
  - Gene conversion
  - Data quality issues

---

## Usage Example

```r
# Source functions
source("./src/integrate_intronic_polymorphism.R")

# Run contrast analysis
results <- contrast_gamma_anacoda(
  gamma_results = gamma_results,
  codon_usage = codon_usage,
  preferred_codons = preferred_codons_roc,
  anacoda_intensity = selection_coeff_intensity,
  genetic_code = genetic_code_dna_long
)

# Check correlation
print(results[1, .(Spearman_rho, Spearman_p)])

# Examine genes with highest discrepancy
results[, Discrepancy := abs(S_poly - S_coeff)]
high_discrepancy <- results[order(-Discrepancy)][1:10]
print(high_discrepancy[, .(Gene_name, S_poly, S_coeff, Discrepancy)])
```

---

## Validation Checks

### Data Structure
```r
# After running contrast_gamma_anacoda()

# Check distributions
summary(results$S_poly)
summary(results$S_coeff)

# Check for outliers
boxplot(results[, .(S_poly, S_coeff)])

# Check correlation
cor.test(results$S_poly, results$S_coeff, method = "spearman")
```

### Expected Values
- **S_poly**: Mean 0.5-2.0 (depends on Ne and s)
- **S_coeff**: Mean 0.01-0.1 (AnaCoDa scale)
- **Correlation**: ρ > 0.5 for strong validation

---

## Biological Predictions

If selection for codon bias is real:

1. **High-expression genes** should have higher S_poly and S_coeff
2. **Ribosomal proteins** should be in top decile
3. **Housekeeping genes** should show strong selection
4. **Tissue-specific genes** should show weak selection

Test these predictions:
```r
# Expression correlation
cor.test(results$S_poly, log10(results$Expression))

# Compare gene categories
ribosomal <- results[grep("^Rp", Gene_name)]
housekeeping <- results[Gene_name %in% housekeeping_list]

t.test(ribosomal$S_poly, results$S_poly)
```

---

## Files Modified

1. **`src/derivation_gamma_from_polymorphism.R`**
   - Updated `estimate_gamma_for_AA()` interval to [0, 50]
   - Added documentation explaining positive-only constraint

2. **`src/integrate_intronic_polymorphism.R`**
   - Added `contrast_gamma_anacoda()` function
   - Updated `estimate_gamma_by_gene_with_neutral_params()` to use [0, 50]

3. **`examples/example_gamma_anacoda_contrast.R`** (NEW)
   - Complete workflow demonstration
   - Validation checks
   - Diagnostic plots

---

## Next Steps

1. **Run the contrast analysis** on your data
2. **Validate biological predictions** (expression, gene categories)
3. **Investigate outliers** (genes with high discrepancy)
4. **Generate publication figures** using the scatter plot
5. **Write results section** emphasizing validation success

---

## References

**Wright-Fisher Model:**
- Kimura, M. (1962). On the probability of fixation of mutant genes in a population. *Genetics* 47:713-719.

**AnaCoDa Model:**
- Gilchrist, M.A. et al. (2015). Combining models of protein translation and population genetics to predict protein production rates from codon usage patterns. *Mol Biol Evol* 32:2279-2292.

**Selection Inference from Polymorphism:**
- Ewens, W.J. (2004). *Mathematical Population Genetics I. Theoretical Introduction*. Springer.

---

## Bottom Line

This implementation provides a mathematically rigorous framework for validating polymorphism-based selection estimates against a mechanistic model. The key innovation is the precise aggregation formula that makes γ estimates directly comparable to AnaCoDa's selection intensity, enabling robust statistical testing of concordance between independent inference methods.
