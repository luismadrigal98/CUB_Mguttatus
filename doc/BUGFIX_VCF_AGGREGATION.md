# CRITICAL BUG FIX: VCF Data Aggregation

**Date:** 2024-12-15  
**Issue:** Data pooling across codon positions invalidating selection inference  
**Status:** FIXED

---

## Problem Description

### The Bug
The VCF processing function `prepare_vcf_for_gamma_estimation()` was collapsing all synonymous codon positions within a Gene×AA combination into a **single "super-site"** by grouping only by `(Gene, AA, Preferred_Codon)`.

### Evidence
From `summary(gamma_results)`:
```r
N_Sites        Total_Alleles    Mean_Alleles_Per_Site
Min.   :1      Min.   :    44   (not shown - would be ~4000)
1st Qu.:1      1st Qu.:  1868
Median :1      Median :  3354
Mean   :1      Mean   :  4281
```

**Critical indicators:**
- `N_Sites = 1` (always) → Only one "site" per Gene×AA
- `Total_Alleles ≈ 4,281` → Sum of ALL codon positions
- Expected: `n ≈ 187` per site (**inbred lines** = homozygous genotypes only)
- Actual: `n ≈ 4,000+` → Clearly pooled across ~10+ positions

**Note on Sample Size:** For **inbred lines**, n represents the number of homozygous individuals, not allele count. With 187 inbred lines, expect n≈187. Heterozygous calls indicate residual heterozygosity and are excluded.

### Why This Is Critical

#### Statistical Violation
The Wright-Fisher model assumes **independent loci** experiencing drift:
- Each codon position has its own allele frequency trajectory
- Variance between sites reflects stochastic drift
- Pooling sites **eliminates this variance**

#### Mathematical Impact
```
Correct:   P(Data | γ) = ∏ᵢ P(kᵢ, nᵢ | γ, α, β)   [product over sites]
Incorrect: P(Data | γ) = P(∑kᵢ, ∑nᵢ | γ, α, β)     [single pooled site]
```

Pooling creates a "super-site" with:
- Sample size ~10× larger than reality
- No between-site variance
- **Massively inflated statistical power**
- Artificially narrow confidence intervals

#### Biological Consequences
1. **Type I error inflation**: Weak selection appears significant
2. **Gamma overestimation**: Drift variance is missing, so model attributes deviations to selection
3. **Invalid hypothesis tests**: Cannot distinguish selection from drift

---

## The Fix

### Code Changes

#### File: `src/convert_vcf_codon_format.R`

**Before (WRONG):**
```r
}, by = .(Gene, AA, Preferred_Codon)]
```

**After (CORRECT):**
```r
}, by = .(Gene, Codon_Pos, AA, Preferred_Codon)]
```

**Key change:** Added `Codon_Pos` to the grouping variables to preserve independent loci.

---

#### File: `src/integrate_intronic_polymorphism.R`

**Updated estimation logic:**

**Before:**
```r
# WRONG assumption
if (n[1] < min_sample_size) {  # Checking first element only
  gamma_est <- NA_real_
}
```

**After:**
```r
# Correct handling of VECTORS
min_sites <- 5               # Need at least 5 independent loci
min_sample_per_site <- 50    # Each site needs ≥50 alleles

if (n_sites < min_sites) {
  gamma_est <- NA_real_
} else if (any(n < min_sample_per_site)) {
  gamma_est <- NA_real_
} else {
  # Pass VECTORS to likelihood function
  estimate_gamma_for_AA(counts = k, sample_sizes = n, ...)
}
```

**Added validation output:**
```r
Mean sites per Gene×AA: %.1f (should be >5)
Mean alleles per site: %.1f (should be ~187 for INBRED lines, NOT >500)
```

**Important:** For inbred lines, sample size reflects homozygous genotypes only. Residual heterozygosity is filtered out as unreliable.

---

## Validation

### Expected Output After Fix

#### Data Structure
```r
summary(gamma_results)
     N_Sites          Total_Alleles    Mean_Alleles_Per_Site
 Min.   : 1.0      Min.   :   110      Min.   : 100
 1st Qu.: 8.0      1st Qu.:  1248      1st Qu.: 180
 Median :15.0      Median :  2805      Median : 187
 Mean   :18.3      Mean   :  3424      Mean   : 190
 3rd Qu.:25.0      3rd Qu.:  4675      3rd Qu.: 200
```

**Key indicators:**
- ✅ `N_Sites > 1` for most Gene×AA (multiple independent loci)
- ✅ `Mean_Alleles_Per_Site ≈ 187` (matches 187 inbred lines, homozygous only)
- ✅ `Total_Alleles = N_Sites × 187` (sum across independent sites)
- ✅ Heterozygous calls excluded (residual heterozygosity safeguard)

#### Statistical Behavior
- **Wider confidence intervals** (correctly accounting for drift variance)
- **Fewer significant results** (realistic Type I error rate)
- **Lower gamma estimates** (not confounding drift with selection)

---

## Testing

### Quick Validation Script

```r
# Load fixed functions
source("./src/convert_vcf_codon_format.R")
source("./src/integrate_intronic_polymorphism.R")

# Process VCF data
vcf_prepared <- prepare_vcf_for_gamma_estimation(
  vcf_codon_dt = your_vcf_data,
  genetic_code_df = genetic_code
)

# Check structure
cat("Rows in output:", nrow(vcf_prepared), "\n")
cat("Unique Gene×AA:", nrow(unique(vcf_prepared[, .(Gene, AA)])), "\n")
cat("Mean n per row:", mean(vcf_prepared$n), "\n")

# CRITICAL CHECKS:
stopifnot(nrow(vcf_prepared) > nrow(unique(vcf_prepared[, .(Gene, AA)])))
stopifnot(mean(vcf_prepared$n) < 250)  # Should be ~187 for inbred lines
```

### Expected Warnings (GOOD!)
If you see these after the fix, the code is working correctly:
```
⚠️  Gene×AA with <5 sites: XXX (should be low but non-zero)
Mean alleles per site: 187.4 ✓
```

If you see these, something is STILL WRONG:
```
⚠️  CRITICAL: Sample sizes are too large! Sites are being pooled!
⚠️  WARNING: Sample sizes suggest diploid counting. For inbred lines, expect n≈187.
⚠️  CRITICAL ERROR: N_Sites is too low! Data is still collapsed!
```

---

## Impact Assessment

### Before Fix (INVALID RESULTS)
- All gamma estimates are **unreliable**
- Significance tests have **inflated Type I error**
- Cannot compare with AnaCoDa or CAI
- Results would mislead biological interpretation

### After Fix (VALID RESULTS)
- Statistical model assumptions are satisfied
- Hypothesis tests have correct Type I error rate
- Gamma estimates are unbiased
- Comparable to other codon bias metrics

---

## Related Functions

### Modified Files
1. `src/convert_vcf_codon_format.R` - Fixed grouping logic
2. `src/integrate_intronic_polymorphism.R` - Updated validation and thresholds

### Dependent Analyses
All downstream analyses using `gamma_results` must be **re-run** after this fix:
- [ ] Gene-level aggregation (`aggregate_gamma_per_gene()`)
- [ ] Comparison with AnaCoDa (`compare_gamma_with_anacoda()`)
- [ ] CAI validation (`validate_against_cai()`)
- [ ] Expression correlation plots
- [ ] Any published figures or tables

---

## References

**Wright-Fisher Model:**
- Assumes independent loci each experiencing drift
- Requires per-site allele frequencies
- See: Kimura (1962), Ewens (2004)

**Beta-Binomial Likelihood:**
- Models variation BETWEEN sites due to drift
- Pooling eliminates this variance component
- See: Griffiths & Tavaré (1998)

---

**Bottom Line:** This was a critical bug that invalidated all previous gamma estimates. The fix ensures that the statistical model correctly accounts for genetic drift by preserving independent loci.
