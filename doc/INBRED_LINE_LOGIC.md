# Inbred Line Genotype Counting Logic

**Date:** 2024-12-16  
**Context:** Sample size calculation for Wright-Fisher selection inference

---

## Biological Context

### Inbred Lines vs. Outbred Populations

**Inbred Lines (F ≈ 1):**
- Near-complete homozygosity at most loci
- Each individual carries essentially one allelic state
- Sample size n = number of homozygous individuals
- Expected n ≈ 187 (for 187 lines)

**Outbred Diploids (F ≈ 0):**
- High heterozygosity expected under Hardy-Weinberg
- Each individual contributes two alleles
- Sample size n = 2 × number of individuals
- Expected n ≈ 374 (for 187 diploids)

### Residual Heterozygosity

**Sources:**
1. Recent outcrossing events
2. Relaxed selection against heterozygosity
3. Genotyping errors
4. Structural variation (duplications/deletions)

**Problem for Selection Inference:**
- Heterozygous calls in inbred lines are **unreliable**
- Could represent true segregation OR technical artifacts
- Including them biases sample size upward and introduces noise

---

## Implementation

### Python VCF Filtering (`filter_vcf_for_introns.py`)

**Genotype Counting Logic:**
```python
c0 = 0  # Ref allele count (homozygotes only)
c1 = 0  # Alt allele count (homozygotes only)
het_count = 0  # Track heterozygous calls

for sample_str in parts[9:]:
    # ... (depth filtering) ...
    
    gt_val = sample_fields[gt_idx]
    
    # INBRED LINES: Count only homozygous calls
    if gt_val == '0/0' or gt_val == '0|0':
        c0 += 1  # One homozygous individual
    elif gt_val == '1/1' or gt_val == '1|1':
        c1 += 1  # One homozygous individual
    elif '/' in gt_val or '|' in gt_val:
        # Heterozygous call - skip due to residual heterozygosity
        het_count += 1
        continue

total_n = c0 + c1  # Sample size = homozygous individuals only
```

**Key Points:**
1. Only homozygous genotypes (0/0, 1/1) contribute to allele counts
2. Heterozygous genotypes (0/1) are **excluded** as unreliable
3. Sample size n = c0 + c1 ≈ 187 (not doubled)
4. Missing data (./.) already filtered by AD=0,0,0 check

---

### R Validation (`src/convert_vcf_codon_format.R`)

**Updated Thresholds:**
```r
# Check sample sizes (should be ~187 for inbred lines)
cat("Sample size distribution (n = homozygous genotypes per site):\n")
cat(sprintf("  Mean: %.1f\n", mean(result$n)))
cat(sprintf("  Median: %.0f\n", median(result$n)))
cat(sprintf("  Range: %d - %d\n\n", min(result$n), max(result$n)))

if (mean(result$n) > 500) {
  warning("⚠️  CRITICAL: Sample sizes are too large! You may be summing across sites!")
}

if (mean(result$n) > 250) {
  warning("⚠️  WARNING: Sample sizes suggest diploid counting. For inbred lines, expect n≈187.")
}
```

**Rationale:**
- n > 250: Likely counting both alleles per individual (diploid logic)
- n > 500: Definitely pooling across sites (critical bug)
- n ≈ 187: Correct for homozygous-only counting

---

## Wright-Fisher Model Implications

### Effective Sample Size

For the Wright-Fisher diffusion model:
```
P(k | n, α, β, S) = Beta-Binomial with selection
```

**Inbred Lines:**
- n = number of independent **chromosomes** (homozygous individuals)
- Each observation is one allelic state
- Variance reflects genetic drift: Var(p) ∝ p(1-p)/n

**Key Insight:** Excluding heterozygotes **reduces** sample size but **increases** data quality. Better to have n=187 reliable observations than n=374 noisy ones.

---

## Expected Output

### Before Fix (WRONG)
```
Sample size distribution:
  Mean: 374.2    # Diploid counting or heterozygotes included
  Median: 374
  Range: 360 - 400
```

### After Fix (CORRECT)
```
Sample size distribution (n = homozygous genotypes per site):
  Mean: 187.4    # Homozygous individuals only
  Median: 187
  Range: 180 - 195
```

**Interpretation:**
- n ≈ 187: All lines successfully genotyped
- Range 180-195: Minor missing data (~5-10% per site)
- No heterozygotes inflating counts

---

## Heterozygosity Rate Monitoring

To assess data quality, track heterozygosity:

```python
# Add to processing loop
total_het_rate = sum(het_counts) / sum(total_genotypes)
print(f"Residual heterozygosity rate: {total_het_rate:.3%}")
```

**Expected:**
- < 1%: High-quality inbred lines
- 1-5%: Moderate inbreeding
- > 5%: Check for outcrossing or mixed samples

**Note:** This diagnostic was not implemented in the current version but could be added for QC.

---

## Summary

| Aspect | Inbred Lines | Outbred Diploids |
|--------|-------------|------------------|
| Genotypes counted | 0/0, 1/1 only | 0/0, 0/1, 1/1 all |
| Sample size formula | n = homozygotes | n = 2 × individuals |
| Expected n (187 samples) | ~187 | ~374 |
| Heterozygote handling | **Excluded** | Included |
| Justification | F≈1, residual het unreliable | F≈0, het expected |

**Bottom line:** For inbred lines, count individuals (not alleles), and exclude heterozygous calls to avoid noise from residual heterozygosity or genotyping errors.
