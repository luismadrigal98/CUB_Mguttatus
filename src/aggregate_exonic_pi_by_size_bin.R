#!/usr/bin/env Rscript
# Aggregate exonic π by size bins
# Input: per-exon π CSV from calculate_exonic_pi.py
# Output: CSV with mean π and SD by size bin

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: Rscript aggregate_exonic_pi_by_size_bin.R <input.csv> <output.csv>\n", file = stderr())
    quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

# Read data
cat("Reading per-exon π from:", input_file, "\n", file = stderr())
df <- read.csv(input_file, stringsAsFactors = FALSE)

# Define size bins: [1-1000), [1000-2000), ..., [10000+]
create_size_bin <- function(length_bp) {
    if (length_bp < 1000) return("[1-1000)")
    if (length_bp < 2000) return("[1000-2000)")
    if (length_bp < 3000) return("[2000-3000)")
    if (length_bp < 4000) return("[3000-4000)")
    if (length_bp < 5000) return("[4000-5000)")
    if (length_bp < 6000) return("[5000-6000)")
    if (length_bp < 7000) return("[6000-7000)")
    if (length_bp < 8000) return("[7000-8000)")
    if (length_bp < 9000) return("[8000-9000)")
    if (length_bp < 10000) return("[9000-10000)")
    return("[10000+]")
}

# Add size bins
df <- df |>
    mutate(
        size_bin = sapply(length_bp, create_size_bin),
        size_bin = factor(size_bin, levels = c(
            "[1-1000)", "[1000-2000)", "[2000-3000)", "[3000-4000)",
            "[4000-5000)", "[5000-6000)", "[6000-7000)", "[7000-8000)",
            "[8000-9000)", "[9000-10000)", "[10000+]"
        ))
    )

# Aggregate by size_bin
cat("Aggregating by size bins...\n", file = stderr())

result <- df |>
    group_by(size_bin) |>
    summarise(
        n_exons = n(),
        n_total_sites = sum(n_sites),
        mean_pi = mean(pi, na.rm = TRUE),
        sd_pi = sd(pi, na.rm = TRUE),
        min_pi = min(pi, na.rm = TRUE),
        max_pi = max(pi, na.rm = TRUE),
        .groups = 'drop'
    ) |>
    arrange(size_bin) |>
    mutate(
        mean_pi = round(mean_pi, 6),
        sd_pi = ifelse(is.na(sd_pi), 0, round(sd_pi, 6)),
        min_pi = round(min_pi, 6),
        max_pi = round(max_pi, 6)
    )

# Write output
cat("Writing aggregated results to:", output_file, "\n", file = stderr())
write.csv(result, output_file, row.names = FALSE)

# Print summary
cat("\n=== Mean nucleotide diversity (π) by exon size bin ===\n", file = stderr())
print(result)

cat("\nDone.\n", file = stderr())
