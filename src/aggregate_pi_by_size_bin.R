#!/usr/bin/env Rscript
# Aggregate per-feature π by size bins
# Input: per-feature CSV from calculate_pi_by_feature.py
# Output: CSV with mean π by size bin and feature type

library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: Rscript aggregate_pi_by_size_bin.R <input.csv> <output.csv>\n", file = stderr())
    quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

# Read data
cat("Reading per-feature π from:", input_file, "\n", file = stderr())
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

# Aggregate by feature_type and size_bin
cat("Aggregating by size bins...\n", file = stderr())

agg_pi_C <- df |>
    group_by(feature_type, size_bin) |>
    summarise(
        count = n(),
        mean_pi_C = mean(pi_C, na.rm = TRUE),
        sd_pi_C = sd(pi_C, na.rm = TRUE),
        mean_q_C = mean(q_C, na.rm = TRUE),
        .groups = 'drop'
    )

agg_pi_G <- df |>
    group_by(feature_type, size_bin) |>
    summarise(
        mean_pi_G = mean(pi_G, na.rm = TRUE),
        sd_pi_G = sd(pi_G, na.rm = TRUE),
        mean_q_G = mean(q_G, na.rm = TRUE),
        .groups = 'drop'
    )

# Combine
result <- agg_pi_C |>
    left_join(agg_pi_G, by = c("feature_type", "size_bin"))

# Sort and format
result <- result |>
    arrange(feature_type, size_bin) |>
    mutate(
        mean_pi_C = round(mean_pi_C, 6),
        sd_pi_C = ifelse(is.na(sd_pi_C), 0, round(sd_pi_C, 6)),
        mean_pi_G = round(mean_pi_G, 6),
        sd_pi_G = ifelse(is.na(sd_pi_G), 0, round(sd_pi_G, 6)),
        mean_q_C = round(mean_q_C, 6),
        mean_q_G = round(mean_q_G, 6)
    )

# Write output
cat("Writing aggregated results to:", output_file, "\n", file = stderr())
write.csv(result, output_file, row.names = FALSE)

# Also print summary
cat("\n=== Summary of nucleotide diversity by size bin ===\n", file = stderr())
print(result)

cat("\nDone.\n", file = stderr())
