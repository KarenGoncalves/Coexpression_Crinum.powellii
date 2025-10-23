#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
DIR <- args[1]           # Directory containing input and output files
in_suffix <- args[2]     # Suffix for input files (e.g., "_counts.txt")
out_suffix <- args[3]    # Suffix for output files (e.g., "_expressed.txt")
species <- args[4]       # Species name (without spaces)

# Load required libraries
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)


# Construct input and output file paths
input <- paste0(DIR, "/", species, in_suffix)
output <- paste0(DIR, "/", species, out_suffix)

# Read input file
counts <- read_delim(input)

# Rename first column to "isoform"
names(counts)[1] <- "isoform"

# Calculate total expression across all samples
counts$Total <- rowSums(counts[, -1])

# Filter isoforms with non-zero expression
filtered <- counts %>% filter(Total > 0)

# Write filtered isoform names to output file
write.table(filtered$isoform,
            file = output,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Print first few rows for verification
head(filtered)

