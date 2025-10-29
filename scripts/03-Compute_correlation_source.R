#============================================================
#   1 Load Packages and Source Function, and input files
#============================================================

# Install if needed: 
# install.packages("bigstatsr")
# install.packages("bigparallelr")
source("scripts/FUNCTIONS.R") # loads packages too
library(bigstatsr)
library(furrr)
library(purrr)
library(dplyr)
library(readr)

r_cutoff <- 0.7

genes_selected <-
    read_delim("output_tables/high_var_genes_20qnt.tsv")$gene_ID %>% unique

Exp_table_long_averaged_z <- 
    list.files(full.names = T, path = "output_tables/",
               pattern = "Exp_table_long_averaged_z_") %>%
    read_delim(delim = "\t", show_col_types = FALSE) %>% 
    select(-group) %>% 
    filter(gene_ID %in% genes_selected)


metadata <- read_delim("inputs/metadata.txt") %>%
    mutate(tissue = gsub("_\\d+$", "", `Library Name`))

Baits <- read_delim("inputs/Baits.tsv") %>% 
    filter(gene_ID %in% genes_selected)

# ============================================================
#   2 Clean Tissue Names and Define Order
# ============================================================

tissue_order <- c(
    "Root", "Basal_plate", "Bulb", "Leaf", 
    "Dark3", "Dark4", "Light3", "Light4"
)

metadata <- metadata %>%
    mutate(TissueOrdered = factor(tissue, levels = tissue_order),
           In_vivo = case_when(tissue %in% 
                                   paste0(rep(c("Dark", "Light"), each = 2), c(3,4)) ~
                                   "In vitro",
                               .default = "In vivo")
    )

# ============================================================
#   ️3 Prepare z-score Matrix
# ============================================================

z_score_wide <- Exp_table_long_averaged_z %>%
    mutate(tissue = Sample_name) %>% 
    select(gene_ID, tissue, z.score) %>%
    pivot_wider(names_from = tissue, values_from = z.score) %>%
    select(-gene_ID) %>% t %>% unique

colnames(z_score_wide) = unique(Exp_table_long_averaged_z$gene_ID)
nreps <- nrow(z_score_wide)

# ============================================================
#   4 Correlation 
# ============================================================

# Use `bigstatsr` for large matrices

# Convert to FBM (Filebacked Big Matrix)
X <- as_FBM(as.matrix(z_score_wide))
cor_matrix_big <- big_cor(X, block.size = 300)

# Convert to regular matrix
cor_matrix <- cor_matrix_big[]
rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(z_score_wide)
rm(cor_matrix_big)

## Extract correlation values
# Make long table of results
plan(multisession, workers = 12)

# Output file path
output_file <- "cor_table.csv"

# Write header once
tibble(Gene1 = character(), 
       Gene2 = character(), 
       r_vals = numeric(), 
       t_vals = numeric(), 
       p_vals = numeric()) %>% 
write_csv(output_file)

# Generate index pairs
index_pairs <- map(2:nrow(cor_matrix), ~ tibble(i = .x, j = 1:(.x - 1))) %>%
    bind_rows()

# Split into chunks to avoid memory overload
chunk_size <- 10000
chunks <- split(index_pairs, ceiling(seq_along(index_pairs$i) / chunk_size))

# Process each chunk and write to disk
walk(chunks, function(chunk) {
    future_pmap_dfr(.progress = T, chunk, function(i, j) {
        tibble(
            Gene1 = colnames(z_score_wide)[i],
            Gene2 = colnames(z_score_wide)[j],
            r_vals = cor_matrix[i, j]
        ) %>%
            mutate(
                t_vals = ifelse(abs(r_vals) == 1, Inf, 
                                r_vals * sqrt(nreps - 2) / sqrt(1 - r_vals^2)),
                p_vals = ifelse(abs(r_vals) == 1, 0, 
                                pt(-abs(t_vals), df = nreps - 2))
            )
    }) %>% 
        write_csv(output_file, append = TRUE)
})
# Adjust p-values using FDR

cor_table <- read_csv("cor_table.csv",
                      num_threads = parallel::detectCores())
    
cor_table$fdr_vals <- p.adjust(cor_table$p_vals, method = "fdr")

ngroups=40
cor_table %>%
    mutate(group = ntile(row_number(), ngroups)) %>%
    group_split(group) %>%
    walk2(1:ngroups, ~ 
              write_tsv(.x, paste0("output_tables/correlation_full_20qnt_", .y, ".tsv"),
                        num_threads = parallel::detectCores())
          )
