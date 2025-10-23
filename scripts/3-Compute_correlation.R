#============================================================
#   1 Load Packages and Source Function, and input files
#============================================================

source("scripts/FUNCTIONS.R") # loads packages too
r_cutoff <- 0.7

Exp_table_long_averaged_z <- 
    list.files(full.names = T, path = "output_tables/",
               pattern = "Exp_table_long_averaged_z_") %>%
    read_delim(delim = "\t", show_col_types = FALSE) %>% 
    select(-group)
genes_expressed <- unique(Exp_table_long_averaged_z$gene_ID)

metadata <- read_delim("inputs/metadata.txt") %>%
    mutate(tissue = gsub("_\\d+$", "", `Library Name`))

Baits <- read_delim("inputs/Genes_of_interest.txt") %>% 
    filter(gene_ID %in% genes_expressed)

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
    as.data.frame()

nreps <- ncol(z_score_wide) - 1

# ============================================================
#   4 Correlation Against Baits
# ============================================================

correlation_against_baits <- map_dfr(Baits$gene_ID, function(x) {
    bait <- z_score_wide %>% filter(gene_ID == x) %>% 
        select(-gene_ID) %>% as.matrix() %>% t()
    other_mat <- z_score_wide %>% filter(gene_ID != x) %>% 
        select(-gene_ID) %>% as.matrix() %>% t()
    
    cor_vals <- cor(bait, other_mat) %>% as.vector()
    r2 <- cor_vals^2
    t_stat <- cor_vals * sqrt((nreps - 2) / (1 - r2))
    p_val <- pt(t_stat, df = nreps - 2, lower.tail = cor_vals <= 0)
    fdr <- p.adjust(p_val, method = "fdr")
    
    tibble(
        Gene1 = x,
        Gene2 = (z_score_wide %>% filter(gene_ID != x))$gene_ID,
        correlation = cor_vals,
        r2 = r2,
        t = t_stat,
        p.value = p_val,
        FDR = fdr,
        significant = FDR < 0.01
    )
})

# ============================================================
#   5 Correlation Distribution Plot
# ============================================================

correlation_against_baits %>%
    ggplot() +
    geom_histogram(aes(x = correlation), fill = "black", bins = 100) +
    geom_vline(xintercept = 0.9, color = "red", linewidth = 0.85) +
    labs(
        title = "Gene correlations - Bait genes",
        x = "Number of correlations"
    ) +
    scale_x_continuous(breaks = seq(-1, 1, 0.4)) +
    theme_classic() +
    theme(
        text = element_text(size = 14),
        axis.text = element_text(color = "black")
    )

# ============================================================
#   6 Filter Correlations
# ============================================================

edge_table <- correlation_against_baits %>%
    filter(correlation > r_cutoff, significant) %>%
    rename(from = Gene1, to = Gene2) %>%
    relocate(from, to, .before = everything())

# ============================================================
#    7 All-vs-All Correlations
# ============================================================

genes_to_keep <- edge_table %>%
    distinct(to) %>%
    pull(to) %>%
    union(unique(edge_table$from))

Expr_filtered_long <- Exp_table_long_averaged_z %>%
    filter(gene_ID %in% genes_to_keep)

Expr_filtered_wide <- Expr_filtered_long %>%
    mutate(tissue = Sample_name) %>% 
    select(gene_ID, tissue, z.score) %>%
    pivot_wider(names_from = tissue, values_from = z.score) %>%
    column_to_rownames("gene_ID")

cor_matrix <- cor(t(Expr_filtered_wide))
lower_tri_indices <- which(lower.tri(cor_matrix), arr.ind = TRUE)

# Extract correlation values
r_vals <- cor_matrix[lower_tri_indices]

# Compute t-statistics
t_stats <- ifelse(abs(r_vals) == 1, 
                  Inf, r_vals * sqrt(nreps - 2) / sqrt(1 - r_vals^2))

# Compute p-values
p_vals <- ifelse(abs(r_vals) == 1, 0, 
                 pt(-abs(t_stats), df = nreps - 2))
# Adjust p-values using FDR
fdr_vals <- p.adjust(p_vals, method = "fdr")

# Create result data frame with variable names
cor_long <- data.frame(
    Gene1 = rownames(cor_matrix)[lower_tri_indices[, 1]],
    Gene2 = colnames(cor_matrix)[lower_tri_indices[, 2]],
    correlation = r_vals,
    r2 = r_vals^2,
    t_stat = t_stats,
    p_value = p_vals,
    FDR = fdr_vals
)
rm(t_stats, r_vals, p_vals, fdr_vals, lower_tri_indices)

cor_long_filtered <- cor_long %>%
    filter(FDR < 0.05, correlation > r_cutoff) #%>% nrow

write_tsv(cor_long_filtered, 
          "output_tables/GeneCoexpression_all_vs_all_filtered.tsv")

