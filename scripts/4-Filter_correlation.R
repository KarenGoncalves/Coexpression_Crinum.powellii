#============================================================
#   1 Load Packages and Source Function, and input files
#============================================================

source("scripts/FUNCTIONS.R") # loads packages too
r_cutoff <- 0.9

Exp_table_long_averaged_z <- 
    list.files(full.names = T, path = "output_tables/",
               pattern = "Exp_table_long_averaged_z_") %>%
    read_delim(delim = "\t", show_col_types = FALSE) %>% 
    select(-group) %>% 
    inner_join(read_delim("output_tables/high_var_genes_20pct.tsv"), 
               by = "gene_ID")

metadata <- read_delim("inputs/metadata.txt") %>%
    mutate(tissue = gsub("_\\d+$", "", `Library Name`))

Baits <- read_delim("inputs/Baits.tsv") %>% 
    filter(gene_ID %in% Exp_table_long_averaged_z$gene_ID)

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
#   ️Plot correlation
# ============================================================

(list.files(full.names = T,
           path="output_tables",
           pattern="correlation_full_20qnt_") %>% 
    read_delim() %>%
    mutate(Significant = fdr_vals < 0.01) %>% 
    ggplot() +
    geom_histogram(aes(x = r_vals, fill = Significant), bins = 100) +
    geom_vline(xintercept = r_cutoff, color = "red", linewidth = 0.85) +
    labs(
        title = "Gene correlations",
        y = "Number of correlations",
        x = "Correlation score (r)",
        caption = "Correlations between top 20% most variable genes"
    ) +
    scale_x_continuous(breaks = seq(-1, 1, 0.4)) +
    theme_classic() +
    theme(
        text = element_text(size = 14),
        axis.text = element_text(color = "black"))) %>%
 ggsave(filename = "plots/Histogram_significant_correlations_20qnt.png")

cor_long_filtered <- list.files(full.names = T,
                                path="output_tables",
                                pattern="correlation_full_20qnt_") %>% 
    read_delim() %>%
    select(-group) %>% 
    mutate(Significant = fdr_vals < 0.01) %>% 
    filter(Significant, r_vals > r_cutoff)

ngroups=5
cor_long_filtered %>%
    mutate(group = ntile(row_number(), ngroups)) %>%
    group_split(group) %>%
    walk2(1:ngroups, ~ write_tsv(
        .x %>% select(-group), paste0("output_tables/GeneCoexpression_all_vs_all_filtered_", .y, ".tsv"),
        num_threads = 12)
        )

cor_long_filtered %>% 
ggplot() +
    geom_histogram(aes(x = r_vals, fill = Significant), bins = 100) +
    labs(
        title = "Gene correlations",
        y = "Number of correlations",
        x = "Correlation score (r)",
        caption = "Correlations between top 20% most variable genes"
    ) +
    # scale_x_continuous(breaks = seq(-1, 1, 0.4)) +
    theme_classic() +
    theme(
        text = element_text(size = 14),
        axis.text = element_text(color = "black"))

perfect_cors = cor_long_filtered %>% filter(r_vals == 1, is.infinite(t_vals))

Exp_table_long_averaged_z %>% 
    filter(gene_ID %in% unique(perfect_cors$Gene1, perfect_cors$Gene2)) %>%
    View
# cor_long_filtered <- 
#     list.files(path="output_tables/", pattern="GeneCoexpression_all_vs_all",
#                full.names = T) %>% 
#     read_delim()
genes <- c(cor_long_filtered$Gene1, cor_long_filtered$Gene2) %>% unique
Baits %>% filter(gene_ID %in% genes)
