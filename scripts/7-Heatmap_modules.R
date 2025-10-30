source("scripts/FUNCTIONS.R") # loads packages too,
r_cutoff = 0.9
resolution = 2
RDATA = list.files("RDATA/MainAnalysis_Cor0.9/", full.names = T)
load(RDATA, verbose = T)
# ============================================================
#   2 Clean Tissue Names and Define Order
# ============================================================

tissue_table <- tibble(Sample_name =  c(
    "Root", "Basal_plate", "Bulb", "Leaves", 
    "Dark3", "Dark4", "Light3", "Light4"
)) %>%
    mutate(
        Culture = if_else(startsWith(Sample_name, "Dark") | 
                              startsWith(Sample_name, "Light"), "In vitro", "In vivo"),
        Treatment = case_when(
            startsWith(Sample_name, "Dark3")  ~ "Dark\n2m/L",
            startsWith(Sample_name, "Dark4")  ~ "Dark\n4m/L",
            startsWith(Sample_name, "Light3") ~ "Light\n2m/L",
            startsWith(Sample_name, "Light4") ~ "Light\n4m/L",
            .default=""
        )
    ) %>% 
    mutate(Tissue = ifelse(Culture == "In vitro", Treatment,
                                Sample_name) %>% 
               factor(., levels = unique(.)))
metadata <- read_delim("inputs/metadata.txt") %>%
    mutate(tissue = gsub("_\\d+$", "", `Library Name`)) %>% 
    select(-`Library Name`, -`Run Accession`) %>% 
    right_join(tissue_table, by = join_by("tissue" == "Sample_name")) %>% 
    rename(Sample_name = "tissue")


# Now the modules in the table module_peak_exp are a factor ordered by the 
# Sample in which they are most expressed, so we join the tables:
# metadata, module_mean_z and module_peak_exp
module_peak_exp_mod <- 
    mutate(module_peak_exp, 
           orderedSamples = sapply(Sample_name, \(x) {
                which(tissue_table$Sample_name == x)
           })) %>%
    dplyr::select(module, Sample_name,
                  orderedSamples)

heatmap_data = 
    full_join(modules_mean_z, module_peak_exp_mod[,-c(2, 3)], 
              by = c('module')) %>%
    mutate(ordered_modules = factor(module,
                                    levels = unique(.$module))) %>% 
    # Join with metadata, removing the column Replicate
    full_join(metadata, by = "Sample_name")

# genes per module
genes_per_module = Expr_averaged_z_high_var_modules %>% 
    group_by(module) %>% 
    summarize(gene_ID = unique(gene_ID)) %>% 
    count()

heatmap_data = full_join(heatmap_data, genes_per_module,
                         by = "module") %>%
    arrange(desc(n))

### Heatmap of module mean expression z-score
(module_heatmap <- heatmap_data %>% 
    ggplot(aes(x = Tissue,
               y = fct_rev(ordered_modules))) +
    geom_tile(aes(fill = mean.z), color = "grey80") +
    scale_fill_gradient2(mid = "white",
                         high = "#67001F",
                         low = "#053061",
                         breaks = c(-.75, 0, 2), 
                         labels = c("< -.75", "0", "> 2")) +
    labs(x = NULL, y = "Module", fill = "z-score",
         caption = paste0("r threshold = ", r_cutoff, 
                          "\nresolution = ", 2)
    ) +
    facet_grid(cols = vars(Culture),
               scales = "free") +
    heatmap_theme + 
        theme(strip.text = element_text(),
              strip.background = element_blank()
              ))

(module_heatmap_nGenes = heatmap_data %>%
    select(ordered_modules, n) %>%
    unique %>%
    ggplot(aes(y = fct_rev(ordered_modules), x = "",
               label = n)) +
    geom_tile(fill = "white", color = "white") +
    geom_text() + labs(x = "Genes\nin\nmodule") +
    annotation_theme)

wrap_plots(module_heatmap, 
           module_heatmap_nGenes, 
           ncol = 2, 
           widths = c(1, 0.08))
plotName = paste0("plots/Heatmap_Modules_",
                  "r", r_cutoff, "_res", resolution, 
                  "highVarTop20qnt.tiff")
ggsave(plotName, height = 4.25, width = 8)

# annotation <- funct_anno %>%
#     inner_join(my_network_modules,
#                by = join_by("geneID" == 'gene_ID')) %>%
#     mutate(Module_name = curModule) %>%
#     arrange(module) %>%
#     write_delim(file = paste0("output_tables/Module", 
#                               curModule,
#                               "_annotations.tsv"),
#                 delim = "\t", col_names = T, na = "", 
#                 quote = "none")
