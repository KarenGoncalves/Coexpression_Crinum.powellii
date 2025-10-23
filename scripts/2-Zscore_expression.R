#########################################
## Create Z-transformed expression table
#########################################


source("scripts/FUNCTIONS.R") # loads packages too

## Input files ##
Exp_table <- read.delim("inputs/Filtered_gene_TPM.tsv",
                      header = T)
metadata <- read_delim("inputs/metadata.txt", 
                       col_names = T, delim = "\t") %>% 
    mutate(Sample_name = gsub("_\\d$", "", `Library Name`))

names(Exp_table)[1] <- "gene_ID"


## Long exp_table ##

Exp_table_long <- Exp_table %>% 
    pivot_longer(cols = !gene_ID, 
                 names_to = "library", 
                 values_to = "tpm") %>% 
    mutate(logTPM = log10(tpm + 1)) %>% unique

Exp_table_log_wide <- Exp_table_long %>% 
    pivot_wider(names_from = library, 
                values_from = logTPM, 
                id_cols = gene_ID)

## PCA ##
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

PCA_coord <- my_pca$x[, 1:10] %>% 
    as.data.frame() %>% 
    mutate(`Run Accession` = row.names(.)) %>% 
    full_join(metadata, by = "Run Accession")

axis_titles = sapply(1:2, \(x) {
    paste("PC", x, " (", 
          pc_importance[x, 2] %>% signif(3)*100, 
          "% of Variance)", sep = "")
})

PCA_by_sample <- PCA_coord %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Sample_name), 
               color = "grey20", 
               shape = 21, 
               size = 3, alpha = 0.8) +
    scale_fill_manual(values = brewer.pal(n = 8, "Accent")) +
    labs(x = axis_titles[1], 
         y = axis_titles[2],
         fill = NULL) +  
    theme_bw() +
    theme(text = element_text(size= 14),
          axis.text = element_text(color = "black"),
          legend.position = "bottom"
    )

PCA_by_sample

dir.create("plots/MainAnalysis")
ggsave("plots/MainAnalysis/PCA.png", height = 5, width = 5, bg = "white")

### Gene co-expression analysis ###
# Average up the reps #
Exp_table_long_averaged <- Exp_table_long %>% 
    full_join(metadata, 
              by = join_by("library" == "Run Accession")) %>% 
    group_by(gene_ID, Sample_name) %>% 
    summarise(mean.logTPM = mean(logTPM)) %>% 
    ungroup()  

head(Exp_table_long_averaged)

# Z-score #

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
    group_by(gene_ID) %>% 
    mutate(z.score = zscore(mean.logTPM)) %>% 
    ungroup()
split_at <- 400000

Exp_table_long_averaged_z %>% 
    mutate(group = (row_number() - 1) %/% !! split_at) %>%
    group_split(group) %>%
    map(.f = ~{
        fileName = paste0("output_tables/Exp_table_long_averaged_z_",
                          unique(.x$group), ".tsv")
        write_delim(.x, fileName,
                    quote = "none", append = F,
                    col_names = T, delim = "\t")
    }) 
