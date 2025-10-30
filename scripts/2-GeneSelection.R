#============================================================
#   1 Load Packages and Source Function, and input files
#============================================================

source("scripts/FUNCTIONS.R") # loads packages too

Exp_table_long_averaged_z <- 
    list.files(full.names = T, path = "output_tables/",
               pattern = "Exp_table_long_averaged_z_") %>%
    read_delim(delim = "\t", show_col_types = FALSE) %>% 
    select(-group) %>% 
    mutate(log10meanTPM = log10(mean.TPM + 1))

expr_gt.1TPM = Exp_table_long_averaged_z %>% 
    group_by(gene_ID) %>% 
    slice_max(order_by = mean.TPM) %>% 
    filter(mean.TPM > 1) %>% 
    select(gene_ID) %>% 
    unique
length(expr_gt.1TPM$gene_ID)

Exp_table_long_averaged_z <-
    Exp_table_long_averaged_z %>% 
    filter(gene_ID %in% expr_gt.1TPM$gene_ID)

ngenes = unique(Exp_table_long_averaged_z$gene_ID) %>% length
genes_expressed <- unique(Exp_table_long_averaged_z$gene_ID)

metadata <- read_delim("inputs/metadata.txt") %>%
    mutate(tissue = gsub("_\\d+$", "", `Library Name`))

Baits <- read_delim("inputs/Genes_of_interes.txt") %>% 
    filter(gene_ID %in% genes_expressed)

## Long exp_table ##
# Gene selection #

all_var_and_ranks <- Exp_table_long_averaged_z %>% 
	group_by(gene_ID) %>% 
	summarise(CV.TPM = sd(mean.TPM)/mean(mean.TPM),
	          CV.logTPM = sd(log10meanTPM)/mean(log10meanTPM)) %>% 
	ungroup() %>% 
	mutate(rank = rank(CV.TPM, ties.method = "average"),
	       rank = rank(CV.logTPM, ties.method = "average")) 
write_delim(all_var_and_ranks,
	    "output_tables/all_var_and_ranks.tsv",
	    col_names = T, quote = "none",
	    delim = "\t", na = '')

high_var_genes_20qnt <- all_var_and_ranks %>%
	filter(CV.logTPM > quantile(CV.logTPM, 0.8))

# high_var_genes_20pct <- high_var_genes %>% 
# 	slice_max(order_by = CV.TPM, 
# 		  n = round(ngenes/5, digits = 0)
# 	)

write_delim(high_var_genes_20qnt,
	    file = "output_tables/high_var_genes_20qnt.tsv",
	    col_names = T, delim = "\t", quote = "none"
)

bait_var <- high_var_genes_20qnt %>% 
	inner_join(Baits, by = "gene_ID") %>% 
	group_by(gene_ID) %>% 
	slice_max(n = 1, order_by = CV.TPM)

Exp_table_long_averaged_z_high_var <- 
	Exp_table_long_averaged_z %>% 
	filter(gene_ID %in% high_var_genes_20qnt$gene_ID)

save(Exp_table_long_averaged_z_high_var,
     high_var_genes_20qnt, 
     file = "RDATA/GeneSelection_objects.RData"
     )
# Check where the bait genes fall along the variance distribution
var_plot = all_var_and_ranks %>% 
	ggplot(aes(x = CV.logTPM, y = rank))  +
	geom_hline(
		data = bait_var, aes(yintercept = rank),
		color = "tomato1", linewidth = 0.1, alpha = 0.5
	) +
	geom_vline(
		data = bait_var, aes(xintercept = CV.logTPM), 
		color = "tomato1", linewidth = 0.1, alpha = 0.5
	) +
	geom_rect( 
		xmax = max(high_var_genes_20qnt$CV.logTPM), 
		xmin = min(high_var_genes_20qnt$CV.logTPM),
		ymax = nrow(all_var_and_ranks),
		ymin = nrow(all_var_and_ranks) - nrow(high_var_genes_20qnt),
		fill = "dodgerblue2", alpha = 0.1
	) + 
	geom_line(linewidth = 0.8) +
	labs(y = "rank",
	     x = "CV(log10(TPM))",
	     caption = "Blue box = top 20 quantile high CV genes.\nRed lines = bait genes.") +
	theme_classic() +
	theme(
		text = element_text(size = 11),
		axis.text = element_text(color = "black"),
		plot.caption = element_text(hjust = 0)
	) 

ggsave(plot = var_plot, 
       filename = "plots/gene_var_distribution.png", 
       height = 5, width = 5)
