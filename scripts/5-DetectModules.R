rm(list = ls())
gc()
source("scripts/FUNCTIONS.R") # loads packages too
Baits = read_delim("inputs/Baits.tsv")
minGenes = 5
r_cutoff = 0.9
# future::plan() - sets parallel session (using multisession())
plan(multisession, workers = detectCores())
# Find gene co-expression modules with Leiden algorithm (graph-based method)
edge_table_select = 
    list.files(path="output_tables/",
               pattern="GeneCoexpression_all_vs_all",
               full.names = T) %>% 
    read_delim()
names(edge_table_select)[1:2] <- 
    c("to", "from")
dim(edge_table_select)
genes_coexpressed = 
    c(edge_table_select$to, edge_table_select$from) %>% 
    unique

Exp_table_long_averaged_z_high_var <- 
    list.files("output_tables/",
               pattern = "Exp_table_long_averaged",
               full.names = T) %>% 
    read_delim() %>% select(-group) %>% 
    filter(gene_ID %in% genes_coexpressed)
    
    
# Load gene annotation table
funct_anno = 
    list.files("inputs", pattern="Filtered_gene_annotation",
               full.names = T) %>% 
    read_delim(na = ".") %>% 
    select(gene_id, sprot_Top_BLASTP_hit, Pfam, 
           gene_ontology_BLASTP, gene_ontology_Pfam, 
           Kegg, EggNM.COG_category, EggNM.Description, 
           EggNM.GOs) %>% 
    filter(gene_id %in% genes_coexpressed) %>%
    mutate(BLASTP_hit = gsub("^([A-Z0-9_]+)\\^\\1.+", "\\1", sprot_Top_BLASTP_hit),
           BLASTP_name = gsub(".+RecName: Full=(.+);\\^.+", "\\1", sprot_Top_BLASTP_hit),
           Pfam_ID = gsub("^(PF\\d+\\.\\d+)\\^.+", "\\1", Pfam),
           Pfam_name = gsub("[PF0-9\\.]+\\^([A-Za-z0-9_\\-]+).+", "\\1", Pfam)) %>% 
    select(-sprot_Top_BLASTP_hit, -Pfam) %>% 
    unique

# Merge annotation and edge tables
node_table <- data.frame(
	gene_ID = genes_coexpressed) %>% 
	left_join(funct_anno %>% 
	              group_by(gene_id) %>% 
	              summarize(Description = paste0(BLASTP_name, collapse = ";")), 
	          by = join_by("gene_ID" == "gene_id")) %>%
    unique()

head(node_table)
dim(node_table)

# Create graph
my_network <- graph_from_data_frame(
	edge_table_select,
	vertices = node_table,
	directed = F
)

#' Too low resolution leads to forcing genes with 
#' different expression patterns into the same module.
#' Too high resolution leads to many genes not contained in any one module.

optimization_results <- future_map(.options = furrr_options(seed = TRUE),
	.x = seq(from = 0.25, to = 5, by = 0.25),
	\(x) optimize_resolution(network = my_network, 
				 resolution = x,
				 minGenes = minGenes)) %>% 
	lapply(\(x) data.frame(V1 = x[1], V2 = x[2])) %>%
	list_rbind %>% 
	mutate(resolution = seq(from = 0.25, to = 5, by = 0.25)) %>% 
 	rename(num_module = V1,
 	       num_contained_gene = V2)

head(optimization_results)

## Plot optimization results ##
transform_factor = 
	with(optimization_results,
	     max(num_contained_gene)/max(num_module)
	)

optimization_results %>% 
	ggplot(aes(x = resolution)) +
	geom_line(aes(y = num_module),
		  linewidth = 1.1, alpha = 0.8, color = "dodgerblue2") +
	geom_point(aes(y = num_module),
		   size = 3, alpha = 0.7) +
	geom_line(aes(y = num_contained_gene/transform_factor),
		  linewidth = 1.1, alpha = 0.8, color = "violetred2") +
	geom_point(aes(y = num_contained_gene/transform_factor),
		   size = 3, alpha = 0.7) +
	scale_y_continuous(
		sec.axis = sec_axis(
			~. * transform_factor,
			name = paste("num. genes in\nmodules w/ >=", minGenes, "genes")
		)) +
	labs(x = "resolution parameter",
	     y = paste("num. modules w/ >=", minGenes, "genes")) +
	theme_optimization


ggsave("plots/Optimize_resolution.png", 
       height = 5, width = 5, bg = "white")

## Select resolution and detect modules
resolution <- 
    (optimization_results %>% 
    slice_max(order_by = num_contained_gene) %>% 
    slice_max(order_by = resolution))$resolution

modules <- cluster_leiden_mod(
    my_network, r = resolution
)

# Merge edge and node tables
my_network_modules <- data.frame(
    gene_ID = names(membership(modules)),
    module = as.vector(membership(modules))) %>% 
    inner_join(node_table, by = "gene_ID")

# Continue only with modules with 5 or more genes
module_minGenes <- my_network_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= minGenes)

my_network_modules <- my_network_modules %>% 
    filter(module %in% module_minGenes$module)

## Peak expression of modules
Expr_averaged_z_high_var_modules <- 
    Exp_table_long_averaged_z_high_var %>% 
    inner_join(my_network_modules, by = "gene_ID")

modules_mean_z <- Expr_averaged_z_high_var_modules %>% 
    group_by(module, Sample_name) %>% 
    summarise(mean.z = mean(z.score)) %>% 
    ungroup()

# This will create a table with the module number and the Sample in which the module is most expressed.
module_peak_exp <- modules_mean_z %>% 
    group_by(module) %>% 
    slice_max(order_by = mean.z, n = 1)

paste0("RDATA/MainAnalysis_Cor", r_cutoff) %>% dir.create()
save(Expr_averaged_z_high_var_modules,
     module_peak_exp, modules_mean_z, my_network_modules,
     Baits, funct_anno,
     file = paste0(
         "RDATA/MainAnalysis_Cor", r_cutoff, "/ModuleData_forPlot_res",
         resolution, ".RData")
)
