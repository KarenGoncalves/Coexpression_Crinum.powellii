pkgs=c("tidyverse", "igraph", "ggraph",
       "readxl", "patchwork", "RColorBrewer",
       "viridis", "parallel", "doParallel",
       "furrr", "future"
)

pkgs.to.install = pkgs[!pkgs %in% installed.packages()]
if (length(pkgs.to.install) != 0) install.packages(pkgs.to.install)
for (pkg in pkgs) library(pkg, character.only = T)
rm(list = ls())

### ggplot theme ###
theme_optimization = 
	theme_classic() +
	theme(text = element_text(size = 14),
	      axis.text = element_text(color = "black"),
	      axis.title.y.right = element_text(color = "violetred2"),
	      axis.title.y.left = element_text(color = "dodgerblue2")
	)

# Z-score function #
zscore <- function(numeric_vector) {
	(numeric_vector - mean(numeric_vector)) / 
		sd(numeric_vector)
}

# optimize_resolution function #
# Method to find best resolution for cluster detection
optimize_resolution <- function(network, resolution, minGenes){
	modules = network %>% 
		cluster_leiden(
			resolution_parameter = resolution,
			objective_function = "modularity")
	
	parsed_modules = data.frame(
		gene_ID = names(membership(modules)),
		module = as.vector(membership(modules)) 
	)
	
	# counts the number of modules found with more than 5 genes
	num_module_minGenes = parsed_modules %>% 
		group_by(module) %>% 
		count() %>% 
		arrange(-n) %>% 
		filter(n >= minGenes) %>% 
		nrow() %>% 
		as.numeric()
	
	# Counts the number of genes contained in all modules with more than 5 genes
	num_genes_contained = parsed_modules %>% 
		group_by(module) %>% 
		count() %>% 
		arrange(-n) %>% 
		filter(n >= minGenes) %>% 
		ungroup() %>% 
		summarise(sum = sum(n)) %>% 
		as.numeric()
	
	return(c(num_module_minGenes, num_genes_contained))
	
} 

get_optimization_table <- function(network, 
				   min_res, 
				   max_res, 
				   step,
				   minGenes) {
	future_map(.options = furrr_options(seed = TRUE),
		   .x = seq(from = min_res, to = max_res, by = step),
		   \(x) optimize_resolution(network = network, 
		   			 resolution = x,
		   			 minGenes)) %>% 
		lapply(\(x) data.frame(V1 = x[1], V2 = x[2])) %>%
		list_rbind %>% 
		mutate(resolution = seq(from = min_res, to = max_res, by = step)) %>% 
		rename(num_module = V1,
		       num_contained_gene = V2)
}

select = dplyr::select


cluster_leiden_mod <- function(network_data, r) {
    cl <- cluster_leiden(network_data, 
                   resolution_parameter = r, 
                   objective_function = "modularity"
    )
    membership <- membership(cl)
    
    # Convert to tibble
    module_df <- tibble(
        gene_ID = names(membership),
        old_module = membership
    )
    
    # Reorder module labels by size
    module_df <- module_df %>%
        count(old_module, name = "size") %>%
        arrange(desc(size)) %>%
        mutate(new_module = row_number()) %>% 
        right_join(module_df, by = "old_module")
    
    # Update the igraph object with new module labels
    cl$membership <- sapply(1:length(cl$membership), \(rowID) {
        gene = cl$names[rowID]
        module_df %>% filter(gene_ID == gene) %>% 
            pull(new_module)
    })
    return(cl)
}

get_T_value <- function(r, nreps) {
	t = r*sqrt( (nreps-2) / (1-r^2) )
	return(t)
}
get_p_value <- function(t, nreps) {
	pvalue = case_when(
		t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
		t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
	)
	return(pvalue)
}

heatmap_theme <- theme_classic() +
	theme(text = element_text(size = 14),
	      axis.text.y = element_text(color = "black"),
	      #axis.text.x = element_blank(),
	      legend.position = "top",
	      strip.text = element_blank(),
	      panel.spacing = unit(0, "lines")
	)
annotation_theme <- theme_void() +
	theme(axis.title.x = element_text(vjust = 8))

##### Heatmap of all genes in the module #####
heatmap_all_genes <- function(cor_matrix, 
			      Exp_table_long_averaged_z,
			      caption, title) {
	gene_clustering <- as.dist(1 - cor_matrix) %>%
		hclust(method = "complete")
	Exp_table_long_averaged_z$gene_clustering_order <- 
		factor(Exp_table_long_averaged_z$gene_ID, 
		       levels = rownames(cor_matrix)[
		       	gene_clustering$order]
		)
	breaks = with(Exp_table_long_averaged_z,
		      c(min(z.score), 0, max(z.score))) 
	heatmap_genes <-
		Exp_table_long_averaged_z %>%
		ggplot(aes(x = Sample_name, y = gene_clustering_order)) +
		geom_tile(aes(fill = z.score), color = NA) +
		scale_fill_gradient2(mid = "white",
				     high = "#67001F",
				     low = "#053061",
				     breaks = breaks,
				     labels = breaks %>% round(2)
				     ) +
		labs(x = NULL, y = "Genes", fill = "z-score",
		     title = title,
		     caption = paste0(caption, "\nGenes in module: ",
		     		 length(gene_clustering$labels))
		)
		return(heatmap_genes)
		       
}
