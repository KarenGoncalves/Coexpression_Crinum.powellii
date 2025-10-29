### Improve baits file

library(tidyverse)

genes_interest = read_delim("RAW_DATA/Genes_of_interest.txt") %>% 
    group_by(gene_ID) %>% 
    summarize(Description = paste0(Name, collapse = ",")) %>% 
    write_tsv("inputs/Genes_of_interes.txt")
 