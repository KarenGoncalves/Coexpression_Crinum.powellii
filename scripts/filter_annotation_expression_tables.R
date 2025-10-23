library(tidyverse)

report = read_delim("RAW_DATA/Crinum_powellii.tsv")
protein_coding = report %>% filter(prot_id != ".")
names(protein_coding)[1] = "gene_id"

expressed_iso = read_delim("RAW_DATA/Crinum_powellii_isoform_expressed.txt", 
                           delim = "\t", col_names = "transcript_id")
expressed_gene = read_delim("RAW_DATA/Crinum_powellii_gene_expressed.txt", 
                            delim = "\t", col_names = "gene_id")


inner_join(expressed_gene, protein_coding, by = "gene_id") %>% nrow
TPM_iso = read_delim("RAW_DATA/Crinum_powellii_kallisto.isoform.TPM.not_cross_norm")
names(TPM_iso)[1] = "transcript_id"

TPM_gene = read_delim("RAW_DATA/Crinum_powellii_kallisto.gene.TPM.not_cross_norm")
names(TPM_gene)[1] = "gene_id"

Filtered_TPM_gene_prot = TPM_gene %>% inner_join(expressed_gene, by = "gene_id") %>% 
    inner_join(protein_coding %>% select(gene_id), by = "gene_id")

Filtered_TPM_iso_prot = TPM_iso %>% inner_join(expressed_iso, by = "transcript_id") %>% 
    inner_join(protein_coding %>% select(transcript_id), by = "transcript_id")

max_lines = 50000
Filtered_TPM_gene_prot %>% select(gene_id) %>% 
    inner_join(protein_coding, by = "gene_id") %>% 
    mutate(file_id = (row_number() - 1) %/% max_lines + 1) %>% 
    group_by(file_id) %>%
    group_walk(~ write_tsv(.x, paste0("inputs/Filtered_gene_annotation_", 
                                      .y$file_id, ".csv")))
Filtered_TPM_gene_prot %>% 
    write.table("inputs/Filtered_gene_TPM.tsv", 
                append = F, quote = F, sep = "\t", row.names = F)


Filtered_TPM_iso_prot %>% 
    write.table("inputs/Filtered_iso_TPM.tsv", 
                append = F, quote = F, sep = "\t", row.names = F)
