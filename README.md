Data from _de novo_ assembly

- Used kallisto to map and count reads;
- Used TPM table and removed genes with total expression < 1 TPM
- Kept only protein-encoding genes
- Used coefficient of variation of log10(TPM) as measure of variability
- Select genes in the 80th> quantile of variability for construction of network
- Used bigstatsr::bigcor for computing correlation
- Computed FDR with the whole vector of p-values, instead of one by one in a mutate function

