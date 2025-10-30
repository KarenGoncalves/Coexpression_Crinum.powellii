#### Calculate minimum expression threshold
#### Install and load packages ####
pkgs = c("tidyverse", "DESeq2", "devtools", 
	 "GenomicFeatures", "ggpubr")
pkgs.To.Install = ! pkgs %in% installed.packages()

if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])

# CustomSelection has to be installed from github, so it cannot be grouped with the others yet
if (! "CustomSelection" %in% installed.packages() ) {
	install_github("KarenGoncalves/CustomSelection")
}
pkgs = c(pkgs, "CustomSelection")

for (curPkg in pkgs) library(curPkg, character.only = T) 


#### CustomSelection ####
Exp_table <- read.csv("Data/kallisto_TPM_TMM_normalized.csv",
		      header = T, row.names = 1)

cutv = DAFS(tpm = Exp_table)
data.frame(Replicates = names(cutv),
	   thresholds = unname(cutv)) %>%
	write.table(file = "Data/cutv.tsv",
		    append = F, quote = F, sep = "\t", 
		    row.names = F, col.names = T)
