library(tidyverse)
library(data.table)
library(clusterProfiler)
library(viridis)

DEGs <- fread("output/dexseq/dtu_exposed_analysis/sig_degs.csv")

###############
## GO annots ##
###############

##GO annot tables
term_to_gene <- fread("output/dexseq/enrichment/GO_to_geneID.csv")
term_to_name<- fread("output/dexseq/enrichment/GO_to_GOname.csv")
go_annot_table <- fread("output/dexseq/enrichment/GO_annots.csv")
go_annot_table <- go_annot_table[,c(1,2)]

results <- enricher(DEGs$groupID, TERM2GENE=term_to_gene, TERM2NAME = term_to_name)
##seems legit enough?
results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)
##merge with annots table
plot_results_dt <- merge(results_dt, go_annot_table, by.x="ID", by.y="pathway")
plot_results_dt$total_DEG_GOs <- tstrsplit(plot_results_dt$GeneRatio, "/", keep=c(2))
plot_results_dt$total_DEG_GOs <- as.numeric(plot_results_dt$total_DEG_GOs)
plot_results_dt$GeneRatio <- (plot_results_dt$Count/plot_results_dt$total_DEG_GOs)
fwrite(plot_results_dt, "output/dexseq/dtu_exposed_analysis/exposed_GO_enrichment.csv")

## GeneRatio - number of DEGs with that GO annot/total DEGs with GO annots - DEG list
## BgRatio - total no. genes with that GO annot/total number of genes with GO annots - whole transcriptome
## Q value?

#################
## Pfam annots ##
#################

pfam_to_gene <- fread("output/dexseq/enrichment/pfam_to_geneID.csv")
pfam_to_name<- distinct(fread("output/dexseq/enrichment/pfam_to_domain_name.csv"))
pfamannot_table <- fread("output/dexseq/enrichment/pfam_annots.csv")
pfamannot_table <- pfamannot_table[,c(1,2)]
pfamannot_table <- distinct(pfamannot_table)

results <- enricher(DEGs$groupID, TERM2GENE=pfam_to_gene, TERM2NAME = pfam_to_name)
##seems legit enough?
pfam_results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)
##merge with annots table
pfam_plot_results_dt <- merge(pfam_results_dt, pfamannot_table, by.x="ID", by.y="Pfam_ID", all.x=TRUE)
pfam_plot_results_dt$total_DEG_PFAMs <- tstrsplit(pfam_plot_results_dt$GeneRatio, "/", keep=c(2))
pfam_plot_results_dt$total_DEG_PFAMs <- as.numeric(pfam_plot_results_dt$total_DEG_PFAMs)
pfam_plot_results_dt$GeneRatio <- (pfam_plot_results_dt$Count/pfam_plot_results_dt$total_DEG_PFAMs)
fwrite(pfam_plot_results_dt, "output/dexseq/dtu_exposed_analysis/exposed_Pfam_enrichment.csv")
