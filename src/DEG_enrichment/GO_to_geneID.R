library(data.table)
library(tidyverse)

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
##fairly sure we always used Pfam GO annots for this
go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))
fwrite(go_annot_table, "output/dexseq/enrichment/GO_annots.csv")

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

##extract GO terms for isoforms of genes
go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
##table of GO term to gene name
term_to_gene <- go_term_table[,c(2,1)]
term_to_name <- go_annot_table[,c(1,3)]
fwrite(term_to_gene, "output/dexseq/enrichment/GO_to_geneID.csv")
fwrite(term_to_name, "output/dexseq/enrichment/GO_to_GOname.csv")
