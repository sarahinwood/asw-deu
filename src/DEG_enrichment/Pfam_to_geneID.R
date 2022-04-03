library(data.table)
library(tidyverse)

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
##fairly sure we always used Pfam GO annots for this
pfam_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(Pfam, "`")))])
pfam_annot_table <- pfam_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
pfam_annot_table <- pfam_annot_table[,c(1,2,3)]
pfam_annot_table<-setnames(pfam_annot_table, old=c("V1", "V2", "V3"), new=c("Pfam_ID", "domain", "name"))
fwrite(distinct(pfam_annot_table), "output/dexseq/enrichment/pfam_annots.csv")

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_pfam_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

##extract GO terms for isoforms of genes
gene_ids <- trinotate_report[!is.na(Pfam), unique(`#gene_id`)]
pfam_term_list <- lapply(gene_ids, EXTRACT_pfam_TERMS, trinotate_report=trinotate_report)
pfam_term_table <- rbindlist(pfam_term_list)
##table of GO term to gene name
term_to_gene <- pfam_term_table[,c(2,1)]
term_to_name <- pfam_annot_table[,c(1,3)]
fwrite(distinct(term_to_gene), "output/dexseq/enrichment/pfam_to_geneID.csv")
fwrite(distinct(term_to_name), "output/dexseq/enrichment/pfam_to_domain_name.csv")
