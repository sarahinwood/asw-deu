library(data.table)
library(fgsea)
library(ggplot2)
library(viridis)
library(stringr)

set.seed(10)

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
##fairly sure we always used Pfam GO annots for this
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
res_group <- fread("output/dexseq/dtu_exposed_analysis/all_res.csv")

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
setorder(res_group, stat)
ranks <- res_group[!is.na(stat), stat]
names(ranks) <- res_group[!is.na(stat), groupID]

##old fgsea needed nperm=, but appears new version does not
fgsea_res <- fgsea(pathways, ranks)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sorted_fgsea_res_no_na <- sorted_fgsea_res[!is.na(padj)]
sum(sorted_fgsea_res_no_na$padj<0.05)

sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
##need number of genes in leading edge for lollipop plot
annot_sig_fgsea$leadingEdge_size <- str_count(annot_sig_fgsea$leadingEdge, "TRINITY_DN")
fwrite(annot_sig_fgsea, "output/deseq2/stage_WT/sig_GO_enrichment.csv")

#####################
## all in one plot ##
#####################
###swap _ for space in pathway kind
annot_sig_fgsea$pathway_kind <- gsub("_", " ", annot_sig_fgsea$pathway_kind)
annot_sig_fgsea$pathway_name <- tstrsplit(annot_sig_fgsea$pathway_name, ", ", keep=c(1))
##reorder - sorts by pathway_kind reverse alphabetically but can't figure out how to do any better
annot_sig_fgsea$pathway_name <- factor(annot_sig_fgsea$pathway_name, levels=annot_sig_fgsea$pathway_name[order(annot_sig_fgsea$pathway_kind, annot_sig_fgsea$NES)])
##plot
ggplot(annot_sig_fgsea, aes(pathway_name, NES)) +
  geom_segment(aes(y=0, yend=annot_sig_fgsea$NES, x=pathway_name, xend=pathway_name), alpha=0.4)+
  geom_point(aes(colour=pathway_kind, size=leadingEdge_size)) + 
  labs(x="Gene ontology terms", y="FGSEA normalized enrichment score",
       colour="GO domain", size="Leading\nedge size") +
  ylab("FGSEA normalized enrichment score")+
  coord_flip() +
  scale_colour_viridis(discrete=TRUE)+
  theme_bw()

##leading edge
memb_res <- annot_sig_fgsea[annot_sig_fgsea$pathway == "GO:0016021",]
memb_res_leading_edge <- data.frame(memb_res$leadingEdge)
setnames(memb_res_leading_edge, old=c("c..TRINITY_DN21376_c1_g1....TRINITY_DN2656_c0_g1....TRINITY_DN15880_c0_g1..."), new=c("gene_id"))
memb_leading_annots <- merge(memb_res_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(memb_leading_annots, "output/deseq2/asw/interaction/interaction/fgsea/membrane_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0016021"]], ranks) + labs(title="integral component of membrane")
