library(data.table)
library(VennDiagram)

##full dtu results table
exposed_dtu <- read.table("output/dtu_exposed/dtu_exposed_.dexseq.results.dat", sep="	", fill=TRUE, row.names=NULL)
##trinotate annotations
trinotate <- fread("data/trinotate_annots.csv", na.strings = "")

##subset for only sig genes
sig_exposed_dtu <- subset(exposed_dtu, padj<0.05)
#how many genes have DE exons? - 806
length(unique(sig_exposed_dtu$groupID))

##reduce results table
sig_dtu <- sig_exposed_dtu[,c(1:8,13:15,45)]
##merge with trinotate annotations
sig_dtu_annots <- merge(sig_dtu, trinotate, by.x="groupID", by.y="#gene_id")

##what number have blastx - 477
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
##what number have blastp - 459
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
##what number have pfam GO - 328
sum(!is.na(sig_dtu_annots$gene_ontology_Pfam))

##overlap with exposed DESeq analysis
exposed <- fread("data/exposed_analysis/output/deseq2/asw/all_exposed/all_sig_annots.csv")
location <- fread("data/exposed_analysis/output/deseq2/asw/all_location/all_sig_annots.csv")
interaction <- fread("data/exposed_analysis/output/deseq2/asw/interaction/all_degs_tables/all_annot_degs.csv")

#Draw Venn Diagram
vd <- venn.diagram(x = list("interaction"=interaction$rn, "dtu"=sig_dtu_annots$groupID), filename=NULL,
                   fill=c("#440154FF", "#21908CFF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd)
##exposed overlap - 0
##location overlap - 3
##interaction overlap - 0

##dtu analysis is detecting a different set of genes than full gene DE analysis was
