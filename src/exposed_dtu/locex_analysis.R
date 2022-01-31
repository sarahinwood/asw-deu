library(data.table)
library(DEXSeq)

##dxd object
dxd <- readRDS("output/dexseq/dtu_exposed/dtu_exposed.dds")
##trinotate annotations
trinotate <- fread("data/trinotate_annots.csv", na.strings = "")

##factors and design
dxd$PC1_sign <- factor(paste(dxd$PC1_sign))
##control location PC1, test exposure
design(dxd) <- ~PC1_sign+sample+exon+condition:exon

##run dexseq
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )

##save full list for FGSEA
dxr1.sorted = dxr1[order(dxr1$padj),]
fwrite(dxr1.sorted, "output/dexseq/dtu_exposed/loc_ex/all_res.csv")

##sig degs
sig_exposed_dtu <- subset(dxr1.sorted, padj<0.05)
frwite(sig_exposed_dtu, "output/dexseq/dtu_exposed/loc_ex/sig_degs.csv")

#how many genes have DE exons? - 806
length(unique(sig_exposed_dtu$groupID))

##reduce results table
sig_dtu <- sig_exposed_dtu[,c(1:8,13:15,45)]
##merge with trinotate annotations
sig_dtu_annots <- merge(sig_dtu, trinotate, by.x="groupID", by.y="#gene_id")
frwite(sig_dtu_annots, "output/dexseq/dtu_exposed/loc_ex/sig_degs_annots.csv")

##what number have blastx - 477
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
##what number have blastp - 459
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
##what number have pfam GO - 328
sum(!is.na(sig_dtu_annots$gene_ontology_Pfam))

##plot dexseq plots for top 50
pdf("output/dexseq/dtu_exposed/loc_ex/dtu_exposed_.dexseq.pdf")
top_genes = unique(dxr1.sorted$groupID[dxr1.sorted$padj < 0.1 & ! is.na(dxr1.sorted$padj)])
top_genes = top_genes[1:min(50, length(top_genes))]
message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
for (gene in top_genes) { 
  plotDEXSeq( dxr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
}
