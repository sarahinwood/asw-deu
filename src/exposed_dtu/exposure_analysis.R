#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(DEXSeq)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds"]]
trinotate_file <- snakemake@input[["trinotate"]]

########
# MAIN #
########

##dxd object
dxd <- readRDS(dds_file)
##trinotate annotations
trinotate <- fread(trinotate_file, na.strings = "")

##factors and design
dxd$Location <- factor(paste(dxd$Location))
dxd$PC1_sign <- factor(paste(dxd$PC1_sign))
dxd$Treatment <- factor(paste(dxd$Treatment))
##control location PC1, test exposure
design(dxd) <- ~Location+PC1_sign+sample+exon+Treatment:exon

##run dexseq on 6/8 cores
dxdr1 <- DEXSeq(dxd, fullModel=design(dxd), BPPARAM=MulticoreParam(workers=1), fitExpToVar="Treatment")
saveRDS(dxdr1, snakemake@output[["dds_res"]])

##run dexseq
#dxd = estimateSizeFactors( dxd )
#dxd = estimateDispersions( dxd )
#plotDispEsts( dxd )
#dxd = testForDEU( dxd )
#dxd = estimateExonFoldChanges( dxd, fitExpToVar="Treatment")
#dxr1 = DEXSeqResults( dxd )

##save full list for FGSEA
dxr1.sorted = dxr1[order(dxr1$padj),]
fwrite(dxr1.sorted, "output/dexseq/dtu_exposed/exposure/all_res.csv")

##sig degs
sig_exposed_dtu <- subset(dxr1.sorted, padj<0.05)
frwite(sig_exposed_dtu, snakemake@output[["sig_res"]])

#how many genes have DE exons? - 806
length(unique(sig_exposed_dtu$groupID))

##reduce results table
sig_dtu <- sig_exposed_dtu[,c(1:8,13:15,45)]
##merge with trinotate annotations
sig_dtu_annots <- merge(sig_dtu, trinotate, by.x="groupID", by.y="#gene_id")
frwite(sig_dtu_annots, snakemake@output[["sig_res_annots"]])

##what number have blastx - 477
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
##what number have blastp - 459
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
##what number have pfam GO - 328
sum(!is.na(sig_dtu_annots$gene_ontology_Pfam))

##plot dexseq plots for top 50
pdf("output/dexseq/dtu_exposed/exposure/dtu_exposed_.dexseq.pdf")
top_genes = unique(dxr1.sorted$groupID[dxr1.sorted$padj < 0.1 & ! is.na(dxr1.sorted$padj)])
top_genes = top_genes[1:min(50, length(top_genes))]
message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
for (gene in top_genes) { 
  plotDEXSeq( dxr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
}

# write log
sessionInfo()
