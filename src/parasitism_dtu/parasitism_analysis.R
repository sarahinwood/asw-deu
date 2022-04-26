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
library(BiocParallel)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]
trinotate_file <- snakemake@input[["trinotate file"]]

threads <- snakemake@params[["bp_threads"]]
# for faster runtime
BPPARAM = MulticoreParam(threads, progressbar=T)

########
# MAIN #
########

##dxd object
dxd <- readRDS(dds_file)

# filter lowly expressed exons to keep exons with counts above 10 in more than 25 samples
ToFilter <- apply(counts(dxd), 1, function(x) sum(x > 10)) >= 25
table(ToFilter)
dxd <- dxd[ToFilter,]


##factors and design
dxd$parasitism <- factor(paste(dxd$parasitism))
design(dxd) <- ~sample+exon+parasitism:exon

##run dexseq - already run size factors
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM, quiet=F)
dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="parasitism", BPPARAM=BPPARAM)
saveRDS(dxd, snakemake@output[["dds_res"]])

dxdr1_res = DEXSeqResults(dxd)
##save full list for FGSEA
dxdr1.sorted = as.data.table(dxdr1_res[order(dxdr1_res$padj),])
fwrite(dxdr1.sorted, snakemake@output[["all_res"]])

##sig dtus
sig_location_dtu <- subset(dxdr1.sorted, padj<0.05)

#how many genes have DE exons?
length(unique(sig_location_dtu$groupID))

##merge sig dtus with trinotate annotations
trinotate <- fread(trinotate_file)
sig_dtu_annots <- merge(sig_location_dtu, trinotate, by.x="groupID", by.y="#gene_id")
fwrite(sig_dtu_annots, snakemake@output[["sig_dtu_annots"]])

##plot dexseq plots for top 50
pdf(snakemake@output[["dtu_plots"]])
top_genes = unique(dxdr1.sorted$groupID[dxdr1.sorted$padj < 0.1 & ! is.na(dxdr1.sorted$padj)])
top_genes = top_genes[1:min(50, length(top_genes))]
message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
for (gene in top_genes) { 
  plotDEXSeq( dxdr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
}
dev.off()

# write log
sessionInfo()
