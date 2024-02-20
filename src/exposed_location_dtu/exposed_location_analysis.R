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
trinotate_file <- snakemake@input[["trinotate_file"]]
asw_location <- snakemake@params[["location"]]

threads <- snakemake@params[["bp_threads"]]
# for faster runtime
BPPARAM = MulticoreParam(threads, progressbar=T)

########
# MAIN #
########

#trinotate annots
trinotate <- fread(trinotate_file, na.strings="")
# dxd object
dxd <- readRDS(dds_file)

ToFilter <- apply(counts(dxd), 1, function(x) sum(x > 10)) >=3
table(ToFilter)
dxd <- dxd[ToFilter,]
paste("filtered exons")

# subset for location (snakemake should run through both)
dxd_location <- dxd[,dxd$Location==asw_location]
# factor and design
dxd_location$condition <- factor(dxd_location$condition)
dxd_location$PC1_sign <- factor(paste(dxd_location$PC1_sign))
design(dxd_location) <- ~sample+exon+PC1_sign:exon+condition:exon

# run dexseq - already run size factors
dxd_location <- estimateDispersions(dxd_location, BPPARAM=BPPARAM, quiet=F)
paste("est disp done")
dxd_location <- testForDEU(dxd_location, BPPARAM=BPPARAM, reducedModel =~sample+exon+PC1_sign:exon)
paste("test DEU done")
dxd_location <- estimateExonFoldChanges(dxd_location, fitExpToVar="condition", BPPARAM=BPPARAM)
paste("est FCs done")
saveRDS(dxd_location, snakemake@output[["dds_res"]])
paste("dexseq analysis done")

dxdr1_res = DEXSeqResults(dxd_location)
##save full list for FGSEA
dxdr1.sorted = as.data.table(dxdr1_res[order(dxdr1_res$padj),])
fwrite(dxdr1.sorted, snakemake@output[["all_res"]])
paste("dexseq results written")

##sig dtus
sig_location_dtu <- subset(dxdr1.sorted, padj<0.05)

#how many genes have DE exons?
length(unique(sig_location_dtu$groupID))

##merge sig dtus with trinotate annotations
sig_dtu_annots <- merge(sig_location_dtu, trinotate, by.x="groupID", by.y="#gene_id")
fwrite(sig_dtu_annots, snakemake@output[["sig_dtu_annots"]])
paste("sig results written")

##plot dexseq plots for top 50
pdf(snakemake@output[["dtu_plots"]])
top_genes = unique(dxdr1.sorted$groupID[dxdr1.sorted$padj < 0.1 & ! is.na(dxdr1.sorted$padj)])
top_genes = top_genes[1:min(50, length(top_genes))]
message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
for (gene in top_genes) { 
  plotDEXSeq( dxdr1_res , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
}
dev.off()
paste("plots made")

# write log
sessionInfo()
