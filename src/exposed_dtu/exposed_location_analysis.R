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
library(BatchJobs)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]
trinotate_file <- snakemake@input[["trinotate_file"]]
threads <- snakemake@params[["threads"]]

# for faster runtime
BPPARAM=BatchJobsParam(threads)

########
# MAIN #
########

##dxd object
dxd <- readRDS(dds_file)
##trinotate annotations
trinotate <- fread(trinotate_file, na.strings = "")
paste("trinotate and dxd read in")

##factors and design
dxd$group <- factor(paste(dxd$Location, dxd$Treatment, sep="_"))
##control location PC1, test exposure - because we control for sample, do we even need to control for location/PC1? - don't think so
##https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#4_Standard_analysis_workflow
##all observations from same sample are also from same condition e.g. location - in previous analyses haven't included sample in design
design(dxd) <- ~sample+exon+group:exon
##should control for location-sp dtu - do we want it to though?
paste("dxd design")

##run dexseq - already run size factors
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM, quiet=F) ##long step
paste("estdisp complete")
plotDispEsts(dxd)
dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
paste("testforDEU complete")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="group", BPPARAM=BPPARAM)
paste("estimateExonFoldChanges complete")

##pipeline command - but dxd file had no results
#dxd <- DEXSeq(dxd, fullModel=design(dxd), BPPARAM=MulticoreParam(workers=1), fitExpToVar="Treatment", quiet=F)

saveRDS(dxd, snakemake@output[["dds_res"]])

dxdr1_res = DEXSeqResults(dxd)
##save full list for FGSEA
dxdr1.sorted = as.data.table(dxdr1_res[order(dxdr1_res$padj),])
fwrite(dxdr1.sorted, "output/dexseq/dtu_exposed_location_analysis/all_res.csv")

##sig degs
sig_exposed_dtu <- subset(dxdr1.sorted, padj<0.05)
fwrite(sig_exposed_dtu, "output/dexseq/dtu_exposed_location_analysis/sig_degs.csv")

#how many genes have DE exons?
length(unique(sig_exposed_dtu$groupID))

##merge with trinotate annotations
sig_dtu_annots <- merge(sig_exposed_dtu, trinotate, by.x="groupID", by.y="#gene_id")
fwrite(sig_dtu_annots, "output/dexseq/dtu_exposed_location_analysis/sig_degs_annots.csv")

##what number have blastx?
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
##what number have blastp?
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
##what number have pfam GO?
sum(!is.na(sig_dtu_annots$gene_ontology_Pfam))

##plot dexseq plots for top 50
#pdf("output/dexseq/dtu_exposed_analysis/dtu_exposed_dexseq.pdf")
#top_genes = unique(dxdr1.sorted$groupID[dxdr1.sorted$padj < 0.1 & ! is.na(dxdr1.sorted$padj)])
#top_genes = top_genes[1:min(50, length(top_genes))]
#message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
#for (gene in top_genes) { 
#  plotDEXSeq( dxdr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
#}

# write log
sessionInfo()
