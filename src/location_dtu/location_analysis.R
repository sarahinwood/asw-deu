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

# for faster runtime
BPPARAM=BatchJobsParam(1)

########
# MAIN #
########

##dxd object
dxd <- readRDS(dds_file)
##trinotate annotations
trinotate <- fread(trinotate_file, na.strings = "")
paste("read in dds and trinotate")

##filter down to less files to see if error remains
dxd <- dxd[,dxd$PC1_sign=="positive"]
paste("subset to positive PC only")

##factors and design
dxd$Location <- factor(paste(dxd$Location))
#dxd$sample <- factor(paste(dxd$sample))
#dxd$exon <- factor(paste(dxd$exon))
paste("dxd factors")

#filter out low counts
#keep <- rowSums(counts(dxd)) >= 30
#dxd <- dxd[keep,]
#paste("filter out low counts")

##control location PC1, test exposure - because we control for sample, do we even need to control for location/PC1? - don't think so
##https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#4_Standard_analysis_workflow
##all observations from same sample are also from same condition e.g. location - in previous analyses haven't included sample in design
design(dxd) <- ~sample+exon+Location:exon

##can try design = ~sample+exon+Location:exon+Attack:exon in future analysis with reduced = ~sample+exon+Location:exon - does this look for differences I care about?

paste("design")

##run dexseq - already run size factors
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM, quiet=F) ##long step - crashing with atomic vector error
paste("estimateDispersions complete")
saveRDS(dxd, "output/dexseq/dtu_location_analysis/asw_estdisp.dds")

plotDispEsts(dxd)
paste("plotDispEsts complete")

dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
paste("testForDEU complete")
saveRDS(dxd, "output/dexseq/dtu_location_analysis/asw_testdeu.dds")

dxd <- estimateExonFoldChanges(dxd, fitExpToVar="Location", BPPARAM=BPPARAM)
paste("estimateExonFoldChanges complete")
saveRDS(dxd, snakemake@output[["dds_res"]])
saveRDS(dxd, "output/dexseq/dtu_location_analysis/asw_foldchanges.dds")

dxdr1_res = DEXSeqResults(dxd)
##save full list for FGSEA
dxdr1.sorted = as.data.table(dxdr1_res[order(dxdr1_res$padj),])
fwrite(dxdr1.sorted, "output/dexseq/dtu_location_analysis/all_res.csv")

##sig degs
sig_location_dtu <- subset(dxdr1.sorted, padj<0.05)
fwrite(sig_location_dtu, "output/dexseq/dtu_location_analysis/sig_degs.csv")

#how many genes have DE exons?
length(unique(sig_location_dtu$groupID))

##merge with trinotate annotations
sig_dtu_annots <- merge(sig_location_dtu, trinotate, by.x="groupID", by.y="#gene_id")
fwrite(sig_dtu_annots, "output/dexseq/dtu_location_analysis/sig_degs_annots.csv")

##what number have blastx?
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
##what number have blastp?
sum(!is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
##what number have pfam GO?
sum(!is.na(sig_dtu_annots$gene_ontology_Pfam))

##plot dexseq plots for top 50
#pdf("output/dexseq/dtu_location_analysis/dtu_location_dexseq.pdf")
#top_genes = unique(dxdr1.sorted$groupID[dxdr1.sorted$padj < 0.1 & ! is.na(dxdr1.sorted$padj)])
#top_genes = top_genes[1:min(50, length(top_genes))]
#message("Top 50 genes: (", paste(top_genes, collapse=','), ")")
#for (gene in top_genes) { 
#  plotDEXSeq( dxdr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
#}

# write log
sessionInfo()
