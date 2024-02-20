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
library(pheatmap)
library(viridis)
library(tibble)

###########
# GLOBALS #
###########

trinotate_file <- snakemake@input[["trinotate_file"]]
dds_file <- snakemake@input[["dds_file"]]

##########
## DEGs ##
##########

trinotate <- fread(trinotate_file, na.strings = "")

# dexseq results
dxd <- readRDS(dds_file)
dxr <- DEXSeqResults(dxd)

# save full list for FGSEA
dxr.sorted = as.data.table(dxr[order(dxr$padj),])
fwrite(dxr.sorted, snakemake@output[["all_res"]])

# sig degs
sig_exposed_dtu <- subset(dxr.sorted, padj<0.05)
length(unique(sig_exposed_dtu$groupID)) # 154 genes with DEU
sig_exposed_dtu$gene_exon <- paste(sig_exposed_dtu$groupID, sig_exposed_dtu$featureID, sep=":")


# merge with trinotate annotations
sig_dtu_annots <- merge(sig_exposed_dtu, trinotate, by.x="groupID", by.y="#gene_id")
fwrite(sig_dtu_annots, snakemake@output[["sig_dtu_annots"]])

sig_annots_bx <- subset(sig_dtu_annots, !is.na(sprot_Top_BLASTX_hit))
sig_annots_bp <- subset(sig_dtu_annots, !is.na(sprot_Top_BLASTP_hit))
sig_blast_annots <- full_join(sig_annots_bx, sig_annots_bp)
length(unique(sig_blast_annots$groupID)) # 53 blast annotated genes with DEU

##blast for unann DEGs?!

plotDEXSeq(dxr, "ASW_TRINITY_DN4192_c0_g1", legend=TRUE)

heatmap <- subset(sig_exposed_dtu, gene_exon=="ASW_TRINITY_DN4192_c0_g1:E006")

#############
## heatmap ## first half of vst table = counts for exon (rest = counts for other exons of sample gene - this exon or other)
#############

##vst transform
asw_vst <- varianceStabilizingTransformation(dxd, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% sig_dtu_annots$gene_exon)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##second half of columns are reads aligning to rest of gene - don't need that
asw_vst_degs <- asw_vst_degs[,-c(13:24)]

##sample no.s and row no.s don't add up

##reorder for plot
#asw_vst_degs <- asw_vst_degs[,c(13:18, 1:6, 19:24, 7:12)]

##get treatment label info - careful here - need to check, can't just copy paste as rownames!=sample names
sample_to_treatment <- distinct(data.table(data.frame(colData(dxd)[,c("Treatment", "sample")])))
sample_to_treatment$sample <- seq.int(nrow(sample_to_treatment))
sample_to_treatment <- sample_to_treatment %>% remove_rownames %>% column_to_rownames(var="sample")
treatment_colours <- list(Treatment=c(exposed="#3B0F70FF", ex_control="#FCFDBFFF"))

##plot
# not clustered by sample
pdf(snakemake@output[["dtu_heatmap_unclustered"]])
pheatmap(asw_vst_degs, cluster_rows=T, cluster_cols=F, show_rownames=F,
         annotation_col=sample_to_treatment, annotation_colors=treatment_colours, annotation_names_col=F,
         show_colnames = F, border_color=NA, color=viridis(50))
dev.off()

# clustered by sample
pdf(snakemake@output[["dtu_heatmap_clustered"]])
pheatmap(asw_vst_degs, cluster_rows=T, cluster_cols=T, show_rownames=T,
         annotation_col=sample_to_treatment, annotation_colors=treatment_colours, annotation_names_col=F,
         show_colnames = F, border_color=NA, color=viridis(50))
dev.off()

# write log
sessionInfo()

