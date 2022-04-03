library(data.table)
library(dplyr)
library(DEXSeq)
library(pheatmap)
library(viridis)
library(tibble)

sig_dtu <- fread("output/dexseq/dtu_exposed_location_analysis/sig_degs_annots.csv", na.strings="")
sig_dtu$gene_exon <- paste(sig_dtu$groupID, sig_dtu$featureID, sep=":")

length(unique(sig_dtu$groupID))

sum(!is.na(sig_dtu$sprot_Top_BLASTX_hit))

annot_degs_bx <- subset(sig_dtu, !is.na(sig_dtu$sprot_Top_BLASTX_hit))
annot_degs_bp <- subset(sig_dtu, !is.na(sig_dtu$sprot_Top_BLASTP_hit))
annot_degs <- full_join(annot_degs_bx, annot_degs_bp)
length(unique(annot_degs$groupID))

fwrite(annot_degs, "output/dexseq/dtu_exposed_location_analysis/annotated_degs.csv")

# plots
dxd <- readRDS("output/dexseq/dtu_exposed_location_analysis/asw_exposed_location_res.dds")
dxdr1_res = DEXSeqResults(dxd)
plotDEXSeq(dxdr1_res, "ASW_TRINITY_DN5643_c0_g1", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)



#############
## heatmap ##
#############
##vst transform
asw_vst <- varianceStabilizingTransformation(dxd, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% sig_dtu$gene_exon)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##second half of columns are reads aligning ti rest of gene - don't need that
asw_vst_degs <- asw_vst_degs[,-c(25:48)]

##reorder for plot
asw_vst_degs <- asw_vst_degs[,c(13:18, 1:6, 19:24, 7:12)]

##get tissue label info
sample_to_location <- data.table(data.frame(colData(dxd)[,c("Treatment", "Location", "sample", "exon")]))
sample_to_location <- subset(sample_to_location, exon=="this")
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample")
sample_to_location <- sample_to_location[,-c(3)]


location_colours <- list(Treatment=c(exposed="#3B0F70FF", ex_control="#FCFDBFFF"), Location=c(Dunedin="#DE4968FF", Ruakura="#FE9F6DFF"))
##plot
# not clustered by sample
pheatmap(asw_vst_degs, cluster_rows=T, cluster_cols=F, show_rownames=F,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=F,
         show_colnames = F, border_color=NA, color=viridis(50))

# clustered by sample
pheatmap(asw_vst_degs, cluster_rows=T, cluster_cols=T, show_rownames=F,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=F,
         show_colnames = F, border_color=NA, color=viridis(50))
