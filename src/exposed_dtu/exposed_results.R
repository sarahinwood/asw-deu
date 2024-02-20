library(data.table)
library(dplyr)

sig_dtu_annots <- fread("output/dexseq/dtu_exposed_analysis_le/sig_degs_annots.csv", na.strings = "")
# no exons DU
length(unique(sig_dtu_annots$featureID))
# no genes DTU exons in
length(unique(sig_dtu_annots$groupID))

exposed_deseq <- fread("data/asw-rnaseq-pc1/output/triple_deseq2/ASW/exposure_samples/exposure/exposure_res_group.csv")

deseq_dexseq <- subset(exposed_deseq, rn %in% sig_dtu_annots$groupID)


sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))

annot_degs_bx <- subset(sig_dtu_annots, !is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
annot_degs_bp <- subset(sig_dtu_annots, !is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
annot_degs <- full_join(annot_degs_bx, annot_degs_bp)
length(unique(annot_degs$groupID))

fwrite(annot_degs, "output/dexseq/dtu_exposed_analysis/annotated_degs.csv")

dxd <- readRDS("output/dexseq/dtu_exposed_analysis/asw_exposed_res.dds")
