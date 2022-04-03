library(data.table)
library(dplyr)

sig_dtu <- fread("output/dexseq/dtu_exposed_analysis/sig_degs.csv", na.strings="")

sig_dtu_annots <- fread("output/dexseq/dtu_exposed_analysis/sig_degs_annots.csv", na.strings = "")
length(unique(sig_dtu_annots$groupID))

sum(!is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))

annot_degs_bx <- subset(sig_dtu_annots, !is.na(sig_dtu_annots$sprot_Top_BLASTX_hit))
annot_degs_bp <- subset(sig_dtu_annots, !is.na(sig_dtu_annots$sprot_Top_BLASTP_hit))
annot_degs <- full_join(annot_degs_bx, annot_degs_bp)
length(unique(annot_degs$groupID))

fwrite(annot_degs, "output/dexseq/dtu_exposed_analysis/annotated_degs.csv")

dxd <- readRDS("output/dexseq/dtu_exposed_analysis/asw_exposed_res.dds")
