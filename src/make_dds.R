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

full_info_file <- snakemake@input[["full_info_file"]]
counts_info_file <- snakemake@input[["counts_info_file"]]
flatfile <- snakemake@input[["flatfile"]]
full_path <- snakemake@params[["full_path"]]

########
# MAIN #
########

# sample name and counts filepath
counts_info <- fread(counts_info_file, header=T)
sample_counts_file <- counts_info[,c(1,3)]
# full table
full_info <- fread(full_info_file)

# merge sample table with counts file-paths
sample_table <- merge(full_info, sample_counts_file, by.y="V1", by.x=ifelse(counts_info_file=="output/dtu_exposed/dtu_exposed_.sample_table", "replicate_name",
                                                                            ifelse(counts_info_file=="output/dtu_location/dtu_location_.sample_table", "location_name", "error")))
# changing names for dexseq
setnames(sample_table, old=c("Treatment", "sample_name"), new=c("condition", "sample"))
# full counts path
sample_table$counts_filename <- paste(full_path, sample_table$counts_filename, sep="")

# read in with dexseq
samples_info = sample_table
paste("made sample table for dexseq")
dxd <- DEXSeqDataSetFromHTSeq(as.vector(samples_info$counts_filename), sampleData=samples_info, flattenedfile=flatfile)
paste("made dxd object")
##estimate size factors on all data before splitting for species
dxd <- estimateSizeFactors(dxd)
paste("estimated size factors")
saveRDS(dxd, snakemake@output[["dual_dxd"]])

############################################
##subset gene table into ASW & Mh & genes ##
############################################

##ASW
asw_dxd <- dxd[grep("ASW", rownames(dxd)), ]
saveRDS(asw_dxd, snakemake@output[["asw_dxd"]])

##Mh
mh_dxd <- dxd[grep("MH", rownames(dxd)), ]
saveRDS(mh_dxd, snakemake@output[["mh_dxd"]])

# write log
sessionInfo()
