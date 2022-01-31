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

sample_info_file <- snakemake@input[["sample_info_file"]]
full_info_file <- snakemake@input[["full_info_file"]]
counts_info_file <- snakemake@input[["counts_info_file"]]
flatfile <- snakemake@input[["flatfile"]]

#sample_info_file = "data/sample_table_exposed.txt"
#full_info_file = "data/full_sample_table.csv"
#counts_info_file = "output/dtu_exposed/dtu_exposed_.sample_table"
#flatfile = "output/supertranscripts/trinity_genes.gtf.dexseq.gff"

########
# MAIN #
########

sample_info <- fread(sample_info_file, header=F)
full_info <- fread(full_info_file)
counts_info <- fread(counts_info_file, header=T)
##merge sample table with extra sample info
full_sample_table <- merge(sample_info, full_info, by.x="V3", by.y="path_to_reads1", all.x=TRUE, all.y=FALSE)
##merge with counts table info
counts_sample_table <- merge(full_sample_table, counts_info, by.x="V2", by.y="V1")
#simplify
sample_table <- counts_sample_table[,c(1,3,17,6,9,14)]
setnames(sample_table, old=c("V1", "V2"), new=c("condition", "sample"))
##full counts path
sample_table$counts_filename <- paste("output/dtu_exposed/", sample_table$counts_filename, sep="")

samples_info = sample_table
dxd = DEXSeqDataSetFromHTSeq(as.vector(samples_info$counts_filename), sampleData=samples_info, flattenedfile=flatfile)
saveRDS(dxd, snakemake@output[["dds"]])

# write log
sessionInfo()
