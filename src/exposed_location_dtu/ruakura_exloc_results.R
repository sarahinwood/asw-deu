library(data.table)
library(DEXSeq)

dxd <- readRDS("output/dexseq/dtu_exposed_location_analysis/asw_exposed_location_res.dds")
dxdr1_res = DEXSeqResults(dxd)

# design was ~sample+exon+group:exon
# where group = location_treatment