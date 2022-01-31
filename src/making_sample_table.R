library(data.table)

table <- fread("data/starting_sample_table.csv")
table$location <- tstrsplit(table$sample_name, "", keep=c(1))
table$condition_name <- paste(table$condition_name, table$location, sep="_")

table$replicate_name1 <- tstrsplit(table$sample_name, "", keep=c(3))
table$replicate_name2 <- tstrsplit(table$sample_name, "", keep=c(4))
table$replicate_name <- paste(table$replicate_name1, table$replicate_name2, sep="")

table$path_to_reads1 <- paste(table$path_to_reads1, table$sample_name, "_r1.fq.gz", sep="")
table$path_to_reads2 <- paste(table$path_to_reads2, table$sample_name, "_r2.fq.gz", sep="")
fwrite(table, "data/new_sample_table.csv")
