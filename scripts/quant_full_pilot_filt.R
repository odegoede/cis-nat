#!/usr/bin/env Rscript

# Usage: Rscript 

args <- commandArgs(trailingOnly=TRUE)
# args[1] = bed file, args[2] = sample name

options(stringsAsFactors=F)

out_name <- paste0("scratch/pilot_indiv_quant/", args[2], ".txt")

ind_quant <- read.table(args[1])
colnames(ind_quant) <- c("chr", "start", "end", "bases", "cat", "strand", "score")
ind_quant$site_name <- paste(ind_quant$chr, ind_quant$start, ind_quant$end, sep = "_")
rownames(ind_quant) <- ind_quant$site_name

# get coverage and edit rows separate
ind_quant$coverage <- as.numeric(unlist(lapply(strsplit(ind_quant$score, ":"), "[[", 2)))
ind_quant$edited_reads <- as.numeric(unlist(lapply(strsplit(ind_quant$score, ":"), "[[", 1)))
ind_quant$editLevel <- ind_quant$edited_reads / ind_quant$coverage
# remove rows with coverage <10
filt_quant <- ind_quant[ind_quant$coverage >= 10, c("chr", "start", "end", "site_name", "coverage", "editLevel")]

# write output table
write.table(filt_quant, file = out_name, quote = F, row.names = F, col.names = T, sep = "\t")

