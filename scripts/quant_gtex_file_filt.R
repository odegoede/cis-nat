#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

.libPaths(c(.libPaths(), "/home/odegoede/R/x86_64-pc-linux-gnu-library/3.5/"))
options(stringsAsFactors=F)

library(data.table)
source("scripts/source_temp_unzip.R")

out_name <- paste0("scratch/filt_indiv_quant/", gsub("data/indiv_quant_files/", "", gsub("\\.txt\\.gz", "", args[1])), ".txt")

ind_quant <- temp_unzip(args[1], fread, data.table = F)
colnames(ind_quant) <- c("chr", "position", "coverage", "editedreads", "editlevel")
ind_quant$orig_site_name <- paste(ind_quant$chr, ind_quant$position, sep = "_")
rownames(ind_quant) <- ind_quant$orig_site_name

to_keep <- read.table(args[2])
to_keep <- as.character(to_keep$V1)

filt_quant <- ind_quant[to_keep, ]

write.table(filt_quant, file = out_name, quote = F, row.names = F, col.names = T, sep = "\t")

