#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 01: Identifying all overlapping transcripts in human genome
## By: Olivia de Goede, 2021
#####

####
## Set options and load required libraries
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
source("source_temp_unzip.R")



####
## Assign arguments
option_list <- list( 
  make_option("--gtf", default = NULL, 
              help = "GTF file [default \"%default\"]"),
  make_option("--outdir", default = "./", 
              help = "Output directory to write transcript overlaps to [default is current working directory]"),
  make_option("--keepanno", default = "neither", 
              help = "Should the gene, transcript, and exon annotations be saved? Files will be saved to outdir. \"text\" for .txt files, \"RData\" for .RData, \"both\" for both .txt and .RData files. Anything else will be interpreted as \"neither\". [default %default]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))



####
## Basic first test of inputs
# Check the required arguments (GTF file) is provided
if (is.null(opt$gtf)) { 
  stop("GTF file not provided, exiting\n") 
}

# Check that the output directory exists
if (!dir.exists(opt$outdir)) {
  stop("Output directory does not exist, exiting\n")
}
# For saving files, create variable OUTDIR that ends in a /
if (endsWith(opt$outdir, "/")) {
  OUTDIR = opt$outdir
} else {
  OUTDIR = paste0(opt$outdir, "/")
}


####
## Read in and organize GTF
if (endsWith(opt$gtf, "gz") | endsWith(opt$gtf, "zip") | endsWith(opt$gtf, "bzip2") | endsWith(opt$gtf, "xz")) {
  gtf_anno <- suppressWarnings(temp_unzip(opt$gtf, fread, data.table = F))
} else {
  gtf_anno <- suppressWarnings(fread(opt$gtf, data.table = F))
}

# Test: check that fields are as expected
# all of column 1 should start with "chr"; column 3 should include c("gene", "transcript", "exon"); 
# columns 4 and 5 should be integers, and col4 should be <= col5; column 7 should all be "+" or "-";
# column 9 contains a bunch of stuff, but a basic check is "gene_id" should always be there
if (!all(startsWith(gtf_anno[,1], "chr") == T)) {
  stop("Error in GTF: not all entries in column 1 (seqname) are chr--")
} else if (!(any(gtf_anno[,3] == "gene") & any(gtf_anno[,3] == "transcript") & any(gtf_anno[,3] == "exon"))) {
  stop("Error in GTF: missing either gene, transcript, or exon level information in column 3 (type)")
} else if (!(is.integer(gtf_anno[,4]) & is.integer(gtf_anno[,5]))) {
  stop("Error in GTF: columns 4 (start) and 5 (end) are not integers")
} else if (!all(gtf_anno[,4] <= gtf_anno[,5])) {
  stop("Error in GTF: column 4 (start) is not always <= column 5 (end)")
} else if (!all(gtf_anno[,7] %in% c("+", "-"))) {
  stop("Error in GTF: column 7 is not defining strand")
} else if (!all(grepl("gene_id", gtf_anno[,9]))) {
  stop("Error in GTF: column 9 (attribute) is missing gene_id, the minimum info")
}

# Set up column names, remove mitochondrial chromosome
gtf_anno <- gtf_anno[,-c(2,6,8)] # these are unnecessary fields
colnames(gtf_anno) <- c("chr", "type", "start", "end", "strand", "attribute")
gtf_anno <- gtf_anno[gtf_anno$chr != "ChrM", ]



####
## Make the gene_anno, transcript_anno, and exon_anno files
# add gene information from attribute field (all 3 anno files need this)
gtf_anno$attribute <- gsub("; ", ";", as.character(gtf_anno$attribute))
gtf_anno$ensgene <- gsub("\"", "", gsub("gene_id ", "", unlist(lapply(strsplit(gtf_anno$attribute, ";"), grep, pattern = "gene_id", value = T))))
gtf_anno$symbol <- gsub("\"", "", gsub("gene_name ", "", unlist(lapply(strsplit(gtf_anno$attribute, ";"), grep, pattern = "gene_name", value = T))))
gtf_anno$biotype <- gsub("\"", "", gsub("gene_type ", "", unlist(lapply(strsplit(gtf_anno$attribute, ";"), grep, pattern = "gene_type", value = T))))
gtf_anno$gencode_level <- gsub("\"", "", gsub(" ", "_", unlist(lapply(strsplit(gtf_anno$attribute, ";"), grep, pattern = "^level ", value = T))))

# separate gene, transcript, and exon
gene_anno <- gtf_anno[gtf_anno$type == "gene", ]
gene_anno <- gene_anno[,-grep("attribute", colnames(gene_anno))]
rownames(gene_anno) <- gene_anno$ensgene
tx_anno <- gtf_anno[gtf_anno$type == "transcript", ]
exon_anno <- gtf_anno[gtf_anno$type == "exon", ]

# for transcript_anno, add transcript information from attribute field
tx_anno$tx_id <- gsub("\"", "", gsub("transcript_id ", "", unlist(lapply(strsplit(tx_anno$attribute, ";"), grep, pattern = "transcript_id", value = T))))
tx_anno$tx_type <- gsub("\"", "", gsub("transcript_type ", "", unlist(lapply(strsplit(tx_anno$attribute, ";"), grep, pattern = "transcript_type", value = T))))
# for transcript support level, not all transcripts have TSL values in their attribute field: set as NA first
tx_anno$tsl <- NA
tx_anno[grepl("transcript_support_level", tx_anno$attribute), ]$tsl <- gsub("\"", "", gsub("transcript_support_level ", "tsl_", unlist(lapply(strsplit(tx_anno[grepl("transcript_support_level", tx_anno$attribute), ]$attribute, ";"), grep, pattern = "transcript_support_level", value = T))))
tx_anno[is.na(tx_anno$tsl), ]$tsl <- "tsl_NA"
tx_anno <- tx_anno[,-grep("attribute", colnames(tx_anno))]
rownames(tx_anno) <- tx_anno$tx_id

# for exon_anno, add transcript information and exon number from attribute field
exon_anno$tx_id <- gsub("\"", "", gsub("transcript_id ", "", unlist(lapply(strsplit(exon_anno$attribute, ";"), grep, pattern = "transcript_id", value = T))))
exon_anno$tx_type <- gsub("\"", "", gsub("transcript_type ", "", unlist(lapply(strsplit(exon_anno$attribute, ";"), grep, pattern = "transcript_type", value = T))))
exon_anno$exon_number <- as.numeric(gsub("exon_number ", "", unlist(lapply(strsplit(exon_anno$attribute, ";"), grep, pattern = "exon_number", value = T))))
exon_anno <- exon_anno[,-grep("attribute", colnames(exon_anno))]

# if the user put something in --keepanno, save these files to outdir
if (opt$keepanno %in% c("text", "both")) {
  write.table(gene_anno, file = paste0(OUTDIR, "gene_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
  write.table(tx_anno, file = paste0(OUTDIR, "transcript_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
  write.table(exon_anno, file = paste0(OUTDIR, "exon_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
}
if (opt$keepanno %in% c("RData", "both")) {
  save(gene_anno, tx_anno, exon_anno, file = paste0(OUTDIR, "gene_tx_exon_anno_files.RData"))
}



####
## Find overlaps! 

## Criteria for overlaps:
# 1. transcripts come from opposite strands
# 2. transcripts either have good transcript support level scores (TSL 1-3), or are the transcripts
#    with the most evidence for their gene
# 3. have at least 100 bp *continuous* overlap


# Filter transcripts based on TSL


# Separate into plus and minus strand
# (overlaps need to be on opposite strands)


# Check for overlaps


# Check for any overlapping transcript pairs that share start&end coordinates in an exon
# (this is rare, but this means that once processed the transcripts could have continuous overlaps 
# across multiple exons, which would need to be summed)


# If any multi-exon overlaps, flag them


# Organize overlaps into a data.frame
# (e.g. remove any of the exact same overlap between different isoforms of same genes)


# Combine entries for the multi-exon overlaps


# Flag duplicate overlaps, filter down to unique overlap regions
# duplicate_overlap field: marks overlaps involving alternative transcripts of the same gene that 
# have the exact same overlapping region


# Save overlap files
