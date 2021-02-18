#!/usr/bin/Rscript

#####
## cis-NAT project
## Script 01: Identifying all overlapping transcripts in human genome
## By: Olivia de Goede, 2021
#####

####
## Set options and load required libraries
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
# suppressPackageStartupMessages(library(GenomicRanges))
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
# gtf_anno <- gtf_anno[,-c(2,6,8)] # these are unnecessary fields
# colnames(gtf_anno) <- c("chr", "type", "start", "end", "strand", "attribute")
# gtf_anno <- gtf_anno[gtf_anno$chr != "ChrM", ]


####
## Make (and possibly save) the gene_anno, transcript_anno, and exon_anno files



####
## Find overlaps! 

## Criteria for overlaps:
# 1. transcripts come from opposite strands
# 2. transcripts either have good transcript support level scores (TSL 1-3), or are the transcripts
#    with the most evidence for their gene
# 3. have at least 100 bp *continuous* overlap

