#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 02: Make BED file of coordinates of transcript overlaps
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/02_save_overlap_locations.R --overlap output/01_putative_cisNAT_uniqueRegionsOnly.RData --outdir output/


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--overlap", default = NULL, 
              help = "RData file of transcript overlaps [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory to write transcript overlaps to [default is current working directory]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## DEFINE FUNCTIONS

# rep_per_exon(): for multi-exon overlaps: separate out genomic location for each exon involved in overlap
rep_per_exon <- function(tx_pair, overlap_dat = bonus_rows) {
  temp <- overlap_dat[tx_pair, ]
  rep_num <- length(unlist(strsplit(temp$longest_overlap_start, ",")))
  out_df <- data.frame(chr = rep(temp$chr, rep_num),
                       longest_overlap_start = unlist(strsplit(temp$longest_overlap_start, ",")), 
                       longest_overlap_end = unlist(strsplit(temp$longest_overlap_end, ",")), 
                       minus_gene = rep(temp$minus_gene, rep_num),
                       plus_gene = rep(temp$plus_gene, rep_num), 
                       row.names = paste(tx_pair, c(1:rep_num), sep = "_"))
  out_df
}


####
## INPUT TESTS
# Check the required arguments (overlaps file) are provided
if (is.null(opt$overlap)) { 
  stop("Overlaps file not provided, exiting\n") 
}

# Create variable out_dir that is consistent (doesn't have trailing "/")
out_dir <- file.path(opt$outdir)
# Check that the directory exists
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist, exiting\n")
}


####
## LOAD SOURCE SCRIPTS
# n/a


####
## READ IN INPUT FILES
load(file.path(opt$overlap)) # object name from script 01: putative_cisnat


####
## Merge putative_cisnat into simplified file of genomic locations of any overlapping regions
# (to identify priority regions for read alignment)
# fields: chr, start, end, Alu_flag

# note: positions are based on GENCODE GTFs (1-based); will be saved as 1-start, fully closed
# these saved coordinates are just based on overlapping regions, not extended to the rest of the transcript

overlap_loc <- putative_cisnat[-grep(",", putative_cisnat$longest_overlap_start), c("chr", "longest_overlap_start", "longest_overlap_end", "minus_gene", "plus_gene")]

# pull out the multi-exon overlaps (where genome coordinates of overlap are a little trickier)
bonus_rows <- putative_cisnat[grep(",", putative_cisnat$longest_overlap_start), c("chr", "longest_overlap_start", "longest_overlap_end", "minus_gene", "plus_gene")]
bonus_dat <- do.call(rbind, lapply(rownames(bonus_rows), rep_per_exon))

# bring them back into the overlap_loc df, make into granges
overlap_loc <- rbind(overlap_loc, bonus_dat)
overlap_ranges <- makeGRangesFromDataFrame(overlap_loc, keep.extra.columns = T, seqnames.field = "chr",
                                           start.field = "longest_overlap_start", end.field = "longest_overlap_end")
length(overlap_ranges) # 15,989

# narrow the ranges down to combine overlaps that intersect (i.e. overlapping overlaps)
length(overlap_ranges_reduced <- reduce(overlap_ranges)) # 7,351

# if overlaps are very close together, might as well just merge - but what range?
between <- gaps(overlap_ranges_reduced)
summary(width(between) < 1000) # 861 are, 6490 aren't
summary(width(between) < 100) # 148 are, 7203 aren't
# going to combine if within 100 bases
length(overlap_ranges_reduced) # 7,351
length(between) # also 7,351

# (a bit of exploratory code with magic numbers):
which(width(between) < 100) # 142 is the first index
overlap_ranges_reduced[141:143] # the short dist is between 141 (ends at 32817689) and 142 (starts at 32817699)
# output from (which(width(between) < 100)) gives the index of the second region
# (of the two that are within 100 bp of each other)

start_row_inds <- which(width(between) < 100) - 1 # where to get start pos from
end_row_inds <- which(width(between) < 100) # where to get end pos from
to_add <- data.frame(chr = seqnames(overlap_ranges_reduced[end_row_inds]),
                     start = start(overlap_ranges_reduced[start_row_inds]),
                     end = end(overlap_ranges_reduced[end_row_inds]))
# get rid of the rows used to make "to_add", put "to_add" into the granges object
to_remove <- c(start_row_inds, end_row_inds)
overlap_ranges_final <- c(overlap_ranges_reduced[-to_remove],
                          makeGRangesFromDataFrame(to_add))
overlap_ranges_final <- sort(overlap_ranges_final)

# double check - are all original granges overlapping with this final set?
countOverlaps(overlap_ranges[1:6], overlap_ranges_final, type = "within") 
# ^ named integers, showing how many ranges in overlap_ranges_final intersect the query 9k from overlap_ranges
# so if there are no zeroes in countOverlaps(overlap_ranges, overlap_ranges_final), that means overlap_ranges_final covers everybody
length(countOverlaps(overlap_ranges, overlap_ranges_final, type = "within")) # 15,989 - good, checking each overlap from initial list
summary(countOverlaps(overlap_ranges, overlap_ranges_final, type = "within"))
# ^ great! all are 1 or 2

# Save these coordinates - will save as overlap_ranges
# coordinates are 1-based, fully-closed
rm(overlap_ranges)
overlap_ranges <- data.frame(chr = seqnames(overlap_ranges_final),
                             start = start(overlap_ranges_final),
                             end = end(overlap_ranges_final))
save(overlap_ranges, file = file.path(out_dir, "02_overlappingTranscript_genomicRanges_1start_fullClose.RData"))
write.table(overlap_ranges, file = file.path(out_dir, "02_overlappingTranscript_genomicRanges_1start_fullClose.txt"),
            quote = F, sep = "\t", col.names = T, row.names = F)

