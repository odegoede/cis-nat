#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 03: Flag overlapping transcripts that also contain IRAlu regions
## By: Olivia de Goede, 2021
#####

# if an overlapping transcript also contains part of an inverted-repeat Alu, we will need to dive into 
# it to figure out if RNA editing is on the cis-NAT dsRNA (from the overlapping transcripts), or the 
# IRAlu dsRNA. A reasonable first check would be, is expression of both tx required for editing?

# Note: the Alu annotation input file was obtained from UCSC RepeatMasker; the "score" is perhaps the Smith-Waterman score of the match

####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--overlap", default = NULL, 
              help = "RData file of transcript overlaps [default \"%default\"]"),
  make_option("--alufile", default = NULL,
              help = "File with genomic annotation of Alu regions [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory to write transcript overlaps to [default \"%default\"]"),
  make_option("--scriptdir", default = "./scripts",
              help = "Scripts directory (to load functions saved in source scripts). [default assumes running from repo, scripts are in \"%default\"]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## DEFINE FUNCTIONS
# extend() will extend the upstream and downstream boundaries of genomic ranges
# adjust the values of upstream and downstream to whatever you want
extend <- function(gr, upstream = 700, downstream = 700){
  # subtract "upstream" from start and add "downstream" to ends
  new_start <- start(gr) - upstream
  new_end <- end(gr) + downstream
  # update the gr object
  ranges(gr) <- IRanges(new_start, new_end)
  # correct out of bounds indices
  trim(gr)
}

# combine_alus() returns alu row number and a comma separated string of alus overlapping query alu
combine_alus <- function(grouped_df, key, alus_df) {
  # grouped_df and key come from the group_map function
  #   - grouped_df = each group separately
  #   - key = key of group
  #   - alus_df = alus genomic ranges as a df
  # grab row numbers of overlapped alus for current query alu
  row_nums <- grouped_df$subjectHits
  # grab the rows from full alu df
  alu_rows <- dplyr::slice(alus_df, row_nums)
  # get the names of the corresponding alus
  alu_names <- alu_rows$name
  # combine all of the names together as comma separated list
  alu_string <- paste(alu_names, collapse = ",")
  # return the alu row index and the string of corresponding probes
  output <- data.frame(row_number = unlist(key), alus = alu_string)
  output
}


####
## INPUT TESTS
# Check the required arguments (overlaps file) are provided
if (is.null(opt$overlap)) { 
  stop("Overlaps file not provided, exiting\n") 
}
if (is.null(opt$alufile)) { 
  stop("Alu file not provided, exiting\n") 
}

# Create variable out_dir and script_dir that are consistent (doesn't have trailing "/")
out_dir <- file.path(opt$outdir)
script_dir <- file.path(opt$scriptdir)
# Check that the directories exist
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist, exiting\n")
}
if (!dir.exists(script_dir)) {
  stop("Script directory does not exist, exiting\n")
}


####
## LOAD SOURCE SCRIPTS
source(file.path(script_dir, "source_temp_unzip.R"))


####
## READ IN INPUT FILES
load(file.path(opt$overlap)) # object name from script 01: putative_cisnat
alu_file <- file.path(opt$alufile)
if (endsWith(alu_file, "gz") | endsWith(alu_file, "zip") | endsWith(alu_file, "bzip2") | endsWith(alu_file, "xz")) {
  alu_locs <- suppressWarnings(temp_unzip(alu_file, fread, data.table = data.table))
} else {
  alu_locs <- suppressWarnings(fread(alu_file, data.table = data.table))
}
colnames(alu_locs) <- c("chr", "start", "end", "name", "score", "strand")


####
## Filter alu_locs to only those that are parts of inverted pairs

# Notes from Qin: potential false negatives:
# 1) two or more Alus not in inverted, but tandem orientations; 
# 2) inverted Alus more than 700 nt apart from each other (unlikely to form stable dsRNAs).

## This code will get sets of Alus 700 nt apart
# will run twice: hits with ignore.strand = T that are NOT hits with ignore.strand = F in
# findOverlaps() are overlaps on opposite strands

# give each Alu a unique name, make into granges
alu_locs$name <- paste(alu_locs$name, c(1:nrow(alu_locs)), sep = "_")
alu_gr <- makeGRangesFromDataFrame(alu_locs, keep.extra.columns = TRUE)
# create an expanded granges object, of the Alu ranges +/-700 
expanded_alu_gr <- extend(alu_gr)

# find overlaps between alu locations and expanded alu ranges
# ignore.strand = T (all strand matches)
alu_overlaps_allStrand <- as.data.frame(findOverlaps(expanded_alu_gr, alu_gr, ignore.strand = T))
# and ignore.strand = F (same strand only)
alu_overlaps_sameStrand <- as.data.frame(findOverlaps(expanded_alu_gr, alu_gr, ignore.strand = F))
# remove rows from allStrand that are also found in sameStrand (we want opposite strand only)
# this will also remove rows that are just the Alu matching with itself
alu_overlaps_allStrand$pairHits <- paste(alu_overlaps_allStrand$queryHits, alu_overlaps_allStrand$subjectHits, sep = "_")
alu_overlaps_sameStrand$pairHits <- paste(alu_overlaps_sameStrand$queryHits, alu_overlaps_sameStrand$subjectHits, sep = "_")
summary(alu_overlaps_allStrand$pairHits %in% alu_overlaps_sameStrand$pairHits) # 615,686 FALSE (on opposite strand, are inverted Alus); 2,222,492 TRUE
alu_overlaps_df <- alu_overlaps_allStrand[!(alu_overlaps_allStrand$pairHits %in% alu_overlaps_sameStrand$pairHits), ]

# call the function to take overlap indices
# first, group by the queryHits, then pass each group into our combine_alus function
# this returns the row number for alu-of-interest and a string of alus that overlap the alu-of-interest's extension
# SLOW CODE ALERT: doesn't scale well, so broke it into a for loop
# (~7,000 queryHits of alu_overlaps_df at a time seems close to optimal)
hit_vals <- unique(alu_overlaps_df$queryHits)
n_at_a_time <- 7000
for (i in c(1:round(length(hit_vals)/n_at_a_time))) {
  start_ind <- (i-1)*n_at_a_time + 1
  end_ind <- i*n_at_a_time
  if (i == round(length(hit_vals)/n_at_a_time)) {
    end_ind <-  length(hit_vals)
  }
  query_vals <- hit_vals[start_ind:end_ind]
  temp <- group_by(alu_overlaps_df[alu_overlaps_df$queryHits %in% query_vals, ], queryHits) %>%
    group_map(combine_alus, alus_df = alu_locs) %>% reduce(rbind)
  if (i > 1) {
    collapsed_alus <- rbind(collapsed_alus, temp)
  }
  else {
    collapsed_alus <- temp
  }
}

# combine the alu data with the overlapping probes
alus_with_overlap <- as.data.frame(alu_gr) %>% dplyr::slice(collapsed_alus$row_number) %>% mutate(adjacent_alus = collapsed_alus$alus)
# make a simple chr name (some of them have more detailed annotation in seqnames, e.g. "chr22_KI270736v1_random")
alus_with_overlap$seqnames <- as.character(alus_with_overlap$seqnames)
alus_with_overlap$chr <- unlist(lapply(strsplit(alus_with_overlap$seqnames, "_"), "[[", 1))
table(alus_with_overlap$chr)
# save alu overlap file
save(alus_with_overlap, file = file.path(out_dir, "03_aluOverlaps_700b_oppStrands.RData"))


####
## Add Alu flag field to overlaps that indicates if an overlap region includes known inverted Alu sequences

# prepare GRanges objects
# multi-exon overlaps: will take the smallest start and the largest end position (widest range, will include introns)
temp <- putative_cisnat
temp[temp$multi_exon_overlap, ]$longest_overlap_start <- as.numeric(unlist(lapply(strsplit(temp[temp$multi_exon_overlap, ]$longest_overlap_start, ","), min)))
temp[temp$multi_exon_overlap, ]$longest_overlap_end <- as.numeric(unlist(lapply(strsplit(temp[temp$multi_exon_overlap, ]$longest_overlap_end, ","), max)))
cisnat_gr <- makeGRangesFromDataFrame(temp, 
                                      keep.extra.columns = F, seqnames.field = "chr", start.field = "longest_overlap_start",
                                      end.field = "longest_overlap_end", strand.field = "*")
iralu_gr <- makeGRangesFromDataFrame(alus_with_overlap, seqnames.field = "chr")

# findOverlaps of putative_cisnat with any Alus that have an overlap
# (strand does not matter, since the overlapping regions include both strands)
cisnat_alu_overlap <- as.data.frame(findOverlaps(cisnat_gr, iralu_gr, ignore.strand = T))
collapsed_overlaps <- group_by(cisnat_alu_overlap, queryHits) %>%
  group_map(combine_alus, alus_df = alus_with_overlap) %>% reduce(rbind)
putative_cisnat$alu_overlap <- NA
putative_cisnat[collapsed_overlaps$row_number, ]$alu_overlap <- collapsed_overlaps$alus

# check multi-exon overlaps - do they overlap with Alu
# (if so, should dive in and check whether the overlap is truly in the exon)
putative_cisnat[putative_cisnat$multi_exon_overlap,] # all NA, no worries

# save
save(putative_cisnat, file = file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData"))
write.table(putative_cisnat, file = file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.txt"),
            quote = F, sep = "\t", col.names = T, row.names = F)

