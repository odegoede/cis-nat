#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 01: Identifying all overlapping transcripts in human genome
## By: Olivia de Goede, 2021
#####

####
## SET OPTIONS, LOAD LIBRARIES
# renv::restore() # had problems installing BioConductor packages
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--gtf", default = NULL, 
              help = "GTF file [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory to write transcript overlaps to [default is current working directory]"),
  make_option("--scriptdir", default = "./scripts",
              help = "Scripts directory (to load functions saved in source scripts). [default assumes running from repo, scripts are in \"%default\"]"),
  make_option("--scratchdir", default = NULL, 
              help = "Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default \"%default\"]"),
  make_option("--keepanno", default = "neither", 
              help = "Should the gene, transcript, and exon annotations be saved? Files will be saved to outdir. \"text\" for .txt files, \"RData\" for .RData, \"both\" for both .txt and .RData files. Anything else will be interpreted as \"neither\". [default \"%default\"]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## DEFINE FUNCTIONS

# boundary_match(): for a tx_pair, see if any of its exon boundaries are shared (potential for 
# an overlap that spans multiple exons)
boundary_match <- function(df, id_field = "tx_pair") {
  bound_matches <- (df %>% 
                      mutate(start_match = minus_start == plus_start) %>%
                      mutate(end_match = minus_end == plus_end) %>%
                      filter(start_match == TRUE | end_match == TRUE))[[id_field]]
  bound_matches <- unique(bound_matches)
  bound_matches
}


# next_exon_check(): labels tx_pairs as TRUE/FALSE for multi-exon overlaps (output is just T/F, 
# not info on where the overlap is; that assessment is tricky enough that we should check manually)
next_exon_check <- function(df, tx_pair) {
  tx1 <- unlist(strsplit(tx_pair, "__"))[1]
  tx2 <- unlist(strsplit(tx_pair, "__"))[2]
  # first, are there multiple overlaps for this pair? (if not, won't overlap across multiple exons)
  if (!nrow(df[df$tx_pair == tx_pair, ]) > 1) {
    return(FALSE)
  }
  # second, make sure there's at least one shared start and one shared end (necessary for multi-exon overlap)
  shared_start <- any(df[df$tx_pair == tx_pair, ]$minus_start == df[df$tx_pair == tx_pair, ]$plus_start)
  shared_end <- any(df[df$tx_pair == tx_pair, ]$minus_end == df[df$tx_pair == tx_pair, ]$plus_end)
  if (!(shared_start & shared_end)) {
    return(FALSE)
  }
  # third, find shared boundary coordinates
  tx1_starts <- exon_anno[exon_anno$tx_id == tx1, ]$start
  tx1_starts <- tx1_starts[order(tx1_starts)]
  tx2_starts <- exon_anno[exon_anno$tx_id == tx2, ]$start
  tx2_starts <- tx2_starts[order(tx2_starts)]
  tx1_ends <- exon_anno[exon_anno$tx_id == tx1, ]$end
  tx1_ends <- tx1_ends[order(tx1_ends)]
  tx2_ends <- exon_anno[exon_anno$tx_id == tx2, ]$end
  tx2_ends <- tx2_ends[order(tx2_ends)]
  tx_1in2_start <- which(tx1_starts %in% tx2_starts)
  tx_1in2_end <- which(tx1_ends %in% tx2_ends)
  tx_2in1_start <- which(tx2_starts %in% tx1_starts)
  tx_2in1_end <- which(tx2_ends %in% tx1_ends)
  # fourth, check that the shared start position is from an exon that is 1 greater than the shared end
  # (may need to do awkward loops in the event of multiple shared start/end positions)
  tx1_pass <- FALSE
  for (s_ind in tx_1in2_start) {
    for (e_ind in tx_1in2_end) {
      if (s_ind - e_ind == 1) {tx1_pass <- TRUE}
    }
  }
  tx2_pass <- FALSE
  for (s_ind in tx_2in1_start) {
    for (e_ind in tx_2in1_end) {
      if (s_ind - e_ind == 1) {tx2_pass <- TRUE}
    }
  }
  if (tx1_pass & tx2_pass) {
    return(TRUE)
  }
}


# combine_dup_rows(): to parse through and combine info from the transcript pairs with multiple
# overlaps (i.e. tx_pairs with more than 1 row in overlap_df)
combine_dup_rows <- function(df = overlap_df, tx_pair) {
  temp <- df[df$tx_pair == tx_pair, ]
  out_df <- data.frame(chr = unique(temp$chr),
                       minus_gene = unique(temp$minus_gene),
                       minus_symbol = NA,
                       minus_gene_type = NA,
                       plus_gene = unique(temp$plus_gene),
                       plus_symbol = NA,
                       plus_gene_type = NA,
                       minus_tx = unique(temp$minus_tx),
                       minus_tx_type = NA,
                       plus_tx = unique(temp$plus_tx),
                       plus_tx_type = NA,
                       minus_exons_in_overlap = paste(temp$minus_exon_number, collapse = ","),
                       plus_exons_in_overlap = paste(temp$plus_exon_number, collapse = ","),
                       all_overlap_widths = paste(temp$overlap_width, collapse = ","),
                       longest_overlap_width = max(temp$overlap_width),
                       longest_overlap_start = temp[which.max(temp$overlap_width), ]$overlap_start,
                       longest_overlap_end = temp[which.max(temp$overlap_width), ]$overlap_end,
                       duplicate_overlap = NA,
                       row.names = tx_pair)
  out_df
}


# choose_which_dup(): for duplicate overlap events (same gene pair, exact same coordinates of overlap), 
# pick the "best" one (based on TSL and transcript types)
choose_which_dup <- function(df = to_filter, overlap_id) {
  temp <- df[df$overlap_id == overlap_id, ]
  # first, check if we can use TSL for priority:
  if (length(unique(temp$minus_tsl)) > 1 | length(unique(temp$plus_tsl)) > 1) {
    minus_ind <- which(temp$minus_tsl == min(temp$minus_tsl)) # rows with the best minus strand TSL
    plus_ind <- which(temp$plus_tsl == min(temp$plus_tsl)) # rows with the best plus strand TSL
    # see if we can just filter to rows with the best minus and plus strand TSLs:
    if (length(intersect(minus_ind, plus_ind)) > 0) {
      temp <- temp[intersect(minus_ind, plus_ind), ]
    }
    # if we can't just filter to rows with the best minus&plus strand TSLs, do a variety of TSL checks:
    else {
      # check if one has a bigger TSL difference
      minus_diff <- nth(unique(as.numeric(temp$minus_tsl)), 2) - min(as.numeric(temp$minus_tsl))
      plus_diff <- nth(unique(as.numeric(temp$plus_tsl)), 2) - min(as.numeric(temp$plus_tsl))
      if (minus_diff != plus_diff) {
        if (minus_diff > plus_diff) {temp <- temp[minus_ind, ]}
        else {temp <- temp[plus_ind, ]}
      }
      # if no diff, take whichever has the worst-best TSL
      else if (min(as.numeric(temp$minus_tsl)) != min(as.numeric(temp$plus_tsl))) {
        if (min(as.numeric(temp$minus_tsl)) > min(as.numeric(temp$plus_tsl))) {temp <- temp[minus_ind, ]}
        else {temp <- temp[plus_ind, ]}
      }
      # if no diff, just take the union and check for tx type diffs
      else {
        temp <- temp[union(minus_ind, plus_ind), ]
      }
    }
  }
  # next, check if we can use transcript type for priority:
  if (length(unique(temp$minus_tx_type)) > 1 | length(unique(temp$plus_tx_type)) > 1) {
    minus_ind <- which(temp$minus_tx_type == min(temp$minus_tx_type))
    plus_ind <- which(temp$plus_tx_type == min(temp$plus_tx_type))
    if (length(intersect(minus_ind, plus_ind)) > 0) {
      temp <- temp[intersect(minus_ind, plus_ind), ]
    }
  }
  # finally, if all TSL and tx types are the same or if there are still multiple rows in temp, just take the first one (very similar at this point)
  if (nrow(temp) > 1) {
    temp <- temp[1, ]
  }
  rownames(temp)
}



####
## INPUT TESTS
# Check the required arguments (GTF file) is provided
if (is.null(opt$gtf)) { 
  stop("GTF file not provided, exiting\n") 
}

# Create variable out_dir and script_dir that is consistent (doesn't have trailing "/")
out_dir <- file.path(opt$outdir)
script_dir <- file.path(opt$scriptdir)
# Check that the directories exist
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist, exiting\n")
}
if (!dir.exists(script_dir)) {
  stop("Script directory does not exist, exiting\n")
}

# If scratch directory is provided, create variable scratch_dir that is consistent
if (!is.null(opt$scratchdir)) {
  scratch_dir <- file.path(opt$scratchdir)
  # Check that scratch directory exists
  if (!dir.exists(scratch_dir)) {
    stop("Scratch directory does not exist, exiting\n")
  }
}



####
## LOAD SOURCE SCRIPTS
source(file.path(script_dir, "source_temp_unzip.R"))
source(file.path(script_dir, "source_gtf_reader.R"))



####
## READ IN INPUT FILES
gtf_anno <- read_gtf(filename = opt$gtf)

# Remove unnecessary columns and mitochondrial chromosome
gtf_anno <- gtf_anno[,-c(2,6,8)] # these are unnecessary fields
gtf_anno <- gtf_anno[gtf_anno$chr != "ChrM", ]



####
## Make the gene_anno, transcript_anno, and exon_anno files
# add gene information from attribute field (all 3 anno files need this)
gtf_anno <- get_basic_attr(gtf_anno)

# filter to gene rows only
gene_anno <- gtf_anno[gtf_anno$type == "gene", ]
gene_anno <- gene_anno[,-grep("attribute", colnames(gene_anno))]
rownames(gene_anno) <- gene_anno$ensgene

# filter gtf_anno to transcript rows, add extra transcript information
tx_anno <- get_tx_attr(gtf_anno)
tx_anno$tsl <- "tsl_NA"
attr_list <- get_attr_list(tx_anno[grepl("transcript_support_level", tx_anno$attribute), ])
tx_anno[grepl("transcript_support_level", tx_anno$attribute), ]$tsl <- paste0("tsl_", pull_attr(attr_list, "transcript_support_level"))
tx_anno <- tx_anno[,-grep("attribute", colnames(tx_anno))]
rownames(tx_anno) <- tx_anno$tx_id

# filter gtf_anno to exon rows, add extra transcript information and exon number
# NOTE: for exons on the minus strand, the GENCODE GTF numbers them in opposite order (exon 1 isn't 
# the first exon, but the last - it's the one with the lowest position values)
exon_anno <- get_exon_attr(gtf_anno)
exon_anno$tsl <- "tsl_NA"
attr_list <- get_attr_list(exon_anno[grepl("transcript_support_level", exon_anno$attribute), ])
exon_anno[grepl("transcript_support_level", exon_anno$attribute), ]$tsl <- paste0("tsl_", pull_attr(attr_list, "transcript_support_level"))
exon_anno <- exon_anno[,-grep("attribute", colnames(exon_anno))]

# if the user put something in --keepanno, save these files to outdir
if (opt$keepanno %in% c("text", "both")) {
  write.table(gene_anno, file = paste0(out_dir, "gene_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
  write.table(tx_anno, file = paste0(out_dir, "transcript_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
  write.table(exon_anno, file = paste0(out_dir, "exon_anno.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
}
if (opt$keepanno %in% c("RData", "both")) {
  save(gene_anno, tx_anno, exon_anno, file = file.path(out_dir, "gene_tx_exon_anno_files.RData"))
}



####
## Find overlaps! 
## Criteria for overlaps:
# 1. transcripts come from opposite strands
# 2. transcripts either have good transcript support level scores (TSL 1-3), or are the transcripts
#    with the most evidence for their gene
# 3. have at least 100 bp *continuous* overlap


## Step 1: Filter transcripts based on TSL
temp_tx <- tx_anno

good_scores <- c("tsl_1", "tsl_2", "tsl_3")
bad_scores <- c("tsl_4", "tsl_5", "tsl_NA")
temp_tx$to_keep <- FALSE
temp_tx[temp_tx$tsl %in% good_scores, ]$to_keep <- TRUE
# which genes have at least one transcript with a good score?:
good_genes <- unique(temp_tx[temp_tx$to_keep, ]$ensgene)
# which genes don't have any transcripts with good scores?:
bad_genes <- unique(temp_tx$ensgene)[-which(unique(temp_tx$ensgene) %in% good_genes)]
# For the bad_genes, get the transcripts with the best possible TSL scores (even though they're
# subthreshold), just so that we have some representation of that gene.
# Loop through each "bad score" (from least-bad to worst), get the set of bad_genes that have
# at least one transcript of that bad-but-not-worst score; mark the transcripts (rows of temp_tx)
# that are in these subthresh_genes and have this bad-but-not-worst score as rows to keep:
for (i in bad_scores) {
  subthresh_genes <- unique(temp_tx[temp_tx$ensgene %in% bad_genes & temp_tx$tsl == i, ]$ensgene)
  temp_tx[temp_tx$ensgene %in% subthresh_genes & temp_tx$tsl == i, ]$to_keep <- TRUE
  bad_genes <- bad_genes[-which(bad_genes %in% subthresh_genes)]
}
table(temp_tx$to_keep) # FALSE, 80314 ; TRUE, 149266
temp_tx <- temp_tx[temp_tx$to_keep, ]
temp_exon <- exon_anno[exon_anno$tx_id %in% temp_tx$tx_id, ]
# save these intermediate files if scratchdir was provided:
if (exists("scratch_dir")) {
  save(temp_tx, temp_exon, file = file.path(scratch_dir, "01_starterFiles_find_tx_overlaps.RData"))
}



## Step 2: Separate into plus and minus strand GRanges objects
# (overlaps need to be on opposite strands)
table(temp_tx$strand) # 73,237 - ; 76,029 +

plus_tx <- temp_tx[temp_tx$strand == "+", ]
minus_tx <- temp_tx[temp_tx$strand == "-", ]
plus_exon <- makeGRangesFromDataFrame(temp_exon[temp_exon$strand == "+", ], keep.extra.columns = T)
minus_exon <- makeGRangesFromDataFrame(temp_exon[temp_exon$strand == "-", ], keep.extra.columns = T)
# give the exons unique identifiers (add the exon number to the end of the tx_id):
names(plus_exon) <- paste(plus_exon$tx_id, plus_exon$exon_number, sep = "_")
names(minus_exon) <- paste(minus_exon$tx_id, minus_exon$exon_number, sep = "_")



## Step 3: Find overlaps
hits <- findOverlaps(plus_exon, minus_exon, select = "all", ignore.strand = T, minoverlap = 100)

minus_ind <- subjectHits(hits) # indices of overlapping regions in minus_exon GRanges object
plus_ind <- queryHits(hits) # indices of overlapping regions in plus_exon GRanges object
length(minus_ind) # 19,971 overlaps 

# make a data.frame summarizing all of these overlaps
overlap_df <- data.frame(minus_ind = minus_ind, plus_ind = plus_ind, 
                         minus_gene = minus_exon[minus_ind]$ensgene, plus_gene = plus_exon[plus_ind]$ensgene,
                         minus_tx = minus_exon[minus_ind]$tx_id, plus_tx = plus_exon[plus_ind]$tx_id,
                         minus_exon = names(minus_exon[minus_ind]), plus_exon = names(plus_exon[plus_ind]),
                         chr = as.character(seqnames(minus_exon[minus_ind])), 
                         minus_start = start(minus_exon[minus_ind]), minus_end = end(minus_exon[minus_ind]),
                         plus_start = start(plus_exon[plus_ind]), plus_end = end(plus_exon[plus_ind]))

overlap_df$tx_pair <- paste(overlap_df$minus_tx, overlap_df$plus_tx, sep = "__")
overlap_df$gene_pair <- paste(overlap_df$minus_gene, overlap_df$plus_gene, sep = "__")
length(unique(overlap_df$tx_pair)) # 16,096 unique tx pairs in these overlaps
length(unique(overlap_df$gene_pair)) # 4,837 unique gene pairs in these overlaps

overlap_df$overlap_start <- apply(overlap_df[,c("minus_start", "plus_start")], 1, max)
overlap_df$overlap_end <- apply(overlap_df[,c("minus_end", "plus_end")], 1, min)
# get overlap width (width = end-start+1)
overlap_df$overlap_width <- overlap_df$overlap_end - overlap_df$overlap_start + 1



## Step 4: Check for any overlapping transcript pairs that share start OR end coordinates in an exon
# (this is rare, but this means that once processed the transcripts could have continuous overlaps 
# across multiple exons, which would need to be summed)
length(start_end_matches <- boundary_match(overlap_df)) 
# ^ 217 tx pairs where the exon match, right up to at least one exon boundary (and may extend into next)

# If these boundary match pairs also overlap at the next exon starting right from its boundary, 
# they are multi-exon overlaps
multi_check <- sapply(start_end_matches, next_exon_check, df = overlap_df)
length(multi_overlaps <- names(multi_check[multi_check])) # 5 transcript pairs with multi-exon overlaps



## Step 5: Organize overlaps into a data.frame
# each row is a unique tx_pair
# fields: chr, minus_gene, minus_symbol, minus_gene_type, plus_gene, plus_symbol, plus_gene_type, 
# minus_tx, minus_tx_type, plus_tx, plus_tx_type, minus_exons_in_overlap, plus_exons_in_overlap, 
# all_overlap_widths, longest_overlap_width, longest_overlap_start, longest_overlap_end, duplicate_overlap
# NOTE: the "duplicate_overlap" field marks overlaps involving the same gene pair, and that have the exact 
# same overlap coordinates (occurs because of alternative transcripts of the genes having same overlap)
overlap_df$minus_exon_number <- unlist(lapply(strsplit(overlap_df$minus_exon, "_"), "[[", 2))
overlap_df$plus_exon_number <- unlist(lapply(strsplit(overlap_df$plus_exon, "_"), "[[", 2))

c_names <- c("chr", "minus_gene", "minus_symbol", "minus_gene_type", "plus_gene", "plus_symbol", "plus_gene_type", 
             "minus_tx", "minus_tx_type", "plus_tx", "plus_tx_type", "minus_exons_in_overlap", "plus_exons_in_overlap", 
             "all_overlap_widths", "longest_overlap_width", "longest_overlap_start", "longest_overlap_end", 
             "duplicate_overlap")
tx_pair_dat <- data.frame(matrix(NA, nrow = length(unique(overlap_df$tx_pair)), ncol = length(c_names)),
                          row.names = unique(overlap_df$tx_pair))
colnames(tx_pair_dat) <- c_names

# pull out tx_pairs that have duplicates (or more) in overlap_df; 
# will parse through their overlap info, combine it, and add it into just one row in tx_pair_dat
dupl_tx_pair <- unique(overlap_df[which(duplicated(overlap_df$tx_pair)), ]$tx_pair)
non_dupl_tx_dat <- overlap_df[-which(overlap_df$tx_pair %in% dupl_tx_pair), ] # df of non-duplicate tx_pairs
rownames(non_dupl_tx_dat) <- non_dupl_tx_dat$tx_pair
# can put info from non_dupl_tx_dat directly into tx_pair_dat (don't need to combine info from multiple rows):
tx_pair_dat[,c(1,2,5,8,10,12,13,14,15,16,17)] <- non_dupl_tx_dat[rownames(tx_pair_dat), c("chr", "minus_gene", "plus_gene", "minus_tx", "plus_tx", 
                                                                         "minus_exon_number", "plus_exon_number", "overlap_width", 
                                                                         "overlap_width", "overlap_start", "overlap_end")]
# now parse through the tx_pairs with >1 row in overlap_df:
dupl_dat <- do.call(rbind, lapply(dupl_tx_pair, combine_dup_rows, df = overlap_df))
tx_pair_dat[dupl_tx_pair, ] <- dupl_dat[dupl_tx_pair, ]

tx_pair_dat$minus_symbol <- gene_anno[tx_pair_dat$minus_gene, ]$symbol
tx_pair_dat$minus_gene_type <- gene_anno[tx_pair_dat$minus_gene, ]$biotype
tx_pair_dat$minus_tx_type <- tx_anno[tx_pair_dat$minus_tx, ]$tx_type
tx_pair_dat$plus_symbol <- gene_anno[tx_pair_dat$plus_gene, ]$symbol
tx_pair_dat$plus_gene_type <- gene_anno[tx_pair_dat$plus_gene, ]$biotype
tx_pair_dat$plus_tx_type <- tx_anno[tx_pair_dat$plus_tx, ]$tx_type



## Step 6: Combine entries for the multi-exon overlaps
# for this step I have updated each situation by hand rather than tried to write a function, since there
# are only 5 of these multi-exon overlaps, and each one will be a different situation (would take a lot
# of time and effort to generalize into a function)
tx_pair_dat$multi_exon_overlap <- FALSE

## "ENST00000599229.2__ENST00000382100.8"
overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382100.8", ]
tx_pair_dat["ENST00000599229.2__ENST00000382100.8", ]
exon_anno[exon_anno$tx_id == "ENST00000382100.8", ]
exon_anno[exon_anno$tx_id == "ENST00000599229.2", ]
# all five exons with overlap are shared-boundary ones (are consecutive exon number and the intermediate exons' starts and ends are the same), can just sum their overlap_width
tx_pair_dat["ENST00000599229.2__ENST00000382100.8", ]$multi_exon_overlap <- TRUE
tx_pair_dat["ENST00000599229.2__ENST00000382100.8", ]$longest_overlap_width <- sum(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382100.8", ]$overlap_width)
tx_pair_dat["ENST00000599229.2__ENST00000382100.8", ]$longest_overlap_start <- paste(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382100.8", ]$overlap_start, collapse = ",")
tx_pair_dat["ENST00000599229.2__ENST00000382100.8", ]$longest_overlap_end <- paste(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382100.8", ]$overlap_end, collapse = ",")

## "ENST00000599229.2__ENST00000382099.2"
overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382099.2", ]
tx_pair_dat["ENST00000599229.2__ENST00000382099.2", ]
exon_anno[exon_anno$tx_id == "ENST00000382099.2", ]
exon_anno[exon_anno$tx_id == "ENST00000599229.2", ]
# all five exons with overlap are shared-boundary ones (are consecutive exon number and the intermediate exons' starts and ends are the same), can just sum their overlap_width
tx_pair_dat["ENST00000599229.2__ENST00000382099.2", ]$multi_exon_overlap <- TRUE
tx_pair_dat["ENST00000599229.2__ENST00000382099.2", ]$longest_overlap_width <- sum(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382099.2", ]$overlap_width)
tx_pair_dat["ENST00000599229.2__ENST00000382099.2", ]$longest_overlap_start <- paste(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382099.2", ]$overlap_start, collapse = ",")
tx_pair_dat["ENST00000599229.2__ENST00000382099.2", ]$longest_overlap_end <- paste(overlap_df[overlap_df$tx_pair == "ENST00000599229.2__ENST00000382099.2", ]$overlap_end, collapse = ",")

## "ENST00000319770.7__ENST00000396799.3"
overlap_df[overlap_df$tx_pair == "ENST00000319770.7__ENST00000396799.3", ]
tx_pair_dat["ENST00000319770.7__ENST00000396799.3", ]
exon_anno[exon_anno$tx_id == "ENST00000319770.7", ]
exon_anno[exon_anno$tx_id == "ENST00000396799.3", ]
# the only two exons with overlap are the shared-boundary ones, can just sum their overlap_width
tx_pair_dat["ENST00000319770.7__ENST00000396799.3", ]$multi_exon_overlap <- TRUE
tx_pair_dat["ENST00000319770.7__ENST00000396799.3", ]$longest_overlap_width <- sum(overlap_df[overlap_df$tx_pair == "ENST00000319770.7__ENST00000396799.3", ]$overlap_width)
tx_pair_dat["ENST00000319770.7__ENST00000396799.3", ]$longest_overlap_start <- paste(overlap_df[overlap_df$tx_pair == "ENST00000319770.7__ENST00000396799.3", ]$overlap_start, collapse = ",")
tx_pair_dat["ENST00000319770.7__ENST00000396799.3", ]$longest_overlap_end <- paste(overlap_df[overlap_df$tx_pair == "ENST00000319770.7__ENST00000396799.3", ]$overlap_end, collapse = ",")

## "ENST00000396801.7__ENST00000396799.3"
overlap_df[overlap_df$tx_pair == "ENST00000396801.7__ENST00000396799.3", ]
tx_pair_dat["ENST00000396801.7__ENST00000396799.3", ]
exon_anno[exon_anno$tx_id == "ENST00000396799.3", ]
exon_anno[exon_anno$tx_id == "ENST00000396801.7", ]
# the only two exons with overlap are the shared-boundary ones, can just sum their overlap_width
tx_pair_dat["ENST00000396801.7__ENST00000396799.3", ]$multi_exon_overlap <- TRUE
tx_pair_dat["ENST00000396801.7__ENST00000396799.3", ]$longest_overlap_width <- sum(overlap_df[overlap_df$tx_pair == "ENST00000396801.7__ENST00000396799.3", ]$overlap_width)
tx_pair_dat["ENST00000396801.7__ENST00000396799.3", ]$longest_overlap_start <- paste(overlap_df[overlap_df$tx_pair == "ENST00000396801.7__ENST00000396799.3", ]$overlap_start, collapse = ",")
tx_pair_dat["ENST00000396801.7__ENST00000396799.3", ]$longest_overlap_end <- paste(overlap_df[overlap_df$tx_pair == "ENST00000396801.7__ENST00000396799.3", ]$overlap_end, collapse = ",")

## "ENST00000355772.8__ENST00000396799.3"
overlap_df[overlap_df$tx_pair == "ENST00000355772.8__ENST00000396799.3", ]
tx_pair_dat["ENST00000355772.8__ENST00000396799.3", ]
exon_anno[exon_anno$tx_id == "ENST00000396799.3", ]
exon_anno[exon_anno$tx_id == "ENST00000355772.8", ]
# the only two exons with overlap are the shared-boundary ones, can just sum their overlap_width
tx_pair_dat["ENST00000355772.8__ENST00000396799.3", ]$multi_exon_overlap <- TRUE
tx_pair_dat["ENST00000355772.8__ENST00000396799.3", ]$longest_overlap_width <- sum(overlap_df[overlap_df$tx_pair == "ENST00000355772.8__ENST00000396799.3", ]$overlap_width)
tx_pair_dat["ENST00000355772.8__ENST00000396799.3", ]$longest_overlap_start <- paste(overlap_df[overlap_df$tx_pair == "ENST00000355772.8__ENST00000396799.3", ]$overlap_start, collapse = ",")
tx_pair_dat["ENST00000355772.8__ENST00000396799.3", ]$longest_overlap_end <- paste(overlap_df[overlap_df$tx_pair == "ENST00000355772.8__ENST00000396799.3", ]$overlap_end, collapse = ",")



## Step 7: Flag duplicate overlaps
# duplicate_overlap field: marks overlaps involving alternative transcripts of the same gene that 
# have the exact same overlapping region
tx_pair_dat$overlap_id <- paste(paste(tx_pair_dat$minus_gene, tx_pair_dat$plus_gene, sep = "_"),
                                as.character(apply(tx_pair_dat[,c("longest_overlap_start", "longest_overlap_end", "longest_overlap_width")], 1, paste, collapse = ";")), sep = "__")
tx_pair_dat$duplicate_overlap <- FALSE
duplicate_ids <- tx_pair_dat[duplicated(tx_pair_dat$overlap_id), ]$overlap_id
tx_pair_dat[tx_pair_dat$overlap_id %in% duplicate_ids, ]$duplicate_overlap <- TRUE
summary(tx_pair_dat$duplicate_overlap) # 6,174 unique overlaps; 9,922 overlaps represent the same gene pair and overlapping region
# save these intermediate files if scratchdir was provided
if (exists("scratch_dir")) {
  save(tx_pair_dat, file = file.path(scratch_dir, "01_full_transcript_overlap_table.RData"))
}



## Step 8: filter down to unique overlap regions
# some tx_pairs are alternative transcripts of genes that have the exact same overlap
# this filters to just unique overlap regions
# (but note: the same gene pair is still often represented multiple times, e.g. if alternative transcripts
# have slightly different overlap ranges)
unique_overlaps <- tx_pair_dat[tx_pair_dat$duplicate_overlap == FALSE, ]

# how to filter?
# 1. based on stronger TSL score
# 2. based on transcript annotation
tsl_priority <- c("tsl_1", "tsl_2", "tsl_3", "tsl_4", "tsl_5", "tsl_NA")
tx_priority <- c("protein_coding", "lncRNA", "processed_transcript", "retained_intron", "nonsense_mediated_decay",
                 "non_stop_decay", "pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene",
                 "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene",
                 "snRNA", "snoRNA", "misc_RNA", "TEC")

to_filter <- tx_pair_dat[tx_pair_dat$duplicate_overlap == TRUE, ]

to_filter$minus_tsl <- factor(tx_anno[to_filter$minus_tx, ]$tsl, levels = tsl_priority, ordered = T) # min is best (tsl_1), max is worst (tsl_NA)
to_filter$plus_tsl <- factor(tx_anno[to_filter$plus_tx, ]$tsl, levels = tsl_priority, ordered = T)

to_filter$minus_tx_type <- factor(to_filter$minus_tx_type, levels = tx_priority, ordered = T) # min is best (protein-coding), max is worst (TEC)
to_filter$plus_tx_type <- factor(to_filter$plus_tx_type, levels = tx_priority, ordered = T)

dup_ids <- unique(to_filter$overlap_id) # 3,011 rows from these 9,922 rows
keep_ids <- as.character(sapply(dup_ids, choose_which_dup, df = to_filter))
to_keep <- to_filter[keep_ids, colnames(unique_overlaps)]
dim(to_keep) # 3,011 overlaps

putative_cisnat <- rbind(unique_overlaps, to_keep) # 9,185 rows
putative_cisnat <- putative_cisnat[,-c(grep("duplicate_overlap", colnames(putative_cisnat)),
                                       grep("overlap_id", colnames(putative_cisnat)))]



####
## Save overlap annotation file
save(putative_cisnat, file = file.path(out_dir, "01_putative_cisNAT_uniqueRegionsOnly.RData"))
write.table(putative_cisnat, file = file.path(out_dir, "01_putative_cisNAT_uniqueRegionsOnly.txt"),
            quote = F, sep = "\t", col.names = T, row.names = F)


