#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 05: Find overlaps with RNA editing data
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/05_check_editing_sites.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --editanno data/gtex_edit/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.AnnoVar_Genes.Clusered.txt.gz --examplequant data/gtex_edit/indiv_files/GTEX-11ZTS-1026-SM-5LU8O.txt.gz

## Idea:
# narrow down the list of overlapping transcripts based on whether 
# there's at least 1 RNA editing site within the region of overlap
#
# also want to save list of editing sites in these overlaps, to 
# filter & combine individual coverage and editing files

# Important file note:
# not all of the sites quantified are in clusters (so each individual quant file 
# has rows that cluster_dat does not)


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
options(datatable.fread.datatable = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--editanno", default = NULL, 
              help = "File with editing clusters/site locations [default \"%default\"]"),
  make_option("--examplequant", default = NULL, 
              help = "Example of individual RNA editing quantification file [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory [default is \"%default\"]"),
  make_option("--figdir", default = "./figures", 
              help = "Figure directory [default is \"%default\"]"),
  make_option("--scriptdir", default = "./scripts", 
              help = "Script directory (to load functions saved in source scripts). [default is \"%default\"]"),
  make_option("--scratchdir", default = NULL, 
              help = "Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default \"%default\"]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## INPUT TESTS
# Check the required files are provided
if (is.null(opt$editanno)) { 
  stop("Editing annotation file not provided, exiting\n") 
}
if (is.null(opt$examplequant)) { 
  stop("Example quantification file not provided, exiting\n") 
}
# Check that the required files exist
if (!file.exists(opt$editanno)) {
  stop("Editing annotation file does not exist as named, exiting\n")
}
if (!file.exists(opt$examplequant)) {
  stop("Example quantification file does not exist as named, exiting\n")
}

# Create consistent directory variables
out_dir <- file.path(opt$outdir)
script_dir <- file.path(opt$scriptdir)
fig_dir <- file.path(opt$figdir)
# Check that the directories exist
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist, exiting\n")
}
if (!dir.exists(script_dir)) {
  stop("Script directory does not exist, exiting\n")
}
if (!dir.exists(fig_dir)) {
  stop("Figure directory does not exist, exiting\n")
}

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


####
## READ IN INPUT FILES
load(file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData")) # object name from script 03: putative_cisnat
load(file.path(out_dir, "gene_tx_exon_anno_files.RData")) # object name from script 01: gene_anno, tx_anno, and exon_anno

edit_anno_file <- file.path(opt$editanno)
if (endsWith(edit_anno_file, "gz") | endsWith(edit_anno_file, "zip") | endsWith(edit_anno_file, "bzip2") | endsWith(edit_anno_file, "xz")) {
  cluster_dat <- suppressWarnings(temp_unzip(edit_anno_file, fread))
} else {
  cluster_dat <- suppressWarnings(fread(edit_anno_file))
}
colnames(cluster_dat) <- c("chr", "site_start", "site_end", "gene_name", "cluster_name", "strand")

ex_quant_file <- file.path(opt$examplequant)
if (endsWith(ex_quant_file, "gz") | endsWith(ex_quant_file, "zip") | endsWith(ex_quant_file, "bzip2") | endsWith(ex_quant_file, "xz")) {
  indiv_quant <- suppressWarnings(temp_unzip(ex_quant_file, fread))
} else {
  indiv_quant <- suppressWarnings(fread(ex_quant_file))
}
colnames(indiv_quant) <- c("chr", "position", "coverage", "editedreads", "editlevel")


####
## DEFINE FUNCTIONS
# get_pair_type() returns the categories of gene pairs, ordered with some priority level (e.g so you don't get
# both protein_coding__lncRNA and lncRNA__protein_coding as pair types)
get_pair_type <- function(df = plot_dat, pair_type = "tx", tx_priority = tx_priority) {
  if (!pair_type %in% c("tx", "gene")) {stop("pair_type isn't \"tx\" or \"gene\"")}
  pseudo_types <- c("pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", 
                    "transcribed_unprocessed_pseudogene", "processed_pseudogene", "unitary_pseudogene", 
                    "polymorphic_pseudogene", "unprocessed_pseudogene", "IG_V_pseudogene", "rRNA_pseudogene")
  small_types <- c("scaRNA", "miRNA", "snRNA", "snoRNA")
  other_types <- c("misc_RNA", "TEC", "TR_C_gene")
  if (pair_type == "tx") {
    plot_dat[plot_dat$plus_tx_type %in% pseudo_types, ]$plus_tx_type <- "pseudogene"
    plot_dat[plot_dat$plus_tx_type %in% small_types, ]$plus_tx_type <- "small_RNA"
    plot_dat[plot_dat$plus_tx_type %in% other_types, ]$plus_tx_type <- "other"
    plot_dat[plot_dat$minus_tx_type %in% pseudo_types, ]$minus_tx_type <- "pseudogene"
    plot_dat[plot_dat$minus_tx_type %in% small_types, ]$minus_tx_type <- "small_RNA"
    plot_dat[plot_dat$minus_tx_type %in% other_types, ]$minus_tx_type <- "other"
    plot_dat$plus_tx_type <- factor(plot_dat$plus_tx_type, levels = tx_priority, ordered = T)
    plot_dat$minus_tx_type <- factor(plot_dat$minus_tx_type, levels = tx_priority, ordered = T)
    plot_dat$minus_number <- as.numeric(plot_dat$minus_tx_type)
    plot_dat$plus_number <- as.numeric(plot_dat$plus_tx_type)
    plot_dat$plus_first <- plot_dat$plus_number <= plot_dat$minus_number
    plot_dat$xcat_tx <- NA
    plot_dat[plot_dat$plus_first, ]$xcat_tx <- paste(as.character(plot_dat[plot_dat$plus_first, ]$plus_tx_type), 
                                                     as.character(plot_dat[plot_dat$plus_first, ]$minus_tx_type), sep = "__")
    plot_dat[!plot_dat$plus_first, ]$xcat_tx <- paste(as.character(plot_dat[!plot_dat$plus_first, ]$minus_tx_type), 
                                                      as.character(plot_dat[!plot_dat$plus_first, ]$plus_tx_type), sep = "__")
  }
  else if (pair_type == "gene") {
    plot_dat[plot_dat$plus_gene_type %in% pseudo_types, ]$plus_gene_type <- "pseudogene"
    plot_dat[plot_dat$plus_gene_type %in% small_types, ]$plus_gene_type <- "small_RNA"
    plot_dat[plot_dat$plus_gene_type %in% other_types, ]$plus_gene_type <- "other"
    plot_dat[plot_dat$minus_gene_type %in% pseudo_types, ]$minus_gene_type <- "pseudogene"
    plot_dat[plot_dat$minus_gene_type %in% small_types, ]$minus_gene_type <- "small_RNA"
    plot_dat[plot_dat$minus_gene_type %in% other_types, ]$minus_gene_type <- "other"
    plot_dat$plus_gene_type <- factor(plot_dat$plus_gene_type, levels = tx_priority, ordered = T)
    plot_dat$minus_gene_type <- factor(plot_dat$minus_gene_type, levels = tx_priority, ordered = T)
    plot_dat$minus_number <- as.numeric(plot_dat$minus_gene_type)
    plot_dat$plus_number <- as.numeric(plot_dat$plus_gene_type)
    plot_dat$plus_first <- plot_dat$plus_number <= plot_dat$minus_number
    plot_dat$xcat_gene <- NA
    plot_dat[plot_dat$plus_first, ]$xcat_gene <- paste(as.character(plot_dat[plot_dat$plus_first, ]$plus_gene_type), 
                                                       as.character(plot_dat[plot_dat$plus_first, ]$minus_gene_type), sep = "__")
    plot_dat[!plot_dat$plus_first, ]$xcat_gene <- paste(as.character(plot_dat[!plot_dat$plus_first, ]$minus_gene_type), 
                                                        as.character(plot_dat[!plot_dat$plus_first, ]$plus_gene_type), sep = "__")
  }
  plot_dat
}


####
## Save cluster file for faster loading in the future
# also: position info is different between the cluster_dat annotation and the individual
# quantification files; individual quant files have position 1 lower than cluster_dat annotation
cluster_dat$site_name <- paste(cluster_dat$chr, cluster_dat$site_start, sep = "_")
cluster_dat$site_name_minusOne <- paste(cluster_dat$chr, cluster_dat$site_start-1, sep = "_")

indiv_quant$site_name <- paste(indiv_quant$chr, indiv_quant$position, sep = "_")
summary(cluster_dat$site_name %in% indiv_quant$site_name) # only 266633/2.1 million
summary(cluster_dat$site_name_minusOne %in% indiv_quant$site_name) # 2182197/2.1 million (only 3 not matched up)
summary(duplicated(cluster_dat$site_name)) # all unique
if (!dir.exists("data/gtex_edit")) {
  dir.create("data/gtex_edit")
}
save(cluster_dat, file = file.path("data/gtex_edit/All.AG.stranded.annovar.RData"))


####
## Prep for checking tx overlaps for editing sites
# remember that a couple of overlaps span multiple exons, so there are commas in the start/end fields; separate those out
summary(grepl(",", putative_cisnat$longest_overlap_start)) # 2 of them
# since so few, just make a separate grange for each manually (magic numbers alert)
# (otherwise will need to keep track of names, and sum all editing sites falling into any of the 
# ranges associated with that overlap name... would be annoying)
putative_cisnat[grepl(",", putative_cisnat$longest_overlap_start),] # magic numbers source
comp_one <- GRanges(c("chr12", "chr12"), IRanges(start = c(6666730, 6669031), end = c(6668115, 6669184), 
                                                 names = c("ENST00000361959.7__ENST00000396799.3", "ENST00000361959.7__ENST00000396799.3")))
comp_two <- GRanges(rep("chr9", 3), IRanges(start = c(2622100,2635453,2639859), end = c(2622271,2635572,2639981), 
                                            names = rep("ENST00000599229.2__ENST00000382096.5", 3)))
temp <- putative_cisnat[-which(rownames(putative_cisnat) %in% c("ENST00000361959.7__ENST00000396799.3", "ENST00000599229.2__ENST00000382096.5")), ]
put_cisnat_grange <- makeGRangesFromDataFrame(temp, seqnames.field = "chr", start.field = "longest_overlap_start", 
                                              end.field = "longest_overlap_end", keep.extra.columns = T)


####
## Check 1: how many transcript overlap regions contain any editing site?
site_grange <- makeGRangesFromDataFrame(indiv_quant, seqnames.field = "chr", start.field = "position", end.field = "position", keep.extra.columns = T)

put_cisnat_to_site <- findOverlaps(put_cisnat_grange, site_grange, select = "all", ignore.strand = T)
put_cisnat_ind <- queryHits(put_cisnat_to_site)
site_ind <- subjectHits(put_cisnat_to_site)
# random sanity check:
put_cisnat_grange[53]
site_grange[site_ind[put_cisnat_ind == 53]]

length(unique(put_cisnat_ind)) # 1,083 of the overlaps contain at least 1 RNA editing site
# what about the two multi-exon overlaps?:
findOverlaps(comp_one, site_grange, select = "all", ignore.strand = T) # no overlaps
findOverlaps(comp_two, site_grange, select = "all", ignore.strand = T) # no overlaps
# 1,083 / 15,986 contain RNA-editing sites

## Plot histogram of number of editing sites/overlap
# (was originally going to have a 0-overlaps bar, but that will dwarf the plot)
plot_dat <- data.frame(put_cisnat_ind = as.integer(names(table(put_cisnat_ind))),
                       num_sites = as.integer(table(put_cisnat_ind)))
plot_dat$put_cisnat_name <- names(put_cisnat_grange[plot_dat$put_cisnat_ind])
summary(plot_dat$num_sites) # range from 1 to 215
(g <- ggplot(plot_dat, aes(x = num_sites)) + geom_bar() + theme_bw() +
    ggtitle("Number of RNA-editing sites in identified overlaps (1,083 with at least 1 site)"))
ggsave(g, file = file.path(fig_dir, "05_01_numSites_in_overlaps.png"), height = 6, width = 7)
ggsave(g, file = file.path(fig_dir, "05_01_numSites_in_overlaps.pdf"), height = 6, width = 7)

dim(plot_dat[plot_dat$num_sites > 1,]) # 816
(g <- ggplot(plot_dat[plot_dat$num_sites > 1,], aes(x = num_sites)) + geom_bar() + theme_bw() +
    scale_x_continuous(minor_breaks = seq(0,220,10)) +
    ggtitle("Number of RNA-editing sites in identified overlaps (816 with >1 site)"))
ggsave(g, file = file.path(fig_dir, "05_01_numSites_in_overlaps_moreThanOne.png"), height = 6, width = 7)
ggsave(g, file = file.path(fig_dir, "05_01_numSites_in_overlaps_moreThanOne.pdf"), height = 6, width = 7)


## Make data.frame filtered with just these tx overlaps that contain an RNA editing site
cisnat_withSites <- putative_cisnat[plot_dat$put_cisnat_name, ]
head(cisnat_withSites)
# how many unique gene pairs is this?
length(unique(paste(cisnat_withSites$minus_gene, cisnat_withSites$plus_gene, sep = "_"))) # 486 gene pairs with editing site
length(unique(paste(putative_cisnat$minus_gene, putative_cisnat$plus_gene, sep = "_"))) # 6,155 gene pairs in original putative list
# how many overlap a possible IR-Alu?
summary(is.na(cisnat_withSites$alu_overlap)) # TRUE for 902 (no Alu overlap), FALSE for 181 (Alu overlap)
# to keep in mind - some of these might be IRAlu driven instead of cis-NATs
# reflects the bias of the editing site annotation file

# add site info to df
all(plot_dat$put_cisnat_name == rownames(cisnat_withSites)) # TRUE
rownames(plot_dat) <- plot_dat$put_cisnat_name
cisnat_withSites$sites <- NA
for (i in rownames(cisnat_withSites)) {
  cisnat_withSites[i, ]$sites <- paste(site_grange[site_ind[put_cisnat_ind == plot_dat[i, ]$put_cisnat_ind]]$site_name, collapse = ";")
}


####
## Check 2: how many overlap regions include an editing cluster
cluster_dat$cluster_start <- as.integer(unlist(lapply(strsplit(cluster_dat$cluster_name, "_"), "[[", 2)))
cluster_dat$cluster_end <- as.integer(unlist(lapply(strsplit(cluster_dat$cluster_name, "_"), "[[", 3)))
cluster_unique <- cluster_dat[-which(duplicated(cluster_dat$cluster_name)), c("chr", "cluster_name", "cluster_start", "cluster_end")]
rownames(cluster_unique) <- cluster_unique$cluster_name
dim(cluster_unique) # 168,623 clusters
cl_grange <- makeGRangesFromDataFrame(cluster_unique, seqnames.field = "chr", start.field = "cluster_start", 
                                      end.field = "cluster_end", keep.extra.columns = T)

put_cisnat_to_cluster <- findOverlaps(put_cisnat_grange, cl_grange, select = "all", ignore.strand = T, minoverlap = 1)
put_cisnat_cl_ind <- queryHits(put_cisnat_to_cluster)
cluster_ind <- subjectHits(put_cisnat_to_cluster)
# sanity check:
put_cisnat_grange[53]
cl_grange[cluster_ind[put_cisnat_cl_ind == 53]]

length(unique(put_cisnat_cl_ind)) # 687 of the overlaps contain RNA editing clusters

## Plot histogram of number of editing sites/overlap
# (was originally going to have a 0-overlaps bar, but that will dwarf the plot)
plot_dat <- data.frame(put_cisnat_ind = as.integer(names(table(put_cisnat_cl_ind))),
                       num_cluster = as.integer(table(put_cisnat_cl_ind)))
plot_dat$put_cisnat_name <- names(put_cisnat_grange[plot_dat$put_cisnat_ind])
(g <- ggplot(plot_dat, aes(x = num_cluster)) + geom_bar() + theme_bw() +
    ggtitle("Number of RNA-editing clusters in identified overlaps (687 include a cluster)"))
ggsave(g, file = file.path(fig_dir, "05_02_numCluster_in_overlaps.png"), height = 6, width = 6)
ggsave(g, file = file.path(fig_dir, "05_02_numCluster_in_overlaps.pdf"), height = 6, width = 6)

# (filter to just these ones)
cisnat_withCluster <- putative_cisnat[plot_dat$put_cisnat_name, ]
head(cisnat_withCluster)
# how many unique gene pairs is this?
length(unique(paste(cisnat_withCluster$minus_gene, cisnat_withCluster$plus_gene, sep = "_"))) # 298
# ^ compared to 486 gene pairs with at least one editing site, and 6,155 gene pairs in original putative list

# sanity check that the ones which include clusters are also in the cisnat_withSites df (if no site overlap, then defs no cluster):
summary(rownames(cisnat_withCluster) %in% rownames(cisnat_withSites)) # all are - good

# add cluster info
all(plot_dat$put_cisnat_name == rownames(cisnat_withCluster)) # TRUE
rownames(plot_dat) <- plot_dat$put_cisnat_name
cisnat_withCluster$clusters <- NA
for (i in rownames(cisnat_withCluster)) {
  cisnat_withCluster[i, ]$clusters <- paste(names(cl_grange[cluster_ind[put_cisnat_cl_ind == plot_dat[i, ]$put_cisnat_ind]]), collapse = ";")
}



####
## what are the gene types of cisNATs overlapping RNA editing sites?
# First, define tx types for clearer plotting:
tx_priority_full <- c("protein_coding", "lncRNA", "processed_transcript", "retained_intron", "nonsense_mediated_decay", 
                      "non_stop_decay", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", 
                      "transcribed_unprocessed_pseudogene", "processed_pseudogene", "unitary_pseudogene", 
                      "polymorphic_pseudogene", "unprocessed_pseudogene", "IG_V_pseudogene", "TR_C_gene",
                      "scaRNA", "miRNA", "snRNA", "snoRNA", "rRNA_pseudogene", "misc_RNA", "TEC")
tx_priority <- c("protein_coding", "lncRNA", "processed_transcript", "retained_intron", "nonsense_mediated_decay", 
                 "non_stop_decay", "pseudogene", "small_RNA", "other")
## Transcript-level
nrow(plot_dat <- cisnat_withSites) # 1,083 unique tx pairs
# get pair types:
plot_dat <- get_pair_type(df = plot_dat, pair_type = "tx", tx_priority = tx_priority)
length(unique(plot_dat$xcat_tx)) # 27
unique(plot_dat$xcat_tx)
# plot:
xlab_lev <- names(table(plot_dat$xcat_tx)[order(table(plot_dat$xcat_tx), decreasing = T)])
plot_dat$xcat_tx <- factor(plot_dat$xcat_tx, levels = xlab_lev)
labeltable <- as.data.frame(table(plot_dat$xcat_tx))
colnames(labeltable)[1] <- "xcat_tx"
g <- ggplot(plot_dat, aes(x = xcat_tx)) + geom_bar() + 
  geom_text(data = labeltable, aes(y = Freq, label = Freq), vjust = -0.3, size = 3) + 
  ylab(NULL) + xlab(NULL) + theme_bw() + 
  ggtitle(paste("Transcript pair types (number of unique overlap regions =", sum(labeltable$Freq), ")")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(g, file = file.path(fig_dir, "05_03_overlap_regions_withSites_transcript_types.png"), width = 7, height = 6)
ggsave(g, file = file.path(fig_dir, "05_03_overlap_regions_withSites_transcript_types.pdf"), width = 7, height = 6)

## Unique gene pairs, gene type
plot_dat <- cisnat_withSites
# filter to unique gene pairs (putative_cisnat is unique tx_pairs)
plot_dat$gene_pair <- paste(plot_dat$minus_gene, plot_dat$plus_gene, sep = "__")
nrow(plot_dat <- plot_dat[!duplicated(plot_dat$gene_pair), ]) # 486 unique gene pairs
# get pair types:
plot_dat <- get_pair_type(df = plot_dat, pair_type = "gene", tx_priority = tx_priority)
length(unique(plot_dat$xcat_gene)) # 11
unique(plot_dat$xcat_gene)
# plot:
xlab_lev <- names(table(plot_dat$xcat_gene)[order(table(plot_dat$xcat_gene), decreasing = T)])
plot_dat$xcat_gene <- factor(plot_dat$xcat_gene, levels = xlab_lev)
labeltable <- as.data.frame(table(plot_dat$xcat_gene))
colnames(labeltable)[1] <- "xcat_gene"
g <- ggplot(plot_dat, aes(x = xcat_gene)) + geom_bar() + 
  geom_text(data = labeltable, aes(y = Freq, label = Freq), vjust = -0.3, size = 3) + 
  ylab(NULL) + xlab(NULL) + theme_bw() + 
  ggtitle(paste("Gene pair types (number of unique gene pairs with overlap regions =", sum(labeltable$Freq), ")")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(g, file = file.path(fig_dir, "05_04_overlap_regions_withSites_gene_types.png"), width = 6, height = 6)
ggsave(g, file = file.path(fig_dir, "05_04_overlap_regions_withSites_gene_types.pdf"), width = 6, height = 6)



####
## Plot length of overlaps that include editing sites
# Overlap length (in plot, partition length by number of editing sites in region)
cisnat_withSites$n_sites <- unlist(lapply(strsplit(cisnat_withSites$sites, ";"), "length"))
plot_dat <- cisnat_withSites
plot_dat$n_sites_cat <- "1"
plot_dat[plot_dat$n_sites > 1 & plot_dat$n_sites <= 5, ]$n_sites_cat <- "2-5"
plot_dat[plot_dat$n_sites > 5 & plot_dat$n_sites <= 10, ]$n_sites_cat <- "6-10"
plot_dat[plot_dat$n_sites > 10 & plot_dat$n_sites <= 25, ]$n_sites_cat <- "11-25"
plot_dat[plot_dat$n_sites > 25 & plot_dat$n_sites <= 50, ]$n_sites_cat <- "26-50"
plot_dat[plot_dat$n_sites > 50 & plot_dat$n_sites <= 100, ]$n_sites_cat <- "51-100"
plot_dat[plot_dat$n_sites > 100, ]$n_sites_cat <- "101+"
plot_dat$n_sites_cat <- factor(plot_dat$n_sites_cat, levels = rev(c("1", "2-5", "6-10", "11-25", "26-50", "51-100", "101+")))
plot_dat$bp_range <- "100-200"
plot_dat[plot_dat$longest_overlap_width > 200 & plot_dat$longest_overlap_width <= 300, ]$bp_range <- "201-300"
plot_dat[plot_dat$longest_overlap_width > 300 & plot_dat$longest_overlap_width <= 400, ]$bp_range <- "301-400"
plot_dat[plot_dat$longest_overlap_width > 400 & plot_dat$longest_overlap_width <= 500, ]$bp_range <- "401-500"
plot_dat[plot_dat$longest_overlap_width > 500 & plot_dat$longest_overlap_width <= 600, ]$bp_range <- "501-600"
plot_dat[plot_dat$longest_overlap_width > 600 & plot_dat$longest_overlap_width <= 700, ]$bp_range <- "601-700"
plot_dat[plot_dat$longest_overlap_width > 700 & plot_dat$longest_overlap_width <= 800, ]$bp_range <- "701-800"
plot_dat[plot_dat$longest_overlap_width > 800 & plot_dat$longest_overlap_width <= 900, ]$bp_range <- "801-900"
plot_dat[plot_dat$longest_overlap_width > 900 & plot_dat$longest_overlap_width <= 1000, ]$bp_range <- "901-1000"
plot_dat[plot_dat$longest_overlap_width > 1000 & plot_dat$longest_overlap_width <= 1500, ]$bp_range <- "1001-1500"
plot_dat[plot_dat$longest_overlap_width > 1500 & plot_dat$longest_overlap_width <= 2000, ]$bp_range <- "1501-2000"
plot_dat[plot_dat$longest_overlap_width > 2000 & plot_dat$longest_overlap_width <= 2500, ]$bp_range <- "2001-2500"
plot_dat[plot_dat$longest_overlap_width > 2500 & plot_dat$longest_overlap_width <= 3000, ]$bp_range <- "2501-3000"
plot_dat[plot_dat$longest_overlap_width > 3000 & plot_dat$longest_overlap_width <= 4000, ]$bp_range <- "3001-4000"
plot_dat[plot_dat$longest_overlap_width > 4000 & plot_dat$longest_overlap_width <= 5000, ]$bp_range <- "4001-5000"
plot_dat[plot_dat$longest_overlap_width > 5000 & plot_dat$longest_overlap_width <= 6000, ]$bp_range <- "5001-6000"
plot_dat[plot_dat$longest_overlap_width > 6000, ]$bp_range <- "6001+"
plot_dat$bp_range <- factor(plot_dat$bp_range, levels = c("100-200", "201-300", "301-400", "401-500", "501-600",
                                                          "601-700", "701-800", "801-900", "901-1000", "1001-1500",
                                                          "1501-2000", "2001-2500", "2501-3000", "3001-4000",
                                                          "4001-5000", "5001-6000", "6001+"))
table(plot_dat$bp_range, plot_dat$n_sites_cat)
g <- ggplot(plot_dat, aes(x = bp_range, fill = n_sites_cat)) + geom_bar() + 
  theme_bw() + xlab("overlap width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(g, file = file.path(fig_dir, "05_05_number_at_overlapLengthThresholds_withSites.png"), width = 6, height = 5)
ggsave(g, file = file.path(fig_dir, "05_05_number_at_overlapLengthThresholds_withSites.pdf"), width = 6, height = 5)



####
## Save the cis-NATs with RNA-editing overlap
save(cisnat_withSites, cisnat_withCluster, file = file.path(out_dir, "05_cisNAT_with_RNAedit.RData"))


####
## Save informative tables of the 828 overlaps that include an RNA-editing site
cisnat_withSites$num_edit_sites <- unlist(lapply(strsplit(cisnat_withSites$sites, ";"), length))
to_write <- cisnat_withSites[,-c(12,13,18,19)]
to_write$edit_cluster <- cisnat_withCluster[rownames(to_write), ]$clusters
write.table(to_write, file = file.path(out_dir, "05_cisNAT_with_RNAedit.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")


####
## Save text file with all editing sites in at least one transcript overlap

# recall that site_name in cisnat_withSites came from indiv_quant; position is 
# the lower value. this will match what is in each individual quant file.

sites_at_quant_pos <- unique(unlist(strsplit(cisnat_withSites$sites, ";")))
length(sites_at_quant_pos) # 5,604 sites
summary(sites_at_quant_pos %in% indiv_quant$site_name) # sanity check: all are in the indiv_quant file
write(sites_at_quant_pos, file = file.path(out_dir, "05_sites_in_cisnat_quantFilePos.txt"), sep = "\n")

