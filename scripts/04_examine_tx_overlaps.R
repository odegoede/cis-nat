#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 04: Explore putative cis-NATs
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/04_examine_tx_overlaps.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --edittable data/gtex_editing_paper_tableS4.xlsx --coloctable data/gtex_lncrna_paper_hits_summary_table_allGeneTypes.RData


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--outdir", default = "./output", 
              help = "Output directory [default is \"%default\"]"),
  make_option("--figdir", default = "./figures", 
              help = "Figure directory [default is \"%default\"]"),
  make_option("--scriptdir", default = "./scripts",
              help = "Scripts directory (to load functions saved in source scripts). [default is \"%default\"]"),
  make_option("--scratchdir", default = NULL, 
              help = "Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default \"%default\"]"),
  make_option("--edittable", default = NULL,
              help = "File of cis-NATs from RNA-editing paper (optional) [default \"%default\"]"),
  make_option("--coloctable", default = NULL,
              help = "File of QTL colocalization from GTEx lncRNA paper (optional) [default \"%default\"]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## INPUT TESTS
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
# n/a


####
## READ IN INPUT FILES
load(file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData")) # object name from script 03: putative_cisnat
load(file.path(out_dir, "gene_tx_exon_anno_files.RData")) # object name from script 01: gene_anno, tx_anno, and exon_anno


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
## Found missing annotation for some transcripts in putative cisnat; update and rewrite
missing_minus <- unique(putative_cisnat[is.na(putative_cisnat$minus_tx_type),]$minus_tx)
for (t in missing_minus) {
  putative_cisnat[putative_cisnat$minus_tx == t, ]$minus_tx_type <- tx_anno[t, ]$tx_type
}
missing_plus <- unique(putative_cisnat[is.na(putative_cisnat$plus_tx_type),]$plus_tx)
for (t in missing_plus) {
  putative_cisnat[putative_cisnat$plus_tx == t, ]$plus_tx_type <- tx_anno[t, ]$tx_type
}

save(putative_cisnat, file = file.path(opt$overlap))


####
## What type of transcript/gene pair are these overlaps?

## First, define tx types for clearer plotting:
tx_priority_full <- c("protein_coding", "lncRNA", "processed_transcript", "retained_intron", "nonsense_mediated_decay", 
                      "non_stop_decay", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", 
                      "transcribed_unprocessed_pseudogene", "processed_pseudogene", "unitary_pseudogene", 
                      "polymorphic_pseudogene", "unprocessed_pseudogene", "IG_V_pseudogene", "TR_C_gene",
                      "scaRNA", "miRNA", "snRNA", "snoRNA", "rRNA_pseudogene", "misc_RNA", "TEC")
tx_priority <- c("protein_coding", "lncRNA", "processed_transcript", "retained_intron", "nonsense_mediated_decay", 
                 "non_stop_decay", "pseudogene", "small_RNA", "other")

## Unique tx_pairs, transcript type
plot_dat <- putative_cisnat
# get pair types:
plot_dat <- get_pair_type(df = plot_dat, pair_type = "tx", tx_priority = tx_priority)
length(unique(plot_dat$xcat_tx)) # 37
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
ggsave(g, file = file.path(fig_dir, "04_01_overlap_regions_transcript_types.png"), width = 7, height = 6)
ggsave(g, file = file.path(fig_dir, "04_01_overlap_regions_transcript_types.pdf"), width = 7, height = 6)

## Unique gene pairs, gene type
plot_dat <- putative_cisnat
# filter to unique gene pairs (putative_cisnat is unique tx_pairs)
plot_dat$gene_pair <- paste(plot_dat$minus_gene, plot_dat$plus_gene, sep = "__")
plot_dat <- plot_dat[!duplicated(plot_dat$gene_pair), ]
# get pair types:
plot_dat <- get_pair_type(df = plot_dat, pair_type = "gene", tx_priority = tx_priority)
length(unique(plot_dat$xcat_gene)) # 15
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
ggsave(g, file = file.path(fig_dir, "04_02_overlap_regions_gene_types.png"), width = 6, height = 6)
ggsave(g, file = file.path(fig_dir, "04_02_overlap_regions_gene_types.pdf"), width = 6, height = 6)


####
## How long are the overlaps?
# Note: initial filter setting was minimum overlap = 100
plot_dat <- putative_cisnat

g <- ggplot(plot_dat, aes(x = longest_overlap_width)) + geom_histogram(binwidth = 50) +
    geom_vline(xintercept = 200, color = "red") +
    theme_bw() + ggtitle(paste("Length of", nrow(plot_dat), "unique regions, binwidth = 50"))
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_fullX.png"), width = 6, height = 4)
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_fullX.pdf"), width = 6, height = 4)

g <- ggplot(plot_dat, aes(x = longest_overlap_width)) + geom_histogram(binwidth = 10) +
    coord_cartesian(xlim = c(0,2000)) + theme_bw() + geom_vline(xintercept = 200, color = "red") +
    ggtitle(paste("Length of", nrow(plot_dat[plot_dat$longest_overlap_width <= 2000, ]), "unique regions, binwidth = 10, X-axis crop at 2000"))
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_cropX.png"), width = 6, height = 4)
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_cropX.pdf"), width = 6, height = 4)

g <- ggplot(plot_dat, aes(x = longest_overlap_width)) + geom_histogram(binwidth = 10) +
    coord_cartesian(xlim = c(0,1000)) + theme_bw() + geom_vline(xintercept = 200, color = "red") +
    ggtitle(paste("Length of", nrow(plot_dat[plot_dat$longest_overlap_width <= 1000, ]), "unique regions, binwidth = 10, X-axis crop at 1000"))
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_verycropX.png"), width = 6, height = 4)
ggsave(g, file = file.path(fig_dir, "04_03_lengthOfOverlaps_verycropX.pdf"), width = 6, height = 4)


####
## How many overlaps are there at ___ bp threshold?
plot_dat <- data.frame(bp_thresh = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000,
                                     2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000),
                       num_regions = NA)
for (i in c(1:nrow(plot_dat))) {
  plot_dat[i, ]$num_regions <- nrow(putative_cisnat[putative_cisnat$longest_overlap_width >= plot_dat[i,]$bp_thresh, ])
}
plot_dat$prop_of_overlaps <- plot_dat$num_regions/plot_dat[plot_dat$bp_thresh == 100, ]$num_regions
g <- ggplot(plot_dat, aes(x = as.factor(bp_thresh), y = prop_of_overlaps)) + geom_col() + 
  theme_bw() + geom_text(aes(y = prop_of_overlaps, label = num_regions), vjust = -0.3, size = 3) +
  xlab("overlap length (bp) threshold") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(g, file = file.path(fig_dir, "04_04_number_at_overlapLengthThresholds.png"), width = 6, height = 5)
ggsave(g, file = file.path(fig_dir, "04_04_number_at_overlapLengthThresholds.pdf"), width = 6, height = 5)


####
## Chromosomal representation
# compare the total number of genes on each chromosome (background) to the total number of overlaps on each chromosome (hits)
# at the gene-level, so filter overlaps to just unique gene pairs
temp <- putative_cisnat
temp$gene_pair <- paste(temp$plus_gene, temp$minus_gene, sep = "__")
temp <- temp[!duplicated(temp$gene_pair), ]
plot_dat <- data.frame(cat = c(rep("cisnat", length(unique(temp$chr))),
                               rep("genes", length(unique(gene_anno$chr)))),
                       chr = c(names(table(temp$chr)),
                               names(table(gene_anno$chr))),
                       num_regions = c(as.numeric(table(temp$chr)),
                                       as.numeric(table(gene_anno$chr))))
plot_dat$chr <- factor(plot_dat$chr, levels = c(paste0("chr", c(1:22)), "chrX", "chrY", "chrM"))
g <- ggplot(plot_dat, aes(x = chr, y = num_regions, fill = cat)) + geom_col() + 
    facet_wrap(~cat, scales = "free_y") + theme_bw() + 
    geom_text(aes(y = num_regions, label = num_regions), vjust = -0.3, size = 3) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(g, file = file.path(fig_dir, "04_05_number_overlaps_by_chromosome_v1number.png"), width = 10, height = 6)
ggsave(g, file = file.path(fig_dir, "04_05_number_overlaps_by_chromosome_v1number.pdf"), width = 10, height = 6)

plot_dat$prop_regions <- NA
plot_dat[plot_dat$cat == "cisnat", ]$prop_regions <- plot_dat[plot_dat$cat == "cisnat", ]$num_regions/sum(plot_dat[plot_dat$cat == "cisnat", ]$num_regions)
plot_dat[plot_dat$cat == "genes", ]$prop_regions <- plot_dat[plot_dat$cat == "genes", ]$num_regions/sum(plot_dat[plot_dat$cat == "genes", ]$num_regions)
g <- ggplot(plot_dat, aes(x = chr, y = prop_regions, fill = cat)) + geom_col(position = position_dodge()) + 
    theme_bw() + geom_text(aes(y = prop_regions, label = num_regions, color = cat), vjust = -0.3) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(g, file = file.path(fig_dir, "04_05_number_overlaps_by_chromosome_v2proportion.png"), width = 10, height = 6)
ggsave(g, file = file.path(fig_dir, "04_05_number_overlaps_by_chromosome_v2proportion.pdf"), width = 10, height = 6)


####
## Optional checks: How many of these possible cis-NATs overlap with my previously identified trait-relevant lncRNAs? (from GTEx colocalization)
# just summary of how many - will go into detail on these in another script
if (!is.null(opt$coloctable)) {
  load(file.path(opt$coloctable))
  # filter to lncRNAs only
  coloc_lnc <- results_df[results_df$gene_type == "lncrna", ]
  coloc_genes <- unique(c(unique(coloc_lnc[coloc_lnc$hit_type == "eQTL", ]$feature), unique(coloc_lnc[coloc_lnc$hit_type == "sQTL", ]$sqtl_gene)))
  summary(coloc_genes %in% rownames(gene_anno)) # 79 are, 51 aren't
  summary(gsub("\\..*", "", coloc_genes) %in% gsub("\\..*", "", rownames(gene_anno))) # when removing trailing .## part of engene ID: 121 are, 9 aren't
  # after manual search, the 9 absent genes appear to be affiliated with GRCh37. Some have disappeared entirely 
  # (e.g. RP11-876N24.4 / ENSG00000262222.1), others might now be another gene (GTEx's LINC00593 / ENSG00000259703.5 
  # could have become DRAIC / ENSG00000245750.10)
  to_check <- gsub("\\..*", "", coloc_genes) # 146 genes
  length(txs <- rownames(tx_anno[gsub("\\..*", "", tx_anno$ensgene) %in% to_check, ])) # 603 transcripts
  length(txs <- txs[which(txs %in% unique(c(putative_cisnat$minus_tx, putative_cisnat$plus_tx)))])
  # ^ 103 transcripts that are colocalized with some trait and are also in putative cisNATs
  length(unique(tx_anno[txs, ]$ensgene)) # correspond to 49 genes
}
