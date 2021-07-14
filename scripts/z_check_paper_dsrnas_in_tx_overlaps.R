#!/usr/bin/env Rscript

#####
## cis-NAT project
## Exploratory script 01: Examine overlap between annotation-based overlaps and immunogenic dsRNAs from RNA-editing paper
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/exA_check_paper_dsrnas_in_tx_overlaps.R --overlap output/03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData --annofile output/gene_tx_exon_anno_files.RData --edittable data/gtex_editing_paper_tableS4.xlsx --gtex_gtf data/gencode.v26.GRCh38.genes.gtf --outdir output/ --figdir figures/exploratory/ --scriptdir scripts/


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicRanges))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--overlap", default = NULL, 
              help = "RData file of transcript overlaps [default \"%default\"]"),
  make_option("--annofile", default = NULL,
              help = "File with genomic annotation of gene/tx/exon annotation [default \"%default\"]"),
  make_option("--edittable", default = NULL,
              help = "File of cis-NATs from RNA-editing paper (optional) [default \"%default\"]"),
  make_option("--gtex_gtf", default = NULL,
              help = "File with GTEx's GTF [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory to write transcript overlaps to [default is \"%default\"]"),
  make_option("--figdir", default = "./figures", 
              help = "Figure directory to put figures in [default is \"%default\"]"),
  make_option("--scriptdir", default = "./scripts", 
              help = "Script directory to get source code from [default is \"%default\"]")
)

# get command line options; if help option encountered, print help and exit;
# otherwise if options not found on command line, then set defaults.
opt <- parse_args(OptionParser(option_list=option_list))


####
## DEFINE FUNCTIONS

# confirm_no_overlap() checks a candidate gene and makes sure no tx are overlapping on opposite strand
# if no overlap (expected), returns nothing; if overlap, returns gene name
confirm_no_overlap <- function(g, exon_df = exon_anno, verbose = F) {
  gene_of_interest <- makeGRangesFromDataFrame(exon_df[exon_df$symbol == g, ], keep.extra.columns = T)
  rest_of_chr <- makeGRangesFromDataFrame(exon_df[exon_df$symbol != g & exon_df$chr == levels(seqnames(gene_of_interest)), ],
                                          keep.extra.columns = T)
  # look for findOverlaps that are found with ignore.strand = T, but not found with ignore.strand = F
  all_strands <- findOverlaps(gene_of_interest, rest_of_chr, ignore.strand = T)
  same_strand <- findOverlaps(gene_of_interest, rest_of_chr, ignore.strand = F)
  if (verbose) {
    print(paste("gene is:", g))
    print(paste("number of either-strand overlaps:", length(all_strands)))
    print(paste("number of same-strand overlaps:", length(same_strand)))
  }
  if (length(all_strands) - length(same_strand) == 0) {
    print(paste("** Gene has no opp-strand overlaps:", g))
  } else {
    print(paste("** Gene has overlaps!:", g))
    if (verbose) {
      all_matches <- paste(queryHits(all_strands), subjectHits(all_strands), sep = "__")
      same_matches <- paste(queryHits(same_strand), subjectHits(same_strand), sep = "__")
      opp_strand <- all_strands[!all_matches %in% same_matches]
      print(gene_of_interest[unique(queryHits(opp_strand))])
      print(rest_of_chr[unique(subjectHits(opp_strand))])
    }
  }
}



# plot_nearby() makes a really basic plot with boxes for exons; coloring is gene name; y-axis position is unique transcript ID
# (negative y value = minus strand, positive = plus strand)
plot_nearby <- function(g, file_prefix = "exA_nearbyGenes_", extender = NULL) {
  file_png <- file.path(fig_dir, paste0(file_prefix, g, ".png"))
  file_pdf <- file.path(fig_dir, paste0(file_prefix, g, ".pdf"))
  if (is.null(extender)) {
    extender <- 2 * (max(tx_anno[tx_anno$symbol == g, ]$end) - min(tx_anno[tx_anno$symbol == g, ]$start))
  }
  range_start <- min(tx_anno[tx_anno$symbol == g, ]$start) - extender
  range_end <- max(tx_anno[tx_anno$symbol == g, ]$end) + extender
  plot_dat <- exon_anno[exon_anno$chr == unique(exon_anno[exon_anno$symbol == g, ]$chr), ]
  plot_dat <- plot_dat[plot_dat$start %in% c(range_start:range_end) | plot_dat$end %in% c(range_start:range_end), ]
  plot_dat$y_val <- NA
  plot_dat[plot_dat$strand == "-", ]$y_val <- -1 * (as.numeric(factor(plot_dat[plot_dat$strand == "-", ]$tx_id)))
  plot_dat[plot_dat$strand == "+", ]$y_val <- as.numeric(factor(plot_dat[plot_dat$strand == "+", ]$tx_id))
  g <- ggplot(plot_dat) + theme_bw() + ggtitle(paste(g)) + geom_hline(yintercept = 0, col = "black", lwd = 0.75) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = (y_val-0.5), ymax = (y_val+0.5), fill = symbol))
  ggsave(g, filename = file_png, width = 7, height = 7)
  ggsave(g, filename = file_pdf, width = 7, height = 7)
}


####
## INPUT TESTS
# Check the required arguments (overlaps file) are provided
if (is.null(opt$overlap)) { 
  stop("Overlaps file not provided, exiting\n") 
}
if (is.null(opt$annofile)) { 
  stop("Annotation file not provided, exiting\n") 
}
if (is.null(opt$gtex_gtf)) { 
  stop("GTEx GTF file not provided, exiting\n") 
}
if (is.null(opt$edittable)) { 
  stop("RNA editing table not provided, exiting\n") 
}

# Create consistent directory variables
out_dir <- file.path(opt$outdir)
fig_dir <- file.path(opt$figdir)
script_dir <- file.path(opt$scriptdir)
# Check that directory exists
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist, exiting\n")
}
if (!dir.exists(fig_dir)) {
  stop("Figure directory does not exist, exiting\n")
}
if (!dir.exists(script_dir)) {
  stop("Script directory does not exist, exiting\n")
}


####
## LOAD SOURCE SCRIPTS
source(file.path(script_dir, "source_gtf_reader.R"))
source(file.path(script_dir, "source_temp_unzip.R"))


####
## READ IN INPUT FILES
load(file.path(opt$overlap)) # object name from script 03: putative_cisnat
load(file.path(opt$annofile)) # object name from script 01: gene_anno, tx_anno, and exon_anno
coloc_cisnat <- as.data.frame(readxl::read_xlsx(file.path(opt$edittable)))


####
## Summarize number of immune-trait associated dsRNAs from GTEx RNA-editing paper that are in putative_cisnat

# In the supplementary table of the paper, cis-NATs are provided as comma-separated gene symbols
paper_set <- unique(coloc_cisnat$`cis-NATs`)
paper_genes <- unique(unlist(strsplit(paper_set, ",")))

summary(paper_genes %in% unique(gene_anno$symbol))
# ^ 59 immunogenic dsRNA genes are in the gencode v35 annotation, 4 aren't

# manually search for and fill in the paper_genes that don't have the same symbol in gencode v35
# TODO: could do this through online search (get gencode v26, as in GTEx, get ensgene_id, check ensgene_id in gencode v35)
temp <- paper_genes[!(paper_genes %in% unique(gene_anno$symbol))]
# missing: C5orf56, FAM69A, AHSA2, LINC00094
# online search (manual): C5orf56 = IRF1-AS1, FAM69A = DIPK1A, AHSA2 = ENSG00000173209, LINC00094 = BRD3OS
c("IRF1-AS1", "DIPK1A", "BRD3OS") %in% unique(gene_anno$symbol) # all are
gene_anno[grepl("ENSG00000173209", gene_anno$ensgene), ] # AHSA2 = AHSA2P
paper_genes <- gsub("C5orf56", "IRF1-AS1", gsub("FAM69A", "DIPK1A", gsub("LINC00094", "BRD3OS", gsub("AHSA2", "AHSA2P", paper_genes))))
summary(paper_genes %in% gene_anno$symbol) # now all 63 are in the annotation

# how many paper_genes are in the putative_cisnat table?
summary(paper_genes %in% c(unique(putative_cisnat$minus_symbol), unique(putative_cisnat$plus_symbol)))
# 46 paper-identified immunogenic dsRNA genes are present in this set of possible cis-NATs; 17 are not
shared_cisnat <- paper_genes[(paper_genes %in% c(unique(putative_cisnat$minus_symbol), unique(putative_cisnat$plus_symbol)))]
not_found <- paper_genes[!(paper_genes %in% c(unique(putative_cisnat$minus_symbol), unique(putative_cisnat$plus_symbol)))]
out_df <- data.frame(editPaper_dsrna_gene = c(shared_cisnat, not_found), 
                     in_putative_cisnat_table = c(rep(TRUE, length(shared_cisnat)), rep(FALSE, length(not_found))))
write.table(out_df, file = file.path(out_dir, "exA_immunogenic_dsrna_status.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


####
## What's up with these "not_found" genes?
# (all of the cis-NATs in the supplementary table should be detected as tx overlaps from the annotation)
row_ind <- c()
for (i in not_found) {row_ind <- c(row_ind, grep(i, coloc_cisnat$`cis-NATs`))}
not_found_rows <- coloc_cisnat[row_ind, ]
table(not_found_rows$`cis-NATs`) # all single genes, except for ANXA9/MINDY1 pair (ANXA9 is also a standalone cis-NAT)
# (the ones that are doubles would give a hint of what specific overlap to check for - for the singles, it's just 
# checking if any genes are close / could be overlapping)

# loop through each gene, make a gene neighborhood picture
for (gene in not_found) {
  confirm_no_overlap(g = gene)
  plot_nearby(g = gene)
}

# four genes had overlaps: MCAM, TCTA, DAP3, and STAT6
for (gene in c("MCAM", "TCTA", "DAP3", "STAT6")) {
  confirm_no_overlap(g = gene, verbose = T)
}
# looks like 3 were missed due to the TSL thresholds, and 1 was missed due to length of overlap
not_found <- not_found[!not_found %in% c("MCAM", "TCTA", "DAP3", "STAT6")]


####
## TODO: make same checks for not_found, but using GTEx gene annotation
gtex_anno <- read_gtf(filename = opt$gtex_gtf)
# Remove unnecessary columns and mitochondrial chromosome
gtex_anno <- gtex_anno[,-c(2,6,8)] # these are unnecessary fields
gtex_anno <- gtex_anno[gtex_anno$chr != "ChrM", ]
# get the exon information to run "confirm_no_overlap()"
gtex_anno <- get_basic_attr(gtex_anno)
# filter gtex_anno to exon rows, add extra transcript information and exon number
gtex_exon <- get_exon_attr(gtex_anno)
gtex_exon <- gtex_exon[,-grep("attribute", colnames(gtex_exon))]

summary(not_found %in% gtex_anno$symbol)

# as with the GENCODE v35 annotation, loop through each gene double check for overlaps in the GTEx gtf
for (gene in not_found) {
  confirm_no_overlap(g = gene, exon_df = gtex_exon, verbose = T)
}
# genes with overlaps in full anno: CFAP126
# ^ manual checks show that this is just because of coordinates changing between versions
# 12 genes still don't show overlaps: "ANXA9"   "MINDY1"  "EHBP1L1" "EML1"    "IL6R"    "GNA12"   "AMZ1"    "SWAP70"  "IFIT3"   "FUT2"    "QRSL1"   "TYMP"  

# Conclusion: will need to check with Qin how he identified whether a dsRNA was a cis-NAT


