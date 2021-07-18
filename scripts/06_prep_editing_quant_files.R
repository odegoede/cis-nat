#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 06: Prep files to check editing quantification at tx overlaps
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/06_check_editing_quantification.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --editcoverage data/gtex_edit/coverage_matrix.txt.gz --editlevel data/gtex_edit/editLevel_matrix.txt.gz --tpmfile data/cisnat_gene_tpm.gct.gz --gtexanno data/rnaseqc_genes_info.RData

## Idea:
# want to get an idea of samples' coverage and RNA editing levels in these regions
# where RNA editing was assessed.
# (no point in keeping regions where no individuals or tissue types have coverage or variation in RNA editing).



####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
options(datatable.fread.datatable = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--samplefile", default = NULL,
              help = "File with GTEx sample annotation [default \"%default\"]"),
  make_option("--editcoverage", default = NULL, 
              help = "File with coverage for the editing sites of interest"),
  make_option("--editlevel", default = NULL, 
              help = "File with editLevels for the editing sites of interest"),
  make_option("--tpmfile", default = NULL, 
              help = "File with gene expression (TPM) values for genes of interest"),
  make_option("--gtexanno", default = NULL, 
              help = "GTEx gene annotation file"),
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
# Check the required arguments (overlaps file) are provided
if (is.null(opt$editanno)) { 
  stop("Editing annotation file not provided, exiting\n") 
}
if (is.null(opt$editcoverage)) { 
  stop("Editing coverage file not provided, exiting\n") 
}
if (is.null(opt$editlevel)) { 
  stop("Editing levels file not provided, exiting\n") 
}
if (is.null(opt$tpmfile)) { 
  stop("TPM file not provided, exiting\n") 
}
if (is.null(opt$gtexanno)) { 
  stop("GTEx annotation file not provided, exiting\n") 
}
if (is.null(opt$samplefile)) { 
  stop("Sample info file not provided, exiting\n") 
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
load(file.path(opt$samplefile)) # object from previous GTEx project, copied into this project's data directory: sample.info
sample.info <- sample.info[sample.info$ANALYTE_TYPE == "RNA:Total RNA", ]
load(file.path(out_dir, "gene_tx_exon_anno_files.RData")) # object name from script 01: gene_anno, tx_anno, and exon_anno
load(file.path("data/gtex_edit/All.AG.stranded.annovar.RData")) # RData version of editing site annotation, saved in script 05: cluster_dat
load(file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData")) # object name from script 03: putative_cisnat
load(file.path(out_dir, "05_cisNAT_with_RNAedit.RData")) # two objects from script 05, named: cisnat_withSites, cisnat_withCluster
load(file.path(opt$gtexanno))

cov_file <- file.path(opt$editcoverage)
if (endsWith(cov_file, "gz") | endsWith(cov_file, "zip") | endsWith(cov_file, "bzip2") | endsWith(cov_file, "xz")) {
  coverage_dat <- suppressWarnings(temp_unzip(cov_file, fread, data.table = F))
} else {
  coverage_dat <- suppressWarnings(fread(cov_file, data.table = F))
}
rownames(coverage_dat) <- coverage_dat$site_name
coverage_dat <- coverage_dat[,-1]

lev_file <- file.path(opt$editlevel)
if (endsWith(lev_file, "gz") | endsWith(lev_file, "zip") | endsWith(lev_file, "bzip2") | endsWith(lev_file, "xz")) {
  editLevel_dat <- suppressWarnings(temp_unzip(lev_file, fread, data.table = F))
} else {
  editLevel_dat <- suppressWarnings(fread(lev_file, data.table = F))
}
temp <- editLevel_dat$site_name
editLevel_dat <- editLevel_dat[,-1]
editLevel_dat[editLevel_dat == "N/A"] <- NA
editLevel_dat <- sapply(editLevel_dat, as.numeric)
rownames(editLevel_dat) <- temp

tpm_file <- file.path(opt$tpmfile)
tpm_dat <- temp_unzip(tpm_file, fread, data.table = F)
rownames(tpm_dat) <- tpm_dat$Name
tpm_dat <- tpm_dat[,-c(1:2)]
tpm_dat <- tpm_dat[,which(colnames(tpm_dat) %in% sample.info$SAMPID)]


####
## DEFINE FUNCTIONS
get_id_match <- function(query_id, query_df = rnaseqc.genes, subject_df = gene_anno) {
  enchar <- rnaseqc.genes[query_id, ]$ensgene_char
  symb <- query_df[query_id, ]$symbol
  if (enchar %in% subject_df$ensgene_char) {
    return(rownames(subject_df[subject_df$ensgene_char == enchar, ])[1])
  }
  else if (symb %in% subject_df$symbol) {
    return(rownames(subject_df[subject_df$symbol == symb, ])[1])
  }
  else {
    return(NA)
  }
}


####
## Process and save coverage and editlevel files for faster loading in the future
cluster_dat$site_name <- paste(cluster_dat$chr, cluster_dat$site_start, sep = "_")
cluster_dat$quant_site_name <- paste(cluster_dat$chr, cluster_dat$site_start-1, sep = "_")
summary(duplicated(cluster_dat$site_name)) # all unique
summary(colnames(coverage_dat) %in% sample.info$SAMPID) # all are
summary(colnames(editLevel_dat) %in% sample.info$SAMPID) # all are
all(colnames(editLevel_dat) %in% colnames(coverage_dat))
save(coverage_dat, editLevel_dat, file = "data/coverage_and_editLevel.RData")


####
## Create a df matching sample names to tissue types, save for future use
tissue_sample_match <- tibble::tibble(gtex_tissue = names(table(sample.info$SMTSD)[order(table(sample.info$SMTSD), decreasing = T)]),
                                      work_tissue = NA,
                                      n_samples = as.numeric(table(sample.info$SMTSD)[order(table(sample.info$SMTSD), decreasing = T)]),
                                      samples = rep("A", length(table(sample.info$SMTSD))))
tissue_sample_match$work_tissue <- gsub(" ", "_", gsub("\\)", "", gsub(" \\(", "_", gsub(" - ", "_", tissue_sample_match$gtex_tissue))))
for (i in tissue_sample_match$gtex_tissue) {
  x <- sample.info[sample.info$SMTSD == i, ]$SAMPID
  tissue_sample_match[tissue_sample_match$gtex_tissue == i, ]$samples <- paste(x, collapse = ",")
}
save(tissue_sample_match, file = "data/sample_tissue_match.RData")


####
## make gene name match df
# needed because putative cisnats were identified with gencode v35, but gtex is annotated with gencode v26
gene_id_match <- data.frame(gtex = rownames(tpm_dat), cisnat = NA)
gene_anno$ensgene_char <- gsub("\\..*", "", rownames(gene_anno))
gene_id_match$cisnat <- as.character(sapply(gene_id_match$gtex, get_id_match))
save(gene_id_match, file = "data/gtex_match_v35_gene_ids.RData")


####
## Make TPM summary files (median, and expression thresholds)
tpm_exp <- tpm_dat >= 0.1
all(colnames(tpm_exp) == sample.info$SAMPID) # TRUE

tissue_median <- as.data.frame(matrix(nrow = nrow(tpm_dat),
                                      ncol = nrow(tissue_sample_match)),
                               row.names = rownames(tpm_dat))
colnames(tissue_median) <- tissue_sample_match$work_tissue
tissue_01_exp <- as.data.frame(matrix(nrow = nrow(tpm_dat),
                                      ncol = nrow(tissue_sample_match)),
                               row.names = rownames(tpm_dat))
colnames(tissue_01_exp) <- tissue_sample_match$work_tissue
for (i in c(1:ncol(tissue_median))) {
  tissue <- colnames(tissue_median)[i]
  print(tissue)
  samples <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == colnames(tissue_median)[i], ]$samples, ","))
  # expression (T/F)
  exps <- rowSums(tpm_exp[, samples])
  exps <- exps >= length(samples)/5
  stopifnot(all(names(exps) == rownames(tissue_01_exp)))
  tissue_01_exp[,tissue] <- as.logical(exps)
  # medians
  mdns <- apply(tpm_dat[, samples], 1, median, na.rm = T)
  stopifnot(all(names(mdns) == rownames(tissue_median)))
  tissue_median[,tissue] <- as.numeric(mdns)
}
save(tissue_median, tissue_01_exp, file = "data/tissue_gExp_summary.RData")


####
## make edit_df
# each row is a edit_site * tissue combo
sites <- unique(unlist(strsplit(cisnat_withSites$sites, ";")))
total_rows <- length(sites) * nrow(tissue_sample_match)
df_cols <- c("site_id", "tissue", "n_samples", "min_coverage", "med_coverage", "max_coverage", "n_cov_tenPlus",
             "min_edit", "med_edit", "mean_edit", "max_edit", "n_edit_zero", "sd", "coefVar")
edit_df <- data.frame(matrix(NA, nrow = total_rows, ncol = length(df_cols)))
colnames(edit_df) <- df_cols
edit_df$tissue <- rep(tissue_sample_match$work_tissue, length(sites))
edit_df$n_samples <- rep(tissue_sample_match$n_samples, length(sites))
edit_df$site_id <- rep(sites, each = nrow(tissue_sample_match))
summary(duplicated(paste(edit_df$site_id, edit_df$tissue))) # no dups, which is good (each row is a unique site/tissue combo)
summary(sites %in% rownames(coverage_dat)) # site position is same between these dfs (minusOne), no need to adjust for matching

for (tis in unique(edit_df$tissue)) {
  print(tis)
  print(date())
  samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tis, ]$samples, ","))
  # coverage
  coverage_temp <- coverage_dat[, samps]
  cmin <- apply(coverage_temp, 1, min, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$min_coverage <- cmin[sites]
  cmed <- apply(coverage_temp, 1, median, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$med_coverage <- cmed[sites]
  cmax <- apply(coverage_temp, 1, max, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$max_coverage <- cmax[sites]
  ccovten <- rowSums(coverage_temp >= 10)
  edit_df[edit_df$tissue == tis, ]$n_cov_tenPlus <- ccovten[sites]
  # editing
  editLevel_temp <- editLevel_dat[ ,samps]
  editLevel_temp_filt <- editLevel_temp
  editLevel_temp_filt[coverage_temp<10] <- NA # block out edit level values where the coverage is too low
  emin <- apply(editLevel_temp_filt, 1, min, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$min_edit <- emin[sites]
  emed <- apply(editLevel_temp_filt, 1, median, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$med_edit <- emed[sites]
  emean <- apply(editLevel_temp_filt, 1, mean, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$mean_edit <- emean[sites]
  emax <- apply(editLevel_temp_filt, 1, max, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$max_edit <- emax[sites]
  eZero <- rowSums(editLevel_temp_filt == 0, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$n_edit_zero <- eZero[sites]
  esd <- apply(editLevel_temp_filt, 1, sd, na.rm = T)
  edit_df[edit_df$tissue == tis, ]$sd <- esd[sites]
}

edit_df$coefVar <- edit_df$sd/edit_df$mean_edit
edit_df[is.infinite(edit_df$min_edit), ]$min_edit <- NA
edit_df[is.infinite(edit_df$max_edit), ]$max_edit <- NA
summary(edit_df) # 35k consistent NAs across the edits; they're ones where all sample coverage is <10
# ^ (more NAs in the sd and coefVar columns, I think because of 0/0 division)
save(edit_df, file = file.path(out_dir, "06_edit_df.RData"))



####
## make exp_df
# each row is a gene * tissue combo
genes <- rownames(tpm_dat)
total_rows <- length(genes) * nrow(tissue_sample_match)
df_cols <- c("gene_id", "tissue", "n_samples", "n_samples_tpm_detected", "min_tpm", "med_tpm", 
             "mean_tpm", "max_tpm", "sd_all", "sd_detected", "coefVar_all", "coefVar_detected")
exp_df <- data.frame(matrix(NA, nrow = total_rows, ncol = length(df_cols)))
colnames(exp_df) <- df_cols
exp_df$tissue <- rep(tissue_sample_match$work_tissue, length(genes))
exp_df$n_samples <- rep(tissue_sample_match$n_samples, length(genes))
exp_df$gene_id <- rep(genes, each = nrow(tissue_sample_match))
summary(duplicated(paste(exp_df$gene_id, exp_df$tissue))) # no dups, which is good (each row is a unique site/tissue combo)

for (tis in unique(exp_df$tissue)) {
  print(tis)
  print(date())
  samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tis, ]$samples, ","))
  tpm_temp <- tpm_dat[, samps]
  tmin <- apply(tpm_temp, 1, min, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$min_tpm <- tmin[genes]
  tmed <- apply(tpm_temp, 1, median, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$med_tpm <- tmed[genes]
  tmean <- apply(tpm_temp, 1, mean, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$mean_tpm <- tmean[genes]
  tmax <- apply(tpm_temp, 1, max, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$max_tpm <- tmax[genes]
  tdetect <- rowSums(tpm_temp > 0)
  exp_df[exp_df$tissue == tis, ]$n_samples_tpm_detected <- tdetect[genes]
  tsd <- apply(tpm_temp, 1, sd, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$sd_all <- tsd[genes]
  tpm_temp_filt <- tpm_temp
  tpm_temp_filt[tpm_temp == 0] <- NA # block out tpm values of zero
  tsddetect <- apply(tpm_temp_filt, 1, sd, na.rm = T)
  exp_df[exp_df$tissue == tis, ]$sd_detected <- tsddetect[genes]
}

exp_df$coefVar_all <- exp_df$sd_all/exp_df$mean_tpm
exp_df$coefVar_detected <- exp_df$sd_detected/exp_df$mean_tpm
rownames(gene_id_match) <- gene_id_match$gtex
exp_df$cisnat_gene <- gene_id_match[exp_df$gene_id, ]$cisnat
summary(exp_df) # for the NAs in coefVar and sd fields: 7422 rows have no samples with TPM > 0 (NAs in sd_detected, coefVar_all, and coefVar_detected), and 3009 rows have one sample with TPM > 0 (NAs in sd_detected and coefVar_detected)
save(exp_df, file = file.path(out_dir, "06_exp_df.RData"))


