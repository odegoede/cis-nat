#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script pilot_01: filtering and combining the pilot full-genome sprint output
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/pilot_01_filter_SPRINT_output.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
library(gtools)


####
## ASSIGN ARGUMENTS
option_list <- list( 
  make_option("--samplefile", default = NULL,
              help = "File with GTEx sample annotation [default \"%default\"]"),
  make_option("--outdir", default = "./output", 
              help = "Output directory [default is current working directory]"),
  make_option("--scriptdir", default = "./scripts",
              help = "Scripts directory (to load functions saved in source scripts). [default is \"%default\"]"),
  make_option("--scratchdir", default = NULL, 
              help = "Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default \"%default\"]"),
  make_option("--figdir", default = "./figures", 
              help = "Figure directory [default is \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))


####
## INPUT TESTS
if (is.null(opt$samplefile)) { 
  stop("Sample info file not provided, exiting\n") 
}

# Create variables for the directories that are consistent (doesn't have trailing "/")
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
## READ IN INPUT FILES
load(file.path(out_dir, "gene_tx_exon_anno_files.RData")) # object name from script 01: gene_anno, tx_anno, and exon_anno
load(file.path("data/gtex_edit/All.AG.stranded.annovar.RData")) # RData version of editing site annotation, saved in script 05: cluster_dat
load("data/gtex_match_v35_gene_ids.RData") # R object from script 06, named: gene_id_match
load("data/sample_tissue_match.RData") # R object from script 06, named: tissue_sample_match
load(file.path(opt$samplefile))
sample.info <- sample.info[sample.info$ANALYTE_TYPE == "RNA:Total RNA", ]


####
## DEFINE FUNCTIONS
get_total_sites <- function(tissue_quant_file) {
  temp <- read.table(file.path(paste0("data/full_pilot/tissue_quant/", tissue_quant_file)), 
                     header = T, sep = "\t")
  return(unique(temp$site_name))
}

get_sample_names <- function(tissue_quant_file) {
  temp <- read.table(file.path(paste0("data/full_pilot/tissue_quant/", tissue_quant_file)), 
                     header = T, sep = "\t")
  return(unique(temp$sample))
}

enough_sites_edit <- function(in_table, tis, n_sample) {
  temp <- in_table[, grep(tis, colnames(in_table))]
  temp <- temp > 0
  temp_filt <- temp[as.logical(rowSums(temp) >= n_sample), ]
  return(rownames(temp_filt))
}



####
## Get all editing sites across all files
total_sites <- c()
sample_names <- c()

for (tis in tissue_sample_match$work_tissue) {
  print(tis)
  file_name <- grep(tis, list.files("data/full_pilot/tissue_quant/"), value = T)
  total_sites <- union(total_sites, get_total_sites(file_name))
  sample_names <- union(sample_names, get_sample_names(file_name))
}

length(total_sites) # 2,330,470
length(sample_names) # 289



####
## Plot number of samples / tissue
plot_dat <- data.frame(samps = paste0("GTEX", sub(".*GTEX", "", sample_names)),
                       tissue = sub("_GTEX.*", "", sample_names))
g <- ggplot(plot_dat, aes(x = tissue)) + geom_bar() + theme_bw()
# ^ not a very cool plot; nearly all tissues have 6 samples
table(plot_dat$tissue)[table(plot_dat$tissue) != 6]
# Kidney is 4, salivary gland is 4, small intestine is 5

####
## Load and combine each tissue file
total_sites <- mixedsort(total_sites)
pilot_cov_dat <- data.frame(matrix(NA, nrow = length(total_sites), ncol = length(sample_names)),
                            row.names = total_sites)
colnames(pilot_cov_dat) <- sample_names
pilot_edit_dat <- data.frame(matrix(NA, nrow = length(total_sites), ncol = length(sample_names)),
                             row.names = total_sites)
colnames(pilot_edit_dat) <- sample_names

for (tis in tissue_sample_match$work_tissue) {
  print(tis)
  file_name <- grep(tis, list.files("data/full_pilot/tissue_quant/"), value = T)
  temp <- read.table(file.path(paste0("data/full_pilot/tissue_quant/", file_name)), 
                     header = T, sep = "\t")
  print(summary(temp$editLevel == 0))
  for (samp in unique(temp$sample)) {
    temp_samp <- temp[temp$sample == samp, ]
    rownames(temp_samp) <- temp_samp$site_name
    pilot_edit_dat[rownames(temp_samp), samp] <- temp_samp$editLevel
    pilot_cov_dat[rownames(temp_samp), samp] <- temp_samp$coverage
  }
}
dim(pilot_edit_dat)
dim(pilot_cov_dat)
pilot_edit_dat[1:6,1:6]
if (!is.null(opt$scratchdir)) {
  save(pilot_cov_dat, pilot_edit_dat, file = file.path(scratch_dir, "pilot_data_initial_combine.rds"))
}


####
## Filter sites
## remove sites that are only present in one sample
rowSums(!is.na(head(pilot_edit_dat)))
summary(rowSums(!is.na(pilot_edit_dat)) == 1) # 1,100,817 are NAs for all but 1 sample; 1,229,653 have real values for more than 1 sample
summary(rowSums(!is.na(pilot_cov_dat)) == 1) # should be the same as line above, and it is: 1,100,817 are NAs for all but 1 sample; 1,229,653 have real values for more than 1 sample
dim(pilot_edit_filt <- pilot_edit_dat[as.logical(rowSums(!is.na(pilot_edit_dat)) > 1), ]) # 1,229,653 rows, 578 cols
dim(pilot_cov_filt <- pilot_cov_dat[as.logical(rowSums(!is.na(pilot_cov_dat)) > 1), ]) # 1,229,653, 578 cols

## remove sites where there isn't a single tissue with evidence of editing (>0 edited reads in at least 1 sample)
summary(rowSums(pilot_edit_filt, na.rm = T) == 0) # all have evidence in at least 1 sample
# maybe require that at least 2 samples of one tissue have editing values?
# (pretty strict, since some tissues only have 4 samples; keep both the filt and the strict_filt datasets)
strict_filt_sites <- c()
for (tis in unique(unlist(lapply(strsplit(colnames(pilot_edit_filt), "_"), "[[", 1)))) {
  print(tis)
  strict_filt_sites <- union(strict_filt_sites, enough_sites_edit(in_table = pilot_edit_filt, tis = tis, n_sample = 2))
}
length(strict_filt_sites) # 84,939
head(strict_filt_sites)
strict_filt_sites <- strict_filt_sites[!is.na(strict_filt_sites)] # first one is NA - 84,938 sites
pilot_edit_strict_filt <- pilot_edit_filt[strict_filt_sites, ]

save(pilot_edit_strict_filt, pilot_edit_filt, pilot_cov_filt, file = "data/full_pilot/filtered_pilot_data.rds")

