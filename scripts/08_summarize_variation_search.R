#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 08: summarizing interesting variation in editing across all tissues
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/08_summarize_variation_search.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz

####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
options(datatable.fread.datatable = F)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape2))

option_list <- list( 
  make_option("--tissue", default = NULL, 
              help = "name of tissue"),
  make_option("--samplefile", default = NULL,
              help = "File with GTEx sample annotation [default \"%default\"]"),
  make_option("--tpmfile", default = NULL, 
              help = "File with gene expression (TPM) values for genes of interest"),
  make_option("--outdir", default = "./output", 
              help = "Output directory [default is \"%default\"]"),
  make_option("--figdir", default = "./figures", 
              help = "Figure directory [default is \"%default\"]"),
  make_option("--scriptdir", default = "./scripts", 
              help = "Script directory (to load functions saved in source scripts). [default is \"%default\"]"),
  make_option("--scratchdir", default = NULL, 
              help = "Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default \"%default\"]")
)
opt <- parse_args(OptionParser(option_list=option_list))


####
## ASSIGN ARGUMENTS


####
## INPUT TESTS
if (is.null(opt$tpmfile)) { 
  stop("TPM file not provided, exiting\n") 
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
load(file.path(opt$samplefile))
sample.info <- sample.info[sample.info$ANALYTE_TYPE == "RNA:Total RNA", ]

tpm_file <- file.path(opt$tpmfile)
tpm_dat <- temp_unzip(tpm_file, fread, data.table = F)
rownames(tpm_dat) <- tpm_dat$Name
tpm_dat <- tpm_dat[,-c(1:2)]
tpm_dat <- tpm_dat[,which(colnames(tpm_dat) %in% sample.info$SAMPID)]

load(file.path(out_dir, "gene_tx_exon_anno_files.RData")) # object name from script 01: gene_anno, tx_anno, and exon_anno
load(file.path("data/gtex_edit/All.AG.stranded.annovar.RData")) # RData version of editing site annotation, saved in script 05: cluster_dat
load(file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData")) # object name from script 03: putative_cisnat
load(file.path(out_dir, "05_cisNAT_with_RNAedit.RData")) # two objects from script 05, named: cisnat_withSites, cisnat_withCluster
load("data/gtex_match_v35_gene_ids.RData") # R object from script 06, named: gene_id_match
load("data/sample_tissue_match.RData") # R object from script 06, named: tissue_sample_match
load("data/tissue_gExp_summary.RData") # Two R objects from script 06, named: tissue_median, tissue_01_exp
load("data/coverage_and_editLevel.RData") # Two R objects from script 06, named: coverage_dat, editLevel_dat
load(file.path(out_dir, "06_exp_df.RData")) # R object from script 06, named: exp_df
load(file.path(out_dir, "06_edit_df.RData")) # R object from script 06, named: edit_df


####
## DEFINE FUNCTIONS
dsrna_boxplot <- function(dsrna, tissue) {
  # all tissues plot: including all sites in dsrna, regardless of if variable or not (too hard to pull data conditionally on tissue)
  if (tissue == "all") {
    sites <- unlist(strsplit(cisnat_withSites[dsrna, ]$sites, ";"))
    if (length(sites) > 50) {return("Over 50 sites - exiting plot")}
    plot_dat <- data.frame()
    for (tis in tissue_sample_match$work_tissue) {
      samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tis, ]$samples, ","))
      cov_tis <- coverage_dat[sites, samps]
      edit_tis <- editLevel_dat[sites, samps]
      if (length(unique(sites)) > 50) {return("Over 50 sites - exiting plot")}
      c1 <- t(cov_tis[sites, , drop = F])
      e1 <- t(edit_tis[sites, , drop = F])
      e1_covFilt <- e1
      if (any(c1 < 10)) {
        e1_covFilt[c1 < 10] <- NA
      }
      plot_dat_tis <- reshape2::melt(e1_covFilt)
      colnames(plot_dat_tis) <- c("sample", "site", "edit_level")
      plot_dat_tis$var_site <- "n/a"
      plot_dat <- rbind(plot_dat, plot_dat_tis)
    }
  }
  else {
    var_dsrna_tis <- var_dsrna_allTissues[var_dsrna_allTissues$tissue == tissue, ]
    sites <- unlist(strsplit(cisnat_withSites[dsrna, ]$sites, ";"))
    var_sites <- unlist(strsplit(var_dsrna_tis[var_dsrna_tis$dsrna_id == dsrna, ]$var_sites, ";"))
    if (length(sites) > 50) {return("Over 50 sites - exiting plot")}
    samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tissue, ]$samples, ","))
    cov_tis <- coverage_dat[sites, samps]
    edit_tis <- editLevel_dat[sites, samps]
    c1 <- t(cov_tis[sites, , drop = F])
    e1 <- t(edit_tis[sites, , drop = F])
    e1_covFilt <- e1
    if (any(c1 < 10)) {
      e1_covFilt[c1 < 10] <- NA
    }
    plot_dat <- reshape2::melt(e1_covFilt)
    colnames(plot_dat) <- c("sample", "site", "edit_level")
    plot_dat$var_site <- "not_var_site"
    plot_dat[plot_dat$site %in% var_sites, ]$var_site <- "a_var_site"
  }
  plot_dat$site <- as.character(plot_dat$site)
  # each site is a factor (not to scale)
  plot_dat$site <- factor(plot_dat$site, levels = mixedsort(unique(plot_dat$site)))
  # plotting with flexible edit axis
  g <- ggplot(plot_dat, aes(x = site, y = edit_level, color = var_site)) + geom_boxplot() + xlab(NULL) + theme_bw() +
    ggtitle(paste(tissue, dsrna, cisnat_withSites[dsrna, ]$minus_symbol, cisnat_withSites[dsrna, ]$plus_symbol))
  return(g)
}

tpm_boxplot <- function(dsrna) {
  g1_cn <- cisnat_withSites[dsrna, ]$minus_gene
  g2_cn <- cisnat_withSites[dsrna, ]$plus_gene
  g1 <- gene_id_match[gene_id_match$cisnat == g1_cn, ]$gtex
  g2 <- gene_id_match[gene_id_match$cisnat == g2_cn, ]$gtex
  if (is.na(g1) | is.na(g2)) { return("gene annotation not present in GTEx") }
  tdat <- as.data.frame(t(tpm_dat[c(g1,g2), ]))
  tdat$tissue <- NA
  for (tis in tissue_sample_match$work_tissue) {
    samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tis, ]$samples, ","))
    tdat[samps, ]$tissue <- tis
  }
  p1 <- ggplot(tdat, aes(x = tissue, y = tdat[,g1])) + geom_boxplot() + xlab(NULL) + theme_bw() + ggtitle(paste(g1, gene_anno[g1_cn, ]$symbol))
  p2 <- ggplot(tdat, aes(x = tissue, y = tdat[,g2])) + geom_boxplot() + xlab(NULL) + theme_bw() + ggtitle(paste(g2, gene_anno[g2_cn, ]$symbol))
  return((p1 / p2))
}


####
## Read in and combine each tissue's results from script 7
var_dsrna_allTissues <- data.frame()
var_sites_allTissues <- data.frame()

for (tis in tissue_sample_match$work_tissue) {
  print(tis)
  # dsrna file
  read_file <- file.path(out_dir, paste0("07_perTissue_variation/", tis, "_dsrnasWithVarSites.txt"))
  temp <- read.table(file = read_file, header = T)
  temp$tissue <- tis
  temp <- temp[,c("tissue", "chr", "minus_gene", "minus_symbol", "minus_gene_type", "plus_gene", "plus_symbol", "plus_gene_type", "minus_tx", "minus_tx_type", "plus_tx", "plus_tx_type", "minus_exons_in_overlap", "plus_exons_in_overlap", "all_overlap_widths", "longest_overlap_width", "longest_overlap_start", "longest_overlap_end", "multi_exon_overlap", "alu_overlap", "n_sites", "n_var_sites", "var_types", "site_flags", "var_sites")]
  var_dsrna_allTissues <- rbind(var_dsrna_allTissues, temp)
  # site file
  read_file <- file.path(out_dir, paste0("07_perTissue_variation/", tis, "_sitesWithVariation.txt"))
  temp <- read.table(file = read_file, header = T)
  var_sites_allTissues <- rbind(var_sites_allTissues, temp)
}
save(var_dsrna_allTissues, var_sites_allTissues, file = file.path(out_dir, "08_allTissues_variation_summaryTables.RData"))


####
## How many dsRNAs with variable sites are there per tissue?
# total variable dsRNAs:
nrow(var_dsrna_allTissues) # 33,449
length(unique(paste(var_dsrna_allTissues$minus_tx, var_dsrna_allTissues$plus_tx, sep = "__"))) # 858 tx pairs
length(unique(paste(var_dsrna_allTissues$minus_gene, var_dsrna_allTissues$plus_gene, sep = "__"))) # 323 pairs
table(var_dsrna_allTissues$tissue)[order(table(var_dsrna_allTissues$tissue))] 
# ^ ranges from Kidney_Cortex at 564 to Testis at 767
nrow(temp <- var_dsrna_allTissues[-which(duplicated(var_dsrna_allTissues$dsrna_id)), ])
summary(is.na(temp$alu_overlap)) # 732 do not overlap IRAlu, 126 do

# variable dsRNAs with >1 site:
nrow(var_filt <- var_dsrna_allTissues[var_dsrna_allTissues$n_var_sites > 1, ]) # 26,509
length(unique(paste(var_filt$minus_tx, var_filt$plus_tx, sep = "__"))) # 651 tx pairs
length(unique(paste(var_filt$minus_gene, var_filt$plus_gene, sep = "__"))) # 223 pairs
table(var_filt$tissue)[order(table(var_filt$tissue))] 
# ^ ranges from Whole_Blood at 440 to Thyroid at 602 (Testis dropped down in the rankings quite a bit)
nrow(temp <- var_filt[-which(duplicated(var_filt$dsrna_id)), ])
summary(is.na(temp$alu_overlap)) # 533 do not overlap IRAlu, 118 do

# variable dsRNAs with >=5 sites:
nrow(var_filt <- var_dsrna_allTissues[var_dsrna_allTissues$n_var_sites > 5, ]) # 20,233
length(unique(paste(var_filt$minus_tx, var_filt$plus_tx, sep = "__"))) # 492 tx pairs
length(unique(paste(var_filt$minus_gene, var_filt$plus_gene, sep = "__"))) # 166 pairs
table(var_filt$tissue)[order(table(var_filt$tissue))] 
# ^ ranges from Kidney_Cortex at 335 to Thyroid at 460
nrow(temp <- var_filt[-which(duplicated(var_filt$dsrna_id)), ])
summary(is.na(temp$alu_overlap)) # 394 do not overlap IRAlu, 98 do



####
## Plot: number of variable dsRNAs per tissue, number of variable sites per tissue
# dsrnas (fill is number of VARIABLE sites within the dsrna)
plot_dat <- var_dsrna_allTissues
plot_dat$n_var_sites_cat <- "1"
plot_dat[plot_dat$n_var_sites > 1 & plot_dat$n_var_sites <= 5, ]$n_var_sites_cat <- "2-5"
plot_dat[plot_dat$n_var_sites > 5 & plot_dat$n_var_sites <= 10, ]$n_var_sites_cat <- "6-10"
plot_dat[plot_dat$n_var_sites > 10 & plot_dat$n_var_sites <= 25, ]$n_var_sites_cat <- "11-25"
plot_dat[plot_dat$n_var_sites > 25 & plot_dat$n_var_sites <= 50, ]$n_var_sites_cat <- "26-50"
plot_dat[plot_dat$n_var_sites > 50 & plot_dat$n_var_sites <= 100, ]$n_var_sites_cat <- "51-100"
plot_dat[plot_dat$n_var_sites > 100, ]$n_var_sites_cat <- "101+"
plot_dat$n_var_sites_cat <- factor(plot_dat$n_var_sites_cat, levels = rev(c("1", "2-5", "6-10", "11-25", "26-50", "51-100", "101+")))
g <- ggplot(plot_dat, aes(x = tissue, fill = n_var_sites_cat)) + geom_bar() + 
  theme_bw() + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(g, file = file.path(fig_dir, "08_01_number_at_overlapLengthThresholds_withSites.png"), width = 6, height = 5)
ggsave(g, file = file.path(fig_dir, "08_01_number_at_overlapLengthThresholds_withSites.pdf"), width = 6, height = 5)





####
## Plot: dsRNAs by number of tissues observed in
var_dsrna_allTissues$dsrna_id <- paste(var_dsrna_allTissues$minus_tx, var_dsrna_allTissues$plus_tx, sep = "__")
var_dsrna_allTissues$gene_pair_id <- paste(var_dsrna_allTissues$minus_gene, var_dsrna_allTissues$plus_gene, sep = "__")
# sanity check: no tx pair is represented >49 (the number of tissues)
summary(as.numeric(table(var_dsrna_allTissues$dsrna_id))) # good, max is 49
table(as.numeric(table(var_dsrna_allTissues$dsrna_id))) # of the 858 tx pairs, 427 have a variable site in all tissues; 619 have a variable site in at least 40 tissues

# (unique tx pairs, all)
plot_dat <- data.frame(t(table(var_dsrna_allTissues$dsrna_id)))
colnames(plot_dat) <- c("misc", "dsrna", "num_tissues")
plot_dat$dsrna <- as.character(plot_dat$dsrna)
g <- ggplot(plot_dat, aes(x = num_tissues)) + geom_bar() + 
  scale_x_continuous(breaks = seq(0,50,5), minor_breaks = c(0:49)) + 
  ggtitle("Variable dsRNAs, all, tx_pairs") + theme_bw()
ggsave(g, file = file.path(fig_dir, "08_03_num_tis_vardsrna_all_tx.png"), width = 7, height = 7, units = "in")
ggsave(g, file = file.path(fig_dir, "08_03_num_tis_vardsrna_all_tx.pdf"), width = 7, height = 7, units = "in")

# (unique gene pairs, all)
temp <- var_dsrna_allTissues
temp$gp_tis <- paste(temp$gene_pair_id, temp$tissue, sep = "___")
temp <- temp[-which(duplicated(temp$gp_tis)), ]
plot_dat <- data.frame(t(table(temp$gene_pair_id)))
colnames(plot_dat) <- c("misc", "gene_pair", "num_tissues")
plot_dat$gene_pair <- as.character(plot_dat$gene_pair)
g <- ggplot(plot_dat, aes(x = num_tissues)) + geom_bar() + 
  scale_x_continuous(breaks = seq(0,50,5), minor_breaks = c(0:49)) + 
  ggtitle("Variable dsRNAs, all, gene_pairs") + theme_bw()
ggsave(g, file = file.path(fig_dir, "08_04_num_tis_vardsrna_all_gene.png"), width = 7, height = 7, units = "in")
ggsave(g, file = file.path(fig_dir, "08_04_num_tis_vardsrna_all_gene.pdf"), width = 7, height = 7, units = "in")

# (expression distribution of genes involved in overlaps with RNA editing sites)
# (expression threshold: detected (TPM>0) in at least 25% of samples)
genes <- unique(c(var_dsrna_allTissues$minus_gene, var_dsrna_allTissues$plus_gene))
gtex_genes <- gene_id_match[gene_id_match$cisnat %in% genes, ]$gtex
temp <- exp_df
temp <- temp[temp$gene_id %in% gtex_genes, ]
temp <- temp[temp$med_tpm >= 0.1, ]
plot_dat <- data.frame(gene = names(table(temp$gene_id)),
                       num_tissues = as.numeric(table(temp$gene_id)))
g <- ggplot(plot_dat, aes(x = num_tissues)) + geom_bar() + 
  scale_x_continuous(breaks = seq(0,50,5), minor_breaks = c(0:49)) + 
  ggtitle("Number of tissues expressing genes in variable dsRNAs, all (exp = median TPM >= 0.1)") + theme_bw()
ggsave(g, file = file.path(fig_dir, "08_03_num_tis_expressingGenesIn_vardsrna_all.png"), width = 7, height = 7, units = "in")
ggsave(g, file = file.path(fig_dir, "08_03_num_tis_expressingGenesIn_vardsrna_all.pdf"), width = 7, height = 7, units = "in")



# (unique tx pairs, >1 variable sites)
nrow(var_filt <- var_dsrna_allTissues[var_dsrna_allTissues$n_var_sites > 1, ]) # 26,509
plot_dat <- data.frame(t(table(var_filt$dsrna_id)))
colnames(plot_dat) <- c("misc", "dsrna", "num_tissues")
plot_dat$dsrna <- as.character(plot_dat$dsrna)
g <- ggplot(plot_dat, aes(x = num_tissues)) + geom_bar() + 
  scale_x_continuous(breaks = seq(0,50,5), minor_breaks = c(0:49)) + 
  ggtitle("Variable dsRNAs, more than 1 var site, tx_pairs") + theme_bw()
ggsave(g, file = file.path(fig_dir, "08_05_num_tis_vardsrna_moreThanOneSite_tx.png"), width = 7, height = 7, units = "in")
ggsave(g, file = file.path(fig_dir, "08_05_num_tis_vardsrna_moreThanOneSite_tx.pdf"), width = 7, height = 7, units = "in")

# (unique gene pairs, >1 variable sites)
temp <- var_filt
temp$gp_tis <- paste(temp$gene_pair_id, temp$tissue, sep = "___")
temp <- temp[-which(duplicated(temp$gp_tis)), ]
plot_dat <- data.frame(t(table(temp$gene_pair_id)))
colnames(plot_dat) <- c("misc", "gene_pair", "num_tissues")
plot_dat$gene_pair <- as.character(plot_dat$gene_pair)
g <- ggplot(plot_dat, aes(x = num_tissues)) + geom_bar() + 
  scale_x_continuous(breaks = seq(0,50,5), minor_breaks = c(0:49)) + 
  ggtitle("Variable dsRNAs, more than 1 var site, gene_pairs") + theme_bw()
ggsave(g, file = file.path(fig_dir, "08_06_num_tis_vardsrna_moreThanOneSite_gene.png"), width = 7, height = 7, units = "in")
ggsave(g, file = file.path(fig_dir, "08_06_num_tis_vardsrna_moreThanOneSite_gene.pdf"), width = 7, height = 7, units = "in")



####
## Select example dsRNAs to plot editing in all tissues
# panels: editing at sites in one example tissue; editing at sites across all involved tissues
# tpm across all tissues, gene1; tpm, gene2; 
# wanted to also do edit x tpm across all tissues for each gene, but realized this would be confusing as there are multiple sites / gene

# choose a few dsrnas
nrow(var_filt <- var_dsrna_allTissues[var_dsrna_allTissues$n_var_sites > 1, ]) # 26,509
filt_dat <- data.frame(t(table(var_filt$dsrna_id)))
colnames(filt_dat) <- c("misc", "dsrna", "num_tissues")
filt_dat$dsrna <- as.character(filt_dat$dsrna)
common_dat <- filt_dat[filt_dat$num_tissues == 49, ]
common_dat$num_sites <- unlist(lapply(strsplit(cisnat_withSites[common_dat$dsrna, ]$sites, ";"), length))
common_dat$alu_overlap <- cisnat_withSites[common_dat$dsrna, ]$alu_overlap
table(is.na(common_dat$alu_overlap)) # 285 are not overlapping IRAlus, 64 are
common_dat <- common_dat[is.na(common_dat$alu_overlap), ]

## dsrna: ENST00000229390.8__ENST00000551765.6
ds <- "ENST00000229390.8__ENST00000551765.6"
p1 <- dsrna_boxplot(ds, tissue = "all")
p2 <- dsrna_boxplot(ds, tissue = "Adipose_Subcutaneous")
p3 <- dsrna_boxplot(ds, tissue = "Muscle_Skeletal")
p4 <- dsrna_boxplot(ds, tissue = "Brain_Cerebellum")
p5 <- tpm_boxplot(ds)
gp <- (p1 / p2 / p3 / p4 / p5 + plot_layout(heights = c(1,1,1,1,2)))
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000229390_ENST00000551765.png"), width = 7, height = 10, units = "in")
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000229390_ENST00000551765.pdf"), width = 7, height = 10, units = "in")


## dsrna: ENST00000231948.8__ENST00000434583.5
ds <- "ENST00000231948.8__ENST00000434583.5"
p1 <- dsrna_boxplot(ds, tissue = "all")
p2 <- dsrna_boxplot(ds, tissue = "Adipose_Subcutaneous")
p3 <- dsrna_boxplot(ds, tissue = "Muscle_Skeletal")
p4 <- dsrna_boxplot(ds, tissue = "Brain_Cerebellum")
p5 <- tpm_boxplot(ds)
gp <- (p1 / p2 / p3 / p4 / p5 + plot_layout(heights = c(1,1,1,1,2)))
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000231948_ENST00000434583.png"), width = 7, height = 10, units = "in")
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000231948_ENST00000434583.pdf"), width = 7, height = 10, units = "in")


## dsrna: ENST00000357235.6__ENST00000368926.6
ds <- "ENST00000357235.6__ENST00000368926.6"
p1 <- dsrna_boxplot(ds, tissue = "all")
p2 <- dsrna_boxplot(ds, tissue = "Adipose_Subcutaneous")
p3 <- dsrna_boxplot(ds, tissue = "Muscle_Skeletal")
p4 <- dsrna_boxplot(ds, tissue = "Brain_Cerebellum")
p5 <- tpm_boxplot(ds)
gp <- (p1 / p2 / p3 / p4 / p5 + plot_layout(heights = c(1,1,1,1,2)))
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000357235_ENST00000368926.png"), width = 7, height = 10, units = "in")
ggsave(gp, file = file.path(fig_dir, "08_zexample_dsrna_all49tissues_noIRAluOverlap_ENST00000357235_ENST00000368926.pdf"), width = 7, height = 10, units = "in")


## dsrna that is a lnc combo
pc_lnc <- common_dat
pc_lnc$minus_tx_type <- putative_cisnat[pc_lnc$dsrna, ]$minus_tx_type
pc_lnc$plus_tx_type <- putative_cisnat[pc_lnc$dsrna, ]$plus_tx_type
pc_lnc <- pc_lnc[pc_lnc$minus_tx_type == "lncRNA" | pc_lnc$plus_tx_type == "lncRNA",]
for (i in c(1:nrow(pc_lnc))) {
  ds <- pc_lnc$dsrna[i]
  p1 <- dsrna_boxplot(ds, tissue = "all")
  p2 <- dsrna_boxplot(ds, tissue = "Adipose_Subcutaneous")
  p3 <- dsrna_boxplot(ds, tissue = "Muscle_Skeletal")
  p4 <- dsrna_boxplot(ds, tissue = "Brain_Cerebellum")
  p5 <- tpm_boxplot(ds)
  gp <- (p1 / p2 / p3 / p4 / p5 + plot_layout(heights = c(1,1,1,1,2)))
  ggsave(gp, file = file.path(fig_dir, paste0("08_yexample_dsrna_all49tissues_noIRAluOverlap_", ds, ".png")), width = 7, height = 10, units = "in")
  ggsave(gp, file = file.path(fig_dir, paste0("08_yexample_dsrna_all49tissues_noIRAluOverlap_", ds, ".pdf")), width = 7, height = 10, units = "in")
}
