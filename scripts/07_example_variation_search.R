#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script 07: testing ways to find interesting variation in editing within one tissue
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/07_example_variation_search.R --tissue $tis --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz

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
tis <- opt$tissue
print(tis)
# note: for examples and comments in this code, tis <- "Adipose_Subcutaneous"


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

# check_site: workhorse testing function for sites of interest:
# when check == "low_exp": assesses how a dsRNA site's editing level is matched by its genes' expression levels (e.g. when editing is zero, is it because expression is zero)
# when check == "cov_corr": tests correlation between editing of a site in a dsRNA and coverage at that site
# when check == "pos_corr": tests correlation between editing of a site in a dsRNA and expression levels of the two genes
# when check == "plot": produces several plots summarizing the editing/coverage/expression at site of interest
check_site <- function(site, check) {
  c1 <- t(cov_tis[site, , drop = F])
  e1 <- t(edit_tis[site, , drop = F])
  genes <- unique(c(cisnat_withSites[grep(site, cisnat_withSites$sites), ]$minus_gene, 
                    cisnat_withSites[grep(site, cisnat_withSites$sites), ]$plus_gene))
  gtex_genes <- unique(gene_id_match[gene_id_match$cisnat %in% genes, ]$gtex)
  to_keep <- rownames(c1[as.logical(c1[, site] >= 10), , drop = FALSE])
  c1 <- c1[to_keep, , drop = F]
  colnames(c1) <- "cov"
  e1 <- e1[to_keep, , drop = F]
  colnames(e1) <- "edit"
  t1 <- t(tpm_tis[gtex_genes, to_keep])
  stopifnot(all(rownames(t1) == rownames(c1)))
  stopifnot(all(rownames(t1) == rownames(e1)))
  plot_dat <- as.data.frame(cbind(c1, e1, t1))
  ds <- unique(paste(cisnat_withSites[grep(site, cisnat_withSites$sites), ]$minus_gene, cisnat_withSites[grep(site, cisnat_withSites$sites), ]$plus_gene, sep = "__"))
  if (check == "plot") {
    if (length(ds) > 1) {
      for (i in ds) {
        return(plot_per_cisnat(dsrna = i, plot_dat = plot_dat, site = site))
      }
    }
    else {
      return(plot_per_cisnat(dsrna = ds, plot_dat = plot_dat, site = site))
    }
  }
  if (check == "low_exp") {
    if (length(ds) > 1) {
      rval <- c()
      for (i in ds) {
        rval <- c(rval, check_low_exp(dsrna = i, plot_dat = plot_dat, site = site))
        rval <- factor(rval, levels = c("PASS", "MINI_FAIL", "NO_TEST", "BIG_FAIL"))
        return(as.character(rval[order(rval)][1]))
      }
    }
    else {
      return(check_low_exp(dsrna = ds, plot_dat = plot_dat, site = site))
    }
  }
  if (check == "pos_corr") {
    if (length(ds) > 1) {
      rval <- c()
      for (i in ds) {
        rval <- c(rval, check_pos_corr(dsrna = i, plot_dat = plot_dat, site = site))
        rval <- factor(rval, levels = c("POS", "NONE", "MINI_NEG", "BIG_NEG"))
        return(as.character(rval[order(rval)][1]))
      }
    }
    else {
      return(check_pos_corr(dsrna = ds, plot_dat = plot_dat, site = site))
    }
  }
  if (check == "cov_corr") {
    return(check_cov_corr(dsrna = ds, plot_dat = plot_dat, site = site))
  }
}

# plot_per_cisnat: for check = "plot":
plot_per_cisnat <- function(dsrna, plot_dat, site) {
  g1 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[1], ]$gtex[1]
  g2 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[2], ]$gtex[1]
  plot_dat$site <- site
  p1 <- ggplot(plot_dat, aes(x = site, y = edit)) + geom_boxplot() + coord_cartesian(ylim = c(0,1)) + xlab(NULL) + theme_bw()
  p2 <- ggplot(plot_dat, aes(x = cov, y = edit)) + geom_point(alpha = 0.7) + coord_cartesian(ylim = c(0,1)) + theme_bw()
  p5 <- ggplot(plot_dat, aes(x = site, y = edit)) + geom_boxplot() + xlab(NULL) + theme_bw()
  p6 <- ggplot(plot_dat, aes(x = cov, y = edit)) + geom_point(alpha = 0.7) + theme_bw()
  if (is.na(g1)) {
    p3 <- ggplot()
    p7 <- ggplot()
  }
  else {
    p3 <- ggplot(plot_dat, aes(x = plot_dat[,g1], y = edit)) + geom_point(alpha = 0.7) + coord_cartesian(ylim = c(0,1), xlim = c(0, max(plot_dat[,g1]))) + xlab(g1) + theme_bw()
    p7 <- ggplot(plot_dat, aes(x = plot_dat[,g1], y = edit)) + geom_point(alpha = 0.7) + coord_cartesian(xlim = c(0, max(plot_dat[,g1]))) + xlab(g1) + theme_bw()
  }
  if (is.na(g2)) {
    p4 <- ggplot()
    p8 <- ggplot()
  }
  else {
    p4 <- ggplot(plot_dat, aes(x = plot_dat[,g2], y = edit)) + geom_point(alpha = 0.7) + coord_cartesian(ylim = c(0,1), xlim = c(0, max(plot_dat[,g2]))) + xlab(g2) + theme_bw()
    p8 <- ggplot(plot_dat, aes(x = plot_dat[,g2], y = edit)) + geom_point(alpha = 0.7) + coord_cartesian(xlim = c(0, max(plot_dat[,g2]))) + xlab(g2) + theme_bw()
  }
  ((p1 / p2) | (p3 / p4) | (p5 / p6) | (p7 / p8))
}

# check_low_exp: for check = "low_exp":
check_low_exp <- function(dsrna, plot_dat, site) {
  g1 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[1], ]$gtex[1]
  g2 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[2], ]$gtex[1]
  if (is.na(g1) | is.na(g2)) {
    return("NO_TEST")
  }
  # make sure that at least one of the genes has SOME zeroes, but isn't ALL zeroes for tpm
  # note which (or both) gene(s) has some zeroes
  genes_of_interest <- c()
  if (any(plot_dat[,g1] == 0) & any(plot_dat[,g1] > 0)) {
    genes_of_interest <- g1
  }
  if (any(plot_dat[,g2] == 0) & any(plot_dat[,g2] > 0)) {
    genes_of_interest <- c(genes_of_interest, g2)
  }
  if (length(genes_of_interest) == 0) {
    return("NO_TEST")
  }
  # add columns for each gene to plot_dat, filled with "TBD" (NA throws error in the final status check, would have to incorporate a bunch of if statements)
  plot_dat$gene_one <- "TBD"
  plot_dat$gene_two <- "TBD"
  # only fill the column if the gene is one of the ones that has some zeroes
  # value is MATCH if presence of editing and presence of expression are synced up (either there's both editLevel and expLevel, or both are zeroes)
  # value is EXP_NO_EDIT if there's expression but no editing (not necessarily a deal breaker, expression could be nuclear or too low to form dsRNA)
  # value is EDIT_NO_EXP if there's editing when only one transcript is expressed (could be a deal breaker, suggests only one transcript needed for editing)
  if (g1 %in% genes_of_interest) {
    plot_dat[plot_dat[,"edit"] == 0 & plot_dat[,g1] == 0, "gene_one"] <- "MATCH_ZERO"
    plot_dat[plot_dat[,"edit"] > 0 & plot_dat[,g1] > 0, "gene_one"] <- "MATCH_EXP"
    plot_dat[plot_dat[,"edit"] == 0 & plot_dat[,g1] > 0, "gene_one"] <- "EXP_NO_EDIT"
    plot_dat[plot_dat[,"edit"] > 0 & plot_dat[,g1] == 0, "gene_one"] <- "EDIT_NO_EXP"
  }
  if (g2 %in% genes_of_interest) {
    plot_dat[plot_dat[,"edit"] == 0 & plot_dat[,g2] == 0, "gene_two"] <- "MATCH_ZERO"
    plot_dat[plot_dat[,"edit"] > 0 & plot_dat[,g2] > 0, "gene_two"] <- "MATCH_EXP"
    plot_dat[plot_dat[,"edit"] == 0 & plot_dat[,g2] > 0, "gene_two"] <- "EXP_NO_EDIT"
    plot_dat[plot_dat[,"edit"] > 0 & plot_dat[,g2] == 0, "gene_two"] <- "EDIT_NO_EXP"
  }
  # get a single useful return value that could go in a column
  # simple, permissive start:
  # if any edit_no_exp, the return is fail
  if (any(plot_dat[,"gene_one"] == "EDIT_NO_EXP") | any(plot_dat[,"gene_two"] == "EDIT_NO_EXP")) {
    # if < 5% of samples have edit with no exp, mini fail (won't filter out based on it)
    if (length(union(grep("EDIT_NO_EXP", plot_dat[,"gene_one"]), grep("EDIT_NO_EXP", plot_dat[,"gene_two"]))) < round(0.05*nrow(plot_dat))) {
      rval <- "MINI_FAIL"
    }
    # otherwise, a big fail and will be filtered out
    else {rval <- "BIG_FAIL"}
  }
  # if any match_zero, the return is pass
  else if (any(plot_dat[,"gene_one"] == "MATCH_ZERO") | any(plot_dat[,"gene_two"] == "MATCH_ZERO")) {
    rval <- "PASS"
  }
  # if somewhere in between, it's neither a flag or a mark of interesting variation
  else {
    rval <- "MODERATE"
  }
  return(rval)
}

# check_pos_corr: for check = "pos_corr":
check_pos_corr <- function(dsrna, plot_dat, site) {
  g1 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[1], ]$gtex[1]
  g2 <- gene_id_match[gene_id_match$cisnat == unlist(strsplit(dsrna, "__"))[2], ]$gtex[1]
  rv1 <- ""
  rv2 <- ""
  if (!is.na(g1)) {
    rv1 <- "NONE"
    res <- cor.test(plot_dat[,"edit"], plot_dat[,g1], alternative = "two.sided", method = "pearson")
    if (res$p.value <= 0.05 & as.numeric(res$estimate) < -0.1) {
      return("BIG_NEG") # negative correlation indicates not a dsRNA relationship, want to flag immediately, doesn't matter what other gene is doing
    }
    else if (res$p.value <= 0.05 & as.numeric(res$estimate) < 0) {
      rv1 <- "MINI_NEG"
    }
    else if (res$p.value <= 0.05 & as.numeric(res$estimate) > 0) {
      rv1 <- "POS"
    }
  }
  if (!is.na(g2)) {
    rv2 <- "NONE"
    res <- cor.test(plot_dat[,"edit"], plot_dat[,g2], alternative = "two.sided", method = "pearson")
    if (res$p.value <= 0.05 & as.numeric(res$estimate) < -0.1) {
      return("BIG_NEG") # negative correlation indicates not a dsRNA relationship, want to flag immediately, doesn't matter what other gene is doing
    }
    else if (res$p.value <= 0.05 & as.numeric(res$estimate) < 0) {
      rv2 <- "MINI_NEG"
    }
    else if (res$p.value <= 0.05 & as.numeric(res$estimate) > 0) {
      rv2 <- "POS"
    }
  }
  rval <- paste(rv1, rv2, sep = "__")
  return(rval)
}

# check_cov_corr: for check = "cov_corr":
# correlation of read coverage at an editing site and editing level
check_cov_corr <- function(dsrna, plot_dat, site) {
  rval <- "NONE"
  res <- cor.test(plot_dat[,"edit"], plot_dat[,"cov"], alternative = "two.sided", method = "pearson")
  if (res$p.value <= 0.05 & as.numeric(res$estimate) < -0.1) {
    rval <- "BIG_NEG"
  }
  else if (res$p.value <= 0.05 & as.numeric(res$estimate) < 0) {
    rval <- "MINI_NEG"
  }
  else if (res$p.value <= 0.05 & as.numeric(res$estimate) > 0) {
    rval <- "POS"
  }
  return(rval)
}

# cisnat_site_check: from a larger set of sites possibly with patterns of interest, check and only return ones that are confirmed interesting candidates
cisnat_site_check <- function(sites_set) {
  sites <- unlist(strsplit(sites_set, ";"))
  if (!any(sites %in% cand_edit$site_id)) {
    return(NA)
  }
  else {
    var_sites <- sites[sites %in% cand_edit$site_id]
    var_set <- paste(var_sites, collapse = ";")
    return(var_set)
  }
}

# var_type_check: from a set of confirmed interesting candidates, return the characteristic of interest
var_type_check <- function(var_sites_set) {
  var_sites <- unlist(strsplit(var_sites_set, ";"))
  rval <- "variable"
  if (any(grepl("POS", cand_edit[cand_edit$site_id %in% var_sites, ]$pos_corr_check))) {
    rval <- c("pos_tpm_corr", rval)
  }
  if (any(cand_edit[cand_edit$site_id %in% var_sites, ]$low_exp_check == "PASS")) {
    rval <- c("low_exp", rval)
  }
  rval <- paste(rval, collapse = ";")
  return(rval)
}

# flag_check: from a set of confirmed interesting candidates, return any possible flags/issues
flag_check <- function(var_sites_set) {
  var_sites <- unlist(strsplit(var_sites_set, ";"))
  rval <- c()
  if (any(grepl("NEG", cand_edit[cand_edit$site_id %in% var_sites, ]$pos_corr_check))) {
    rval <- c(rval, "neg_tpm_corr")
  }
  if (any(grepl("NEG", cand_edit[cand_edit$site_id %in% var_sites, ]$cov_corr_check))) {
    rval <- c(rval, "neg_coverage_corr")
  }
  if (length(rval) == 0)
    rval <- NA
  rval <- paste(rval, collapse = ";")
  return(rval)
}

# dsrna_boxplot: plotting function for editing levels across a dsRNA of interest
dsrna_boxplot <- function(dsrna) {
  sites <- unlist(strsplit(cisnat_with_var[dsrna, ]$sites, ";"))
  var_sites <- unlist(strsplit(cisnat_with_var[dsrna, ]$var_sites, ";"))
  if (length(sites) > 50) {return("Over 50 sites - exiting plot")}
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
  plot_dat$site <- as.character(plot_dat$site)
  if (grepl("chrX", plot_dat$site[1])) {
    plot_dat$pos <- as.numeric(gsub("chrX_", "", plot_dat$site))
  }
  else if (grepl("chrY", plot_dat$site[1])) {
    plot_dat$pos <- as.numeric(gsub("chrY_", "", plot_dat$site))
  }
  else {
    plot_dat$pos <- as.numeric(gsub("chr\\d+_", "", plot_dat$site))
  }
  # g1: each site is a factor (not to scale)
  plot_dat$site <- factor(plot_dat$site, levels = mixedsort(unique(plot_dat$site)))
  g1 <- ggplot(plot_dat, aes(x = site, y = edit_level, color = var_site)) + geom_boxplot() + theme_bw()
  g2 <- ggplot(plot_dat, aes(x = site, y = edit_level, color = var_site)) + geom_boxplot() + coord_cartesian(ylim = c(0,1)) + theme_bw()
  # g2: each site is a position (to scale)
  g3 <- ggplot(plot_dat, aes(x = pos, y = edit_level, color = var_site)) + geom_jitter(alpha = 0.8, width = 0.2, height = 0) + theme_bw()
  g <- (g1 / g2 / g3)
  t <- paste(dsrna, cisnat_with_var[dsrna, ]$minus_symbol, cisnat_with_var[dsrna, ]$plus_symbol)
  return(g + plot_annotation(title = t))
}


####
## Filter data down to just one tissue
# tissue of choice: adipose (subcutaneous)
samps <- unlist(strsplit(tissue_sample_match[tissue_sample_match$work_tissue == tis, ]$samples, ","))
samp_thresh <- 60 # was also thinking of round(length(samps)/10), but for eventual qtl testing would want at least 60 useable samples
cov_tis <- coverage_dat[,samps]
edit_tis <- editLevel_dat[,samps]
edit_df_tis <- edit_df[edit_df$tissue == tis, ]
exp_df_tis <- exp_df[exp_df$tissue == tis, ]
tpm_tis <- tpm_dat[,samps]



####
## Filter just to sites that have sufficient coverage and some variation in the tissue
summary(edit_df_tis$n_cov_tenPlus >= samp_thresh) # 4,015 meet this threshold, 1,589 don't
cand_edit <- edit_df_tis[edit_df_tis$n_cov_tenPlus >= samp_thresh, ]
summary(cand_edit$coefVar)
cand_edit <- cand_edit[which(cand_edit$coefVar > 0), ]
cand_edit <- cand_edit[order(cand_edit$med_edit, decreasing = T), ]
cand_edit$n_edit_moreThanZero <- cand_edit$n_cov_tenPlus - cand_edit$n_edit_zero


####
## Get a list of "red flag" dsRNAs: ones where one of the two genes has zero expression or doesn't have a GTEx annotation match
# if a site only maps to red flag dsRNAs, remove
bad_genes <- unique(exp_df_tis[exp_df_tis$n_samples_tpm_detected < samp_thresh, ]$cisnat_gene)
cn_genes <- unique(c(cisnat_withSites$plus_gene, cisnat_withSites$minus_gene))
bad_genes <- unique(c(bad_genes, cn_genes[!cn_genes %in% gene_id_match$cisnat]))
at_risk_df <- cisnat_withSites[cisnat_withSites$minus_gene %in% bad_genes | cisnat_withSites$plus_gene %in% bad_genes, ]
# ^ big dropoff because genes list is from putative_cisnat, but not many putative_cisnat have quantified editing sites
at_risk_sites <- unique(unlist(strsplit(at_risk_df$sites, ";")))
to_remove <- c()
for (i in at_risk_sites) {
  if (all(rownames(cisnat_withSites[grep(i, cisnat_withSites$sites), ]) %in% rownames(at_risk_df))) {
    to_remove <- c(to_remove, i)
  }
}
to_remove <- unique(to_remove)
cand_edit <- cand_edit[!cand_edit$site_id %in% to_remove, ]



####
## Interesting site type 1: rare low edit
# for low-edit examples, find sites where editing is only >0 in samples where the rare transcript has detectable expression
print("TYPE 1 SITES")
check_site(site = cand_edit$site_id[1], check = "low_exp")
sapply(cand_edit$site_id[1:5], check_site, check = "low_exp")
cand_edit$low_exp_check <- as.character(sapply(cand_edit$site_id, check_site, check = "low_exp"))
table(cand_edit$low_exp_check) # 56 big fail, 122 mini-fail (<5% samples checked had edit no exp), 79 pass
# ^ was suspicious that there were no moderate results, but when there are TPM == 0 events the "EDIT_NO_EXP" is common

# remove candidate edit sites with big fails
cand_edit <- cand_edit[!cand_edit$low_exp_check == "BIG_FAIL", ]
# test plot a few that passed
for (i in c(1:3)) {
  print(paste("TYPE 1 PLOT", i))
  gp <- check_site(cand_edit[cand_edit$low_exp_check == "PASS" & cand_edit$n_edit_moreThanZero >= 25,]$site_id[i], check = "plot")
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_01_low_exp_site", i, "f.png")), width = 10, height = 5, units = "in")
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_01_low_exp_site", i, "f.pdf")), width = 10, height = 5, units = "in")
}


####
## Interesting site type 2: positive correlation edit
# editing level is positively correlated with at least one of the two genes
# adopted to simultaneously filter out any sites that have negative correlations in all of their dsRNAs
# (will still need to check again later for specific dsrna/site pairs, though)
print("TYPE 2 SITES")
check_site(site = cand_edit$site_id[1], check = "pos_corr")
sapply(cand_edit$site_id[1:5], check_site, check = "pos_corr")
cand_edit$pos_corr_check <- as.character(sapply(cand_edit$site_id, check_site, check = "pos_corr"))
table(cand_edit$pos_corr_check) 
# BIG_NEG   MINI_NEG__MINI_NEG     MINI_NEG__NONE      MINI_NEG__POS     NONE__MINI_NEG         NONE__NONE 
# 157                  3                 37                  8                 26               2030 
# NONE__POS      POS__MINI_NEG          POS__NONE           POS__POS 
# 258                  7                 201                 144


# to consider: removing candidate edit sites with big negative correlations (didn't do this time)
# bigneg was based on looking at all significant corr coefs in this tissue; 
# threshold is being negative and with coefficient of greater magnitude than 0.1 (top 75%)

# test plot a few that passed
for (i in c(1:3)) {
  print(paste("TYPE 2 PLOT", i))
  gp <- check_site(cand_edit[grepl("POS", cand_edit$pos_corr_check) & cand_edit$n_edit_moreThanZero >= 25,]$site_id[i], check = "plot")
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_02_pos_corr_site", i, "f.png")), width = 10, height = 5, units = "in")
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_02_pos_corr_site", i, "f.pdf")), width = 10, height = 5, units = "in")
}


####
## Interesting site type 3: variation in editing
# this is pretty much everything else that has met the initial candidate variation thresholds
# can't think of a specific function or threshold that wouldn't penalize based on editing level


####
## Flag: sites that have negative correlations between editing and coverage
# (editing decreases as "expression" increases - variation could just reflect shifts in transcript abundance and not anything editing-related)
check_site(site = cand_edit$site_id[1], check = "cov_corr")
sapply(cand_edit$site_id[1:5], check_site, check = "cov_corr")
cand_edit$cov_corr_check <- as.character(sapply(cand_edit$site_id, check_site, check = "cov_corr"))
table(cand_edit$cov_corr_check) # 102 BIG_NEG, 57 MINI_NEG, 118 POS, the rest 2705 NONE



####
## Get df of cisnats that contain at least 1 interesting site
# prioritize by presence of multiple RNA editing sites within same putative cisnat
# also keep track of variation type (note: the combo of low exp and positive correlation is not very compelling)
cisnat_withSites$var_sites <- as.character(sapply(cisnat_withSites$sites, cisnat_site_check))
summary(is.na(cisnat_withSites$var_sites)) # 371 are NAs, 712 are not (have at least one variable site)
cisnat_withSites$n_sites <- unlist(lapply(strsplit(cisnat_withSites$sites, ";"), "length"))
cisnat_withSites$n_var_sites <- unlist(lapply(strsplit(cisnat_withSites$var_sites, ";"), "length"))
cisnat_withSites[is.na(cisnat_withSites$var_sites), ]$n_var_sites <- 0

plot_dat <- as.data.frame(t(table(cisnat_withSites$n_var_sites)))
plot_dat$Var2 <- as.integer(as.character(plot_dat$Var2))
g <- ggplot(plot_dat, aes(x = Var2, y= as.numeric(Freq))) + geom_col() + theme_bw()
ggsave(g, filename = file.path(fig_dir, paste0("07_", tis, "_z01_dsrnas_numVarSites.png")), width = 6, height = 6, units = "in")
ggsave(g, filename = file.path(fig_dir, paste0("07_", tis, "_z01_dsrnas_numVarSites.pdf")), width = 6, height = 6, units = "in")
plot_dat <- plot_dat[-1,]
g <- ggplot(plot_dat, aes(x = Var2, y= as.numeric(Freq))) + geom_col() + theme_bw()
ggsave(g, filename = file.path(fig_dir, paste0("07_", tis, "_z02_dsrnas_numVarSites_noZero.png")), width = 6, height = 6, units = "in")
ggsave(g, filename = file.path(fig_dir, paste0("07_", tis, "_z02_dsrnas_numVarSites_noZero.pdf")), width = 6, height = 6, units = "in")
print("NUMVARSITES-IN-CISNAT PLOTS DONE")


cisnat_with_var <- cisnat_withSites[cisnat_withSites$n_var_sites > 0, ]
cisnat_with_var$var_types <- as.character(sapply(cisnat_with_var$var_sites, var_type_check))
table(cisnat_with_var$var_types)
# low_exp;pos_tpm_corr;variable   low_exp;variable         pos_tpm_corr;variable            variable 
# 22                               19                           322                           349 

cisnat_with_var$site_flags <- as.character(sapply(cisnat_with_var$var_sites, flag_check))
table(cisnat_with_var$site_flags)
# NA                 neg_coverage_corr                neg_tpm_corr      neg_tpm_corr;neg_coverage_corr 
# 454                             60                       99                            99 

cisnat_with_var <- cisnat_with_var[,c("chr", "minus_gene", "minus_symbol", "minus_gene_type", "plus_gene", "plus_symbol", "plus_gene_type", "minus_tx", "minus_tx_type", "plus_tx", "plus_tx_type", "minus_exons_in_overlap", "plus_exons_in_overlap", "all_overlap_widths", "longest_overlap_width", "longest_overlap_start", "longest_overlap_end", "multi_exon_overlap", "alu_overlap", "n_sites", "n_var_sites", "var_types", "site_flags", "sites", "var_sites")]
print("CISNATS HAVE BEEN FILTERED")



####
## Plot a few example dsrnas
to_plot <- c(1, sample(c(2:nrow(cisnat_with_var)), size = 4, replace = F))
for (i in to_plot) {
  print(paste("PLOT DSRNA EXAMPLE", i))
  gp <- dsrna_boxplot(rownames(cisnat_with_var)[i])
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_03_dsrna_examples", i, "f.png")), width = 10, height = 10, units = "in")
  ggsave(gp, filename = file.path(fig_dir, paste0("07_", tis, "_03_dsrna_examples", i, "f.pdf")), width = 10, height = 10, units = "in")
}


####
## Save an output table of interesting editing variation
print("WRITING OUTPUT TABLES")
# One: at the site level
write.table(cand_edit, file = file.path(out_dir, paste0("07_perTissue_variation/", tis, "_sitesWithVariation.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")

# Two: at the dsRNA level
write.table(cisnat_with_var, file = file.path(out_dir, paste0("07_perTissue_variation/", tis, "_dsrnasWithVarSites.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")



