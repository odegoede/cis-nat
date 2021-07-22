#!/usr/bin/env Rscript

#####
## cis-NAT project
## Script pilot_02: examining pilot data of full-genome editing investigation
## By: Olivia de Goede, 2021
#####

# example usage:
# Rscript scripts/pilot_02_check_SPRINT_output.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData


####
## SET OPTIONS, LOAD LIBRARIES
options(stringsAsFactors = F)
library(GenomicRanges)
library(ggplot2)
library(patchwork)


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

load(file.path(out_dir, "03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData")) # object name from script 03: putative_cisnat 
# ^ (note: expanding back out to this file rather than 05, because 05 was filtered based on presence of GTEx editing data and this 
# is a fresh exploration of pilot data)
load(file.path(out_dir, "05_cisNAT_with_RNAedit.RData")) # two objects from script 05, named: cisnat_withSites, cisnat_withCluster
# ^ (note: this 05 file is still included for comparison of what is detected with GTEx data versus pilot data)

load("data/full_pilot/filtered_pilot_data.rds") # R objects from pilot_01 script: _filt is evidence of editing in at least 1 sample; _strict_filt is evidence of editing in at least 2 samples from the same tissue


####
## DEFINE FUNCTIONS
get_total_sites <- function(tissue_quant_file) {
  temp <- read.table(file.path(paste0("data/full_pilot/tissue_quant/", tissue_quant_file)), 
                     header = T, sep = "\t")
  return(unique(temp$site_name))
}



####
## create a cisnat_with_pilotSites df and cisnat_with_strictPilotSites df
pilot_site_anno <- data.frame(site_name = rownames(pilot_edit_filt), 
                              chr = unlist(lapply(strsplit(rownames(pilot_edit_filt), "_"), "[[", 1)), 
                              start = as.integer(unlist(lapply(strsplit(rownames(pilot_edit_filt), "_"), "[[", 2))), 
                              end = as.integer(unlist(lapply(strsplit(rownames(pilot_edit_filt), "_"), "[[", 3))))
strict_pilot_site_anno <- data.frame(site_name = rownames(pilot_edit_strict_filt), 
                                     chr = unlist(lapply(strsplit(rownames(pilot_edit_strict_filt), "_"), "[[", 1)), 
                                     start = as.integer(unlist(lapply(strsplit(rownames(pilot_edit_strict_filt), "_"), "[[", 2))), 
                                     end = as.integer(unlist(lapply(strsplit(rownames(pilot_edit_strict_filt), "_"), "[[", 3))))


pilot_site_grange <- makeGRangesFromDataFrame(pilot_site_anno, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = T)
strict_pilot_site_grange <- makeGRangesFromDataFrame(strict_pilot_site_anno, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = T)

comp_one <- GRanges(c("chr12", "chr12"), IRanges(start = c(6666730, 6669031), end = c(6668115, 6669184), 
                                                 names = c("ENST00000361959.7__ENST00000396799.3", "ENST00000361959.7__ENST00000396799.3")))
comp_two <- GRanges(rep("chr9", 3), IRanges(start = c(2622100,2635453,2639859), end = c(2622271,2635572,2639981), 
                                            names = rep("ENST00000599229.2__ENST00000382096.5", 3)))
temp <- putative_cisnat[-which(rownames(putative_cisnat) %in% c("ENST00000361959.7__ENST00000396799.3", "ENST00000599229.2__ENST00000382096.5")), ]
put_cisnat_grange <- makeGRangesFromDataFrame(temp, seqnames.field = "chr", start.field = "longest_overlap_start", 
                                              end.field = "longest_overlap_end", keep.extra.columns = T)

## pilot
# recall sites (1,229,653 total) had to have real, measured values in more than one sample (even if that measurement's editing level was 0), and
# some non-zero editing level in at least one sample
put_cisnat_to_pilotSite <- findOverlaps(put_cisnat_grange, pilot_site_grange, select = "all", ignore.strand = T)
put_cisnat_ind <- queryHits(put_cisnat_to_pilotSite)
site_ind <- subjectHits(put_cisnat_to_pilotSite)
# random sanity check:
put_cisnat_grange[18]
pilot_site_grange[site_ind[put_cisnat_ind == 18]]

length(unique(put_cisnat_ind)) # only 484 of the overlaps contain at least 1 RNA editing site - with GTEx data it was 1,083
# what about the two multi-exon overlaps?:
findOverlaps(comp_one, pilot_site_grange, select = "all", ignore.strand = T) # no overlaps
findOverlaps(comp_two, pilot_site_grange, select = "all", ignore.strand = T) # no overlaps
# 484 / 15,986 contain pilot RNA-editing sites

to_keep <- names(put_cisnat_grange[unique(put_cisnat_ind)])
cisnat_with_pilotSites <- putative_cisnat[to_keep, ]
head(cisnat_with_pilotSites)
# how many unique gene pairs is this?
length(unique(paste(cisnat_with_pilotSites$minus_gene, cisnat_with_pilotSites$plus_gene, sep = "_"))) # 222 gene pairs with editing site
length(unique(paste(putative_cisnat$minus_gene, putative_cisnat$plus_gene, sep = "_"))) # 6,155 gene pairs in original putative list
# how many overlap a possible IR-Alu?
summary(is.na(cisnat_with_pilotSites$alu_overlap)) # TRUE for 459 (no Alu overlap), FALSE for 25 (Alu overlap)

# add site info to df
cisnat_with_pilotSites$pilot_sites <- NA
for (i in rownames(cisnat_with_pilotSites)) {
  cisnat_with_pilotSites[i, ]$pilot_sites <- paste(pilot_site_grange[site_ind[put_cisnat_ind == which(names(put_cisnat_grange) == i)]]$site_name, collapse = ";")
}
cisnat_with_pilotSites$n_sites <- as.integer(unlist(lapply(strsplit(cisnat_with_pilotSites$pilot_sites, ";"), length)))
# ^ n_sites is number of pilot sites only
length(unique(unlist(strsplit(cisnat_with_pilotSites$pilot_sites, ";")))) # 3,569 pilot sites in transcript overlaps

# how many are also in cisnat_withSites?
summary(rownames(cisnat_with_pilotSites) %in% rownames(cisnat_withSites)) # 76 are shared, 408 of pilot cisnats are not in GTEx (and around 1000 are in GTEx only)
# add this as a column
cisnat_with_pilotSites$has_site_in_gtexAnalysis <- rownames(cisnat_with_pilotSites) %in% rownames(cisnat_withSites)

# save these results
save(cisnat_with_pilotSites, file = file.path(out_dir, "pilot_02_putative_cisnat_with_pilotSites.RData"))
to_write <- cisnat_with_pilotSites[,c(1:19, 22, 21, 20)]
write.table(to_write, file = file.path(out_dir, "pilot_02_putative_cisnat_with_pilotSites.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")


## strict pilot
# recall strict sites (84,938 total) require that at least 2 samples of the same tissue type have editing values (not just at least 2 samples overall)
put_cisnat_to_strictPilotSite <- findOverlaps(put_cisnat_grange, strict_pilot_site_grange, select = "all", ignore.strand = T)
put_cisnat_ind <- queryHits(put_cisnat_to_strictPilotSite)
site_ind <- subjectHits(put_cisnat_to_strictPilotSite)
# random sanity check:
put_cisnat_grange[321]
strict_pilot_site_grange[site_ind[put_cisnat_ind == 321]]

length(unique(put_cisnat_ind)) # only 79 of the overlaps contain at least 1 "strict" RNA editing site

to_keep <- names(put_cisnat_grange[unique(put_cisnat_ind)])
cisnat_with_strictPilotSites <- putative_cisnat[to_keep, ]
head(cisnat_with_strictPilotSites)
# how many unique gene pairs is this?
length(unique(paste(cisnat_with_strictPilotSites$minus_gene, cisnat_with_strictPilotSites$plus_gene, sep = "_"))) # 33 gene pairs with editing site
# how many overlap a possible IR-Alu?
summary(is.na(cisnat_with_strictPilotSites$alu_overlap)) # TRUE for 71 (no Alu overlap), FALSE for 8 (Alu overlap)

# add site info to df
cisnat_with_strictPilotSites$pilot_sites <- NA
for (i in rownames(cisnat_with_strictPilotSites)) {
  cisnat_with_strictPilotSites[i, ]$pilot_sites <- paste(strict_pilot_site_grange[site_ind[put_cisnat_ind == which(names(put_cisnat_grange) == i)]]$site_name, collapse = ";")
}
cisnat_with_strictPilotSites$n_sites <- as.integer(unlist(lapply(strsplit(cisnat_with_strictPilotSites$pilot_sites, ";"), length)))
length(unique(unlist(strsplit(cisnat_with_strictPilotSites$pilot_sites, ";")))) # 215 strict pilot sites in transcript overlaps

# how many are also in cisnat_withSites
summary(rownames(cisnat_with_strictPilotSites) %in% rownames(cisnat_withSites)) # 18 are shared, 61 of strict pilot cisnats are not in GTEx
# add this as a column
cisnat_with_strictPilotSites$has_site_in_gtexAnalysis <- rownames(cisnat_with_strictPilotSites) %in% rownames(cisnat_withSites)

# save these results
save(cisnat_with_strictPilotSites, file = file.path(out_dir, "pilot_02_putative_cisnat_with_strictPilotSites.RData"))
to_write <- cisnat_with_strictPilotSites[,c(1:19, 22, 21, 20)]
write.table(to_write, file = file.path(out_dir, "pilot_02_putative_cisnat_with_strictPilotSites.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")



####
## Panel plot comparing number of original GTEx editing sites in putative tx overlaps to pilot full-genome editing sites
# NOTE about the original GTEx sites: didn't do any filtering by evidence for editing, because they are confirmed/curated editing sites; 
#     some of these sites could be zero edit across all samples

# (top panel = GTEx, middle panel = pilotSites, bottom panel = strictPilotSites)
# (x-axis is category bins of number of sites; 0, 1, 2-5; 6-10, 11-20, 21-50, 51+)
# stacked bars filled by: has editing sites in both; has editing sites in GTEx only; has editing sites in pilot only

cisnat_withSites$n_sites <- as.integer(unlist(lapply(strsplit(cisnat_withSites$sites, ";"), length)))

# unique tx_pair level
plot_dat <- data.frame(tx_pair = c(rownames(cisnat_withSites), rownames(cisnat_with_pilotSites), rownames(cisnat_with_strictPilotSites)),
                       n_sites = c(cisnat_withSites$n_sites, cisnat_with_pilotSites$n_sites, cisnat_with_strictPilotSites$n_sites),
                       grp = c(rep("gtex", nrow(cisnat_withSites)), rep("pilot", nrow(cisnat_with_pilotSites)), rep("strict_pilot", nrow(cisnat_with_strictPilotSites))),
                       n_site_x_cat = NA, found_in = NA)
shared_pairs <- intersect(rownames(cisnat_withSites), rownames(cisnat_with_pilotSites))
plot_dat[plot_dat$tx_pair %in% shared_pairs, ]$found_in <- "both"
plot_dat[grepl("pilot", plot_dat$grp) & is.na(plot_dat$found_in), ]$found_in <- "pilot_only"
plot_dat[plot_dat$grp == "gtex" & is.na(plot_dat$found_in), ]$found_in <- "gtex_only"
# x_cat: 1, 2-5; 6-10, 11-20, 21-50, 51+
plot_dat[plot_dat$n_sites == 1, ]$n_site_x_cat <- "1"
plot_dat[plot_dat$n_sites >= 2 & plot_dat$n_sites <= 5, ]$n_site_x_cat <- "2-5"
plot_dat[plot_dat$n_sites >= 6 & plot_dat$n_sites <= 10, ]$n_site_x_cat <- "6-10"
plot_dat[plot_dat$n_sites >= 11 & plot_dat$n_sites <= 20, ]$n_site_x_cat <- "11-20"
plot_dat[plot_dat$n_sites >= 21 & plot_dat$n_sites <= 50, ]$n_site_x_cat <- "21-50"
plot_dat[plot_dat$n_sites >= 51, ]$n_site_x_cat <- "51+"
plot_dat$n_site_x_cat <- factor(plot_dat$n_site_x_cat, levels = c("1", "2-5", "6-10", "11-20", "21-50", "51+"))
table(plot_dat$n_site_x_cat)
# 1   2-5  6-10 11-20 21-50   51+ 
#   366   388   349   263   228    52 
g <- ggplot(plot_dat, aes(x = n_site_x_cat, fill = found_in)) + geom_bar(position = "stack") + 
    facet_wrap(~grp, nrow = 3) + ylab("N unique tx_pairs") + theme_bw()
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_tx_pair_with_edit_sites.png"), height = 9, width = 6)
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_tx_pair_with_edit_sites.pdf"), height = 9, width = 6)

# no facets:
plot_dat <- data.frame(tx_pair = c(rownames(cisnat_withSites), rownames(cisnat_with_pilotSites)),
                       sites = c(cisnat_withSites$sites, cisnat_with_pilotSites$pilot_sites),
                       grp = c(rep("gtex", nrow(cisnat_withSites)), rep("pilot", nrow(cisnat_with_pilotSites))),
                       n_sites = c(cisnat_withSites$n_sites, cisnat_with_pilotSites$n_sites),
                       n_site_x_cat = NA, found_in = NA)
shared_pairs <- intersect(rownames(cisnat_withSites), rownames(cisnat_with_pilotSites))
plot_dat[plot_dat$tx_pair %in% shared_pairs, ]$found_in <- "both"
plot_dat[grepl("pilot", plot_dat$grp) & is.na(plot_dat$found_in), ]$found_in <- "pilot_only"
plot_dat[plot_dat$grp == "gtex" & is.na(plot_dat$found_in), ]$found_in <- "gtex_only"
for (i in shared_pairs) {
  tmp <- paste(unique(unlist(strsplit(plot_dat[plot_dat$tx_pair == i, ]$sites, ";"))), collapse = ";")
  plot_dat[plot_dat$tx_pair == i, ]$n_sites <- length(unlist(strsplit(tmp, ";")))
  plot_dat[plot_dat$tx_pair == i, ]$sites <- tmp
}
plot_dat <- plot_dat[!(plot_dat$grp == "pilot" & plot_dat$found_in == "both"), ]
plot_dat[plot_dat$n_sites == 1, ]$n_site_x_cat <- "1"
plot_dat[plot_dat$n_sites >= 2 & plot_dat$n_sites <= 5, ]$n_site_x_cat <- "2-5"
plot_dat[plot_dat$n_sites >= 6 & plot_dat$n_sites <= 10, ]$n_site_x_cat <- "6-10"
plot_dat[plot_dat$n_sites >= 11 & plot_dat$n_sites <= 20, ]$n_site_x_cat <- "11-20"
plot_dat[plot_dat$n_sites >= 21 & plot_dat$n_sites <= 50, ]$n_site_x_cat <- "21-50"
plot_dat[plot_dat$n_sites >= 51, ]$n_site_x_cat <- "51+"
plot_dat$n_site_x_cat <- factor(plot_dat$n_site_x_cat, levels = c("1", "2-5", "6-10", "11-20", "21-50", "51+"))
g <- ggplot(plot_dat, aes(x = n_site_x_cat, fill = found_in)) + geom_bar(position = "stack") + 
  ylab("N unique tx_pairs") + xlab("N editing sites") + 
  scale_fill_manual(values = c("#c84cff", "#ff4c4c", "#4c6cff")) + theme_bw()
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_tx_pair_with_edit_sites_noFacet.png"), height = 6, width = 7)
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_tx_pair_with_edit_sites_noFacet.pdf"), height = 6, width = 7)


# unique gene_pair level
putative_cisnat$gene_pair <- paste(putative_cisnat$minus_gene, putative_cisnat$plus_gene, sep = "__")
plot_dat$gene_pair <- putative_cisnat[plot_dat$tx_pair, ]$gene_pair
plot_dat <- plot_dat[order(plot_dat$n_sites), ]
plot_dat$filter_cat <- paste(plot_dat$gene_pair, plot_dat$grp)
summary(duplicated(plot_dat$filter_cat)) # 905 duplicate gene pair/grp rows, 741 unique
plot_dat <- plot_dat[!duplicated(plot_dat$filter_cat), ]
shared_pairs <- intersect(plot_dat[plot_dat$grp == "gtex", ]$gene_pair, plot_dat[plot_dat$grp == "pilot", ]$gene_pair)
plot_dat$found_in <- NA
plot_dat[plot_dat$gene_pair %in% shared_pairs, ]$found_in <- "both"
plot_dat[grepl("pilot", plot_dat$grp) & is.na(plot_dat$found_in), ]$found_in <- "pilot_only"
plot_dat[plot_dat$grp == "gtex" & is.na(plot_dat$found_in), ]$found_in <- "gtex_only"
(g <- ggplot(plot_dat, aes(x = n_site_x_cat, fill = found_in)) + geom_bar(position = "stack") + 
    facet_wrap(~grp, nrow = 3) + ylab("N unique tx_pairs") + theme_bw())
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_gene_pair_with_edit_sites.png"), height = 9, width = 6)
ggsave(g, file = file.path(fig_dir, "pilot_02_unique_gene_pair_with_edit_sites.pdf"), height = 9, width = 6)



####
## Check out ulcerative colitis candidate lncRNAs in full-genome pilot data
# relevant tissues: Colon_Transverse, Colon_Sigmoid, Small_Intestine_Terminal_Ileum, Spleen, Minor_Salivary_gland
genes <- rownames(gene_anno[gene_anno$symbol %in% c("LINC01475", "AL513542.1"), ])
tx_anno[tx_anno$ensgene %in% genes, ]

for (tx in rownames(tx_anno[tx_anno$ensgene %in% genes, ])) {
  print(paste(tx, tx_anno[tx, ]$symbol))
  print(summary(grepl(tx, rownames(cisnat_with_pilotSites))))
}
# nope, not present in pilot data either

# were they in total sites and then filtered out?
total_sites <- c()
for (tis in tissue_sample_match$work_tissue) {
  print(tis)
  file_name <- grep(tis, list.files("data/full_pilot/tissue_quant/"), value = T)
  total_sites <- union(total_sites, get_total_sites(file_name))
}
chrTen_sites <- total_sites[grep("chr10", total_sites)]
chrTen_pos <- as.integer(unlist(lapply(strsplit(chrTen_sites, "_"), "[[", 2)))
summary(chrTen_pos < max(tx_anno[tx_anno$ensgene %in% genes, ]$end) & chrTen_pos > min(tx_anno[tx_anno$ensgene %in% genes, ]$start))
# ^ 3 sites are within the boundaries of these genes

soi <- chrTen_sites[chrTen_pos < max(tx_anno[tx_anno$ensgene %in% genes, ]$end) & chrTen_pos > min(tx_anno[tx_anno$ensgene %in% genes, ]$start)]
soi # "chr10_99529924_99529925" "chr10_99529935_99529936" "chr10_99529939_99529940"
grep(rownames(tx_anno[tx_anno$symbol == "AL513542.1", ]), rownames(putative_cisnat))
putative_cisnat[grep(rownames(tx_anno[tx_anno$symbol == "AL513542.1", ]), rownames(putative_cisnat)), ]
# ^ overlap ranges from 99527994-99528121; the soi are all downstream of the overlap
tx_anno[tx_anno$ensgene %in% genes, ] # soi only overlap with an isoform of linc01475 that is much longer


