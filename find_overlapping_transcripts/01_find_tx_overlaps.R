#####
## cis-NAT project
## Script 01: Identifying all overlapping transcripts in human genome
## By: Olivia de Goede
#####

## Criteria for overlaps:
# 1. transcripts come from opposite strands
# 2. have at least 100 bp *continuous* overlap

####
## Load libraries, set working directory
setwd(path/to/working/directory)
options(stringsAsFactors = F)

library(GenomicRanges)
library(data.table)


####
## Read in GENCODE GTF
# Version 35
