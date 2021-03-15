# cis-nat: Exploration of putative immunogenic double-stranded RNAs arising from cis-NATs #

## Overview
The innate immune receptor MDA5 binds to double-stranded RNAs (dsRNAs) and initiates immune responses upon their detection, as dsRNAs are often viral genomes. However, the human genome also encodes dsRNAs, either through inverted repetitive sequences or through overlapping transcripts (cis-natural antisense transcripts, or cis-NATs). The mechanism of RNA editing serves to mark “self” dsRNA from “non-self”; if human dsRNAs are not sufficiently edited, they may trigger inappropriate inflammation. Recent work from Billy Li's lab has found a disproportionate number of cis-NATs are implicated in immune-related diseases, which was an unexpected discovery given that they represent the minority of RNA editing sites.

Our goal is to identify all possible cis-NATs across the genome, then build a comprehensive list of cis-NATs that also have clusters of RNA editing. Since RNA editing of “self” RNA prevents MDA5 binding and the subsequent innate immune signaling, the presence of RNA editing clusters will mark cis-NATs that could be immunogenic dsRNAs. Variation in editing of these cis-NATs could influence autoimmune disease risk; defining these relationships would add to our understanding of complex immune-related diseases.


## Background 
### Motivations
We still don't know the molecular mechanisms behind many disease-associated variants' effects on disease risk. Variation in RNA editing machinery is known to be connected to inflammation; within the self-encoded double-stranded RNAs (dsRNAs) of the human genome, cis-NATs are understudied and overrepresented in dsRNAs associated with immune-related diseases. We want to study this group of dsRNAs and identify genetic variants that increase susceptibility to immune-related diseases through their effects on RNA editing.

### Other solutions
n/a

### Current Goals
- Annotate all possible overlapping transcripts. Flag ones that contain Alu repeats (if any).
- Filter & combine each RNA editing quantification file (1 per sample, ~17k files) into one workable file.
- Assess overall editing level and variation, both site-level and region-level.
- Search for genetic variants associated with RNA editing (site-level and region-level)
- Test out calculating a "total immunogenic dsRNA score"; consistency within individual across tissues? applicable to external cohorts with immune trait information and RNA-seq data? (TODO: may need to bring in RNA editing quantification from IRAlus as well)

### Non-Goals
Will not incorporate the steps of RNA-seq read mapping, or quantification of RNA editing level per site. Also currently not intending this to be a pipeline that anyone unaffiliated with the project could run.

### Future goals
Connect to Qin's code for read mapping and RNA editing quantification, and build a start-to-finish pipeline (from raw RNA-seq data --> dsRNA-level editing scores --> identifying variation in editing --> testing for association with genetic variation --> testing for association with immune traits --> total immunogenic burden score?)


## Detailed Design
### User requirements
Required inputs:
- Gene annotation file (or a link to one); assuming GTF type
- Sample information file; must have sample ID field that can be matched to RNA editing quantification file names, and at least fields for tissue and individual
- (Path to directory of) RNA editing quantification files (1/sample); columns: chrom, position, coverage, editedreads, editlevel

For current goals, users are assumed to be people affiliated with the project, or the readers of an eventual paper about this project (for code transparency, not intended for their use).

### New/changed data structures?
Individual RNA editing quant files (.txt.gz) will be combined into one large file (.txt.gz). Output files may be saved as .RData or .txt files for easy access or for tables.

Files will be stored on Oak.

### What APIs will you use/change?
n/a

### Throughput/latency/cost/efficiency concerns?
n/a

### Data validation/what are potential error states?
- Gene annotation file needs to have transcripts on opposite strands that overlap
- Sample information file needs to have sample IDs that can be matched to RNA editing quantification file names, and fields for tissue and individual
- RNA editing quantification files should at least have columns for: chrom, position, coverage, and editedreads (editlevel can be easily calculated); ideally, all quant files will have the same number of rows (i.e. be providing data for the same number of sites) - if not, this should be filled in with NAs so we don't get file sorting/filtering errors

### Logging/monitoring/observability
n/a

### Security/Privacy
Intended to be run by project collaborators on Oak - no actual GTEx data will be put online.

### What will you test?
- necessary files were provided
- for a given script, the necessary directories exist
- files behave as expected (e.g. GTF fields: seqname, source, feature, start, end, score, strand, frame, attribute; e.g. quant file fields: chrom, position, coverage, editedreads, editlevel)
- sample annotation file has the colnames for sample name, tissue, and individual specified. if the sample name isn't the same as the quant file names, needs to also have a quant file name field. (I think this is necessary - depends on how the code for combining quant files shakes out, might instead name columns the same as the sample name field)


## Third Party dependencies
- R version 3.6.1 was used in this project
- R libraries for:
  - parsing arguments: optparse
  - plotting: ggplot2, patchwork
  - data input & structure: data.table, readxl
  - data organization: dplyr, magrittr, purrr
  - working with genomic coordinates data: GenomicRanges
- FastQTL


## Work Estimates

## Alternative Approaches

## Related Work?
[Code related to original GTEx RNA editing quantifcation, for known editing sites](https://github.com/vargasliqin/mpileup)


## Progress / Analysis Steps:
1. Identify all transcript overlaps using genome annotation: [01_find_tx_overlaps.R](https://github.com/odegoede/cis-nat/blob/main/scripts/01_find_tx_overlaps.R), [02_save_overlap_locations.R](https://github.com/odegoede/cis-nat/blob/main/scripts/02_save_overlap_locations.R)
1. Flag overlaps that also contain IRAlus, examine overlap characteristics: [03_flag_tx_overlaps_alu.R](https://github.com/odegoede/cis-nat/blob/main/scripts/03_flag_tx_overlaps_alu.R), [04_examine_tx_overlaps.R](https://github.com/odegoede/cis-nat/blob/main/scripts/04_examine_tx_overlaps.R)
1. Qin: re-map clusters of hyper-editing within these transcript overlaps - opportunity to improve approach?
1. Filter and combine each RNA editing quantification file (1 per sample, ~17k files) into one workable file.
1. With updated RNA-editing cluster IDs, quantify overall levels of editing and variation (site-level and region-level)
1. Filter overlaps to candidate dsRNAs that show inter-individual variability in either: RNA-editing level or expression level
1. Test for associations between genetic variation and editing variation (editing quantitative trait loci, edQTLs).
1. Evaluate a "total immunogenic dsRNA score" within GTEx; consistency across tissues?
1. Examine variable dsRNAs and total dsRNA scores in external cohort with immune trait information and RNA-seq data


## Exploratory Code:
<ol type="A">
  <li>Check how these transcript overlaps are annotated in GTEx RNA-editing paper (follows step 1): [exA_check_paper_dsrnas_in_tx_overlaps.R](https://github.com/odegoede/cis-nat/blob/main/scripts/exA_check_paper_dsrnas_in_tx_overlaps.R) </li> 
  <li>Assess RNA-editing coverage/levels with current quantification (follows exScript A)</li>
  <li>Select candidate dsRNAs with current quantification: some edQTLs, some not; some variable in editing, some variable in expression, some both (follows exScript B)</li>
  <li>Explore editing/expression patterns of candidate dsRNAs (follows exScript C)</li>
</ol>

