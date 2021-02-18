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
<ol>
	<li>Annotate all possible overlapping transcripts. Flag ones that contain Alu repeats (if any).</li>
	<li>Filter & combine each RNA editing quantification file (1 per sample, ~17k files) into one workable file.</li>
	<li>Assess overall editing level and variation, both site-level and region-level.</li>
	<li>Search for genetic variants associated with RNA editing (site-level and region-level)</li>
	<li>Test out calculating a "total immunogenic dsRNA score"; consistency within individual across tissues? applicable to external cohorts with immune trait information and RNA-seq data? (TODO: may need to bring in RNA editing quantification from IRAlus as well)</li>
</ol>


### Non-Goals
<ol>
	<li>Does not incorporate the steps of RNA-seq read mapping, or quantification of RNA editing level per site.</li>
	<li>Currently not intending this to be a pipeline that anyone unaffiliated with the project could run.</li>
</ol>

### Future goals
Connect to Qin's code for read mapping and RNA editing quantification, and build a start-to-finish pipeline (from raw RNA-seq data --> dsRNA-level editing scores --> identifying variation in editing --> testing for association with genetic variation --> testing for association with immune traits --> total immunogenic burden score?)


## Detailed Design
### User requirements
Required inputs:
<ol>
	<li>Gene annotation file (or a link to one); assuming GTF type</li>
	<li>Sample information file; must have sample ID field that can be matched to RNA editing quantification file names, and at least fields for tissue and individual</li>
	<li>(Path to directory of) RNA editing quantification files (1/sample); columns: chrom, position, coverage, editedreads, editlevel</li>
</ol>
For current goals, users are assumed to be people affiliated with the project (myself & Qin), or the readers of an eventual paper about this project (for code transparency, not intended for their use).

### New/changed data structures?
Individual RNA editing quant files (.txt.gz) will be combined into one large file (.txt.gz). Output files may be saved as .RData or .txt files for easy access or for tables.

Files will be stored on Oak.

### What APIs will you use/change?
n/a

### Throughput/latency/cost/efficiency concerns?
n/a

### Data validation/what are potential error states?
<ol>
	<li>Gene annotation file needs to have transcripts on opposite strands that overlap</li>
	<li>Sample information file needs to have sample IDs that can be matched to RNA editing quantification file names, and fields for tissue and individual</li>
	<li>RNA editing quantification files should at least have columns for: chrom, position, coverage, and editedreads (editlevel can be easily calculated); ideally, all quant files will have the same number of rows (i.e. be providing data for the same number of sites) - if not, this should be filled in with NAs so we don't get file sorting/filtering errors</li>
</ol>

### Logging/monitoring/observability
n/a

### Security/Privacy
Intended to be run by project collaborators on Oak - no actual GTEx data will be put online.

### What will you test?
1. Is there variation in expression or editing level of cis-NATs?
2. Is this editing variation connected to genetic variants?
3. Can we combine this expression/editing data across all dsRNAs to create a meaningful total immunogenic dsRNA score?

## Third Party dependencies
<ol>
	<li>R libraries for plotting (ggplot2, patchwork); reading in data (data.table + [temp_unzip code](https://gist.github.com/xhdong-umd/6429e7f96735142fa467f3b1daa91a2c)); </li>
	<li>FastQTL</li>


## Work Estimates

## Alternative Approaches

## Related Work?
[Code related to original GTEx RNA editing quantifcation, for known editing sites](https://github.com/vargasliqin/mpileup)



## Progress / Analysis Steps:
1. Identify all transcript overlaps using genome annotation: [01_find_tx_overlaps.R](https://github.com/odegoede/cis-nat/blob/main/find_overlapping_transcripts/01_find_tx_overlaps.R)
2. Qin: re-map clusters of hyper-editing within these transcript overlaps - opportunity to improve approach
3. With updated RNA-editing cluster IDs, quantify overall levels of editing and variation (site-level and region-level)
4. Filter overlaps to candidate dsRNAs that show inter-individual variability in either: RNA-editing level or expression level
5. Test for associations between genetic variation and editing variation (editing quantitative trait loci, edQTLs).
6. Evaluate a "total immunogenic dsRNA score" within GTEx; consistency across tissues?
7. Examine variable dsRNAs and total dsRNA scores in external cohort with immune trait information and RNA-seq data

## Exploratory Code:
<ol type="A">
  <li>Explore overlap attributes (follows step 1)</li>
  <li>Assess RNA-editing coverage/levels with current quantification (follows exScript A)</li>
  <li>Select candidate dsRNAs with current quantification: some edQTLs, some not; some variable in editing, some variable in expression, some both (follows exScript B)</li>
  <li>Explore editing/expression patterns of candidate dsRNAs (follows exScript C)</li>
</ol>



