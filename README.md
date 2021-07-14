# cis-nat: Exploration of putative immunogenic double-stranded RNAs arising from cis-NATs #

Our goal is to identify all possible overlapping transcripts (cis-natural antisense transcripts, or cis-NATs) across the genome, then build a comprehensive list of cis-NATs that also have clusters of RNA editing. 

Since RNA editing of “self” RNA prevents MDA5 binding and the subsequent innate immune signaling, the presence of RNA editing clusters will mark cis-NATs that could be immunogenic dsRNAs. Variation in editing of these cis-NATs could influence autoimmune disease risk; defining these relationships would add to our understanding of complex immune-related diseases.

The project's design document (includes dependencies) is [here](https://github.com/odegoede/cis-nat/blob/main/design.md).

## Suggested working directory set-up
```
.
├── data                     # Preliminary data files (annotation files, sample info, quant files)
├── figures                  # To write figures to
│   └── exploratory            # (subdir for figures from exploratory scripts)
├── output                   # To write output files to
├── scratch                  # Temporary directory for intermediate files
└── scripts                  # Scripts
```

## Run all R scripts with one command

This uses the script cisnat_script_launcher.sh. It runs scripts 01-08, which examine the GTEx editing data, and does not include the scripts examining the pilot data evaluating all detectable editing sites in the genome. 

To run from working directory:

```sh
scripts/cisnat_script_launcher.sh path/to/work_dir > cisnat_log 2>&1
```


## Walkthrough for code so far:

### Part 1: Identify all possible transcript overlaps using genome annotation

First, find overlap locations.

```sh
$ Rscript scripts/01_find_tx_overlaps.R --gtf data/gencode.v35.annotation.gtf.gz --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --keepanno both
```
```
Usage: scripts/01_find_tx_overlaps.R [options]

Options:
        --gtf=GTF
                GTF file [default "NULL"]

        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]

        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --keepanno=KEEPANNO
                Should the gene, transcript, and exon annotations be saved? Files will be saved to outdir. "RData" for only .RData, "both" for both .txt and .RData files. [default "both"]

        -h, --help
                Show this help message and exit
```
Required input is a GTF file.

Key outputs are:
- annotation files for genes, transcript, and exon (if --keepanno is not "neither")
- text file and R object of all possible transcript overlaps (01_putative_cisNAT_uniqueRegionsOnly.txt and 01_putative_cisNAT_uniqueRegionsOnly.RData)


Second, save overlap locations in BED file format.

```sh
$ Rscript scripts/02_save_overlap_locations.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/
```
```
Usage: scripts/02_save_overlap_locations.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        -h, --help
                Show this help message and exit
```
Key output is a file of genomic regions that cover the transcript overlaps, with minimum BED file columns (02_overlappingTranscript_genomicRanges_1start_fullClose.txt)


### Part 2: Flag overlaps that also contain IRAlus, and examine characteristics of the transcript overlaps

First, flag transcript overlaps that contain inverted repeat Alus (since the IRAlu could be the actual driver of a dsRNA, so these overlaps would need to be examined in-depth).

```sh
$ Rscript scripts/03_flag_tx_overlaps_alu.R --alufile data/hg38_repeats.Alu.bed.gz --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/
```
```
Usage: scripts/03_flag_tx_overlaps_alu.R [options]

Options:
        --alufile=ALUFILE
                File with genomic annotation of Alu regions [default "NULL"]

        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        -h, --help
                Show this help message and exit
```
Required input is a file with genomic locations of Alu sequences. Key output is initial list of putative cis-NATs, with flags for overlaps with IR-Alu sequences (03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData).


Second, examine overlap characteristics (like length of overlap, which gene types are involved, etc.).

```sh
$ Rscript scripts/04_examine_tx_overlaps.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --edittable data/gtex_editing_paper_tableS4.xlsx --coloctable data/gtex_lncrna_paper_hits_summary_table_allGeneTypes.RData
```
```
Usage: scripts/04_examine_tx_overlaps.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --edittable=EDITTABLE
                File of cis-NATs from RNA-editing paper (optional) [default "NULL"]

        --coloctable=COLOCTABLE
                File of QTL colocalization from GTEx lncRNA paper (optional) [default "NULL"]

        -h, --help
                Show this help message and exit
```
Required inputs: a table of the identified cis-NATs in the GTEx editing paper (an Excel file), and a file with the colocalization results from the GTEx lncRNA paper (an RData file). Key outputs are mostly figures.


### Part 3: Incorporate the editing measurements, find overlaps with some indication of RNA editing.

First, narrow identified overlaps down to those that contain at least one RNA editing site.

```sh
$ Rscript scripts/05_check_editing_sites.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --editanno data/gtex_edit/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.AnnoVar_Genes.Clusered.txt.gz --examplequant data/gtex_edit/indiv_files/GTEX-11ZTS-1026-SM-5LU8O.txt.gz 
```
```
Usage: scripts/05_check_editing_sites.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --editanno=EDITANNO
                Annotation of GTEx-identified editing clusters.

        --examplequant=EXAMPLEQUANT
                Example of an editing quantification file for one GTEx sample.

        -h, --help
                Show this help message and exit
```
Required inputs: annotation of editing clusters (as identified by GTEx); one editing quantification file for one GTEx sample (as an example). Key outputs: updated list of putative cis-NATs, to just ones with some evidence of RNA editing (05_cisNAT_with_RNAedit.RData). Also a list of the editing sites in these putative cis-NATs, which was important for filtering the individual GTEx editing quantification files down to just the sites of interest.


Second, prepare files on read coverage and editing level, to more easily query sites of interest.

```sh
$ Rscript scripts/06_check_editing_quantification.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --editcoverage data/gtex_edit/coverage_matrix.txt.gz --editlevel data/gtex_edit/editLevel_matrix.txt.gz --tpmfile data/cisnat_gene_tpm.gct.gz --gtexanno data/rnaseqc_genes_info.RData
```
```
Usage: scripts/06_check_editing_quantification.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --samplefile=SAMPLEFILE
                File with GTEx sample information

        --editcoverage=EDITCOVERAGE
                File with read coverage for the editing sites of interest

        --editlevel=EDITLEVEL
                File with the editing levels for the editing sites of interest

        --tpmfile=TPMFILE
                Gene expression values (TPM) for the genes of interest

        --gtexanno=GTEXANNO
                GTEx gene annotation file (uses GENCODE v26, whereas putative cis-NAT identification was done with GENCODE v35)

        -h, --help
                Show this help message and exit
```
Required inputs: Editing data (both read coverage and editing level), gene expression levels, GTEx sample information, and GTEx gene annotation file. 

Key outputs: 
- big data frames with coverage and editing level info (rows are editing sites of interest, columns are GTEx samples)
- data frame matching GTEx sample names to tissue types
- data frame matching GTEx gene IDs to gene IDs used in putative cis-NAT identification
- summaries of gene expression across tissues (median, and meeting expression threshold of TPM >0.1 in >20% of samples)
- edit_df, a data frame summarizing the most important editing level and coverage metrics - each row is an edit site * tissue combo
- exp_df, a data frame summarizing key expression metrics (expression avg and stdev, number of samples meeting expression thresholds, etc.) - each row is a gene * tissue combo


Third, explored ways to identify sites with interesting variation in editing levels, and to flag sites where the data looks a little bit unreliable.

```sh
$ Rscript scripts/07_example_variation_search.R --tissue $tis --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz
```
```
Usage: scripts/07_example_variation_search.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --samplefile=SAMPLEFILE
                File with GTEx sample information

        --tpmfile=TPMFILE
                Gene expression values (TPM) for the genes of interest

        -h, --help
                Show this help message and exit
```
Required inputs: gene expression data and GTEx sample information. Key outputs: per tissue, a table of sites with interesting editing variation and dsRNAs with interesting editing variation.


Finally, the interesting editing sites per tissue (identified in script 07) are summarized across all tissues, and a few examples are pulled out to plot.

```sh
$ Rscript scripts/08_summarize_variation_search.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz
```
```
Usage: scripts/08_summarize_variation_search.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --samplefile=SAMPLEFILE
                File with GTEx sample information

        --tpmfile=TPMFILE
                Gene expression values (TPM) for the genes of interest

        -h, --help
                Show this help message and exit
```
Required inputs: gene expression data and GTEx sample information. Key outputs: figures summarizing sites or dsRNAs with editing variation across all tissues, plus figures of examples of interest.


## Pilot data examining RNA editing across the whole genome

The pilot data consists of 289 samples where all "A" sites across the genome are assessed for RNA editing, not just the curated list of likely editing sites.


### Filter individual editing quantification files by coverage

To make the individual quantification files a manageable size, filtered individual quant files were generated that only included sites with coverage >10 in that sample.

Here's how to loop through each file name (with the list of file names saved in data/full_pilot/sample_names.txt in this example), from the working directory:

```sh
mkdir scratch/pilot_indiv_quant
PILOT_DIR=path/to/pilot/quant/files
for file in $(cat data/full_pilot/sample_names.txt) ; do echo $file ; Rscript --vanilla scripts/quant_full_pilot_filt.R ${PILOT_DIR}/${file}/*bed $file ; done
```
Required inputs: directory of BED files for each sample. Key outputs: new directory of filtered BED files for each sample.


### Explore filtered pilot data

I combined the filtered BED files for each sample by tissue type, so now each tissue has a file of editing info, with each row being a site_name * individual combination. Columns: chr, start, end, site_name, coverage, editLevel, sample

First, combine all of the pilot data and filter it further. This script removes sites that are only present in one sample, and sites where there isn't a single tissue with evidence of editing (where evidence of editing = >0 edited reads in at least 1 sample).

```sh
$ Rscript scripts/pilot_01_filter_SPRINT_output.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData
```
```
Usage: scripts/pilot_01_filter_SPRINT_output.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --samplefile=SAMPLEFILE
                File with GTEx sample annotation.

        -h, --help
                Show this help message and exit
```
Required inputs: tissue-level combined and coverage-filtered editing quantification BED files. Key outputs: completely combined and completely filtered pilot data (filtered_pilot_data.rds).


Second, examine the information gained through the full genome pilot data over and above the targeted GTEx analysis.

```sh
$ Rscript scripts/pilot_02_check_SPRINT_output.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData
```
```
Usage: scripts/pilot_02_check_SPRINT_output.R [options]

Options:
        --outdir=OUTDIR
                Output directory to write output RData and txt files [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --figdir=FIGDIR
                Figure directory. [default "./figures"]
                
        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "./scratch"]

        --samplefile=SAMPLEFILE
                File with GTEx sample annotation.

        -h, --help
                Show this help message and exit
```
Key output is mostly figures summarizing numbers of informative sites.


## Miscellaneous code

### Filtering GTEx quantification files to just editing sites in putative cis-NATs

This needs a file output from script 05 to run (05_sites_in_cisnat_quantFilePos.txt). The code loops through each sample file name, and filtersthe editing quantification file to just the sites within putative cis-NATs of interest. Since GTEx has around 17000 samples, I split my file listing all GTEx sample names into a bunch of smaller lists (sample_names_A.txt, sample_names_B.txt, etc.), and submitted each of those as a separate job - not the most elegant solution, but it was fine.

Here's an example of if you had 2 files of lists of sample names:

```sh
$ for ind in A B; do 
        echo $ind ; 
        job_name=job_$ind ; 
        file_name=sample_file_names_${ind}.txt ; 
        sbatch --job-name=$job_name scripts/quant_file_filt_wrapper.sh path/to/work_dir $file_name ; 
done
```

This script then launches scripts/quant_file_filt.R, which is the actual workhorse for reading in the sample's editing quantification file, filtering it, and writing out a filtered file for each sample.
