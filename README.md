# cis-nat: Exploration of putative immunogenic double-stranded RNAs arising from cis-NATs #

Our goal is to identify all possible overlapping transcripts (cis-natural antisense transcripts, or cis-NATs) across the genome, then build a comprehensive list of cis-NATs that also have clusters of RNA editing. 

Since RNA editing of “self” RNA prevents MDA5 binding and the subsequent innate immune signaling, the presence of RNA editing clusters will mark cis-NATs that could be immunogenic dsRNAs. Variation in editing of these cis-NATs could influence autoimmune disease risk; defining these relationships would add to our understanding of complex immune-related diseases.

The project's design document (includes dependencies) is [here](https://github.com/odegoede/cis-nat/blob/main/design.md).

## Suggested directory set-up
```
.
├── data                     # Preliminary data files (annotation files, sample info, quant files)
├── figures                  # To write figures to
│   └── exploratory            # (subdir for figures from exploratory scripts)
├── output                   # To write output files to
├── scratch                  # Temporary directory for intermediate files
└── scripts                  # Scripts
```

## Walkthrough for code so far:

### Part 1: Identify all possible transcript overlaps using genome annotation

First, find overlap locations.

```sh
$ Rscript scripts/01_find_tx_overlaps.R --gtf data/gencode.v35.annotation.gtf.gz --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --keepanno both
```
```
Usage: scripts/01_find_tx_overlaps.R [options]

Options:
        --gtf=GTF
                GTF file [default "NULL"]

        --outdir=OUTDIR
                Output directory to write transcript overlaps to [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "NULL"]

        --keepanno=KEEPANNO
                Should the gene, transcript, and exon annotations be saved? Files will be saved to outdir. "text" for .txt files, "RData" for .RData, "both" for both .txt and .RData files. Anything else will be interpreted as "neither". [default "neither"]

        -h, --help
                Show this help message and exit
```
Required input is a GTF file.

Key outputs are:
- annotation files for genes, transcript, and exon (if --keepanno is not "neither")
- text file and R object of all possible transcript overlaps (01_putative_cisNAT_uniqueRegionsOnly.txt and 01_putative_cisNAT_uniqueRegionsOnly.RData)
- bed file of genomic regions that cover the transcript overlaps (TBD)


Second, save overlap locations as BED file.

```sh
$ Rscript scripts/02_save_overlap_locations.R --overlap output/01_putative_cisNAT_uniqueRegionsOnly.RData --outdir output/
```
```
Usage: scripts/02_save_overlap_locations.R [options]

Options:
        --overlap=OVERLAP
                RData file of transcript overlaps [default "NULL"]

        --outdir=OUTDIR
                Output directory to write transcript overlaps to [default is "./output"]

        -h, --help
                Show this help message and exit
```


### Part 2: Flag overlaps that also contain IRAlus, and examine characteristics of the transcript overlaps

First, flag transcript overlaps that contain inverted repeat Alus (since the IRAlu could be the actual driver of a dsRNA, so these overlaps would need to be examined in-depth).

```sh
$ Rscript scripts/03_flag_tx_overlaps_alu.R --overlap output/01_putative_cisNAT_uniqueRegionsOnly.RData --alufile data/hg38_repeats.Alu.bed.gz --outdir output/ --scriptdir scripts/
```
```
Usage: scripts/03_flag_tx_overlaps_alu.R [options]

Options:
        --overlap=OVERLAP
                RData file of transcript overlaps [default "NULL"]

        --alufile=ALUFILE
                File with genomic annotation of Alu regions [default "NULL"]

        --outdir=OUTDIR
                Output directory to write transcript overlaps to [default is "./output"]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default is "./scripts"]

        -h, --help
                Show this help message and exit
```

Second, examine overlap characteristics (like length of overlap, which gene types are involved, etc.).

```sh
$ Rscript scripts/04_examine_tx_overlaps.R --overlap output/03_putative_cisNAT_uniqueRegionsOnly_withAluFlag.RData --annofile output/gene_tx_exon_anno_files.RData --outdir output/ --figdir figures/ --edittable data/gtex_editing_paper_tableS4.xlsx --coloctable data/gtex_lncrna_paper_hits_summary_table_allGeneTypes.RData
```
```
Usage: scripts/04_examine_tx_overlaps.R [options]

Options:
        --overlap=OVERLAP
                RData file of transcript overlaps [default "NULL"]

        --annofile=ANNOFILE
                File with genomic annotation of gene/tx/exon annotation [default "NULL"]

        --outdir=OUTDIR
                Output directory to write transcript overlaps to [default "./output"]

        --figdir=FIGDIR
                Figure directory to put figures in [default "./figures"]

        --edittable=EDITTABLE
                File of cis-NATs from RNA-editing paper (optional) [default "NULL"]

        --coloctable=COLOCTABLE
                File of QTL colocalization from GTEx lncRNA paper (optional) [default "NULL"]

        -h, --help
                Show this help message and exit
```


### Part 3: Combine the individual RNA-editing quantification files into workable summary files.

TBC

