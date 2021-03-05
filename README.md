# cis-nat: Exploration of putative immunogenic double-stranded RNAs arising from cis-NATs #

Our goal is to identify all possible overlapping transcripts (cis-natural antisense transcripts, or cis-NATs) across the genome, then build a comprehensive list of cis-NATs that also have clusters of RNA editing. Since RNA editing of “self” RNA prevents MDA5 binding and the subsequent innate immune signaling, the presence of RNA editing clusters will mark cis-NATs that could be immunogenic dsRNAs. Variation in editing of these cis-NATs could influence autoimmune disease risk; defining these relationships would add to our understanding of complex immune-related diseases.

The project's design document is [here]().

## Walkthrough for code so far:

### Step 1: Get candidate list of all transcript overlaps

```sh
$ Rscript scripts/01_find_tx_overlaps.R --gtf ../../testdir/gencode.v35.annotation.gtf.gz --outdir test_out --scriptdir scripts
```

```
Usage: scripts/01_find_tx_overlaps.R [options]

Options:
        --gtf=GTF
                GTF file [default "NULL"]

        --outdir=OUTDIR
                Output directory to write transcript overlaps to [default is current working directory]

        --scriptdir=SCRIPTDIR
                Scripts directory (to load functions saved in source scripts). [default assumes running from repo, scripts are in "./scripts"]

        --scratchdir=SCRATCHDIR
                Where should intermediate files be saved? If kept as NULL, intermediate files won't be saved. [default "NULL"]

        --keepanno=KEEPANNO
                Should the gene, transcript, and exon annotations be saved? Files will be saved to outdir. "text" for .txt files, "RData" for .RData, "both" for both .txt and .RData files. Anything else will be interpreted as "neither". [default "neither"]

        -h, --help
                Show this help message and exit
```

Key outputs are:
- annotation files for genes, transcript, and exon (if --keepanno is not "neither")
- text file and R object of all possible transcript overlaps (01_putative_cisNAT_uniqueRegionsOnly.txt and 01_putative_cisNAT_uniqueRegionsOnly.RData)
- bed file of genomic regions that cover the transcript overlaps (TBD)
