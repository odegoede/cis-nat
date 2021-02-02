# cis-nat: Exploration of putative immunogenic double-stranded RNAs arising from cis-NATs #

## Introduction ##

The innate immune receptor MDA5 binds to double-stranded RNAs (dsRNAs) and initiates immune responses upon their detection, as dsRNAs are often viral genomes. However, the human genome also encodes dsRNAs, either through inverted repetitive sequences or through overlapping transcripts. The mechanism of RNA-editing serves to mark “self” dsRNA from “non-self”; if human dsRNAs are not sufficiently edited, they may trigger inappropriate inflammation.

In the Li lab's RNA-editing GTEx companion paper, 198 possible immunogenic dsRNAs were identified through colocalization of RNA-editing QTLs with loci from 17 immune-related disease GWASes. Despite inverted repeat Alus containing >99% of RNA-editing sites, they only represented ~67% of the putative immunogenic dsRNAs (for the 13/17 immune-related diseases with >=6 colocalized loci). The other ~33% were cis-natural antisense transcripts (cis-NATs), which are overlapping genes transcribed in opposite directions. This enrichment of cis-NATs in immune-colocalized dsRNAs compared to the overall RNA-editing background suggests that they may be “better” ligands for MDA5 and lead to stronger inflammatory responses.

Many cis-NATs arise from the thousands of antisense lncRNA genes in the genome, which are genes that overlap with a protein-coding gene on the opposite strand. Several of these have been identified as trait-associated through our eQTL-GWAS colocalization analysis. It is possible that for some of these colocalized genes, dsRNA formation is the biology behind their trait association. 

Our goal is to identify all possible cis-NATs across the genome, create a comprehensive list of cis-NATs that also have clusters of RNA editing. Since RNA editing of “self” RNA prevents MDA5 binding and the subsequent innate immune signaling, the presence of RNA editing clusters would mark the cis-NATs that could be immunogenic dsRNAs. Variation in editing of these cis-NATs could influence autoimmune disease risk; defining these relationships would add to our understanding of complex immune-related diseases.

## Analysis Steps: ##
1. Identify all transcript overlaps using genome annotation [01_find_tx_overlaps.R](https://github.com/odegoede/cis-nat/blob/main/find_overlapping_transcripts/01_find_tx_overlaps.R)
2. Qin: re-map clusters of hyper-editing within these transcript overlaps - opportunity to improve approach
3. With updated RNA-editing cluster IDs, quantify levels of editing (at cluster-level instead of site-level?)
4. Filter overlaps to candidate dsRNAs that show inter-individual variability in either: RNA-editing level or expression level

## Exploratory/Experimental Code: ##
<ol type="A">
  <li>Explore overlap attributes (follows script 1)</li>
  <li>Assess RNA-editing coverage/levels with current quantification (follows exScript A)</li>
  <li>Select candidate dsRNAs with current quantification: some edQTLs, some not; some variable in editing, some variable in expression, some both (follows exScript B)</li>
  <li>Explore editing/expression patterns of candidate dsRNAs (follows exScript C)</li>
</ol>

 
