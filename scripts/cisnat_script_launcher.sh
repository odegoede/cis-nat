#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=24:00:00
#SBATCH --mem=8G

# to run from a the command line with a log file: scripts/cisnat_script_launcher.sh path/to/work_dir > cisnat_log 2>&1

WORK_DIR=$1

cd $WORK_DIR
module load r/3.6

echo "Starting Script 01" && date &&
Rscript scripts/01_find_tx_overlaps.R --gtf data/gencode.v35.annotation.gtf.gz --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --keepanno both &&
echo "Starting Script 02" && date &&
Rscript scripts/02_save_overlap_locations.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ &&
echo "Starting Script 03" && date &&
Rscript scripts/03_flag_tx_overlaps_alu.R --alufile data/hg38_repeats.Alu.bed.gz --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ &&
echo "Starting Script 04" && date &&
Rscript scripts/04_examine_tx_overlaps.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --edittable data/gtex_editing_paper_tableS4.xlsx --coloctable data/gtex_lncrna_paper_hits_summary_table_allGeneTypes.RData &&
echo "Starting Script 05" && date &&
Rscript scripts/05_check_editing_sites.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --editanno data/gtex_edit/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.AnnoVar_Genes.Clusered.txt.gz --examplequant data/gtex_edit/indiv_files/GTEX-11ZTS-1026-SM-5LU8O.txt.gz &&
echo "Starting Script 06" && date &&
Rscript scripts/06_check_editing_quantification.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --editcoverage data/gtex_edit/coverage_matrix.txt.gz --editlevel data/gtex_edit/editLevel_matrix.txt.gz --tpmfile data/cisnat_gene_tpm.gct.gz --gtexanno data/rnaseqc_genes_info.RData &&
echo "Starting Script 07" && date &&
Rscript scripts/07_example_variation_search.R --tissue $tis --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz &&
echo "Starting Script 08" && date &&
Rscript scripts/08_summarize_variation_search.R --outdir output/ --scriptdir scripts/ --scratchdir scratch/ --figdir figures/ --samplefile data/sampleInfo.RData --tpmfile data/cisnat_gene_tpm.gct.gz &&
echo "Scripts done" && date
