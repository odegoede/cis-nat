#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# usage: scripts/quant_file_filt_wrapper.sh path/to/work_dir path/to/sample_names_A.txt

WORK_DIR=$1
SAMPLE_NAMES_FILE=$2

cd $WORK_DIR
module load r/3.6

for nom in $(cat $SAMPLE_NAMES_FILE) ; do 
	file="data/indiv_quant_files/${nom}" ;
	Rscript --vanilla scripts/quant_file_filt.R $file data/05_sites_in_cisnat_quantFilePos.txt;
	gzip scratch/filt_indiv_quant/${nom%.gz}
done
