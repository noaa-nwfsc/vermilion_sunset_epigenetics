#!/bin/bash
# 
#SBATCH --time=12:00:00
#SBATCH --ntasks=10 
#SBATCH --job-name=bismark_genome_prep
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH --output=bismark_genome_prep_output.txt
#SBATCH -p medmem

## This script will bisulfite convert a genome. This is required so bowtie (within bismark)
## can align bisulfite converted sequences to a genome


start=`date +%s`  
##### ENVIRONMENT SETUP ##########
GENOME_DIR=/share/nwfsc/refgenomes/vermilion_rockfish

## Load required tools
module load bio/bismark
module load aligners/bowtie2
##############################################################################

REPCORES=10

## 1. Bisulfite convert genome using Bismark
bismark_genome_preparation \
--verbose --path_to_aligner /opt/bioinformatics/aligners/bowtie2/bowtie2-2.5.0/ \
$GENOME_DIR

## 2. Output a run time (helpful for future runs)
end=`date +%s`
runtime=$((end - start ))
echo "runtime: $runtime"


