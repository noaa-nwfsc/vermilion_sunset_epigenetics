#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --job-name=bismarksummary
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH --output=bismark_summary_report_output.txt

#This script will output a summary report of all reads in a directory processed with
#bismark. 


start=`date +%s`

##### ENVIRONMENT SETUP ##########
baseP=/home/awray/epigenetic
outP=${baseP}/03-Bismark
cwd=/home/awray/scripts

## Load required tools
module load bio/bismark
module load bio/samtools
module load aligners/bowtie2
##############################################################################

## 1. Move to directory containing different Bismark alignment, deduplication and methylation 
## extraction (splitting) reports
cd ${outP}


## 2. Run summary generator
bismark2summary --title "all_sample_bismark_summary"

## 3. Output a run time (helpful for future runs)
end=`date +%s`
runtime=$((end - start ))
echo "runtime: $runtime seconds"

