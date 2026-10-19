#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --job-name=intersect_bed2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o filt_samps_%A.txt

# load modules needed

cd /scratch2/awray/epigenetic/post_nextflow_analyses/honeycomb_alignment

module load bio/bedtools/2.31.1

for i in *5x_sorted.tab
do
  intersectBed \
  -a ${i} \
  -b CpG.90filt.all.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment_90filt.bed
done
