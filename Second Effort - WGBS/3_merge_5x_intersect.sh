#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=intersect_bed
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o merge_individuals_%A.txt

#set wd
WORKINGDIR=/scratch2/awray/epigenetic/post_nextflow_analyses/honeycomb_alignment

cd $WORKINGDIR

# load modules needed

module load bio/bedtools/2.31.1

mkdir removed_samples

cat low_reads_low_conversion_samples.txt | while IFS= read -r file; do mv "$file" removed_samples/; done

echo "Removed low quality individuals"

multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

echo "Combined individuals into one file"

awk '$4 == 81' CpG.all.samps.5x_sorted.bed > CpG.95filt.all.samps.5x_sorted.bed
awk '$4 == 77' CpG.all.samps.5x_sorted.bed > CpG.90filt.all.samps.5x_sorted.bed
awk '$4 == 85' CpG.all.samps.5x_sorted.bed > CpG.100filt.all.samps.5x_sorted.bed
