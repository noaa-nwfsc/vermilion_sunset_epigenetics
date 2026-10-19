#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=read_depth
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o read_depth_%A.txt

#set wd
WORKINGDIR=/scratch2/awray/epigenetic/methylseq/vermilion_genome_output/bismark/deduplicated/

cd $WORKINGDIR

module load bio/samtools

 for i in *.bam; do \
	 echo $i
	 filename=$(basename $i) samtools depth -a $i | awk '{c++;s+=$3}END{print s/c}' >> read_depth.txt 
done



