#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=5 # Number of cores
#SBATCH --job-name=trimgalore
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH --output=trimgalore_output.txt

source ~/.bashrc

start=`date +%s`


##### ENVIRONMENT SETUP ##########
baseP=/home/awray/epigenetic
seqP=${baseP}/RAWDATA
outP=${baseP}/CleanedData
SUFFIX1=R1_001.fastq
SUFFIX2=R2_001.fastq

## Load required tools
module load bio/samtools
module load aligners/bowtie2
mamba activate trim-galore-0.6.10

##############################################################################

## 1. Create directories (if not already present)
if [ ! -d "${outP}" ]; then
   mkdir ${outP}
fi

## 2. Move to directory containing raw fastq files (direct from MiSeq)
cd $seqP



## 3. Loop through all of the samples (read in by reading all forward read fastq files)
## run trim_galore which will remove illumina adapters, trim for quality, and clip 15bp 
## from each side of both F and R reads. This will eliminate the non-bisulfite converted
## DNA and improve genome alignment with Bismark

for MYFILE in /home/awray/epigenetic/RAWDATA'/'*$SUFFIX1; do
 MYBASE=`basename --suffix=R1_001.fastq $MYFILE`
 #echo ${MYBASE}$SUFFIX1
 trim_galore -j ${SLURM_NTASKS} --paired --phred33 --illumina \
       --clip_R1 15 --clip_R2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15 \
       --output_dir $outP ${MYBASE}$SUFFIX1 ${MYBASE}$SUFFIX2
done


## 4. Output a run time (helpful for future runs)
end=`date +%s`
runtime=$((end - start ))
echo "runtime: $runtime"

mamba deactivate 
