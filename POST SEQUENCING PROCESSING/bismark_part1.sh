#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --job-name=bismark
#SBATCH --output=bismark_alignment_meth_call_output.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -p medmem


## This script runs 3 Bismark commands for each file/sample within a directory:
##   1. Align (trimmed) sequences to a bisulfite genome using bowtie2
##   2. extracts the methylation call for every single C analysed. The position of every 
##      single C will be written out to a new output file, depending on its context (CpG, 
##      CHG or CHH), whereby methylated Cs will be labelled as forward reads (+), 
##      non-methylated Cs as reverse reads (-).
##   3. Create a HTML containing the alignment and M-bias report

start=`date +%s`

##### ENVIRONMENT SETUP ##########
export baseP=/home/awray/epigenetic
export bamP=${baseP}/03-Bismark
export outP=${baseP}/04-Methylation
export seqP=${baseP}/CleanedData
export refP=/share/nwfsc/refgenomes/vermilion_rockfish/
export cwd=/home/awray/scripts
export tmpd=$cwd/tmp

## Load required tools
module load bio/bismark
module load bio/samtools
module load aligners/bowtie2

#########################################################################################

## 1. Create directories (if not already present)
if [ ! -d "${tmpd}" ]
then
  mkdir -p ${tmpd}
fi


if [ ! -d "${outP}" ]
then
   mkdir -p ${outP}
fi

if [ ! -d "${bamP}" ]
then
   mkdir -p ${bamP}
fi


## 2, Loop through each sample in the 'cleaned data' folder and runs the three commands as
## explained above
for MYFILE in /home/awray/epigenetic/CleanedData'/'*R1_001_val_1.fq; do
 name=`basename --suffix=R1_001_val_1.fq $MYFILE`
 echo "now analyzing $name"
 cd $seqP
# bismark alignment
bismark --bowtie2 -p 4  --output_dir ${bamP}  \
	 --non_directional --temp_dir ${tmpd} --quiet  --genome_folder $refP \
	-1 ${name}R1_001_val_1.fq -2 ${name}R2_001_val_2.fq > ./${name}.bismark.out

cd $bamP

# bismark M bias
bismark_methylation_extractor --multicore 4  --bedGraph --gzip  --output ${outP} --genome_folder $refP \
	${bamP}/${name}R1_001_val_1_bismark_bt2_pe.bam

# generate bismark reports
bismark2report --dir $outP --output ${name}bismark.report.html \
	--alignment_report $bamP/${name}R1_001_val_1_bismark_bt2_PE_report.txt \
	--mbias_report $outP/${name}R1_001_val_1_bismark_bt2_pe.M-bias.txt

done

## 3. Output a run time
end=`date +%s`
runtime=$((end - start ))
echo "runtime: $runtime seconds"

