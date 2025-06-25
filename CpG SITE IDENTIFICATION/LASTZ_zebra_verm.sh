#!/bin/bash
#SBATCH --job-name=LASTZ
#SBATCH -c 1 
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov

##### ENVIRONMENT SETUP ##########
DATADIR=/home/awray/zebrafish_GENOME # The directory containing all of your fastq files
ZEBRAGENOME=/share/nwfsc/refgenomes/zebrafish/GCF_000002035.5_GRCz10_genomic.fna #The name of the zebrafish genome
VERMGENOME=/share/nwfsc/refgenomes/vermilion_rockfish/ver2/GCA_916701275.1_S-miniatus_SEB-74_genomic.fna #The name of the vermilion genome
OUTDIR=/home/awray/zebrafish_GENOME/lastz_alignment # where to store output files

##############################################################################
module load aligners/lastz

## 1.  Make directory for output files
##mkdir $OUTDIR

## 2. Run LASTZ against the single file with zebrafish sequences
lastz $ZEBRAGENOME[multiple] $VERMGENOME[multiple] \
 --step=20 --nogapped\
 --notransition \
 --rdotplot=zebra_vs_verm_aussie_genome_v2.rdp \
 --format=general:score,name1,zstart1,end1,strand1,size1,length1,name2,zstart2,end2,strand2,size2,length1,nmismatch,nmatch,ngap,id%,blastid%,text1,text2,diff,number  > zebra_vs_verm_aussie_method_genome_v2.txt

## 3. Move output files to a set directory
mv zebra_vs_verm.* $OUTDIR

