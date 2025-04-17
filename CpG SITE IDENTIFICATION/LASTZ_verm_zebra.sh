#!/bin/bash
#SBATCH --job-name=LASTZ
#SBATCH -p himem
#SBATCH -c 24
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov

##### ENVIRONMENT SETUP ##########
DATADIR=/home/awray/zebrafish_GENOME # The directory containing all of your fastq files
ZEBRAGENOME=/share/nwfsc/refgenomes/zebrafish/GCF_000002035.5_GRCz10_genomic.fna #The name of the zebrafish genome
VERMGENOME= /share/nwfsc/refgenomes/vermilion_rockfish/GCA_916701115.1.fasta #The name of the vermilion genome
OUTDIR=/home/awray/zebrafish_GENOME/lastz_alignment # where to store output files

## Load required tools
module load aligners/lastz
##############################################################################


## 1.  Make directory for output files
mkdir $OUTDIR

## 2. Run LASTZ against the single file with zebrafish sequences
lastz $ZEBRAGENOME $VERMGENOME \
 --format=general --ambiguous=iupac \
 --step=100 --nogapped \
 --notransition \
 --rdotplot= verm_vs_zebra.rdp \
 --format=maf > verm_vs_zebra.maf

## 3. Move output files to a set directory
mv verm_vs_zebra.* $OUTDIR