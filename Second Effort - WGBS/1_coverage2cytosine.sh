#!/bin/bash

#SBATCH --time=3-00
#SBATCH --nodes=1
#SBATCH --ntasks=10 # Number of cores
#SBATCH -p medmem
#SBATCH --job-name=cov_to_cyto
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o cov_to_cyto__%A.txt

# load modules needed
module load bio/bismark


# run coverage2cytosine merge of strands
OUTDIR=/scratch2/awray/epigenetic/methylseq/output/bismark
GENOMEFOLDER=/share/nwfsc/refgenomes/honeycomb_rockfish/bisulfite_genome

cd $OUTDIR/methylation_calls/methylation_coverage/

 find $OUTDIR/methylation_calls/methylation_coverage/*deduplicated.bismark.cov.gz \
 | xargs basename -s _1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
 | xargs -I{} coverage2cytosine \
 --genome_folder $GENOMEFOLDER \
 -o {} \
 --merge_CpG \
 --zero_based \
$OUTDIR/methylation_calls/methylation_coverage/{}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz

