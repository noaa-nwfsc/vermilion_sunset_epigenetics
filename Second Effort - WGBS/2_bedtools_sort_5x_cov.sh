#!/bin/bash

#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --job-name=5xcovfilt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov


# load modules needed

module load bio/bedtools
WORKINGDIRECTORY=/scratch2/awray/epigenetic/methylseq/output/bismark/methylation_calls/methylation_coverage/

cd $WORKINGDIRECTORY

### Sorting the .cov files so all lines are in the same order
for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  echo $STEM
  bedtools sort -i "${f}" \
  > "${STEM}"_sorted.cov
done

echo "finished bedtools sort"

### Filtering for CpG for 5x coverage. To change the coverage, replace X with your desired coverage in ($5+6 >= X)

for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  echo $STEM
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x_sorted.tab
done

echo "finished 5x filtering"

mv *merged_CpG_evidence.cov /scratch2/awray/epigenetic/post_nextflow_analyses/honeycomb_alignment
mv *_5x_sorted.tab /scratch2/awray/epigenetic/post_nextflow_analyses/honeycomb_alignment
mv *_sorted.cov /scratch2/awray/epigenetic/post_nextflow_analyses/honeycomb_alignment
