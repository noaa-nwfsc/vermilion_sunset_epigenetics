#!/bin/bash

#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=2 # Number of cores
#SBATCH --job-name=alignmentstats
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o alignment_stats_%A.txt

# Define the output file name
OUTPUT_FILE="bismark_summary.tsv"
WORKING_DIR=/scratch2/awray/epigenetic/methylseq/output/bismark/alignments/logs

cd $WORKING_DIR

# Print the header columns directly into the output file (overwriting any existing file with the same name)
echo -e "File_Name\tMapping_Efficiency\tTotal_Sequence_Pairs\tUnique_Best_Hit\tNo_Alignments" > "$OUTPUT_FILE"

# Loop through all files matching the Bismark report pattern in the current directory
# Adjust the extension (*_bismark_bt2_PE_report.txt, etc.) if your filenames differ
for file in *report.txt; do
    # Ensure it's a file, not a directory
    [ -f "$file" ] || continue

    # 1. Extract File Name (grabs the first instance found in the text)
    file_name=$(grep -m 1 "Bismark report for:" "$file" | sed 's/Bismark report for: //')

    # 2. Extract Mapping Efficiency
    mapping_eff=$(grep "Mapping efficiency:" "$file" | awk -F'\t' '{print $2}' | tr -d '[:space:]')

    # 3. Extract Sequence pairs analysed in total
    total_pairs=$(grep "Sequence pairs analysed in total:" "$file" | awk -F'\t' '{print $2}' | tr -d '[:space:]')

    # 4. Extract Number of paired-end alignments with a unique best hit
    unique_hit=$(grep "Number of paired-end alignments with a unique best hit:" "$file" | awk -F'\t' '{print $2}' | tr -d '[:space:]')

    # 5. Extract Sequence pairs with no alignments under any condition
    no_align=$(grep "Sequence pairs with no alignments under any condition:" "$file" | awk -F'\t' '{print $2}' | tr -d '[:space:]')

    # Append the extracted data for this file as a tab-separated row into the output file
    echo -e "${file_name}\t${mapping_eff}\t${total_pairs}\t${unique_hit}\t${no_align}" >> "$OUTPUT_FILE"
done

echo "Parsing complete! Data has been saved to: $OUTPUT_FILE"


