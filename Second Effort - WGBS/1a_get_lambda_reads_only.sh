#!/bin/bash

#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=2 # Number of cores
#SBATCH --job-name=get_lambda
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o get_lambda_%A.txt

cd  /scratch2/awray/epigenetic/methylseq/output/bismark/methylation_calls/methylation_coverage/

# Define the output file
OUTPUT_FILE="combined_phage_lambda_data.txt"

# Clear the output file if it already exists
> "$OUTPUT_FILE"

# Loop over every .csv file in the current directory
for file in *report.txt; do
    
    # Check if files exist to avoid running on a literal '*.csv' string
    [ -e "$file" ] || continue
    
    # Extract the individual's name by removing the .csv extension
    individual_name="${file%.txt}"
    
    # Use awk to process the file:
    # -F, sets the delimiter to a comma (change to -F'\t' if your files are tab-separated)
    # It checks if the 1st column ($1) is exactly 'phage_lambda'
    # If it matches, it prints the individual name followed by the rest of the line
    awk  -v ind="$individual_name" '$1 == "phage_lambda" {print ind "\t" $0}' "$file" >> "$OUTPUT_FILE"

done

echo "Done! Combined data saved to $OUTPUT_FILE"
