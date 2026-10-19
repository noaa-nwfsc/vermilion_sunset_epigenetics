#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --cpus-per-task=4
#SBATCH --job-name=mice_impute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o mice_impute_%j.txt
#SBATCH --mem=400G
#SBATCH --partition=himem

# load modules needed
module load R/4.4.1

Rscript impute_methylation_perc.R

