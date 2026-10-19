#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --cpus-per-task=20
#SBATCH --job-name=model_epi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anita.wray@noaa.gov
#SBATCH -o model_%j.txt
#######--partition=medmem

# load modules needed
module load R/4.6.1

Rscript epigenetic_model_sedna.R
