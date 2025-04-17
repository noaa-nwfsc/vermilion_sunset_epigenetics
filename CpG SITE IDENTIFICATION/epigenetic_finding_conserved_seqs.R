##Epigenetic Aging Comparison

library(tidyverse)

setwd("~/Desktop/Epigenetic/")
#####
##Read in the LASTZ output
alignment_big <- read.delim(file = "Primer Design/LASTZ Outputs/zebra_vs_verm_aussie_method_genome_v2.txt")

#Read in the file from the Mayne paper which outlines which CpG sites were informative of age in zebrafish
aussie_supl_table <- readxl::read_excel("Supplemental_methods/men13440-sup-0001-tables1-s3_updated.xlsx",
                                        skip = 2)

#Add in Chromosome names so it aligns with the LASTZ output
chrom_names <- readxl::read_excel("DanRer_chromosome_names.xlsx")
chrom_names$aussie <- paste("chr", chrom_names$Chromosome, sep = "")

aussie_supl_table2 <- merge(aussie_supl_table, chrom_names, 
                            by.x = 'chromosome', by.y = 'aussie') %>%
  select(-c('GenBank','Chromosome','Size (bp)','GC content (%)','Action','Unlocalized count'))


aussie_supl_table2 <- as.matrix(aussie_supl_table2)    # faster to work with matrices
alignment_big <- as.matrix(alignment_big)

## Make a handful of blank matrices for the subsequent for loop
df <- matrix(NA,nrow=10000,ncol=6)                     # have all the rows there to start with so you don't have to continually build the matrix
df2 <- matrix(NA,nrow=10000,ncol=22)
df_unique_alignments <- matrix(NA,nrow=10000,ncol=22)
df_unique_CpG <- matrix(NA,nrow=10000,ncol=6)

#Set up the blank matrices to have the column names 
colnames(df_unique_CpG) <- colnames(aussie_supl_table2)
colnames(df_unique_alignments) <- colnames(alignment_big)
colnames(df) <- colnames(aussie_supl_table2)
colnames(df2) <- colnames(alignment_big)

ticker <- 1

#This for loop will cycle through the alignments outputted from LASTZ and:
# - Output only unique alignments that overlap the CpG site from the Mayne paper 
#    (i.e. alignments that aligned only once in both the vermilion and zebrafish genome)
# - Output only unique CpG sites 
#    (i.e. output CpG sites that alligned only once in the vermilion genome)

date() #time checker
for (i in 1:nrow(aussie_supl_table2)) {
  #print(i)
  if(sum(rowSums(cbind(aussie_supl_table2[i,6] == alignment_big[,2] & aussie_supl_table2[i,2] >= alignment_big[,3] & aussie_supl_table2[i,2] <= alignment_big[,4] ))) > 0){     # faster to use matrix operations instead of a loop over a vector
    overlaps <- which(rowSums(cbind(aussie_supl_table2[i,6] == alignment_big[,2] & aussie_supl_table2[i,2] >= alignment_big[,3] & aussie_supl_table2[i,2] <= alignment_big[,4] )) > 0)
    for(j in 1:length(overlaps)){
      df2[ticker,] <- alignment_big[overlaps[j],1:22]
      df[ticker,] <- aussie_supl_table2[i,1:6]
    }
  }
  if(sum(rowSums(cbind(aussie_supl_table2[i,6] == alignment_big[,2] & aussie_supl_table2[i,2] >= alignment_big[,3] & aussie_supl_table2[i,2] <= alignment_big[,4] ))) > 0 & sum(rowSums(cbind(aussie_supl_table2[i,6] == alignment_big[,2] & aussie_supl_table2[i,2] >= alignment_big[,3] & aussie_supl_table2[i,2] <= alignment_big[,4] ))) < 2){     # faster to use matrix operations instead of a loop over a vector
    overlaps <- which(rowSums(cbind(aussie_supl_table2[i,6] == alignment_big[,2] & aussie_supl_table2[i,2] >= alignment_big[,3] & aussie_supl_table2[i,2] <= alignment_big[,4] )) > 0)
    for(j in 1:length(overlaps)){
      df_unique_alignments[ticker,] <- alignment_big[overlaps[j],1:22]
      df_unique_CpG[ticker,] <- aussie_supl_table2[i,1:6]
    }
  }
  ticker <- ticker + 1
}
date()

#Remove all of the NA columns and rows
df <- df[which(is.na(df[,1]) == FALSE),]
df2 <- df2[which(is.na(df2[,1]) == FALSE),]
df_unique_alignments <- df_unique_alignments[which(is.na(df_unique_alignments[,1]) == FALSE),]
df_unique_CpG <- df_unique_CpG[which(is.na(df_unique_CpG[,1]) == FALSE),]

#Output the overall alignments (short vs long just determines what info is outputted)
#as well as the 'unique' files which are the sites we are interested in amplifying through PCR
write.csv(df, file = 'Primer Design/LASTZ Outputs/LASTZ_alignment_short_vermilion_genome_v2.csv')
write.csv(df2, file = 'Primer Design/LASTZ Outputs/LASTZ_alignment_long_vermilion_genome_v2.csv')
write.csv(df_unique_alignments, file = "Primer Design/LASTZ Outputs/LASTZ_alignment_uniquegenome_v2.csv")
write.csv(df_unique_CpG, file = "Primer Design/LASTZ Outputs/LASTZ_alignment_unique_CpGgenome_v2.csv")