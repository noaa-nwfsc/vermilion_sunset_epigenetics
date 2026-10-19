# Vermilion Sunset Epigenetic Project

Primary Contacts: Anita Wray,
[anita.wray\@noaa.gov](mailto:anita.wray@noaa.gov)

## Objective:

We leveraged findings from our ongoing research to improve estimates of
age for both vermilion and sunset rockfish. Specifically, we will use
epigenetic analyses of fin clips to estimate the age of individuals in
both species using a subset of our existing 30,000 vermilion and sunset
rockfish samples collected from 2004-2023 during two
industry-collaborative research surveys. Epigenetic ageing provides a
non-lethal method of determining age and lifespan using the degree of
methylation at genomic sites.

### First effort (04/2025):

We used 96 samples from a wide age distribution and good DNA quality
from the fin clip. These samples were both Vermilion and Sunset rockfish
from both the H&L and Trawl surveys. We used the previously published 
genome of vermilion rockfish to identify the location of informative 
methylation markers used zebrafish, and developed primers to target those 
specific regions of the genome. Our results suggest that although our loci 
are somewhat informative of age we do not have enough loci to have a confident 
age call for individuals. In order to increase confidence, we would need to 
locate, amplify, and sequence more CpG sites.
[SEE FIGURE HERE](https://github.com/anita-wray/vermilion_sunset_epigenetics/blob/main/LINEAR%20MODELING/predAge_vs_realAge_glmnet_locusPThresh0.1_MinDep5_4March2025.pdf)
All scripts for this first effort are included in the "First Effort - Mayne Method".

### Current Work (10/2026):

We used whole genome bisulfite sequencing on 104 samples to get a complete look 
at the genomic patterns of methylation. In short, we aligned the sequences to the 
Vermilion Rockfish genome, identified and removed SNPs, then ran through the 
[methylseq nextflow pipeline](https://nf-co.re/methylseq/4.2.0/). Much of this work 
was replicated from [Emma Strand's repo](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2021-04-08-HoloInt-WGBS-Sample-Processing.md). We developed a strong epigenetic clock 
using 95 CpG sites. [SEE FIGURE HERE]()

All scripts for this second effort are included in the "Second Effort - WGBS".

## Disclaimer:

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an 'as is' basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.
