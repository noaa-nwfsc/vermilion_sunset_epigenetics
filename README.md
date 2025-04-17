# Vermilion Sunset Epigenetic Project

Primary Contacts: Anita Wray,
[anita.wray\@noaa.gov](mailto:anita.wray@noaa.gov){.email}, Marty
Kardos, [martin.kardos\@noaa.gov](mailto:martin.kardos@noaa.gov){.email}

## Objective:

We leveraged findings from our ongoing research to improve estimates of
age for both vermilion and sunset rockfish. Specifically, we will use
epigenetic analyses of fin clips to estimate the age of individuals in
both species using a subset of our existing 30,000 vermilion and sunset
rockfish samples collected from 2004-2023 during two
industry-collaborative research surveys. Epigenetic ageing provides a
non-lethal method of determining age and lifespan using the degree of
methylation at genomic sites. These sites are preselected from previous
studies.

### First effort (last update 04/2025):

We used 96 samples from a wide age distribution and good DNA quality
from the fin clip. These samples were both Vermilion and Sunset rockfish
from both the H&L and Trawl surveys.

## Methods:

#### 1. [Find CpG sites in Vermilion Rockfish](https://github.com/noaa-nwfsc/vermilion_sunset_epigenetics/tree/main/CpG%20SITE%20IDENTIFICATION)

Used the previously published genome of vermilion rockfish to identify
the location of informative methylation markers used in European Bass,
Australian lungfish, zebrafish, Murray cod, and Mary River cod
(Anastasiadi & Piferrer, 2019, Mayne et al. 2021, Mayne et al. 2022). We
did this with two methods: 1. Use LASTZ to align the vermilion and
zebrafish genomes then filtered to only unique sites (see bash and R
script) 2. Subsetted the zebrafish genome to only 300bp surrounding the
CpG sites (600bp total, fasta file included in repo) then used NCBI
BLAST against the vermilion genome

#### 2. Develop primers

Developed primers to target those specific regions of the genome. We
used [methprimer](https://www.methprimer.com/) to develop primers based
on the LASTZ and BLAST outputs. This resulted in the .csv file which
gave us primers, annealing temps, and amplicon lengths.

#### 3. Lab Work

Followed protocols outlined in Mayne et al. 2021, Mayne et al. 2022, and
Anastasiadi & Piferrer, 2019 to extract, bisulfite treat the DNA and
amplify, barcode and sequence specific regions in the genome. (NOT
INCLUDED IN THIS REPO)

#### 4. [Post Sequencing Processing](<https://github.com/noaa-nwfsc/vermilion_sunset_epigenetics/tree/main/POST%20SEQUENCING%20PROCESSING>)
Demultiplex, align, and call methylation percentage for each locus
amplified. Most of the scripts are taken from [this
site](https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/WGBS/WGBS),
which was super helpful. In order, you trim, then align, then call the
methylation percentage. The file that we end up using for step 5 is the
.cov file from each individual, produced from the bismark_part1.sh
script. \#### 5. Linear modeling (MARTY ADD HERE)

## Current Results (04/2025):

Our results suggest that although our loci are somewhat informative of
age we do not have enough loci to have a confident age call for
individuals. In order to increase confidence, we would need to locate,
amplify, and sequence more CpG sites.

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
