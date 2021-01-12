MAAPER: Model-based analysis of alternative polyadenylation using 3' end-linked reads
================
Wei Vivian Li, Bin Tian
2021-01-12

<!-- README.md is generated from README.Rmd. Please edit that file -->
<img src="https://github.com/Vivianstats/data-pkg/raw/main/img/MAAPER.png" height="200" align="right" />

## Latest News

> 2020/01/02:

-   Version 1.0.0 released!

## Introduction

MAAPER is a computational method for model-based analysis of alternative polyadenylation using 3' end-linked reads. It uses a probabilistic model to predict polydenylation sites (PASs) for nearSite reads with high accuracy and sensitivity, and examines different types of alternative polyadenylation (APA) events, including those in 3'UTRs and introns, using carefully designed statistics.

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/MAAPER/issues). For suggestions and comments on the method, please contact Vivian (<vivian.li@rutgers.edu>).

## Installation

You can install `MAAPER` from github with:

``` r
# install.packages("devtools")
devtools::install_github("Vivianstats/MAAPER")
```

## Quick start

`maaper` requires three input files:

-   The GTF file of the reference genome;
-   The BAM files of the 3' sequencing data (nearSite reads). The BAM file should be sorted and the index BAI file should be present in the same directory as the BAM file;
-   The PAS annotation file whose version matches the reference genome. We have prepared [PolyA\_DB](https://exon.apps.wistar.org/PolyA_DB/v3/) annotation files for MAAPER, and they can be downloaded from [this page](https://github.com/Vivianstats/data-pkg/tree/main/MAAPER/PolyA_DB).

The final output of `mapper` are two text files named "gene.txt" and "pas.txt", which contain the predicted PASs and APA results.

Below is a basic example which shows how to use the `maaper` function. The bam and gtf files used in this example can be downloaded [here](https://github.com/Vivianstats/data-pkg/tree/main/MAAPER). To save computation time, we are providing a toy example dataset of chr19. In real data application, we do not recommend dividing the files into subsets by chromosomes.

``` r
library(MAAPER)

pas_annotation = readRDS("./mouse.PAS.mm9.rds")
gtf = "./gencode.mm9.chr19.gtf"
# bam file of condition 1 (could be a vector if there are multiple samples)
bam_c1 = "./NT_chr19_example.bam"
# bam file of condition 2 (could be a vector if there are multiple samples)
bam_c2 = "./AS_4h_chr19_example.bam"

maaper(gtf, # full path of the GTF file
       pas_annotation, # PAS annotation
       output_dir = "./", # output directory
       bam_c1, bam_c2, # full path of the BAM files
       read_len = 76, # read length
       ncores = 12,  # number of cores used for parallel computation 
      )
```

Please refer to the package [manual](https://github.com/Vivianstats/MAAPER/blob/master/inst/docs/) for a full list of arguments and detailed usage.
