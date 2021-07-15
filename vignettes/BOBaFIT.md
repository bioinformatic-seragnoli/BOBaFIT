---
title: "BOBaFIT"
author: "Vignette Author"
date: "2021-05-06"
output: BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---

```{=html}
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
```
```{=html}
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>
```
# Overview

The R package BoBafit is composed of three functions which allow the refit and the recalibration of copy number profile of tumor sample. The principal and refitting function was named `DRrefit`, which - throughout a chromosome clustering method and a list of unaltered chromosomes (chromosome list) - recalibrates the copy number values. BoBafit also contains two secondary functions, the `ComputeNormalChromosome`, which generates the chromosome list, and the `PlotCluster`.

# Data

The package checks the diploid region assessment working on pre-estimated segment information, as the copy number and their position. We included a data set `segments` where are showed all the information necessary. The data correspond to segments about 100 breast tumors samples obtained by the project TCGA-BRCA [@Tomczak2015].

| chr |    start |      end | arm |        CN | ID       |    width |
|----:|---------:|---------:|:----|----------:|:---------|---------:|
|   1 |    62920 | 15823420 | 1p  | 2.0676615 | SAMPLE_1 | 15760501 |
|   1 | 15827002 | 15827430 | 1p  | 0.1767767 | SAMPLE_1 |      429 |
|   1 | 15827706 | 16542868 | 1p  | 2.0567460 | SAMPLE_1 |   715163 |
|   1 | 16544783 | 16617312 | 1p  | 0.6907243 | SAMPLE_1 |    72530 |
|   1 | 16617327 | 16864367 | 1p  | 1.5530356 | SAMPLE_1 |   247041 |
|   1 | 16868660 | 16898730 | 1p  | 0.5849412 | SAMPLE_1 |    30071 |

# BOBaFIT workflow

To start the analysis of the diploid region made by BOBaFIT, .tsv files , containing the sample copy number profile, and the chromosome list are needed as inputs. The chromosome list is the list of chromosomes which are the least affected by SCNAs in the tumor and can be manually created or by using the function *`ComputeNormalChromosome`.* We suggest these two sequential steps to allow the right refit and recalibration of sample's diploid region:

1.  `ComputeNormalChromosome()`

2.  `DRrefit()`

Here we performed this analysis workflow on the dataset `segments` described above.

## ComputeNormalChromosome

The chromosome list is a vector that contains the chromosomal arms considered "normal" in the cohort of samples tested. Chromosomes included in the list should be selected when their CN values are subject to minimal fluctuation and tend to remain within the diploid range in the analyzed tumor. *ComputeNormalChromosome* allows to set the chromosomal alteration rate (`perc`). For a more stringent analysis, we suggest an alteration rate of 5% (0.5) ; on the contrary, for a more permissive analysis, we suggest as maximum rate 20-25% (0.20-0.25) . The function input is a sample cohort with their segments.

Here we performed the function in the data set `segments`, using an alteration rate of 15%.

``` {.r}
library(BOBaFIT)
chr_list <- computeNormalChromosomes(segments = segments, perc = 0.15)
```

``` {.r}
chr_list
```

[1] "10q" "11p" "12p" "12q" "13q" "14q" "15q" "18p" "18q" "19p" "19q" "1p" [13] "21q" "22q" "2p" "2q" "3p" "4p" "4q" "6p" "6q" "7p"

Storing the result in the variable `chr_list`, it will be a vector containing the chromosomal arms which present an alteration rate under the indicated `perc` value.

### Figure

Following, the function shows in the Viewer a plot where is possible observe the chromosomal alteration rate of each chromosomal arms and which one have been selected in the chromosome list (Blue bars).

![](figure/Chromosome_list.png "Chromosome_list")

## DRrefit

To create a tumor-specific method that refit and recalibrate the tumor copy number profile, we developed the function `DRrefit`. It uses as input the sample's segments (.tsv) - cohort or single sample-, and the chromosome list. More in details, this function to estimate the right diploid region uses the clustering function `NbClust` [@charrad2014] and the method of clustering can be sets by the user with the argument `clust_method`. The option of this argument are: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid" and "kmeans".

In this example, the `segment` data table and the `chr_list` previously generated are used. The default value of `maxCN` (6) and `clust_method` (ward.d2) are used.

``` {.r}
library(BOBaFIT)

Results <- DRrefit (segments_chort = segments, 
        chrlist = chr_list, 
        plot_output = TRUE, 
        plot_path = path)
```

Setting the argument `plot_output` TRUE, the function will automatically save the plot in the folder indicated by the variable `path`. Therefore, the `plot_path` argument must be given. Instead the `plot_output` is set as FALSE, the `plot_path` can not be indicated. We suggest to store the output of `DRrefit` in a variable (in this example we use `Results`) to view and possibly save the two data frames generated. The plot and the data frames are explained in the next sections.

### Tables

The two output data framesare :

-   The data frame `segments corrected`, containing the segment CN adjusted by the correction factor (CR) - value estimated by the function to correct the diploid region- of all samples analyzed.

``` {.r}
Results$segments_corrected
```

| chr |    start |      end | arm |        CN | ID       |    width | CN_corrected |
|----:|---------:|---------:|:----|----------:|:---------|---------:|-------------:|
|   1 |    62920 | 15823420 | 1p  | 2.0676615 | SAMPLE_1 | 15760500 |    2.0155736 |
|   1 | 15827002 | 15827430 | 1p  | 0.1767767 | SAMPLE_1 |      428 |    0.1246888 |
|   1 | 15827706 | 16542868 | 1p  | 2.0567460 | SAMPLE_1 |   715162 |    2.0046581 |
|   1 | 16544783 | 16617312 | 1p  | 0.6907243 | SAMPLE_1 |    72529 |    0.6386364 |
|   1 | 16617327 | 16864367 | 1p  | 1.5530356 | SAMPLE_1 |   247040 |    1.5009477 |
|   1 | 16868660 | 16898730 | 1p  | 0.5849412 | SAMPLE_1 |    30070 |    0.5328534 |

It is similar to the input one and report the new CN value of each segment calculated by DRrefit (`CN_corrected`).

-   The data frame `report`, which contains all the information about the clustering as the outcome,the number of clusters found in that sample, the chromosome list and the CR used for the adjstment of the diploid region.

``` {.r}
Results$report
```

+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| sample   | clustering | ref_clust_chr                                                                                                                                    | num_clust | correction_factor |
+:=========+:===========+:=================================================================================================================================================+==========:+==================:+
| SAMPLE_1 | SUCCEDED   | 10p, 10q, 11p, 11q, 12p, 12q, 13q, 14q, 15q, 17p, 17q, 18p, 18q, 19p, 19q, 20p, 20q, 21q, 2q, 3p, 3q, 4p, 5p, 5q, 6p, 6q, 7p, 7q, 8p, 8q, 9p, 9q | 5         | -0.0520879        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| SAMPLE_2 | SUCCEDED   | 10p, 11p, 11q, 12q, 13q, 14q, 15q, 16p, 17q, 18p, 18q, 1q, 20q, 21q, 3p, 5q, 6p, 6q, 8q                                                          | 5         | -0.0490066        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| SAMPLE_3 | SUCCEDED   | 12q, 15q, 16p, 16q, 17p, 17q, 18p, 18q, 19p, 19q, 1p, 1q, 20q, 22q, 2p, 2q, 3p, 3q, 6p, 6q, 7p, 7q, 9q                                           | 3         | -0.4587184        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| SAMPLE_4 | SUCCEDED   | 10p, 10q, 11p, 11q, 12p, 12q, 13q, 14q, 15q, 16q, 17q, 18p, 18q, 1q, 20p, 20q, 21q, 2p, 2q, 3p, 3q, 4p, 4q, 5p, 5q, 6p, 6q, 7p, 7q, 8p, 8q       | 3         | -0.0196744        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| SAMPLE_5 | SUCCEDED   | 10q, 11p, 11q, 12q, 13q, 14q, 18p, 18q, 1q, 21q, 2p, 2q, 3p, 3q, 4p, 4q, 5p, 5q, 6p, 6q, 7p, 7q, 8q, 9p                                          | 2         | -0.0206308        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+
| SAMPLE_6 | SUCCEDED   | 11q, 12q, 16q, 17p, 17q, 18p, 18q, 1p, 1q, 21q, 3q, 5p, 5q, 6p, 6q, 7p                                                                           | 6         | -0.0455488        |
+----------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------+-----------+-------------------+

When the column `clustering` reports FAIL, it indicates that , NbClust fails the chromosome clustering in that sample. In this case, the sample will not present clusters, so the input chromosome list will be kept as reference . When the column `clustering` reports SUCCED, NbClust succeeds and and the new chromosome list is chosen. The chromosome list used for each sample are all reported in the column `ref_clust_chr`.

### Figure

Another output of DRrefit is a plot for each sample with the old and new segments positions, where we can appreciate the new copy number profile. The x-axes represent the chromosomes with their genomic position, and the y-axes the copy number value. The dashed line shows the diploid region. Above the plot are reported the sample name, the CR and the chromosomal arm excluded by the new chromosome list.

Based on the CR value two plots can be displayed:

-   CR ≤ 0.1: the new segment and the old segments are orange and red colored, respectively;

![](figure/SAMPLE_1_ward.D2.png "CR < 0.1"){width="704" height="222"}

-   CR \> 0.1: the new segment and the old segments are green and red colored, respectively;

![](figure/SAMPLE_13_ward.D2.png "CR > 0.1"){width="904"}

# PlotCluster

The accessory function is `PlotCluster`. It can be used to visualize the chromosomal cluster in a single sample. It can cluster either a sample cohort or a single sample and the input data is always a .tsv file as the data frame `segments`. The option of `clust_method` argument are the same of `DRrefit`("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid" and "kmeans").

``` {.r}
library(BOBaFIT)
Cluster <- PlotCluster(segs = segments,
                       clust_method = "ward.D2", 
                       plot_path = path)
```

We suggest to store the output on a variable (in this example we use `Cluster`) to view and possibly save the data frame generated. The `PlotCuster` will automatically save the plot in the folder indicated by the variable `path` of the argument `plot_path`. All the outputs are explained and shown in the next section.

### Figure

In the `PlotCluster` plot, the chromosomal arms are labeled and colored according to the cluster they belong to. The y-axis reports the arm CN.

![](figure/SAMPLE_100_PlotCluster.png)

### Tables

The output `report` is a data frame which reports for each sample the clustering outcome (fail or succeeded) and the number of clusters for each sample analyzed.

``` {.r}
Cluster$report
```

| sample   | clustering | num_clust |
|:---------|:-----------|----------:|
| SAMPLE_1 | SUCCEDED   |        NA |
| SAMPLE_1 | SUCCEDED   |         5 |
| SAMPLE_2 | SUCCEDED   |         5 |
| SAMPLE_1 | SUCCEDED   |         5 |

The second output that will be stored in the Cluster variable is a list of dataframes, one for each sample analyzed, where are reported the information about the plot.

``` {.r}
Cluster$plot_tables$SAMPLE_1
```

| chr | cluster  |       CN |
|:----|:---------|---------:|
| 1p  | cluster4 | 1.821790 |
| 1q  | cluster5 | 4.702378 |
| 2p  | cluster4 | 1.844377 |
| 2q  | cluster1 | 2.054458 |
| 3p  | cluster1 | 2.051508 |
| 3q  | cluster1 | 2.045048 |

# Session info {.unnumbered}

    ## \begin{itemize}\raggedright
    ##   \item R version 4.0.5 (2021-03-31), \verb|x86_64-apple-darwin17.0|
    ##   \item Locale: \verb|C/it_IT.UTF-8/it_IT.UTF-8/C/it_IT.UTF-8/it_IT.UTF-8|
    ##   \item Running under: \verb|macOS Catalina 10.15.7|
    ##   \item Matrix products: default
    ##   \item BLAS:   \verb|/Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib|
    ##   \item LAPACK: \verb|/Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib|
    ##   \item Base packages: base, datasets, grDevices, graphics, methods,
    ##     stats, utils
    ##   \item Other packages: BOBaFIT~0.1.0, BiocStyle~2.18.1
    ##   \item Loaded via a namespace (and not attached):
    ##     AnnotationDbi~1.52.0, AnnotationFilter~1.14.0, BSgenome~1.58.0,
    ##     Biobase~2.50.0, BiocFileCache~1.14.0, BiocGenerics~0.36.1,
    ##     BiocManager~1.30.12, BiocParallel~1.24.1, Biostrings~2.58.0,
    ##     DBI~1.1.1, DelayedArray~0.16.3, Formula~1.2-4, GGally~2.1.1,
    ##     GenomeInfoDb~1.26.7, GenomeInfoDbData~1.2.4,
    ##     GenomicAlignments~1.26.0, GenomicFeatures~1.42.3,
    ##     GenomicRanges~1.42.0, Hmisc~4.5-0, IRanges~2.24.1, MASS~7.3-53.1,
    ##     Matrix~1.3-2, MatrixGenerics~1.2.1, NbClust~3.0,
    ##     OrganismDbi~1.32.0, ProtGenerics~1.22.0, R6~2.5.0, RBGL~1.66.0,
    ##     RColorBrewer~1.1-2, RCurl~1.98-1.3, RSQLite~2.2.7, Rcpp~1.0.6,
    ##     Rsamtools~2.6.0, S4Vectors~0.28.1, SummarizedExperiment~1.20.0,
    ##     VariantAnnotation~1.36.0, XML~3.99-0.6, XVector~0.30.0,
    ##     askpass~1.1, assertthat~0.2.1, backports~1.2.1, base64enc~0.1-3,
    ##     biomaRt~2.46.3, biovizBase~1.38.0, bit~4.0.4, bit64~4.0.5,
    ##     bitops~1.0-7, blob~1.2.1, cachem~1.0.4, callr~3.7.0,
    ##     checkmate~2.0.0, cli~2.5.0, cluster~2.1.2, colorspace~2.0-0,
    ##     compiler~4.0.5, crayon~1.4.1, curl~4.3.1, data.table~1.14.0,
    ##     dbplyr~2.1.1, desc~1.3.0, devtools~2.4.0, dichromat~2.0-0,
    ##     digest~0.6.27, dplyr~1.0.5, ellipsis~0.3.2, ensembldb~2.14.1,
    ##     evaluate~0.14, fansi~0.4.2, farver~2.1.0, fastmap~1.1.0,
    ##     foreign~0.8-81, fs~1.5.0, generics~0.1.0, ggbio~1.38.0,
    ##     ggforce~0.3.3, ggplot2~3.3.3, glue~1.4.2, graph~1.68.0, grid~4.0.5,
    ##     gridExtra~2.3, gtable~0.3.0, highr~0.9, hms~1.0.0, htmlTable~2.1.0,
    ##     htmltools~0.5.1.1, htmlwidgets~1.5.3, httr~1.4.2, jpeg~0.1-8.1,
    ##     knitr~1.33, labeling~0.4.2, lattice~0.20-41, latticeExtra~0.6-29,
    ##     lazyeval~0.2.2, lifecycle~1.0.0, magrittr~2.0.1,
    ##     matrixStats~0.58.0, memoise~2.0.0, munsell~0.5.0, nnet~7.3-15,
    ##     openssl~1.4.4, parallel~4.0.5, pillar~1.6.0, pkgbuild~1.2.0,
    ##     pkgconfig~2.0.3, pkgload~1.2.1, plyr~1.8.6, png~0.1-7,
    ##     polyclip~1.10-0, prettyunits~1.1.1, processx~3.5.2, progress~1.2.2,
    ##     ps~1.6.0, purrr~0.3.4, rappdirs~0.3.3, remotes~2.3.0,
    ##     reshape~0.8.8, reshape2~1.4.4, rlang~0.4.11, rmarkdown~2.7,
    ##     roxygen2~7.1.1, rpart~4.1-15, rprojroot~2.0.2, rstudioapi~0.13,
    ##     rtracklayer~1.50.0, scales~1.1.1, sessioninfo~1.1.1, splines~4.0.5,
    ##     stats4~4.0.5, stringi~1.5.3, stringr~1.4.0, survival~3.2-11,
    ##     testthat~3.0.2, tibble~3.1.1, tidyr~1.1.3, tidyselect~1.1.1,
    ##     tools~4.0.5, tweenr~1.0.2, usethis~2.0.1, utf8~1.2.1, vctrs~0.3.8,
    ##     withr~2.4.2, xfun~0.22, xml2~1.3.2, yaml~2.2.1, zlibbioc~1.36.0
    ## \end{itemize}

# Reference {.unnumbered}
