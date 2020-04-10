# Summary

This R package spatialHeatmap is designed to intuitively visualise values from a data matrix onto an SVG image provided by users and promote hypothesis generation.  

The core feature Spatial Heatmap is to map expression profile of a target gene under different conditions to different tissues on a formatted SVG image, where different tissues are predefined. After mapping, the expression profile is represented as different colours across tissues in the image, which are called spatial heatmaps. This feature allows to input multiple target genes. If so, the spatial heatmaps of different genes are generated sequentially on the same page. There is also an option to display these spatial heatmaps by genes or by conditions, which makes it flexible for users to compare expression profiles of the same gene across conditions or different genes across the same condition. In the multiple-layer tissues, if some tissues are covered by front tissues and thus not visible in the spatial heatmaps, then front tissues can be set transparent to exhibit convered ones.  

The data matrix is in the form of “SummarizedExperiment” (Morgan et al. 2018). The “assays” slot contains the gene expression matrix with columns and rows being tissue/conditions and genes respectively, and the “colData” slot contains replicates of samples and conditions. The tissues are pre-defined in an SVG image. To visualise the data successfully, the data matrix and SVG image should meet the format requirements, which includes tissue ids in the SVG image must be identical to tissue names in the “colData”, same tissues should be grouped in the SVG image, etc.  

In addition, the accessory features of Matrix Heatmap and Network display the target gene in the context of corresponding gene module, which make this package more informative. All the utilities can be implemented step-by-step or lanuched as an interactive Shiny App (Chang et al., n.d.; Chang and Borges Ribeiro 2018).  

# Application Scope

The application is not limited to biological data. It is applicable as long as a pair of formatted data matrix in “SummariaedExperiment” and SVG image are provided. Other exmaples include population data collected in different years across different cities, health data of an individual under different conditions, ect.


# Installation 

To install the package, start R (version "3.5") and use the BiocManager::install command:

```{r, eval=FALSE, echo=TRUE, warnings=FALSE} 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("spatialHeatmap")

```
To obtain the most recent updates immediately, one can install it directly from github as follow:
                                                                                                                                                                 
```{r, eval=FALSE, echo=TRUE, warnings=FALSE}

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("jianhaizhang/spatialHeatmap", build_vignettes=TRUE, dependencies=TRUE)

```


Morgan, Martin, Valerie Obenchain, Jim Hester, and Hervé Pagès. 2018. SummarizedExperiment: SummarizedExperiment Container.  
Chang, Winston, and Barbara Borges Ribeiro. 2018. Shinydashboard: Create Dashboards with ’Shiny’. https://CRAN.R-project.org/package=shinydashboard.  
Chang, Winston, Joe Cheng, JJ Allaire, Yihui Xie, and Jonathan McPherson. n.d. Shiny: Web Application Framework for R. https://CRAN.R-project.org/package=shiny.


