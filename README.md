# Summary

This R package spatialHeatmap is for intuitive visualisation of configured data matrix and SVG image provided by users and promote hypothesis generation. Examples include but not limit to RNA-seq, microarray, qPCR, subcellular localisation of proteins, etc., regardless of species.

The core feature “Spatial Heatmap” is to map expression profile of a target gene under different conditions to different cells/tissues/organs (samples) on an associated SVG image, where different samples are predefined. After mapping, the expression profile is represented as different colours across samples in the image, which are called spatial heatmaps. In addition, the accessory features of “Matrix Heatmap” and “Network” display the target gene in the context of corresponding gene module, which make this package more versatile. All the utilities implemented in this R package are also combined as a web-browser based Shiny app, which can be used locally or online on the Shiny server.


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

# Usage

The application of this package is not limited to gene expression data. It can be used as long as a data matrix and a configured SVG image are provided. Other exmaples include population data collected in different years across different cities, health data of an individual under different conditions, ect.


