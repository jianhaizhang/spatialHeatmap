# Overview

The `spatialHeatmap` package provides functionalities for visualizing cell-, tissue- and organ-specific data of biological assays by coloring the corresponding spatial features defined in anatomical images according to a numeric color key. The color scheme used to represent the assay values can be customized by the user. This core functionality of the package is called a spatial heatmap (SHM) plot. It is enhanced with visualization tools for groups of measured items (e.g. gene modules) sharing related abundance profiles, including matrix heatmaps combined with hierarchical clustering dendrograms and network representations. The functionalities of `spatialHeatmap` can be used either in a command-driven mode from within R or a graphical user interface (GUI) provided by a Shiny App that is also part of this package. While the R-based mode provides flexibility to customize and automate analysis routines, the Shiny App includes a variety of convenience features that will appeal to experimentalists and other users less familiar with R. Moreover, the Shiny App can be used on both local computers as well as centralized server-based deployments (e.g. cloud-based or custom servers) that can be accessed remotely as a public web service for using `spatialHeatmap’s` functionalities with community and/or private data. The functionalities of the `spatialHeatmap` package are illustrated in Figure 1.


<center><img title="Overview of spatialHeatmap" src="https://github.com/jianhaizhang/spatialHeatmap/blob/master/vignettes/img/spatialHeatmap_Design.jpg" ></center>

**Figure 1: Overview of `spatialHeatmap`**

# Installation 

To install the package, start R (version "3.6") and use the BiocManager::install command:

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

See user vignette [here](https://github.com/jianhaizhang/spatialHeatmap/tree/master/vignettes)




