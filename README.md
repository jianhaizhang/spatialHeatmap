# Summary

The spatialHeatmap is an R/Bioconductor package and is designed for visualising gene expression data (sequencing, microarray, etc.). The core feature is to map expression values of target genes under different conditions to different tissues, which are called spatial heatmaps. To make it more versatile, the matrix heatmap and network features were also developed, which display the target gene in the context of corresponding gene module. All the utilities implemented in this R package are also combined as a web-browser based Shiny app. The application of this package is not limited to gene expression data. It can be used as long as a data matrix and an associated SVG image are provided, such as population data generated in different years across different cities.

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

The R functions and their usage is given in the main vignette "vignette.html", and the example data matrix and svg image are available in this package. How to properly format and associate a custom SVG image with a data matrix is provided in an independent vignette "SVG_tutorial.html". For a quick test, users can launch the combined Shiny app then either select the top 7 examples from the left menu or download the example data matrix and SVG image in the insturction page. The integrated Shiny app is also available here: https://tgirke.shinyapps.io/spatialHeatmap/


