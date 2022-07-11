\name{spatialHeatmap-package}
\alias{spatialHeatmap-package}
\alias{spatialHeatmap}
\docType{package}
\title{
\packageTitle{spatialHeatmap}
 Spatial Heatmap, Matrix Heatmap, Network
}
\description{
\packageDescription{spatialHeatmap}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spatialHeatmap}
\packageIndices{spatialHeatmap}

The spatialHeatmap package provides functionalities for visualizing cell-, tissue- and organ-specific data of biological assays by coloring the corresponding spatial features defined in anatomical images according to a numeric color key. The color scheme used to represent the assay values can be customized by the user. This core functionality is called a spatial heatmap plot. It is enhanced with nearest neighbor visualization tools for groups of measured items (e.g. gene modules) sharing related abundance profiles, including matrix heatmaps combined with hierarchical clustering dendrograms and network representations. The functionalities of spatialHeatmap can be used either in a command-driven mode from within R or a graphical user interface (GUI) provided by a Shiny App that is also part of this package. While the R-based mode provides flexibility to customize and automate analysis routines, the Shiny App includes a variety of convenience features that will appeal to many biologists. Moreover, the Shiny App has been designed to work on both local computers as well as server-based deployments (e.g. cloud-based or custom servers) that can be accessed remotely as a centralized web service for using spatialHeatmap's functionalities with community and/or private data.

As anatomical images the package supports both tissue maps from public repositories and custom images provided by the user. In general any type of image can be used as long as it can be provided in SVG (Scalable Vector Graphics) format, where the corresponding spatial features have been defined (see aSVG below). The numeric values plotted onto a spatial heatmap are usually quantitative measurements from a wide range of profiling technologies, such as microarrays, next generation sequencing (e.g. RNA-Seq and scRNA-Seq), proteomics, metabolomics, or many other small- or large-scale experiments. For convenience, several preprocessing and normalization methods for the most common use cases are included that support raw and/or preprocessed data. Currently, the main application domains of the spatialHeatmap package are numeric data sets and spatially mapped images from biological and biomedical areas. Moreover, the package has been designed to also work with many other spatial data types, such a population data plotted onto geographic maps. This high level of flexibility is one of the unique features of spatialHeatmap. Related software tools for biological applications in this field are largely based on pure web applications (Winter et al. 2007; Waese et al. 2017) or local tools (Maag 2018; Muschelli, Sweeney, and Crainiceanu 2014) that typically lack customization functionalities. These restrictions limit users to utilizing pre-existing expression data and/or fixed sets of anatomical image collections. To close this gap for biological use cases, we have developed spatialHeatmap as a generic R/Bioconductor package for plotting quantitative values onto any type of spatially mapped images in a programmable environment and/or in an intuitive to use GUI application.

}
\author{
\packageAuthor{spatialHeatmap}
Author: \packageAuthor{spatialHeatmap}
  Jianhai Zhang (PhD candidate at Genetics, Genomics and Bioinformatics, University of California, Riverside), Dr. Thomas Girke (Professor at Department of Botany and Plant Sciences, University of California, Riverside)
Maintainer: \packageMaintainer{spatialHeatmap}
  Jianhai Zhang <jzhan067@@ucr.edu>.
}
\references{
  https://www.w3schools.com/graphics/svg_intro.asp 

  https://shiny.rstudio.com/tutorial/  

  https://shiny.rstudio.com/articles/datatables.html  

  https://rstudio.github.io/DT/010-style.html  

  https://plot.ly/r/heatmaps/  

  https://www.gimp.org/tutorials/ 
 
  https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html  

  http://www.microugly.com/inkscape-quickguide/  

  https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html 

  https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 

  Winter, Debbie, Ben Vinegar, Hardeep Nahal, Ron Ammar, Greg V Wilson, and Nicholas J Provart. 2007. 
  "An 'Electronic Fluorescent Pictograph' Browser for Exploring and Analyzing Large-Scale Biological Data Sets." PLoS One 2 (8): e718

  Waese, Jamie, Jim Fan, Asher Pasha, Hans Yu, Geoffrey Fucile, Ruian Shi, Matthew Cumming, et al. 2017. 
  "EPlant: Visualizing and Exploring Multiple Levels of Data for Hypothesis Generation in Plant Biology." Plant Cell 29 (8): 1806–21

  Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angelica Liechti, et al. 2019.   
  "Gene Expression Across Mammalian Organ Development." Nature 571 (7766): 505–9

  Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas

  Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8

  McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97

  Maag, Jesper L V. 2018. "Gganatogram: An R Package for Modular Visualisation of Anatograms and Tissues Based on Ggplot2." F1000Res. 7 (September): 1576

  Muschelli, John, Elizabeth Sweeney, and Ciprian Crainiceanu. 2014. "BrainR: Interactive 3 and 4D Images of High Resolution Neuroimage Data." R J. 6 (1): 41–48

  Morgan, Martin, Valerie Obenchain, Jim Hester, and Hervé Pagès. 2018. SummarizedExperiment: SummarizedExperiment Container

  Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017).
  shiny: Web Application Framework for R. R package version 1.0.3.
  https://CRAN.R-project.org/package=shiny

  Winston Chang and Barbara Borges Ribeiro (2017). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.6.1.
  https://CRAN.R-project.org/package=shinydashboard

  Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R.
  Journal of Statistical Software, 30(4), 1-37. URL
  http://www.jstatsoft.org/v30/i04/.

  Jeroen Ooms (2017). rsvg: Render SVG Images into PDF, PNG, PostScript, or
  Bitmap Arrays. R package version 1.1. https://CRAN.R-project.org/package=rsvg

  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
  York, 2016.

  Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R
  package version 0.2. https://CRAN.R-project.org/package=DT

  Baptiste Auguie (2016). gridExtra: Miscellaneous Functions for "Grid" Graphics.
  R package version 2.2.1. https://CRAN.R-project.org/package=gridExtra

  Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and
  Tree Diagrams Using 'ggplot2'. R package version 0.1-20.
  https://CRAN.R-project.org/package=ggdendro

  Langfelder P and Horvath S, WGCNA: an R package for weighted correlation
  network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical Software,
  46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

  Simon Urbanek and Jeffrey Horner (2015). Cairo: R graphics device using cairo
  graphics library for creating high-quality bitmap (PNG, JPEG, TIFF), vector
  (PDF, SVG, PostScript) and display (X11 and Win32) output. R package version
  1.5-9. https://CRAN.R-project.org/package=Cairo

  R Core Team (2017). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL
  https://www.R-project.org/.

  Duncan Temple Lang and the CRAN Team (2017). XML: Tools for Parsing and
  Generating XML Within R and S-Plus. R package version 3.98-1.9.
  https://CRAN.R-project.org/package=XML

  Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram,
  Marianne Corvellec and Pedro Despouy (NA). plotly: Create Interactive Web
  Graphics via 'plotly.js'. https://plot.ly/r,
  https://cpsievert.github.io/plotly_book/, https://github.com/ropensci/plotly.

  Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R
  package version 1.10.4. https://CRAN.R-project.org/package=data.table

  R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter:
  methods for filtering genes from high-throughput experiments. R package version
  1.58.1.

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical Software,
  46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

  Almende B.V., Benoit Thieurmel and Titouan Robert (2017). visNetwork: Network
  Visualization using 'vis.js' Library. R package version 2.0.1.
  https://CRAN.R-project.org/package=visNetwork

  Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across 
  Sessions. R package version 1.12.1.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
\keyword{ spatial heatmap}
\seealso{
 \code{\link{norm_data}}, \code{\link{aggr_rep}}, \code{\link{filter_data}}, \code{\link{spatial_hm}}, \code{\link{submatrix}}, \code{\link{adj_mod}}, \code{\link{matrix_hm}}, \code{\link{network}}, \code{\link{return_feature}}, \code{\link{update_feature}}, \code{\link{shiny_shm}}, \code{\link{custom_shiny}}
}
\examples{

## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7
## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are
## included in this package. The complete raw count data are downloaded using the R package 
## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a 
## "data frame" input to exemplify data of simple samples/conditions, while toy data2 as 
## "SummarizedExperiment" to illustrate data involving complex samples/conditions.  

## Set up toy data.
     
# Access toy data1.
cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', 
package='spatialHeatmap')
df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
# Columns follow the namig scheme "sample__condition", where "sample" and "condition" stands 
# for organs and time points respectively.
df.chk[1:3, ]
     
# A column of gene annotation can be appended to the data frame, but is not required.  
ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
df.chk <- cbind(df.chk, ann=ann); df.chk[1:3, ]
     
# Access toy data2. 
cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
count.chk[1:3, 1:5]
     
# A targets file describing samples and conditions is required for toy data2. It should be made
# based on the experiment design, which is accessible through the accession number "E-MTAB-6769"
# in the R package ExpressionAtlas. An example targets file is included in this package and 
# accessed below. 
# Access the example targets file. 
tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
# Every column in toy data2 corresponds with a row in targets file. 
target.chk[1:5, ]
# Store toy data2 in "SummarizedExperiment".

library(SummarizedExperiment)
se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
# The "rowData" slot can store a data frame of gene annotation, but not required.
rowData(se.chk) <- DataFrame(ann=ann)
     
## As conventions, raw sequencing count data should be normalized, aggregated, and filtered to
## reduce noise.
     
# Normalize count data.
# The normalizing function "calcNormFactors" (McCarthy et al. 2012) with default settings is used.  
df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', log2.trans=TRUE)
se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)
# Aggregate count data.
# Aggregate "sample__condition" replicates in toy data1.
df.aggr.chk <- aggr_rep(data=df.nor.chk, aggr='mean')
df.aggr.chk[1:3, ]
# Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part" and
# "condition" is "age". 
se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
assay(se.aggr.chk)[1:3, 1:3]
# Filter out genes with low counts and low variance. Genes with counts over 5 (log2 unit) in at
# least 1\% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
# Filter toy data1.
df.fil.chk <- filter_data(data=df.aggr.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
# Filter toy data2.
se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age', 
pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)

## Spatial heatmaps.

# To make spatial heatmaps, a pair of formatted data and pre-annotated SVG (aSVG) file are 
# required. If the data is a "data frame", the formatting is to use the naming scheme 
# "sample__condition" in column names. If "SummarizedExperiment", the "sample" and "condition"
# replicates should be defined in the "colData" slot. In the aSVG, each spatial feature has a
# unique identifier. The numeric values are mapped to spatial features and translated into 
# colors according to their identifiers programatically. The mapped images are called spatial 
# heatmaps.

# The following shows how to download the corresponding pre-annotated aSVG file from the EBI
# SVG repository based on above tissues and species involved, i.e. c('heart', 'brain') and 
# c('gallus') respectively. See the function "return_feature" for details. An empty directory
# is recommended so as to avoid overwriting existing SVG files. Here "tmp.dir" is used. 

# To meet the package building requirements, the code of querying aSVG remotely is not evaluated. 
# The matching aSVG "gallus_gallus.svg" is included in this package and accessed.

\donttest{
# Make an empty directory "tmp.dir" if not exist.
tmp.dir <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm')
# Remote aSVG repos.
data(aSVG.remote.repo)
tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
tmp.dir.ebi <- paste0(tmp, '/ebi.zip')
tmp.dir.shm <- paste0(tmp, '/shm.zip')
# Download the remote aSVG repos as zip files. According to Bioconductor's 
# requirements, downloadings are not allowed inside functions, so the repos are 
# downloaded before calling "return_feature".  
download.file(aSVG.remote.repo$ebi, tmp.dir.ebi)
download.file(aSVG.remote.repo$shm, tmp.dir.shm)
remote <- list(tmp.dir.ebi, tmp.dir.shm)

# Query aSVGs from remote repos.
feature.df <- return_feature(feature=c('heart', 'brain'), species=c('gallus'), dir=tmp.dir, 
match.only=FALSE, remote=remote)
feature.df
# The path of matching aSVG.
svg.chk <- paste0(tmp.dir, '/gallus_gallus.svg')
}

# Get the matching aSVG path from the package.
svg.chk <- system.file("extdata/shinyApp/example", "gallus_gallus.svg", 
package="spatialHeatmap")

# Plot spatial heatmaps on gene "ENSGALG00000019846". In the middle are spatial heatmaps. Only
# aSVG features with matching countparts in data are colored. On the right is the legend plot,
# only the matching features are labeled.
# Toy data1. 
spatial_hm(svg.path=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4, 
legend.r=1.9, sub.title.size=7, ncol=3)
# Save spaital heatmaps as HTML and video files by assigning "tmp.dir" to "out.dir". 
\donttest{    
tmp.dir <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm')
spatial_hm(svg.path=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4, legend.r=1.9,
sub.title.size=7, ncol=3, out.dir=tmp.dir)
}
# Toy data2.
spatial_hm(svg.path=svg.chk, data=se.fil.chk, ID='ENSGALG00000019846', legend.r=1.9, 
legend.nrow=2, sub.title.size=7, ncol=3)

# When plot spatial heatmaps, the data can also come as as a simple vector. The following 
# gives an example on a vector of 3 random values. 
# Random values.
vec <- sample(1:100, 3)
# Name the vector slots. The last name is assumed as a random sample without a matching 
# feature in aSVG.
names(vec) <- c('brain', 'heart', 'notMapped')
vec
# Plot.
spatial_hm(svg.path=svg.chk, data=vec, ID='geneX', height=0.6, legend.r=1.5, ncol=1)

# Plot spatial heatmaps on aSVGs of two Arabidopsis thaliana development stages.

# Make up a random numeric data frame.
df.test <- data.frame(matrix(sample(x=1:100, size=50, replace=TRUE), nrow=10))
colnames(df.test) <- c('shoot_totalA__condition1', 'shoot_totalA__condition2', 
'shoot_totalB__condition1', 'shoot_totalB__condition2', 'notMapped')
rownames(df.test) <- paste0('gene', 1:10) # Assign row names 
df.test[1:3, ]

# aSVG of development stage 1.
svg1 <- system.file("extdata/shinyApp/example", "arabidopsis.thaliana_organ_shm1.svg",
package="spatialHeatmap")
# aSVG of development stage 2.
svg2 <- system.file("extdata/shinyApp/example", "arabidopsis.thaliana_organ_shm2.svg",
package="spatialHeatmap")
# Spatial heatmaps.
spatial_hm(svg.path=c(svg1, svg2), data=df.test, ID=c('gene1'), height=0.8, legend.r=1.6, 
preserve.scale=TRUE) 

## If users want to use custom identifiers for spatial features in the aSVG file, the function
# "update_feature" should be used. For illustration purpose, the aSVG "gallus_gallus.svg" in 
# this package is copied to 'tmp.dir' as example.

\donttest{
# Make an empty directory "tmp.dir" if not exist.
tmp.dir <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm')
# Make a copy of "gallus_gallus.svg".
file.copy(from=svg.chk, to=tmp.dir, overwrite=FALSE)
# Remote aSVG repos.
data(aSVG.remote.repo)
tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
tmp.dir.ebi <- paste0(tmp, '/ebi.zip')
tmp.dir.shm <- paste0(tmp, '/shm.zip')
# Download the remote aSVG repos as zip files. According to Bioconductor's 
# requirements, downloadings are not allowed inside functions, so the repos are 
# downloaded before calling "return_feature".  
download.file(aSVG.remote.repo$ebi, tmp.dir.ebi)
download.file(aSVG.remote.repo$shm, tmp.dir.shm)
remote <- list(tmp.dir.ebi, tmp.dir.shm)

# Query "gallus_gallus.svg" in remote repos.
feature.df <- return_feature(feature=c('heart', 'brain'), species=c('gallus'), dir=tmp.dir, 
match.only=TRUE, remote=remote)
feature.df

# New features.
ft.new <- c('BRAIN', 'HEART')
# Add new features to the first column.
feature.df.new <- cbind(featureNew=ft.new, feature.df)
feature.df.new
# Update features.
update_feature(df.new=feature.df.new, dir=tmp.dir)
}


## Matrix heatmap

# The matrix heatmap and following network are supplements to the core feature of spatial 
# heatmap. First, nearest neighbors are selected for each target gene according to correlation
# (default) or distance measure independently. There are three alternative parameters used for
# the selection: "p" is the proportion of top nearest neighbors, "n" is the number of top 
# nearest neighbors, and "v" is a specific cutoff value for correlation or distance. Then 
# target genes and their nearest neighbors are hierarchically clustered and visualized in 
# static or interactive matrix heatmap, where target genes are labeled by black lines. If the 
# data is "SummarizedExperiment", the argument "ann" is the column name of gene annotation in
# "rowData" slot. It is only relevant if users want to see annotation when mousing over a node
# in the interactive network below, so it is optional. Here "ann='ann'" is set and the 
# corresponding annotation is appended to selected nearest neighbors.

# Select nearest neighbors for target genes 'ENSGALG00000019846' and 'ENSGALG00000000112'.
df.sub.mat <- submatrix(data=df.fil.chk, ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), p=0.1)
se.sub.mat <- submatrix(data=se.fil.chk, ann='ann', ID=c('ENSGALG00000019846', 
'ENSGALG00000000112'), p=0.1) 

# In the following, "df.sub.mat" and "se.sub.mat" is used in the same way, so only 
# "se.sub.mat" illustrated.

# The subsetted matrix is partially shown below.
se.sub.mat[c('ENSGALG00000019846', 'ENSGALG00000000112'), c(1:2, 63)]

# Static matrix heatmap.
matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, angleCol=80, 
angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=TRUE, 
arg.lis1=list(offsetRow=0.01, offsetCol=0.01))

# Interactive matrix heatmap.
\donttest{ matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=FALSE, 
arg.lis1=list(offsetRow=0.01, offsetCol=0.01)) }

## Network 

# Network analysis with WGCNA (Langfelder and Horvath 2008) is applied on the subsetted matix
# visualized in the matrix heatmap. The gene module containing a specifc target gene is 
# visualized in static and interactive network graphs. Briefly, a correlation matrix or 
# distance matrix is computed on all genes in matrix heatmap, and transformed to an adjacency
# matrix and topological overlap matrix (TOM) sequentially, which are advanced measures to 
# quantify coexpression similarity. Then network modules are identified by hierarchinally 
# clustering the TOM-transformed dissimilarity matrix 1-TOM, which are clusters of genes with 
# highly similar coexpression profiles. The module containing a target gene is finally 
# displayed as network graphs. Refer to function "adj_mod" for details.

# Adjacency matrix and module identification

# The modules are identified by "adj_mod". It returns a list containing an adjacency matrix and
# a data frame of module assignment. 
adj.mod <- adj_mod(data=se.sub.mat)

# The adjacency matrix is a measure of co-expression similarity between genes, where larger 
# value denotes more similarity.
adj.mod[['adj']][1:3, 1:3]

# The modules are identified at two alternative sensitivity levels (ds=2 or 3). From 2 to 3, 
# more modules are identified but module sizes are smaller. The two sets of module assignment
# are returned in a data frame. The first column is ds=2 while the second is ds=3. The numbers
# in each column are module labels, where "0" indicates genes not assigned to any module.
adj.mod[['mod']][1:3, ]

# Static network. In the graph, nodes are genes and edges are adjacencies between genes. The 
# thicker edge denotes higher adjacency (co-expression similarity) while larger node indicates 
# higher gene connectivity (sum of a gene's adjacency with all its direct neighbors). The target
# gene is labeled by "_target". The node connectivity increases from "turquoise" to "violet", 
# and the adjacency increases from "yellow" to "blue".
network(ID="ENSGALG00000019846", data=se.sub.mat, adj.mod=adj.mod, adj.min=0.7, 
vertex.label.cex=1.5, vertex.cex=4, static=TRUE)

# Interactive network. Same with static mode, the target gene ID is appended "_target".  
\donttest{ network(ID="ENSGALG00000019846", data=se.sub.mat, adj.mod=adj.mod, static=FALSE) }

## Shiny App

# In additon to generating spatial heatmaps and corresponding gene context plots from R, 
# spatialHeatmap includes a Shiny App (https://shiny.rstudio.com/) that provides access to the
# same functionalities from an intuitive-to-use web browser interface. Apart from being very 
# user-friendly, this App conveniently organizes the results of the entire visualization 
# workflow in a single browser window with options to adjust the parameters of the individual
# components interactively. This app is launched by the function "shiny_shm" without any 
# parameters. Upon launched, the app automatically displays a pre-formatted example. 
\donttest{shiny_shm()}

# The gene expression data and aSVG image files are uploaded to the Shiny App as tabular 
# text (e.g. in CSV or TSV format) and SVG file, respectively. To also allow users to upload
# gene expression data stored in "SummarizedExperiment" objects, one can export them from R 
# to a tabular file with the "filter_data" function. In this function call, the user sets a 
# desired directory path under "dir" (see below). Within this directory the tabular file will
# be written to "customData.txt" in TSV format. The column names in the exported tabular file
# preserve the experimental design information from the "colData" slot by concatenating the 
# corresponding sample and condition information separated by double underscores. An example 
# of this format is shown in below.

# To interactively view functional descriptions by moving the cursor over network nodes, the 
# corresponding annotation column needs to be present in the "rowData" slot and its column
# name assigned to the "ann" argument. In the exported tabular file the extra annotation
# column is appended to the expression matrix.
\donttest{ se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', 
con.factor='age', pOA=c(0.01, 5), CV=c(0.2, 100), dir='./'); assay(se.fil.chk)[1:3, 1:3] }

# The Shiny app can be customized by including user-provided default examples and default 
# parameters. See the fucntion "custom_shiny" for details.


}
