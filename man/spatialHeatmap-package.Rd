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

The rapid advance in high-throughput biotechnology (e.g. next-generation sequencing, microarray) has generated a large amount of biological data. Accordingly, tools are developed over time to visualise these data, such as eFP (Winter et al. 2007), ePlant (Waese et al. 2017), gganatogram (Maag 2018), brainR (Muschelli, Sweeney, and Crainiceanu 2014), etc. These tools map the gene expression data onto a pre-defined tissue image where the data are measured. The great feature is that they display the data intuitively and interactively and therefore could promote hypothesis generation. However, the visualisation is based on data and image pre-configured by developers. Therefore, users have difficuty in displaying their custom data. To address this issue, we developed the R/Bioconductor package spatialHeatmap.

The spatialHeatmap is designed to allow users to visualise custom data on tissue images in an intuitive way, and therefore promote hypothesis generation. The core feature Spatial Heatmap is to map expression profile of a target gene under different conditions to different tissues on a formatted SVG image, where different tissues are predefined. After mapping, the expression profile is represented as different colours across tissues in the image, which are called spatial heatmaps. This feature allows to input multiple target genes. If so, the spatial heatmaps of different genes are generated sequentially on the same page. There is also an option to display these spatial heatmaps by genes or by conditions, which makes it flexible for users to compare expression profiles of the same gene across conditions or different genes across the same condition. In the multiple-layer tissues, if some tissues are covered by front tissues and thus not visible in the spatial heatmaps, then front tissues can be set transparent to exhibit convered ones, which is shown in the mouse example.

The data is in the form of "SummarizedExperiment" (Morgan et al. 2018). The "assays" slot contains the gene expression matrix with columns and rows being tissue/conditions and genes respectively, and the "colData" slot contains replicates of samples and conditions. The tissues are pre-defined in an SVG image. To visualise the data successfully, the data matrix and SVG image should meet the format requirements, which includes tissue ids in the SVG image must be identical to tissue names in the "colData", same tissues should be grouped in the SVG image, etc.

In addition, the accessory features of Matrix Heatmap and Network display the target gene in the context of corresponding gene module, which make this package more informative. All the utilities can be implemented step-by-step or lanuched as an interactive Shiny App (Chang et al., n.d.; Chang and Borges Ribeiro 2018).

The application is not limited to biological data. It is applicable as long as a pair of formatted data matrix in "SummariaedExperiment" and SVG image are provided. Other exmaples include population data collected in different years across different cities, health data of an individual under different conditions, ect.

}
\author{
\packageAuthor{spatialHeatmap}
Author: \packageAuthor{spatialHeatmap}
  Jianhai Zhang (PhD candidate Genetics, Genomics and Bioinformatics, University of California, Riverside), Dr. Thomas Girke (Professor at Department of Botany and Plant Sciences, University of California, Riverside)
Maintainer: \packageMaintainer{spatialHeatmap}
  Jianhai Zhang <jzhan067@ucr.edu; zhang.jianhai@hotmail.com>.
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

  Winter, Debbie, Ben Vinegar, Hardeep Nahal, Ron Ammar, Greg V Wilson, and Nicholas J Provart. 2007. 
  "An 'Electronic Fluorescent Pictograph' Browser for Exploring and Analyzing Large-Scale Biological Data Sets." PLoS One 2 (8): e718

  Waese, Jamie, Jim Fan, Asher Pasha, Hans Yu, Geoffrey Fucile, Ruian Shi, Matthew Cumming, et al. 2017. 
  "EPlant: Visualizing and Exploring Multiple Levels of Data for Hypothesis Generation in Plant Biology." Plant Cell 29 (8): 1806–21

  Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, 
  Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82

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
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
\keyword{ spatial heatmap}
\seealso{
 \code{\link{filter_data}}, \code{\link{adj_mod}}, \code{\link{shiny_all}}, \code{\link{spatial_hm}}, \code{\link{matrix_hm}}, \code{\link{network}}
}
\examples{
# The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
library(ExpressionAtlas)
rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]

# A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
# The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
target.hum[c(1:3, 41:42), 4:5]
# Place the targets file into "colData" slot as a DataFrame class. 
colData(rse.hum) <- DataFrame(target.hum)
 
# For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data matrix from the brain data.
# Borrow a data matrix.
df <- assay(rse.hum); df[1:2, 1:3]
# Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
 
# The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
se.nor.hum <- norm_data(se=rse.hum, method.norm='ratio', data.trans='log2')

# Average replicates of concatenated sample__condition.
se.aggr.hum <- aggr_rep(se=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.
 
# Genes with low expression level and low variantion are always filtered. 
se.fil.hum <- filter_data(se=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100), dir=NULL)

# Formatted SVG image.
svg.hum <- system.file("extdata/shinyApp/example", "homo_sapiens.brain.svg", package="spatialHeatmap")
# Plot spatial heatmaps of gene ENSG00000008196.
spatial_hm(svg=svg.hum, se=se.fil.hum, ID='ENSG00000008196', col.com=c("yellow", "blue", "purple"), width=1, height=0.5, sub.title.size=11, layout="gene", ncol=2, tis.trans=NULL, legend.position=c(0.5, -0.15), legend.nrow=1)

# Detect modules. 
adj.mod <- adj_mod(se=se.fil.hum, type="signed", minSize=15, dir=NULL)
# The first column is ds=2 while the second is ds=3. The numbers in each column are module labels with "0" meaning genes not assigned to any modules.
adj.mod[['mod']][1:3, ]

# Plot matrix heatmap on gene ENSG00000008196 with ds='3'. Set "static=TRUE" to launch the interactive mode. 
matrix_hm(geneID="ENSG00000008196", se=se.fil.hum, adj.mod=adj.mod, ds="3", scale="no", angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=TRUE, arg.lis1=list(offsetRow=0.1, offsetCol=0.1))

# Plot network on gene ENSG00000008196 with ds='3'. Set "static=TRUE" to launch the interactive mode. 
network(geneID="ENSG00000008196", se=se.fil.hum, ann=NULL, adj.mod=adj.mod, ds="3", adj.min=0.999, con.min=0, vertex.label.cex=1, vertex.cex=0.1, static=TRUE)

# Launch the integrated Shiny app.
\donttest{ shiny_all() }

}
