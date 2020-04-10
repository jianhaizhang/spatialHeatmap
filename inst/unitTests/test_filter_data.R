
test_filter_data <- function() {
# The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
library(ExpressionAtlas); library(SummarizedExperiment) 
rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]

# A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
# The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
target.hum[c(1:3, 41:42), 4:5]
# Place the targets file into "colData" slot as a DataFrame class. 
colData(rse.hum) <- DataFrame(target.hum)

# The count matrix is normalised with estimateSizeFactors (type="ratio").
se.nor.hum <- norm_data(se=rse.hum, method.norm='ratio', data.trans='log2')
# Average replicates of concatenated sample__condition.
se.aggr.hum <- aggr_rep(se=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.

# Genes with low expression level and low variantion are always filtered. 
se.fil.hum <- filter_data(se=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100), dir=NULL)

checkTrue(is(se.fil.hum, "SummarizedExperiment"))

}

