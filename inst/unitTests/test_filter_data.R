
test_filter_data <- function() {
# The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
library(ExpressionAtlas); library(SummarizedExperiment) 

cache.pa <- '~/.cache/shm' # The path of cache.
rse.hum <- read_cache(cache.pa, 'rse.hum') # Read data from cache.
if (is.null(rse.hum)) { # Save downloaded data to cache if it is not cached.
  rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.hum)
}
assay(rse.hum)[1:3, 1:3]
# A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
tar.hum <- system.file('extdata/shinyApp/example/target_human.txt', package='spatialHeatmap')
target.hum <- read.table(tar.hum, header=TRUE, row.names=1, sep='\t')
# The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
target.hum[c(1:3, 41:42), 4:5]
# Place the targets file into "colData" slot as a DataFrame class. 
colData(rse.hum) <- DataFrame(target.hum)

# The count matrix is normalised with estimateSizeFactors (type="ratio").
se.nor.hum <- norm_data(data=rse.hum, norm.fun='ESF', data.trans='log2')
# Average replicates of concatenated sample__condition.
se.aggr.hum <- aggr_rep(data=se.nor.hum, sam.factor='organism_part', con.factor='disease')
assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.

# Genes with low expression level and low variantion are always filtered. 
se.fil.hum <- filter_data(data=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100))

checkTrue(is(se.fil.hum, "SummarizedExperiment"))

}

