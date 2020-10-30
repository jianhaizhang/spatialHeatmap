## Make target file-human brain example.
# Select "E-GEOD-67196", which is RNA-seq data of cerebellum and frontal cortex in brain.
library(ExpressionAtlas)
cache.pa <- '~/.cache/shm' # The path of cache.
rse.hum <- read_cache(cache.pa, 'rse.hum') # Read data from cache.
if (is.null(rse.hum)) { # Save downloaded data to cache if it is not cached.
  rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.hum)
}
# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.hum)); df.con[1:3, ]

# The matching tissues in the SVG template are "cerebellum" and "frontal_cortex" while in the expression data are "cerebellum" and "frontal cortex". The tissue names must be exactly same between the SVG and the expression matrix, so change "frontal cortex" to "frontal_cortex".
# df.con$organism_part <- sub(' ', '_', df.con$organism_part)
# Use abbreviations of the conditions ("individual", "disease", "genotype").
df.con$individual <- sub('individual', 'indi', df.con$individual)
df.con$individual <- sub('Individual ', 'indi', df.con$individual)
df.con$disease <- sub('amyotrophic lateral sclerosis', 'ALS', df.con$disease)
df.con$genotype <- sub('presence of a C9orf72 repeat expansion', 'pre.C9orf72', df.con$genotype)
df.con$genotype <- sub('absence of a C9orf72 repeat expansion', 'ab.C9orf72', df.con$genotype)
df.con$genotype <- sub('not applicable', 'NA', df.con$genotype)

write.table(df.con, 'target_human.txt', col.names=TRUE, row.names=TRUE, sep='\t')

## Make target file-mouse organ example.

# Select "E-MTAB-2801".
rse.mus <- read_cache(cache.pa, 'rse.mus') # Read data from cache.
if (is.null(rse.mus)) { # Save downloaded data to cache if it is not cached.
  rse.mus <- getAtlasData('E-MTAB-2801')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.mus)
}
# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.mus)); df.con[1:3, ]
df.con$organism_part <- gsub('skeletal muscle tissue', 'skeletal muscle', df.con$organism_part)
df.con$strain <- make.names(df.con$strain)

write.table(df.con, 'target_mouse.txt', col.names=TRUE, row.names=TRUE, sep='\t')

## Make target file-chicken organ example.

# Select E-MTAB-6769.
rse.chk <- read_cache(cache.pa, 'rse.chk') # Read data from cache.
if (is.null(rse.chk)) { # Save downloaded data to cache if it is not cached.
  rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.chk)
}
# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.chk)); df.con[1:3, ]
df.con$age <- gsub('(.*) (.*)', '\\2\\1', df.con$age)

write.table(df.con, 'target_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')

## Make targets file of the Arabidopsis shoot example.
library(GEOquery)
# Download GEO dataset GSE14502 (normalised by RMA).
gset <- read_cache(cache.pa, 'gset') # Retrieve data from cache.
if (is.null(gset)) { # Save downloaded data to cache if it is not cached.
  gset <- getGEO("GSE14502", GSEMatrix=TRUE, getGPL=TRUE)[[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, gset)
}

# Assign tissue/condition names.
se <- as(gset, "SummarizedExperiment")
# Targets file.
col.name <- colnames(se)
titles <- make.names(colData(se)$title)
title.factor <- sub('_rep\\d+$', '', titles)
sam <- gsub('(root|shoot)(_)(control|hypoxia)(_.*)', '\\1\\4', title.factor) 
con <- gsub('(root|shoot)(_)(control|hypoxia)(_.*)', '\\3', title.factor) 
target.geo <- data.frame(col.name=col.name, titles=titles, row.names=2, samples=sam, conditions=con, stringsAsFactors=FALSE) 
write.table(target.geo, 'target_geo.txt', col.names=TRUE, row.names=TRUE, sep='\t')

library(ath1121501.db); library(genefilter)
# Gene annotation.
ann <- data.frame(SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=";"), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=";"), stringsAsFactors=FALSE)
ann <- subset(ann, !grepl("AFFX|s_at|a_at|x_at|r_at", rownames(ann)) & !duplicated(SYMBOL) & SYMBOL!="NA")

# Optional: data frame for gene metadata.
mat <- assay(se); colnames(mat) <- rownames(target.geo)
int <- intersect(rownames(mat), rownames(ann))
mat <- mat[int, ]; rownames(mat) <- make.names(ann[int, 'SYMBOL'])
ffun <- filterfun(pOverA(0.03, 6), cv(0.25, 100))
filtered <- genefilter(mat, ffun); mat <- mat[filtered, ]
write.table(mat, 'target_arab.txt', col.names=TRUE, row.names=TRUE, sep='\t')











