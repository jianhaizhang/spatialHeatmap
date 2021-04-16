library(spatialHeatmap); library(ExpressionAtlas)
## Make target file/Shiny app data-human brain example.
# Select "E-GEOD-67196", which is RNA-seq data of cerebellum and frontal cortex in brain.
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

target.hum <- df.con
write.table(target.hum, 'target_human.txt', col.names=TRUE, row.names=TRUE, sep='\t')

colData(rse.hum) <- DataFrame(target.hum)
# Normalise.
se.nor.hum <- norm_data(data=rse.hum, norm.fun='ESF', log2.trans = FALSE)

# Filter.
se.fil.hum <- filter_data(data=se.nor.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.7, 50), CV=c(0.5, 100), dir=NULL); se.fil.hum
# Data matrix.
expr.hum <- as.data.frame(assay(se.fil.hum))

# Row metadata.
library(org.Hs.eg.db)
columns(org.Hs.eg.db); keytypes(org.Hs.eg.db)
row.met <- select(org.Hs.eg.db, keys=rownames(expr.hum), columns=c('ENSEMBL', 'SYMBOL', 'GENENAME'), keytype="ENSEMBL")
row.met <- row.met[!duplicated(row.met$ENSEMBL), ]
row.met$metadata <- paste0(row.met$SYMBOL, ': ', row.met$GENENAME)
# Append row metadata to data matrix.
expr.hum <- expr.hum[row.met$ENSEMBL, ]
expr.hum <- cbind(expr.hum, row.met[, 'metadata', drop = FALSE])

# colnames(expr.hum) <- gsub("_", ".", colnames(expr.hum))
write.table(expr.hum, 'expr_human.txt', col.names=TRUE, row.names=TRUE, sep='\t')





## Make target file/Shiny app data-mouse organ example.

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

target.mus <- df.con
write.table(target.mus, 'target_mouse.txt', col.names=TRUE, row.names=TRUE, sep='\t')

colData(rse.mus) <- DataFrame(target.mus)

se.nor.mus <- norm_data(data=rse.mus, norm.fun='ESF', log2.trans = FALSE)
se.fil.mus <- filter_data(data=se.nor.mus, sam.factor='organism_part', con.factor='strain', pOA=c(0.3, 30), CV=c(1.5, 100), dir=NULL); se.fil.mus
# Data matrix.
expr.mus <- as.data.frame(assay(se.fil.mus))

# Row metadata.
library(org.Mm.eg.db)
columns(org.Mm.eg.db); keytypes(org.Mm.eg.db)
row.met <- select(org.Mm.eg.db, keys=rownames(expr.mus), columns=c('ENSEMBL', 'SYMBOL', 'GENENAME'), keytype="ENSEMBL")
row.met <- row.met[!duplicated(row.met$ENSEMBL), ]
row.met$metadata <- paste0(row.met$SYMBOL, ': ', row.met$GENENAME)
# Append row metadata to data matrix.
expr.mus <- expr.mus[row.met$ENSEMBL, ]
expr.mus <- cbind(expr.mus, row.met[, 'metadata', drop = FALSE])

# colnames(expr.mus) <- gsub("_", ".", colnames(expr.mus))
write.table(expr.mus, 'expr_mouse.txt', col.names=TRUE, row.names=TRUE, sep='\t')



## Make target file/Shiny app data-chicken organ example.

# Select E-MTAB-6769.
rse.chk <- read_cache(cache.pa, 'rse.chk') # Read data from cache.
if (is.null(rse.chk)) { # Save downloaded data to cache if it is not cached.
  rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.chk)
}
# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.chk)); df.con[1:3, ]
df.con$age <- gsub('(.*) (.*)', '\\2\\1', df.con$age)

target.chk <- df.con
write.table(target.chk, 'target_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')

colData(rse.chk) <- DataFrame(target.chk)
se.nor.chk <- norm_data(data=rse.chk, norm.fun='ESF')
# Average the tissue-condition replicates.
se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.05, 5), CV=c(1.5, 100), dir=NULL)
# Data matrix.
expr.chk <- assay(se.fil.chk)
# colnames(expr.chk) <- gsub("_", ".", colnames(expr.chk))
write.table(expr.chk, 'expr_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')



colData(rse.chk) <- DataFrame(target.chk)
se.nor.chk <- norm_data(data=rse.chk, norm.fun='ESF', log2.trans = FALSE)
se.fil.chk <- filter_data(data=se.nor.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.3, 700), CV=c(0.5, 100), dir=NULL); se.fil.chk
# Data matrix.
expr.chk <- assay(se.fil.chk)
# colnames(expr.chk) <- gsub("_", ".", colnames(expr.chk))
write.table(expr.chk, 'expr_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')





## Make targets file/Shiny app data for the Arabidopsis shoot example.
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
target.sh <- mat
write.table(target.sh, 'target_arab.txt', col.names=TRUE, row.names=TRUE, sep='\t')

colData(se.sh) <- DataFrame(target.sh)
# Average the sample-condition replicates.
se.aggr.sh <- aggr_rep(data=se.sh, sam.factor='samples', con.factor='conditions', aggr='mean')
se.fil.sh <- filter_data(data=se.aggr.sh, sam.factor='samples', con.factor='conditions', pOA=c(0.03, 6), CV=c(0.29, 100), dir=NULL)
# Data matrix.
expr.sh <- assay(se.fil.sh)
# colnames(expr.sh) <- gsub("_", ".", colnames(expr.sh))
expr.sh <- cbind.data.frame(expr.sh, ann=rowData(se.fil.sh)[, 'Target.Description'], stringsAsFactors=FALSE)
expr.sh <- expr.sh <- expr.sh[c(9, 1:8, 10:372), ]
write.table(expr.sh, 'expr_arab.txt', col.names=TRUE, row.names=TRUE, sep='\t')


## Make target file/Shiny app data-coleoptile example.
# Retrieve the data. 
rse.clp <- read_cache(cache.pa, 'rse.clp') # Retrieve data from cache.
if (is.null(rse.clp)) { # Save downloaded data to cache if it is not cached.
  rse.clp <- getAtlasData('E-GEOD-115371')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.clp)
}

# Experimental design.
cdat <- colData(rse.clp)

# Edit condition entries.
old <- c('none', 'anaerobic environment', '72 hour anaerobic environment; 24 hour aerobic environment')
new <- c('aerobic', 'anaerobic', 'NA')
tar <- edit_tar(cdat, 'stimulus', old, new)
tar <- as.data.frame(tar)
unique(tar$stimulus)

# Copy stimulus as a condition column.
tar$con <- tar$stimulus
# Edit contidion column.
tar <- edit_tar(tar, 'con', c('aerobic', 'anaerobic'), c('A', 'N'))
unique(tar$con)

# Edit time entries.
tar <- edit_tar(tar, 'age', unique(tar$age), c('0h', '1h', '3h', '12h', '24h', '48h', '72h', '96h'))
tar <- edit_tar(tar, 'age', '96h', '72N24A', tar$stimulus=='NA')
unique(tar$age)

# Edit spatial feature entries.
tar <- edit_tar(tar, 'organism_part', unique(tar$organism_part), c('embryo', 'embryoColeoptile', 'coleoptile'))
unique(tar$organism_part)

# Export targets file.
write.table(tar, 'target_coleoptile.txt', sep='\t')

# Create tissue-time-condition factor.
rse.clp <- com_factor(rse.clp, tar, factors2com=c('organism_part', 'age', 'con'), factor.new='samTimeCon')
colData(rse.clp)[1:3, ]

# Normalize data.
se.nor.clp <- norm_data(data=rse.clp, norm.fun='ESF')
# Aggregate data.
se.aggr.clp <- aggr_rep(data=se.nor.clp, sam.factor='samTimeCon', con.factor=NULL, aggr='mean')
# Filter genes with low counts and low variance.
se.fil.clp <- filter_data(data=se.aggr.clp, sam.factor='samTimeCon', con.factor=NULL, pOA=c(0.07, 7), CV=c(0.7, 100), dir=NULL)
# Export gene expression matrix for use in the Shiny app.
write.table(assay(se.fil.clp), 'expr_coleoptile_samTimeCon.txt', sep='\t')


## Toy data in help files.
# Download the chicken data.
rse.chk <- read_cache(cache.pa, 'rse.chk') # Read data from cache.
if (is.null(rse.chk)) { # Save downloaded data to cache if it is not cached.
  rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.chk)
}
# Targets file.
chk.tar <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
target.chk <- read.table(chk.tar, header=TRUE, row.names=1, sep='\t')
target.chk[1:3, ]
colData(rse.chk) <- DataFrame(target.chk)
cna <- colnames(rse.chk)
# Filter raw counts.
se.chk <- filter_data(data=rse.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.05, 50), CV=c(3, 100), dir=NULL)
count.chk <- assay(se.chk)
# Toy data1: subset 2 tissues, each under 3 times, each tissue-time has 2 replicates.
count.chk.simple <- count.chk[, grepl('brain|heart', colnames(count.chk)) & grepl('day0|day70|day155',colnames(count.chk))]
count.chk.simple <- count.chk.simple[, sort(colnames(count.chk.simple))]
count.chk.simple <- count.chk.simple[, seq(1, ncol(count.chk.simple), by=2)]
count.chk.simple[1:2, ]
write.table(count.chk.simple, 'count_chicken_simple.txt', col.names=TRUE, row.names=TRUE, sep='\t')

# Toy data2.
colnames(count.chk) <- cna
write.table(count.chk, 'count_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')


# US map
# Data source "National, State, and Puerto Rico Commonwealth Totals Datasets: Population, population change, and estimated components of population change: April 1, 2010 to July 1, 2017" (nst-est2017-alldata.csv): https://www.census.gov/content/census/en/data/tables/2017/demo/popest/state-total.html; "Population, population change, and estimated components of population change: April 1, 2010 to July 1, 2018 (NST-EST2018-alldata)" (nst-est2018-alldata.csv): https://www.census.gov/data/datasets/time-series/demo/popest/2010s-state-total.html.

df <- read.table("nst-est2018-alldata1.txt", header=TRUE, row.names=NULL, sep="\t", fill=TRUE); rownames(df) <- df$NAME

year <- function(year) {

  df1 <- df[, -c(1:7)]; df1 <- df1[, grep(year, colnames(df1))]
  df1 <- df1[, c(1, 3, 4, 6, 7)] 
  cna <- c("Population_estimate", "Births", "Deaths", "International_migration", "Domestic_migration"); colnames(df1) <- cna
  df2 <- df1[2:5, ]; rownames(df2) <- paste0(gsub(" ", "_", rownames(df2)), "__", year); return(t(df2))

}

year <- function(year) {

  df1 <- df[, -c(1:7)]
  col.idx <- !grepl("^R", colnames(df1)) & grepl(paste0(year, "$"), colnames(df1)); df1 <- df1[, col.idx]; colnames(df1) <- sub(year, "", colnames(df1)); colnames(df1) <- sub("_", "", colnames(df1))
  df2 <- df1[2:5, ]; rownames(df2) <- paste0(gsub(" ", "_", rownames(df2)), "__", year); return(t(df2)[-1, ])

}

y2010 <- year(2010)
y2011 <- year(2011)
y2012 <- year(2012)

all <- data.frame(y2010, y2011, y2012, Annotation=c("Numeric change in resident total population", "Births", "Deaths", "Natural increase", "Net international migration", "Net domestic migration", "Net migration"), stringsAsFactors=FALSE)
write.table(all, "us_population2018.txt", col.names=TRUE, row.names=TRUE, sep="\t")




