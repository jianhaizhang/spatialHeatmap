
## Shiny app data.
# Human brain data. 
# Access data.
library(ExpressionAtlas)
rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]
# Targets file.
brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
colData(rse.hum) <- DataFrame(target.hum)
# Normalise.
capture.output(se.nor.hum <- norm_data(se=rse.hum, method.norm='ratio', data.trans='log2'), file=tempfile())
# Aggregate.
se.aggr.hum <- aggr_rep(se=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
# Filter.
se.fil.hum <- filter_data(se=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.55, 100), dir=NULL)
# Data matrix.
expr.hum <- assay(se.fil.hum)
# colnames(expr.hum) <- gsub("_", ".", colnames(expr.hum))
write.table(expr.hum, 'expr_human.txt', col.names=TRUE, row.names=TRUE, sep='\t')

# Mouse organ.
rse.mus <- getAtlasData('E-MTAB-2801')[[1]][[1]]
pa.mus <- system.file('extdata/shinyApp/example/target_mouse.txt', package='spatialHeatmap')
# target.mus <- read.table(pa.mus, header=TRUE, row.names=1, sep='\t')
colData(rse.mus) <- DataFrame(target.mus)

capture.output(se.nor.mus <- norm_data(se=rse.mus, method.norm='ratio', data.trans='log2'), file=tempfile())
# Average the tissue-condition replicates.
se.aggr.mus <- aggr_rep(se=se.nor.mus, sam.factor='organism_part', con.factor='strain', aggr='mean')
se.fil.mus <- filter_data(se=se.aggr.mus, sam.factor='organism_part', con.factor='strain', pOA=c(0.15, 5), CV=c(1.7, 100), dir=NULL)
# Data matrix.
expr.mus <- assay(se.fil.mus)
# colnames(expr.mus) <- gsub("_", ".", colnames(expr.mus))
write.table(expr.mus, 'expr_mouse.txt', col.names=TRUE, row.names=TRUE, sep='\t')

# Chicken organ.
rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
pa.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
target.chk <- read.table(pa.chk, header=TRUE, row.names=1, sep='\t')
colData(rse.chk) <- DataFrame(target.chk)
capture.output(se.nor.chk <- norm_data(se=rse.chk, method.norm='ratio', data.trans='log2'), file=tempfile())
# Average the tissue-condition replicates.
se.aggr.chk <- aggr_rep(se=se.nor.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
se.fil.chk <- filter_data(se=se.aggr.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.05, 5), CV=c(1.5, 100), dir=NULL)
# Data matrix.
expr.chk <- assay(se.fil.chk)
# colnames(expr.chk) <- gsub("_", ".", colnames(expr.chk))
write.table(expr.chk, 'expr_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')

# Arabidopsis shoot. 
library(GEOquery)
gset <- getGEO("GSE14502", GSEMatrix=TRUE, getGPL=TRUE)[[1]]
se.sh <- as(gset, "SummarizedExperiment")
rownames(se.sh) <- make.names(rowData(se.sh)[, 'Gene.Symbol'])
pa.sh <- system.file('extdata/shinyApp/example/target_shoot.txt', package='spatialHeatmap')
target.sh <- read.table(pa.sh, header=TRUE, row.names=1, sep='\t')
colData(se.sh) <- DataFrame(target.sh)
# Average the sample-condition replicates.
se.aggr.sh <- aggr_rep(se=se.sh, sam.factor='samples', con.factor='conditions', aggr='mean')
se.fil.sh <- filter_data(se=se.aggr.sh, sam.factor='samples', con.factor='conditions', pOA=c(0.03, 6), CV=c(0.29, 100), dir=NULL)
# Data matrix.
expr.sh <- assay(se.fil.sh)
# colnames(expr.sh) <- gsub("_", ".", colnames(expr.sh))
expr.sh <- cbind.data.frame(expr.sh, ann=rowData(se.fil.sh)[, 'Target.Description'], stringsAsFactors=FALSE)
expr.sh <- expr.sh <- expr.sh[c(9, 1:8, 10:372), ]
write.table(expr.sh, 'expr_arab.txt', col.names=TRUE, row.names=TRUE, sep='\t')


## Toy data.
# Download the chicken data.
library(ExpressionAtlas)
rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
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


