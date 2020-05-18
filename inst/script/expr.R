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


