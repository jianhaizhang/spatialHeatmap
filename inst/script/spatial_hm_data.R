library(affy); library(genefilter)

# Brain data.
library(hgu133a2.db)

# Data source GSE13162: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
data <- ReadAffy(celfile.path="GSE13162_RAW")
mas5 <- mas5(data); exp0 <- exprs(mas5); pma0 <- mas5calls(data)
exp <- log2(exp0); pma <- exprs(pma0)

# Sample title.
sam.title <- system("zcat GSE13162_series_matrix.txt.gz | grep Sample_title", intern=TRUE)
sam.title.sep <- unlist(strsplit(sam.title, "\"\t\""))
sam.title.sep[1] <- "N1A Normal-frontal"; sam.title.sep[56] <- "S10B Sporadic-cerebellum"
sam.title.sep[6] <- "N5B Normal-cerebellum"
colnames(pma) <- colnames(exp) <- gsub("(.*)( )(.*)", "\\3", sam.title.sep)

# Annotation.
keys <- keys(hgu133a2.db)
ann <- select(hgu133a2.db, keys=keys, columns = c("SYMBOL", "GENENAME"))
ann1 <- subset(ann, !grepl("AFFX|s_at|a_at|x_at|r_at", PROBEID))
# Remove duplicated probes and genes.
dup <- ann1[, 1][which(duplicated(ann1[, 1]))]
ann1 <- subset(ann1, !(ann1[, 1] %in% dup))
dup1 <- ann1[, 2][which(duplicated(ann1[, 2]))]
ann1 <- subset(ann1, !(ann1[, 2] %in% dup1))
rownames(ann1) <- ann1[, 1]
# Only keep probes with valid gene symbols.
ann1 <- subset(ann1, !is.na(SYMBOL))

# Only keep probes with annotations.
int <- intersect(rownames(exp), ann1[, 1])
# Use gene symbols as rownames.
rna <- make.names(ann1[int, ][, 2])
exp <- exp[int, ];  pma <- pma[int, ]; rownames(exp) <- rownames(pma) <- rna

# Keep genes having "P" in at least one sample.
pma[pma=="P"] <- 2; pma[pma=="M"] <- 1; pma[pma=="A"] <- 0
pma <- apply(pma, 2, as.numeric) 
pma.ag <- aggregate(t(pma), by=list(colnames(pma)), FUN=mean)
pma.ag <- t(pma.ag); colnames(pma.ag) <- pma.ag[1, ]; pma.ag <- pma.ag[-1, ]
pma.ag <- apply(pma.ag, 2, as.numeric); rownames(pma.ag) <- rna
ffun <- filterfun(pOverA(0.11, 1)); filtered <- genefilter(pma.ag, ffun)
pma.ag <- pma.ag[filtered, ]

# Filter genes according to CV.
exp1 <- exp[rownames(pma.ag), ]
exp.ag <- aggregate(t(exp1), by=list(colnames(exp1)), FUN=mean)
exp.ag <- t(exp.ag); colnames(exp.ag) <- exp.ag[1, ]
exp.ag <- exp.ag[-1, ]; exp.ag <- apply(exp.ag, 2, as.numeric)
rownames(exp.ag) <- rownames(exp1)
ffun <- filterfun(cv(0.15, 10000)); filtered <- genefilter(exp.ag, ffun)
exp.ag <- exp.ag[filtered, ]; dim(exp.ag)
cna<- sub("-", "__", colnames(exp.ag))
colnames(exp.ag) <- gsub("(.*)(__)(.*)", "\\3\\2\\1", cna)

# Append annotation.
rownames(ann1) <- make.names(ann1$SYMBOL)
int1 <- intersect(rownames(ann1), rownames(exp.ag))
ann1 <- ann1[int1, ]; exp.ag <- exp.ag[int1, ]
exp.ag <- as.data.frame(exp.ag)
exp.ag1 <- cbind(exp.ag, Annotation=ann1$GENENAME, stringsAsFactors=FALSE)
write.table(exp.ag1, "brain_expr_ann_row_gen.txt", col.names=TRUE, row.names=TRUE, sep="\t")

# Arabidopsis thaliana root data.
library(ath1121501.db)

# Data source GSE46205: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

data <- ReadAffy(celfile.path="GSE46205_RAW")
mas5 <- mas5(data); exp0 <- exprs(mas5); pma0 <- mas5calls(data)
exp <- log2(exp0); pma <- exprs(pma0)

# Sample title.
sam.title <- system("zcat GSE46205_series_matrix.txt.gz | grep Sample_title", intern=TRUE)
sam.title.sep <- unlist(strsplit(sam.title, "\"\t\""))
sam.title.sep[1] <- "Arabidopsis, root cells, epidermis, standard conditions for 1 hour, replicate 1"; sam.title.sep[96] <- "Arabidopsis, root cells, epidermis, standard conditions for 1 hour, replicate 2"
sam.title.sep <- gsub("\\,", "", sam.title.sep)
sam.title.sep <- sub("Arabidopsis root cells ", "", sam.title.sep)
sam.title.sep <- gsub(" conditions for | NaCl for ", "\\_", sam.title.sep)
sam.title.sep <- gsub(" hour replicate [0-9]| hours replicate [0-9]", "h", sam.title.sep)
sam.title.sep <- gsub(" ", "__", sam.title.sep)
colnames(pma) <- colnames(exp) <- sam.title.sep
idx <- grep("__standard_48h$", sam.title.sep)
pma <- pma[, -idx]; exp <- exp[, -idx]

# Annotation.
ann <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "), stringsAsFactors=FALSE)
rm1 <- grep("AFFX|s_at|a_at|x_at|r_at", rownames(ann)); rm2 <- which(ann[, "SYMBOL"]=="NA"); rm3 <- which(ann$DESC=="NA")
ann1 <- ann[-c(rm1, rm2, rm3), ]
# Remove duplicated probes and genes.
dup <- rownames(ann1)[which(duplicated(rownames(ann1)))]
ann1 <- subset(ann1, !(rownames(ann1) %in% dup))
dup1 <- ann1[, 2][which(duplicated(ann1[, 2]))]
ann1 <- subset(ann1, !(ann1[, 2] %in% dup1))

# Only keep probes with annotations.
int <- intersect(rownames(exp), rownames(ann1))
# Use gene symbols as rownames.
rna <- make.names(ann1[int, ][, 2]); cna <- unique(colnames(exp))
exp <- exp[int, ]; pma <- pma[int, ]; rownames(exp) <- rownames(pma) <- rna

# Keep genes having "P" in at least one sample.
pma[pma=="P"] <- 2; pma[pma=="M"] <- 1; pma[pma=="A"] <- 0
pma <- apply(pma, 2, as.numeric) 
pma.ag <- aggregate(t(pma), by=list(colnames(pma)), FUN=mean)
pma.ag <- t(pma.ag); colnames(pma.ag) <- pma.ag[1, ]; pma.ag <- pma.ag[-1, ]
pma.ag <- apply(pma.ag, 2, as.numeric); rownames(pma.ag) <- rna
ffun <- filterfun(pOverA(0.03, 1)); filtered <- genefilter(pma.ag, ffun)
pma.ag <- pma.ag[filtered, ]

# Filter genes according to CV.
exp1 <- exp[rownames(pma.ag), ]
exp.ag <- aggregate(t(exp1), by=list(colnames(exp1)), FUN=mean)
exp.ag <- t(exp.ag); colnames(exp.ag) <- exp.ag[1, ]; exp.ag <- exp.ag[, cna]
exp.ag <- exp.ag[-1, ]; exp.ag <- apply(exp.ag, 2, as.numeric)
rownames(exp.ag) <- rownames(exp1)
ffun <- filterfun(cv(0.23, 10000)); filtered <- genefilter(exp.ag, ffun)
exp.ag <- exp.ag[filtered, ]; dim(exp.ag)

# Append annotation.
rownames(ann1) <- make.names(ann1$SYMBOL)
int1 <- intersect(rownames(ann1), rownames(exp.ag))
ann1 <- ann1[int1, ]; exp.ag <- exp.ag[int1, ]
exp.ag <- as.data.frame(exp.ag)
exp.ag1 <- cbind(exp.ag, Annotation=ann1$DESC, stringsAsFactors=FALSE)
write.table(exp.ag1, "root_expr_ann_row_gen.txt", col.names=TRUE, row.names=TRUE, sep="\t")

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



