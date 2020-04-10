
## Make targets file of the Arabidopsis shoot example.
  library(GEOquery)
  # Download GEO dataset GSE14502 (normalised by RMA).
  gset <- getGEO("GSE14502", GSEMatrix =TRUE, getGPL=FALSE)[[1]]
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


## Make target file-human brain example.
# Search for human samples related to brain.
all.res <- searchAtlasExperiments(properties="brain", species="Homo sapiens")
# Select "E-GEOD-67196", which is RNA-seq data of cerebellum and frontal cortex in brain.
rse <- getAtlasData('E-GEOD-67196')[[1]][[1]]

# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse)); df.con[1:3, ]

# The matching tissues in the SVG template are "cerebellum" and "frontal_cortex" while in the expression data are "cerebellum" and "frontal cortex". The tissue names must be exactly same between the SVG and the expression matrix, so change "frontal cortex" to "frontal_cortex".
df.con$organism_part <- sub(' ', '_', df.con$organism_part)
# Use abbreviations of the conditions ("individual", "disease", "genotype").
df.con$individual <- sub('individual', 'indi', df.con$individual)
df.con$individual <- sub('Individual ', 'indi', df.con$individual)
df.con$disease <- sub('amyotrophic lateral sclerosis', 'ALS', df.con$disease)
df.con$genotype <- sub('presence of a C9orf72 repeat expansion', 'pre.C9orf72', df.con$genotype)
df.con$genotype <- sub('absence of a C9orf72 repeat expansion', 'ab.C9orf72', df.con$genotype)
df.con$genotype <- sub('not applicable', 'NA', df.con$genotype)

write.table(df.con, 'target_brain.txt', col.names=TRUE, row.names=TRUE, sep='\t')

## Make target file-mouse organ example.

# Search for mouse samples related to heart.
all.mus <- searchAtlasExperiments(properties="heart", species="Mus musculus")
# Select "E-MTAB-2801".
all.mus[7:8, ]
rse.mus <- getAtlasData('E-MTAB-2801')[[1]][[1]]

# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.mus)); df.con[1:3, ]
df.con$organism_part <- gsub('skeletal muscle tissue', 'skeletal.muscle', df.con$organism_part)
df.con$strain <- make.names(df.con$strain)

write.table(df.con, 'target_mouse.txt', col.names=TRUE, row.names=TRUE, sep='\t')

## Make target file-chicken organ example.

# Search for chicken samples related to gallus.
all.chk <- searchAtlasExperiments(properties="", species="gallus")

# Select E-MTAB-6769.
all.chk[11:12, ]
rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]

# Tissue/condition metadata.
df.con <- as.data.frame(colData(rse.chk)); df.con[1:3, ]
df.con$age <- gsub('(.*) (.*)', '\\2\\1', df.con$age)

write.table(df.con, 'target_chicken.txt', col.names=TRUE, row.names=TRUE, sep='\t')













