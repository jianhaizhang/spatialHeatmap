## Example data for manual matching.
library(spatialHeatmap); library(scRNAseq); library(SummarizedExperiment)
sce.mus <- MarquesBrainData()
# Filter cells and genes to obtain example data.
sce.mus.fil <- filter_cell(sce=sce.mus, cutoff=1, p.in.cell=0.85, p.in.gen=0.5)

# Edit colData.
cdat <- colData(sce.mus.fil)[, -1]; colnames(cdat) <- make.names(colnames(cdat))
cdat$source_name <- make.names(cdat$source_name)
colnames(cdat)[colnames(cdat)=='source_name'] <- 'label'
# "variable" is a reserved colname to indicate experiment variables.
colnames(cdat)[colnames(cdat)=='treatment'] <- 'variable'
rownames(cdat) <- make.names(rownames(cdat))

cdat <- edit_tar(cdat, 'variable', '6hr post acute stress', '6h.post.stress')
cdat <- edit_tar(cdat, 'variable', 'No', 'control')
cdat[1:2, ]
colData(sce.mus.fil) <- cdat
sce.mus.fil <- subset(sce.mus.fil, , variable=='control')
# Save the example data.
saveRDS(sce.mus.fil, file='./cell_mouse_brain.rds')

# Quality control, normalization, dimensionality reduction.
sce.dimred <- process_cell_meta(sce.mus.fil, qc.metric=list(subsets=list(Mt=rowData(sce.mus.fil)$featureType=='mito'), threshold=1)) 
# Clustering.
sce.clus <- cluster_cell(sce=sce.dimred, graph.meth='knn', dimred='PCA') 
# Manual cluster labels.
df.clus.mus.sc <- data.frame(cell=rownames(colData(sce.clus)), cluster=colData(sce.clus)$cluster)
write.table(df.clus.mus.sc, 'manual_cluster_mouse_brain.txt', col.names=TRUE, sep='\t')


## Example data for auto-matching/coclustering.

library(spatialHeatmap); source('fun_cocluster.R')
library(SummarizedExperiment); library(SingleCellExperiment)
# Download mouse brain bulk data (bulk_mouse_brain.xls) at https://github.com/jianhaizhang/cocluster_data/tree/master/validate/mouse_brain. (Male mice, WT 30-day old).

# Obtain example data by harsh filtering.
blk.mus.brain <- read.table('bulk_mouse_brain.xls', header=TRUE, sep='\t', row.names=1) 
blk.mus.brain <- filter_data(data=blk.mus.brain, pOA=c(0.3, 6), CV=c(0.55, 100)); dim(blk.mus.brain)
 
# Import single cell data of mouse brain according to instructions at https://github.com/jianhaizhang/cocluster_data. (Male mice strain C57/BL6, 9 weeks old). 
sc.mus.brain <- sc_dat_mus_brain(sc.pa='GSE147747_expr_raw_counts_table.tsv', meta.pa='GSE147747_meta_table.tsv')

# Obtain example data by harsh filtering, and Take overlap genes between bulk and single cells. 
mus.brain <- filter_cell(sce=sc.mus.brain, bulk=blk.mus.brain, gen.rm=NULL, cutoff=1, p.in.cell=0.93, p.in.gen=0.5) 

# Example bulk and single cell data. 
blk.mus <- mus.brain$bulk; sc.mus <- mus.brain$cell 

# Mouse brain aSVG.
svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", package="spatialHeatmap")
feature.df <- return_feature(svg.path=svg.mus.brain.pa)
df.match.mus.brain <- df_match_mus533()
# Bulk tissues are named with aSVG features.
cvt_vecter <- spatialHeatmap:::cvt_vector
colnames(blk.mus) <- cvt_vecter(df.match.mus.brain$trueBulk, df.match.mus.brain$SVGBulk, colnames(blk.mus))
blk.mus[1:3, ]
blk.mus$tissue <- colnames(blk.mus)
saveRDS(blk.mus, file='bulk_mouse_cocluster.rds') 

sc.mus[1:3, 1:5]
intersect(colnames(sc.mus), feature.df$feature)
intersect(colnames(sc.mus), colnames(blk.mus))
sc.mus$cell <- colnames(sc.mus)
saveRDS(sc.mus, file='cell_mouse_cocluster.rds') 

# Example data for auto-matching in Shiny app. 
library(SingleCellExperiment) 
blk.mus$label <- colnames(blk.mus) 
blk.mus$bulkCell <- 'bulk'; blk.mus$tissue <- NULL 
sc.mus$label <- colnames(sc.mus) 
sc.mus$bulkCell <- 'cell'; sc.mus$cell <- NULL 
sce.all <- cbind(blk.mus, sc.mus) 

# Secondary label of clusters.  
sce.nor <- norm_cell(sce.all, count.kp=TRUE) 
sce.dim <- reduce_dim(sce.nor) 
sce <- cluster_cell(sce.dim) 
colnames(colData(sce))[1] <- 'label1' 
cdat <- colData(sce) 
cdat.na <- colnames(cdat)
cdat.na <- c(c('label', 'label1'), setdiff(cdat.na, c('label', 'label1'))) 
colData(sce) <- cdat[, cdat.na] 
sce$sizeFactor <- NULL
assays(sce)$logcounts <- NULL 
sce$variable <- 'control' 
# Bulk tissue labels should always be aSVG features. 
blk.idx <- sce$bulkCell %in% 'bulk' 
sce$label1[blk.idx] <- sce$label[blk.idx] 


library(org.Mm.eg.db)
columns(org.Mm.eg.db); keytypes(org.Mm.eg.db)
row.met <- select(org.Mm.eg.db, keys=rownames(sce), columns=c('SYMBOL', 'GENENAME'), keytype="SYMBOL")
row.met <- row.met[!duplicated(row.met$SYMBOL), ]
row.met$metadata <- row.met$GENENAME
all(rownames(sce)==row.met$SYMBOL)
row.met[, c('SYMBOL', 'GENENAME')] <- NULL
rowData(sce) <- row.met

# sce <- readRDS('../extdata/shinyApp/data/shiny_covis_bulk_cell_mouse_brain.rds')

sce <- cvt_id(db='org.Mm.eg.db', data=sce, from.id='SYMBOL', to.id='ENSEMBL', desc=TRUE)
rdat <- rowData(sce)
# rdat <- subset(rdat, !is.na(SYMBOL) & !is.na(ENSEMBL) & !is.na(GENENAME))
rdat$link <- paste0('https://useast.ensembl.org/Mus_musculus/Gene/Summary?g=', rdat$ENSEMBL)
rdat$metadata <- rdat$desc
rdat$to.id <- rdat$GENENAME <- NULL
rowData(sce) <- rdat; rownames(sce) <- rdat$SYMBOL

df.meta <- data.frame(name=c('assayBulk', 'assayCell', 'aSVG'), link=c('https://www.ncbi.nlm.nih.gov/bioproject/PRJNA725533', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147747', 'https://raw.githubusercontent.com/ebi-gene-expression-group/anatomogram/master/src/svg/mus_musculus.male.svg'), species=c('Mus musculus', 'Mus musculus', 'Mus musculus'), technology=c('RNA-seq', 'RNA-seq', 'NA'), database=c('BioProject', 'GEO', 'EBI anatomogram'), id=c('PRJNA725533', 'GSE147747', 'NA'), description=c('Long-term impact of placental allopregnanolone insufficiency on the mouse transcriptome of cerebral cortex, hippocampus, hypothalamus and cerebellum', 'Hybridized 75 coronal sections from one brain hemisphere of 3 adult mice with Spatial Transcriptomics and defined a molecular atlas using clustering algorithms.', 'An annotated SVG (aSVG) file.'))
metadata(sce)$df.meta <- df.meta
saveRDS(sce, file='shiny_covis_bulk_cell_mouse_brain.rds')

# Tabular files. 
write.table(assay(sce), 'shiny_covis_bulk_cell_mouse_brain_assay.txt', col.names=TRUE, sep='\t')
write.table(colData(sce), 'shiny_covis_bulk_cell_mouse_brain_targets.txt', col.names=TRUE, sep='\t')
write.table(rowData(sce), 'shiny_covis_bulk_cell_mouse_brain_rowdata.txt', col.names=TRUE, sep='\t')
write.table(metadata(sce)$df.meta, 'shiny_covis_bulk_cell_mouse_brain_metadata.txt', col.names=TRUE, sep='\t')

## Example data of Arabidopsis thaliana root for coclustering optimization, downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766. 

# Access Arabidopsis root bulk data according to instructions at https://github.com/jianhaizhang/cocluster_data. 
blk.arab.rt <- bulk_dat(rds='blk_data.rds', cnt.path='RNASeq_counts.csv', meta.path='ici_metadata.csv', marker.blk.path='cell_type_marker_bulk.txt')
blk.arab.rt[1:3, 1:5]

# Matching table between bulk and single cells.
df.match <- df_match()

# Single cell data sc10.
sc.arab.rt10 <- scell_dat('GSM4625995_sc_10_at_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt10[1:3, 1:5])
# Single cell data sc11.
sc.arab.rt11 <- scell_dat('GSM4625996_sc_11_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt11[1:3, 1:5])

blk.arab.rt <- SingleCellExperiment(assays=list(counts=as(as.matrix(blk.arab.rt), 'dgCMatrix')))
sc.arab.rt10 <- SingleCellExperiment(assays=list(counts=as(as.matrix(sc.arab.rt10), 'dgCMatrix')))
sc.arab.rt11 <- SingleCellExperiment(assays=list(counts=as(as.matrix(sc.arab.rt11), 'dgCMatrix')))

# Inital filtering before normalization. 
blk <- filter_data(data=blk.arab.rt, pOA=c(0.2, 15), CV=c(1.5, 100)); dim(blk)
saveRDS(blk, file='bulk_cocluster.rds')

fil.init <- filter_cell(lis=list(sc10=sc.arab.rt10, sc11=sc.arab.rt11), bulk=blk, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.4, p.in.gen=0.1)
saveRDS(fil.init$sc10, file='sc10_cocluster.rds')
saveRDS(fil.init$sc11, file='sc11_cocluster.rds')

