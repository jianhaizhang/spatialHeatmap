## Example data for manual matching.
library(spatialHeatmap); library(scRNAseq)
sce.manual <- MarquesBrainData()
# Filter cells and genes to obtain example data.
sce.manual <- filter_cell(lis=list(sce=sce.manual), min.cnt=1, p.in.cell=0.9, p.in.gen=0.5)[[1]]

# Edit colData.
cdat <- colData(sce.manual)[, -1]; colnames(cdat) <- make.names(colnames(cdat))
cdat$source_name <- make.names(cdat$source_name)
colnames(cdat)[colnames(cdat)=='source_name'] <- 'label'
# "expVar" is a reserved colname to indicate experiment variables.
colnames(cdat)[colnames(cdat)=='treatment'] <- 'expVar'
rownames(cdat) <- make.names(rownames(cdat))

cdat <- edit_tar(cdat, 'expVar', '6hr post acute stress', '6h.post.stress')
cdat <- edit_tar(cdat, 'expVar', 'No', 'control')
cdat[1:2, ]
colData(sce.manual) <- cdat
# Save the example data.
saveRDS(sce.manual, file='./sce_manual_mouse.rds')

# Quality control, normalization, dimensionality reduction.
sce.dimred <- process_cell_meta(sce.manual, qc.metric=list(subsets=list(Mt=rowData(sce.manual)$featureType=='mito'), threshold=1)) 
# Clustering.
sce.clus <- cluster_cell(sce=sce.dimred, graph.meth='knn', dimred='PCA') 
# Manual cluster labels.
df.clus.mus.sc <- data.frame(cell=rownames(colData(sce.clus)), cluster=colData(sce.clus)$cluster)
write.table(df.clus.mus.sc, 'manual_cluster_mouse_brain.txt', col.names=TRUE, sep='\t')


## Example data for auto-matching/coclustering.

library(spatialHeatmap); source('fun_cocluster.R')
# Download mouse brain bulk data (bulk_mouse_brain.xls) at https://github.com/jianhaizhang/cocluster_data/tree/master/validate/mouse_brain.

# Obtain example data by harsh filtering.
blk.mus.brain <- read.table('bulk_mouse_brain.xls', header=TRUE, sep='\t', row.names=1) 
blk.mus.brain <- filter_data(data=blk.mus.brain, pOA=c(0.3, 6), CV=c(0.55, 100)); dim(blk.mus.brain)
 
# Import single cell data of mouse brain according to instructions at https://github.com/jianhaizhang/cocluster_data. 
sc.mus.brain <- sc_dat_mus_brain(sc.pa='GSE147747_expr_raw_counts_table.tsv', meta.pa='GSE147747_meta_table.tsv')

# Obtain example data by harsh filtering, and Take overlap genes between bulk and single cells. 
sc.mus.brain <- filter_cell(lis=list(sc.mus.brain=sc.mus.brain), bulk=blk.mus.brain, gen.rm=NULL, min.cnt=1, p.in.cell=0.93, p.in. 
gen=0.5) 

# Example bulk and single cell data. 
blk.mus <- mus.brain$bulk; sc.mus <- mus.brain$sc.mus.brain 
# blk.mus.pa <- system.file("extdata/shinyApp/example", "bulk_mouse_cocluster.txt", package="spatialHeatmap") 
# blk.mus <- as.matrix(read.table(blk.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
# blk.mus[1:3, ]
# Mouse brain aSVG.
svg.mus.brain.pa <- system.file("extdata/shinyApp/example", "mus_musculus.brain.svg", package="spatialHeatmap")
feature.df <- return_feature(svg.path=svg.mus.brain.pa)
# Matching table between bulk tissues and aSVG features.
match.mus.brain.pa <- system.file("extdata/shinyApp/example", "match_mouse_brain_cocluster.txt", package="spatialHeatmap")
df.match.mus.brain <- read.table(match.mus.brain.pa, header=TRUE, row.names=1, sep='\t')
df.match.mus.brain
# Bulk tissues are named with aSVG features.
cvt_vecter <- spatialHeatmap:::cvt_vector
colnames(blk.mus) <- cvt_vecter(df.match.mus.brain$dataBulk, df.match.mus.brain$SVGBulk, colnames(blk.mus))
blk.mus[1:3, ]
write.table(blk.mus, 'bulk_mouse_cocluster.txt', col.names=TRUE, row.names=TRUE, sep='\t') 

# sc.mus.pa <- system.file("extdata/shinyApp/example", "cell_mouse_cocluster.txt", package="spatialHeatmap") 
# sc.mus <- as.matrix(read.table(sc.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
sc.mus[1:3, 1:5]
intersect(colnames(sc.mus), feature.df$feature)
intersect(colnames(sc.mus), colnames(blk.mus))
write.table(sc.mus, 'cell_mouse_cocluster.txt', col.names=TRUE, row.names=TRUE, sep='\t') 

# Example data for auto-matching in Shiny app.  
library(SingleCellExperiment) 
blk.sc <- as.matrix(cbind(blk.mus, sc.mus)) 
cdat <- DataFrame(bulkCell=c(rep('bulk', ncol(blk.mus)), rep('cell', ncol(sc.mus)))) 
sce.all <- SingleCellExperiment(assays=list(counts=blk.sc), colData=cdat)
 
df.match.mus.pa <- system.file("extdata/shinyApp/example", "match_mouse_brain_cocluster.txt", package="spatialHeatmap") 
df.match <- read.table(df.match.mus.pa, header=TRUE, row.names=1, sep='\t')  
metadata(sce.all)$df.match <- df.match
saveRDS(sce.all, file='sce_auto_bulk_cell_mouse_brain.rds') 

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

