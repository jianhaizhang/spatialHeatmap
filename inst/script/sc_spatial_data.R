library(Seurat); library(SeuratData); library(SummarizedExperiment); library(spatialHeatmap)

# Bulk data.
blk.mus.pa <- system.file("extdata/shinyApp/data", "bulk_mouse_cocluster.rds", package="spatialHeatmap")
blk.mus <- readRDS(blk.mus.pa)
assay(blk.mus)[1:3, 1:5]

# Spatial single-cell data.
if (!'stxBrain' %in% InstalledData()$Dataset) InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
# as.SingleCellExperiment(brain, assay='SCT')
# p <- SpatialFeaturePlot(brain, features = c("Hpca"))

# Joint normalization.
nor.lis <- norm_srsc(cell=brain, assay='Spatial', bulk=blk.mus)
# Separate bulk and cell data.
srt.sc <- nor.lis$cell; bulk <- nor.lis$bulk
blk.aggr <- aggr_rep(data=bulk, assay.na='logcounts', sam.factor='sample', aggr='mean')
saveRDS(blk.aggr, file='bulk_sp.rds')

# Dimension reduction.
srt.sc <- RunPCA(srt.sc, assay = "SCT", verbose = FALSE)
srt.sc <- RunUMAP(srt.sc, assay = "SCT", dims = 1:5)
srt.sc <- RunTSNE(srt.sc, assay = "SCT", reduction = "pca", dims = 1:5)
# Clustering.
srt.sc <- FindNeighbors(srt.sc, reduction = "pca", dims = 1:30)
srt.sc <- FindClusters(srt.sc, verbose = FALSE)
srt.sc$seurat_clusters <- paste0('clus', srt.sc$seurat_clusters)
saveRDS(srt.sc, file='srt_sc.rds')


