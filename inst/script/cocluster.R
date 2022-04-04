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


## Example data for auto-matching/coclustering.

library(spatialHeatmap); source('fun_cocluster.R')
# Download mouse brain bulk data (blk.mus.brain) at https://github.com/jianhaizhang/cocluster_data.                                
                                                                                                                                   
# Obtain example data by hash filtering.                                                                                           
blk.mus.brain <- filter_data(data=, pOA=c(0.3, 6), CV=c(0.6, 100)); dim(blk.mus.brain)                                             
                                                                                                                                   
# Import single cell data of mouse brain according to instructions at https://github.com/jianhaizhang/cocluster_data.              
sc.mus.brain <- sc_dat_mus_brain(sc.pa=file.path(sc.mus.brain.pa, 'GSE147747_expr_raw_counts_table.tsv'), meta.pa=file.path(sc.mus.
brain.meta.pa, 'GSE147747_meta_table.tsv'))                                                                                        
                                                                                                                                   
# Obtain example data by hash filtering.                                                                                           
sc.mus.brain <- filter_cell(lis=list(sc.mus.brain=sc.mus.brain), bulk=NULL, gen.rm=NULL, min.cnt=1, p.in.cell=0.85, p.in.gen=0.5)  
                                                                                                                                   
# Take overlap genes between bulk and single cells.                                                                                
mus.brain <- filter_cell(lis=list(sc.mus.brain=sc.mus.brain$sc.mus.brain), bulk=blk.mus.brain, gen.rm=NULL, min.cnt=0, p.in.cell=0.
5, p.in.gen=0.1)                                                                                                                   
                                                                                                                                   
# Example bulk and single cell data.                                                                                               
blk.mus <- mus.brain$bulk; sc.mus <- mus.brain$sc.mus.brain                                                                        
write.table(blk.mus, 'bulk_mouse_cocluster.txt', col.names=TRUE, row.names=TRUE, sep='\t')                                         
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

