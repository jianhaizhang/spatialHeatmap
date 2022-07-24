# Read bulk tissue data of Arabidopsis root.
bulk_dat <- function(rds, cnt.path, meta.path, marker.blk.path, dev.zone=FALSE) {
  library(readr)
  if (file.exists(rds)) { bulk <- readRDS(rds); return(bulk) } 
  met <- read.table(meta.path, header=TRUE, sep=',')
  # Metadata of bulk data.
  met.li <- subset(met, Reference=='Li et al. 2016')
  w <- which(met.li$Filename %in% c('SRR3664374_counts.txt', 'SRR3664375_counts.txt', 'SRR3664437_counts.txt'))
  met.li[w, 'Marker'] <- 'COR'
  marker <- met.li$Marker; names(marker) <- met.li$Filename
  marker <- marker[!marker %in% c('Elongation', 'Maturation', 'Meristematic', 'Whole_Root')]
  marker <- marker[sort(names(marker))]

  # Bulk data
  cnt <- as.data.frame(read_csv(cnt.path))
  rownames(cnt) <- cnt$Locus
  cnt.li <- cnt[, colnames(cnt) %in% met.li$Filename]
  cnt.li <- cnt.li[, names(marker)]
  colnames(cnt.li) <- marker

  # Bulk data sample names.
  marker.bulk <- read.table(marker.blk.path, header=TRUE, sep='\t')
  colnames(marker.bulk) <- make.names(marker.bulk[1, ])
  marker.bulk <- marker.bulk[-1, ]

  short.na <- marker.bulk$short.name
  names(short.na) <- marker.bulk$marker.name
  short.na <- short.na[marker]
  # Assign sample names to bulk data.
  colnames(cnt.li) <- short.na[colnames(cnt.li)]
  # Should not re-order columns by sorted column names, since duplicated values are created even though use as.matix to retain identical column names.
  # cnt.li <- as.matrix(cnt.li)[, sort(colnames(cnt.li))]
  cna.all <- colnames(cnt.li) 
  if (dev.zone==FALSE) cna.all <- sub('^MAT_|^DEV_', '', cna.all)
  vec1 <- c('COLUMELLA', 'PERICYCLE'); vec2 <- c('COLU', 'PERI')
  for (i in seq_along(vec1)) cna.all <- sub(vec1[i], vec2[2], cna.all)
  colnames(cnt.li) <- cna.all
  saveRDS(cnt.li, file=rds); return(cnt.li)
}

# Read single cell data.
scell_dat <- function(rds, df.match, dev.zone=FALSE) {
  library(Matrix); rds <- readRDS(rds)
  # Single-cell metadata from processed data.
  met.sc <- rds@meta.data
  spl.sct <- rds@assays$spliced_SCT
  spl.sct.cnt <- spl.sct@counts
  # Abbreviate cell anno.
  met.sc$cell.anno <- make.names(met.sc$celltype.anno)
  cell.na <- c('lat.rt.cap', 'procam', 'endo', 'tricho', 'xylem.po.per', 'phlo.po.per', 'atricho', 'colu', 'protoxylem', 'cortex', 'metaphlo.comp.cell', 'metaxylem', 'quies.cent', 'protophlo', 'per', 'phlo', 'xylem', 'put.quies.cent', 'stem.niche')

  names(cell.na) <- c('Lateral.Root.Cap', 'Procambium', 'Endodermis', 'Trichoblast', 'Xylem.Pole.Pericycle', 'Phloem.Pole.Pericycle', 'Atrichoblast', 'Columella', 'Protoxylem', 'Cortex', 'Metaphloem...Companion.Cell', 'Metaxylem', 'Quiescent.Center', 'Protophloem', 'Pericycle', 'Phloem', 'Xylem', 'Putative.Quiescent.Center', 'Stem.Cell.Niche')
  cell.ann.uni <- unique(met.sc$cell.anno)
  id0 <- cell.ann.uni[which(!cell.ann.uni %in% names(cell.na))]
  if (length(id0)>0) { print(id0); stop('New cell names are detected!') }
  for (i in seq_along(cell.na)) {
    cell.anno <- met.sc$cell.anno
    met.sc$cell.anno[cell.anno==names(cell.na[i])] <- cell.na[i]
  }
  colnames(spl.sct.cnt) <- met.sc$cell.anno

  # Abbreviate zones.
  met.sc$zone <- make.names(met.sc$time.anno)
  zone <- c('matur', 'elong', 'meri', 'proxi.lrc', 'proxi.colu', 'dist.lrc', 'dist.colu')
  names(zone) <- make.names(c('Maturation', 'Elongation', 'Meristem', 'Proximal Lateral Root Cap', 'Proximal Columella', 'Distal Lateral Root Cap', 'Distal Columella'))
  for (i in seq_along(zone)) {
    zone.all <- met.sc$zone
    met.sc$zone[zone.all==names(zone[i])] <- zone[i]
  }
  met.sc <- met.sc[, c('orig.ident', 'celltype.anno', 'cell.anno', 'time.anno', 'zone')]
  met.sc$cell.zone <- paste0(make.names(met.sc$cell.anno), '.', met.sc$zone)

  # Transfer cell ids from processed to raw data.
  # Leave out zones.
  # cell.ids <- make.names(met.sc$cell.anno)
  # Assign zones to cells.
  cell.ids <- make.names(met.sc$cell.zone)
  names(cell.ids) <- rownames(met.sc)

  sc.all <- rds@assays$spliced_RNA@counts+rds@assays$unspliced_RNA@counts
  colnames(sc.all) <- sub('_1.*', '', colnames(sc.all))
  colnames(sc.all) <- cell.ids[colnames(sc.all)]
  sc.all <- sc.all[, !is.na(colnames(sc.all))]
  if (dev.zone==FALSE) colnames(sc.all) <- sub('\\.elong$|\\.matur$|\\.meri$', '', colnames(sc.all)) 
  # dgCMatrix: columns should not be re-ordered by sorted column names. 
  # sc.all <- sc.all[, sort(colnames(sc.all))]
  na.uni <- unique(colnames(sc.all))
  no <- na.uni[!na.uni %in% df.match$cell]
  if (length(no) > 0) stop(paste0("These cells are not in the 'df.match': ", paste0(no, collapse=', '), "!"))
  return(sc.all)
}
# Matching table of Arabidopsis root bulk and single cells.
df_match <- function() {
  SVGBulk <- c('NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO', 'LRC', 'LRC', 'LRC', 'LRC', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'QC', 'QC', 'QC', 'HAIR', 'STELE', 'STELE', 'STELE', 'STELE')
  
  cell <- c('atricho', 'colu.dist.colu', 'colu.dist.lrc', 'colu.proxi.colu', 'colu.proxi.lrc', 'colu', 'cortex', 'cortex.dist.lrc', 'endo', 'lat.rt.cap.dist.colu', 'lat.rt.cap.dist.lrc', 'lat.rt.cap.proxi.lrc', 'lat.rt.cap', 'metaphlo.comp.cell', 'metaxylem', 'phlo.po.per', 'procam', 'protophlo', 'protoxylem', 'protoxylem.dist.lrc', 'quies.cent', 'put.quies.cent', 'stem.niche', 'tricho', 'xylem.po.per', 'xylem', 'per', 'phlo')

  dataBulk <-c('NONHAIR,LRC_NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO,ENDO_QC', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'ENDO_QC,QC', 'ENDO_QC,QC', 'ENDO_QC,QC', 'HAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE')

  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, dataBulk=dataBulk)
  return(df.match)
}

# Import mouse brain bulk raw count data generated by systemsPipeR. 
blk_dat_mus <- function(pa=NULL) {
  library(data.table)
  blk.mus <- fread(pa); rna <- blk.mus$V1
  blk.mus <- as(blk.mus[, -1], 'matrix')
  rownames(blk.mus) <- rna; cna <- toupper(colnames(blk.mus))
  if (any(grepl('_\\d+$', cna))) colnames(blk.mus) <- sub('_\\d+$', '', cna) else if (any(grepl('\\d+$', cna))) colnames(blk.mus) <- sub('\\d+$', '', cna) 
  return(blk.mus)
}
# Mouse brain single cell data.
sc_dat_mus_brain <- function(sc.pa=NULL, meta.pa=NULL) {
  library(data.table)
  sc.brain <- fread(sc.pa); sc.brain.met <- fread(meta.pa)
  # Check matching between counts and metadata.
  setkey(sc.brain, V1); setkey(sc.brain.met, V1)
  all(sc.brain$V1==sc.brain.met$V1)
  sc.brain <- sc.brain[sc.brain.met$passed_QC==TRUE]
  sc.brain.met <- sc.brain.met[passed_QC==TRUE]
  all(sc.brain$V1==sc.brain.met$V1)
  sc.brain$V1 <- make.names(sc.brain$V1)
  sc.brain.met$V1 <- make.names(sc.brain.met$V1)

  # Transpose count table.
  rna <- sc.brain$V1; sc.brain <- as(sc.brain[, -1], 'matrix')
  rownames(sc.brain) <- rna; sc.brain <- t(sc.brain)

  # Abbreviate cell names in metadata.
  sc.brain.met$ABA_parent <- make.names(sc.brain.met$ABA_parent)
  cells <- c(Cerebellum='cere', fiber.tracts='fiber.tr', Hippocampal.region='hipp', Isocortex='isocort', Olfactory.areas='olfa', Retrohippocampal.region='retrohipp', Thalamus='thalamus', ventricular.systems='ventri', Cortical.subplate='corti.sub', Hindbrain='hindbrain', Hypothalamus='hypotha', Midbrain='midbrain', Pallidum='pallidum', Striatum='striatum', Undefined.areas='undefined')

  sc.brain.met$cell <- cells[sc.brain.met$ABA_parent]
  # Assign cells names in metadata to count table.
  all(colnames(sc.brain)==sc.brain.met$V1)
  colnames(sc.brain) <- sc.brain.met$cell; return(sc.brain)
}
