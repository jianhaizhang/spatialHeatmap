df_match <- function() {
  SVGBulk <- c('NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO', 'LRC', 'LRC', 'LRC', 'LRC', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'QC', 'QC', 'QC', 'HAIR', 'STELE', 'STELE', 'STELE', 'STELE')
  
  cell <- c('atricho', 'colu.dist.colu', 'colu.dist.lrc', 'colu.proxi.colu', 'colu.proxi.lrc', 'colu', 'cortex', 'cortex.dist.lrc', 'endo', 'lat.rt.cap.dist.colu', 'lat.rt.cap.dist.lrc', 'lat.rt.cap.proxi.lrc', 'lat.rt.cap', 'metaphlo.comp.cell', 'metaxylem', 'phlo.po.per', 'procam', 'protophlo', 'protoxylem', 'protoxylem.dist.lrc', 'quies.cent', 'put.quies.cent', 'stem.niche', 'tricho', 'xylem.po.per', 'xylem', 'per', 'phlo')

  trueBulk <-c('NONHAIR,LRC_NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO,ENDO_QC', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'ENDO_QC,QC', 'ENDO_QC,QC', 'ENDO_QC,QC', 'HAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE')

  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}
df_match <- function() {
  SVGBulk <- c('NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO', 'LRC', 'LRC', 'LRC', 'LRC', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'QC', 'QC', 'QC', 'HAIR', 'STELE', 'STELE', 'STELE', 'STELE')
  
  cell <- c('atricho', 'colu.dist.colu', 'colu.dist.lrc', 'colu.proxi.colu', 'colu.proxi.lrc', 'colu', 'cortex', 'cortex.dist.lrc', 'endo', 'lat.rt.cap.dist.colu', 'lat.rt.cap.dist.lrc', 'lat.rt.cap.proxi.lrc', 'lat.rt.cap', 'metaphlo.comp.cell', 'metaxylem', 'phlo.po.per', 'procam', 'protophlo', 'protoxylem', 'protoxylem.dist.lrc', 'quies.cent', 'put.quies.cent', 'stem.niche', 'tricho', 'xylem.po.per', 'xylem', 'per', 'phlo')

  trueBulk <-c('NONHAIR,LRC_NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO,ENDO_QC', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'ENDO_QC,QC', 'ENDO_QC,QC', 'ENDO_QC,QC', 'HAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE')

  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}
