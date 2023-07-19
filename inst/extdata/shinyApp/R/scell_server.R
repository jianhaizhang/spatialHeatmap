# Module for co-visualization.
scell_server <- function(id, tab, upl.mod.lis, shm.mod.lis, lis.url, parent, session) {
  moduleServer(id, function(input, output, session) {
  cat('Module scell_server ... \n')
  ns <- session$ns;
  # Apply to switching from sce of bulk+cell to sce of cell alone.
  sce.upl <- reactiveValues()
  observeEvent(upl.mod.lis$sce$val, {
    library(Matrix)
    sce.all <- upl.mod.lis$sce$val
    reducedDim(sce.all, 'PCA') <- reducedDim(sce.all, 'UMAP') <- reducedDim(sce.all, 'TSNE') <- NULL
    rowData(sce.all) <- link_dat(rowData(sce.all), link.only=FALSE) 
    df.meta <- metadata(sce.all)$df.meta 
    if (!is.null(df.meta)) { 
      if (ncol(data.frame(df.meta))<2) {  
        lgc.mt <- FALSE; msg <- 'The "df.meta" in the "metadata" slot should be a "data.frame" with at least two columns!' 
        if (!lgc.mt) showModal(modal(msg = msg)); req(!lgc.mt)
      } else {
        metadata(sce.all)$df.meta <- link_dat(df.meta, link.only=FALSE)
      }  
    }
    # save(sce.all, file='sce.all')
    if (!is.null(sce.all)) {
      assay.na <- assayNames(sce.all)
      if ('logcounts' %in% assay.na) logcounts(sce.all) <- as(logcounts(sce.all), 'dgCMatrix')
      if ('count' %in% assay.na) assays(sce.all)$count <- as(assays(sce.all)$count, 'dgCMatrix') 
      if ('counts' %in% assay.na) assays(sce.all)$counts <- as(assays(sce.all)$counts, 'dgCMatrix') 
    }
    sce.upl$cell <- sce.all; if (!is.null(sce.all)) {
      cdat <- colData(sce.all); dat.all <- assay(sce.all)
      int <- all(round(dat.all)==dat.all)
      if (!int) sce.upl$cell <- NULL
      if (!int) showModal(modal(msg='In co-visualization, raw count data are needed !'))
      validate(need(int, '')); sce.upl$method <- 'manual'
      if ('bulkCell' %in% colnames(cdat)) {
        blk.cell.uni <- unique(cdat$bulkCell)
        if (all(c('bulk', 'cell') %in% blk.cell.uni)) {
          cdat.na <- colnames(cdat)
          lab.na <- grep('^label$|^label\\d+', cdat.na, value=TRUE)
          if (length(lab.na)==0) {
            show_mod(length(lab.na)>0, msg='At least one label column in "colData" slot is required, such as "label", "label1", "label2", ...')
          }; validate(need(length(lab.na) > 0, ''))
          sce.upl$method <- 'both'
          bulk <- subset(sce.all, , bulkCell=='bulk')
          # Bulk data are aggregated. 
          bulk <- aggr_rep(data=bulk, sam.factor=lab.na[1], aggr='mean')
          bulk$spFeature <- NULL
          assay(bulk) <- round(assay(bulk)); sce.upl$bulk <- bulk
          sce.upl$cell <- subset(sce.all, , bulkCell=='cell')
        }
      } else {
        updateSelectInput(session, 'methCovis', choices=c('Annotation/Manual'='man'), selected='man')
      }
    }
  })
  # observeEvent(sce.upl$cell, { 
  # output$methCovis <- renderUI({
  #  method <- sce.upl$method
  #  if(is.null(sce.upl$cell) | is.null(method)) return()
  #  ns <- session$ns; sel <- 'auto'
  #  cho <- c('Annotation/manual'='man', 'Automatic'='auto')
  #  if (method=='manual') {
  #    cho <- c('Annotation/manual'='man'); sel <- 'man'
  #  }
  #  selectInput(ns('methCovis'), label='Methods', choices=cho, selected=sel)
  # })
  #output$direc <- renderUI({
  #  method <- sce.upl$method; if(is.null(method)) return() 
  #  ns <- session$ns; sel <- 'toBulk'
  #  cho <- c('Cell2tissue'='toBulk', 'Tissue2cell'='toCell')
  #  if (method=='manual') {
  #    cho <- c('Cell2tissue'='toBulk'); sel <- 'toBulk'
  #  }
  #  selectInput(ns('direc'), label='Mapping direction', choices=cho, selected=sel)
  #})
  #})
  # Initiate three reactive values in the same line can cause problems, since covisGrp slots can be assigned to covis.auto if the former have slots while the latter has no slots.
  # covis.man <- covis.auto <- covisGrp <- reactiveValues()
  covisGrp <- reactiveValues()
  covis.man <- reactiveValues(); covis.auto <- reactiveValues(con.na=FALSE)
  observe({
    meth.covis <- input$methCovis; direc <- input$direc
    cell <- sce.upl$cell; covis.type <- NULL
    if (is.null(meth.covis)|is.null(direc)|is.null(cell)) return()
    if ('man' %in% meth.covis) {
      if (direc=='toBulk') covis.type <- 'toBulk'
      if (direc=='toCell') covis.type <- 'toCell'
    }
    if (meth.covis=='auto' & direc=='toBulk') covis.type <- 'toBulkAuto'
    if (meth.covis=='auto' & direc=='toCell') covis.type <- 'toCellAuto'
    sce.upl$method <- ifelse('man' %in% meth.covis, 'man', 'auto')
    sce.upl$covis.type <- covis.type
    if ('man' %in% sce.upl$method) { 
      showElement('covisMan'); hideElement('covisAuto')
    } else if ('auto' %in% sce.upl$method) {
      showElement('covisAuto'); hideElement('covisMan') 
    }
  })
  observeEvent(sce.upl$method, {
    method <- sce.upl$method; cell <- sce.upl$cell
    if (is.null(method)|is.null(cell)) return()
    # 1) if condition then x <- covis_man_server(arg1=reactiveValue1, ...): server is like a function that returns a value, 2) if condition then covis_man_server(arg1=reactiveValue1, ...): the condition is useless, since reactiveValue1 can cause execution on inner code of covis_man_server regardless of the condition. 
    if ('man' %in% method & is.null(covis.man$dimred)) {
      covis_man_server('covisMan', sce.upl=sce.upl, upl.mod.lis=upl.mod.lis, shm.mod.lis=shm.mod.lis, tab=tab, covis.man=covis.man, lis.url=lis.url, parent=parent, session)
    } else if ('auto' %in% method & is.null(covis.auto$res)) {
      covis_auto_server('covisAuto', sce.upl=sce.upl, upl.mod.lis=upl.mod.lis, shm.mod.lis=shm.mod.lis, tab=tab, covis.auto=covis.auto, lis.url=lis.url, parent=parent, session)
    } 
  })
  # Erase value from previous session.
  observeEvent(sce.upl$covis.type, { covisGrp$val <- NULL })  
  # "covisGrp" from scell section, passed to SHM page.
  observe({
    covis.type <- sce.upl$covis.type
    if (is.null(covis.man$covisGrp) & is.null(covis.auto$covisGrp)) return()
    if (is.null(covis.type)) return()
    if (any(c('toCell', 'toBulk') %in% covis.type)) { 
      covisGrp$val <- covis.man$covisGrp 
    } else if (any(c('toBulkAuto', 'toCellAuto') %in% covis.type)) {
      covisGrp$val <- covis.auto$covisGrp
    }
  })
  sce.res <- reactive({
    cat('Single cell: aggregating cells ... \n')
    if (!grepl(na.sgl, upl.mod.lis$ipt$fileIn)) return()
    covis.type <- sce.upl$covis.type; method <- sce.upl$method
    if (is.null(covis.type)|is.null(method)) return()
    sce.shm <- NULL
    if (covis.type %in% 'toBulk') sce.shm <- covis.man$dimred
    if (covis.type %in% 'toCell') sce.shm <- covis.man$bulk
    if (covis.type %in% 'toCellAuto') { 
      if (is.null(covis.auto$res)) return() 
      sce.shm <- subset(covis.auto$res, , bulkCell=='bulk')
      covisGrp$val <- 'sample'
    }
    if (covis.type %in% 'toBulkAuto') {
      if (is.null(covis.auto$res)) return() 
      sce.shm <- subset(covis.auto$res, , bulkCell=='cell')
      covisGrp$val <- 'assignedBulk'
    }
    if (is.null(sce.shm)) return()
    df.rep <- assay(sce.shm)

    # Column names: sp.ft, exp.var.
    sp.ft <- covisGrp$val; if (is.null(sp.ft)) return()
    cdat <- colData(sce.shm)
    if (!sp.ft %in% colnames(cdat)) return()
    # Auto-formed clusters should not be combined with variable, since all cells under variable are clustered. The cell labels are defined by users under control and condition independently.
    if ('variable' %in% colnames(cdat) & 'man' %in% method) colnames(df.rep) <- paste0(cdat[, sp.ft], '__', cdat[, 'variable']) else { colnames(df.rep) <- cdat[, sp.ft] }
    withProgress(message="Covis: ", value=0, {
      incProgress(0.3, detail="aggregating cells or bulk ...")
      sce.lis <- fread_df(cbind(as(df.rep, 'matrix'), rowData(sce.shm)), rep.aggr='mean')
      incProgress(0.3, detail="aggregating cells or bulk ...")
    })
    lis <- list(sce.lis=sce.lis, covisGrp=sp.ft)
    # save(lis, file='lis')
    cat('Done! \n'); return(lis)
  })
  cat('Module scell_server done! \n')
  onBookmark(function(state) { state })
  return(list(sce.upl=sce.upl, sce.res=sce.res, covis.man=covis.man, covis.auto=covis.auto, input=input))

})}
