assay_dat <- function(input, r1=1, r2=500, c1=1, c2=20) {
  if (r1 < 1) r1 <- 1; if (c1 < 1) c1 <- 1
  lgc.r <- r2 > r1; if (!lgc.r) {
    show_mod(lgc.r, msg='Row End should > Row Start!')
  }; req(lgc.r)  
  lgc.c <- c2 > c1; if (!lgc.c) {
    show_mod(lgc.c, msg='Column End should > Column Start!')
  }; req(lgc.c) 
  dat <- as.matrix(round(assay(input), 2))
  if (nrow(dat) < r2) r2 <- nrow(dat)
  if (ncol(dat) < c2) c2 <- ncol(dat)
  rmet <- rowData(input)
  if ('link' %in% colnames(rmet)) { 
    dat <- cbind.data.frame(rmet[, 'link', drop=FALSE], dat,
stringsAsFactors=FALSE)
  }
  col1 <- NULL; if ('metadata' %in% colnames(rmet)) { 
    col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))                               
    dat <- cbind.data.frame(rmet[, 'metadata', drop=FALSE], dat,
stringsAsFactors=FALSE)
  }
  colnames(dat) <- sub('__', '_', colnames(dat))
  datatable(dat[seq(r1, r2, 1), seq(c1, c2, 1), drop=FALSE], selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
  options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE,  scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, class='cell-border strip hover', columnDefs=col1, fixedColumns=list(leftColumns=1),
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}') 
  )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer') 
}

# Present assay data in tables.
dat_all_ui <- function(id) { 
  ns <- NS(id)
    navbarPage('', id=ns('dat'), selected='expDsg',
      tabPanel("Overview", value='over', 
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", DTOutput(ns("over")), "")),
      ), 
      tabPanel('Experiment design', value='expDsg', 
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("expDsg")), "")), 
      ),
      tabPanel("Assay data", value='assay',
      actionButton(ns("subdat"), "Run", icon=icon("sync"), style = run.col),
      div(id=ns('submsg'), style='width:400px', 
      fluidRow(splitLayout(cellWidths=c('12px', '90px', '1px', '90px', '1px', '95px', '1px', '90px'), '',
        numericInput(ns('r1'), label='Row start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('r2'), label='Row end', value=500, min=2, max=Inf, step=1), '',
        numericInput(ns('c1'), label='Column start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('c2'), label='Column end', value=50, min=2, max=Inf, step=1)
      ))),
      bsTooltip(id=ns('submsg'), title="Subsetting the data matrix for display only, not for downstream analysis.", placement =    "top", trigger = "hover"),
      dataTableOutput(ns("datall"))
      ) # tabPanel
  ) 
}

dat_all_server <- function(id, dat, r2=NULL, c2=NULL) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns
  observeEvent(dat(), { 
    dat <- dat(); if (!check_obj(list(dat))) return()
    if (is(dat, 'list')) {
      blk <- dat$bulk; cell <- dat$cell
      if (!check_obj(list(cell))) return()
      if (!check_obj(list(blk))) dat <- cbind(blk, cell) else dat <- cell
    }
    if (is.null(r2)) r2 <- nrow(dat)
    if (is.null(c2)) c2 <- ncol(dat)
    updateNumericInput(session, inputId="r2", value=r2)
    updateNumericInput(session, inputId="c2", value=c2) 
  }) 
  subdat <- reactiveValues(r1=1, r2=500, c1=1, c2=20) 
  observeEvent(list(input$subdat+1, dat()), ignoreInit=FALSE, { 
    subdat$r1 <- r1 <- input$r1; subdat$r2 <- r2 <- input$r2 
    subdat$c1 <- c1 <- input$c1; subdat$c2 <- c2 <- input$c2 
    dat <- dat() 
    if (!check_obj(list(r1, r2, c1, c2, dat))) return()
    if (is(dat, 'list')) {
      blk <- dat$bulk; cell <- dat$cell
      if (!check_obj(list(cell))) return()
      if (!check_obj(list(blk))) dat <- cbind(blk, cell) else dat <- cell
    }
    if (0 %in% input$subdat) {
      subdat$r2 <- nrow(dat)
      if (ncol(dat) >= 30) subdat$c2 <- 30 else subdat$c2 <- ncol(dat)
      return()
    }   
  })
  observeEvent(list(input$subdat+1, dat()), ignoreInit=FALSE, { 
    withProgress(message="Data table in co-visualization: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    dat <- dat(); r1 <- subdat$r1; r2 <- subdat$r2 
    c1 <- subdat$c1; c2 <- subdat$c2
    if (!check_obj(list(r1, r2, c1, c2, dat))) return()
    if (is(dat, 'list')) {
      blk <- dat$bulk; cell <- dat$cell
      if (!check_obj(list(cell))) return()
      if (!check_obj(list(blk))) dat <- cbind(blk, cell) else dat <- cell
    }
    output$datall <- renderDataTable({
      assay_dat(dat, r1=r1, r2=r2, c1=c1, c2=c2) 
    }); incProgress(0.3, detail="please wait ...") 
   }) 
  }) 
  output$expDsg <- renderDataTable({
    cat('Preparing sample metadata ... \n')
    dat <- dat(); if (!check_obj(list(dat))) return() 
    if (is(dat, 'list')) {
      blk <- dat$bulk; cell <- dat$cell
      if (!check_obj(list(cell))) return()
      if (!check_obj(list(blk))) dat <- cbind(blk, cell) else dat <- cell
    }; cdat <- colData(dat)
    withProgress(message="Sample metadata: ", value = 0, {
      incProgress(0.5, detail="please wait ...") 
      # Tooltip on metadata.
      col1 <- list(list(targets = seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(50, false)")))
      rownames(cdat) <- sub('__', '_', rownames(cdat))
      dtab <- datatable(data.frame(cdat), selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, class='cell-border strip hover', columnDefs=col1, fixedColumns = list(leftColumns=3), 
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
   )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
   # formatRound(colnames(assay.sel), deci)
   incProgress(0.4, detail="please wait ...")
    cat('Done! \n'); dtab
  })
  })
  output$over <- renderDT({
    cat('Preparing assay metadata ... \n')
    dat <- dat(); if (!check_obj(list(dat))) return() 
    if (is(dat, 'list')) {
      cell <- dat$cell 
      if (!check_obj(list(cell))) return() 
      meta <- metadata(cell)$df.meta
    } else meta <- metadata(dat)$df.meta
    if (is.null(meta)) return()
    withProgress(message="Assay/image metadata: ", value = 0, {
      incProgress(0.5, detail="please wait ...") 
      # Tooltip on metadata.
      col1 <- list(list(targets = seq_len(ncol(meta)), render = DT::JS("$.fn.dataTable.render.ellipsis(50, false)")))
      dtab <- datatable(data.frame(meta), selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, class='cell-border strip hover', columnDefs=NULL, fixedColumns = list(leftColumns=3), 
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
   )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
   # formatRound(colnames(assay.sel), deci)
   incProgress(0.4, detail="please wait ...")
    cat('Done! \n'); dtab
  })
  })
})}
