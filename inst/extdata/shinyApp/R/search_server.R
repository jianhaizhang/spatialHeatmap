# Module for searching gene ids.

search_server <- function(id, ids, lis.url, upl.mod.lis, dat.mod.lis, session) {
  moduleServer(id, function(input, output, session) {
    cat('ID search box ... \n'); ns <- session$ns 
    ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
    datIn <- reactiveValues()
    observe({
      se.scl <- check_exp(dat.mod.lis$se.scl())
      if (is(se.scl, 'character')) req('')
      if (!check_obj(list(se.scl))) req('')
      datIn$v <- se.scl
    })
    pre.id <- reactiveValues(id=NULL)
    observe({ # Pre-selected ids in config file.
      if (is.null(datIn$v)) return()
      lis.par <- cfg$lis.par; req(check_obj(lis.par))
      rna <- rownames(datIn$v)
      id <- lis.par$data.matrix['selected.id', 'default']
      if (length(id)==0) { pre.id$id <- rna[1]; return() }
      id <- make.names(strsplit(id, ',')[[1]])
      id <- id[id %in% rna]
      if (length(id)==0) id <- rna[1]; pre.id$id <- id
    })
    sch.sgl <- reactiveValues(cho=NULL, sel=NULL)
    observeEvent(list(datIn$v, pre.id$id, lis.url$par$ids), { # Accounts for IDs in URL for single search box.
      se.scl <- datIn$v; if (!check_obj(list(se.scl))) return()
      rna <- rownames(se.scl); df.met <- rowData(se.scl)
      cho <- paste0(rna, ' ', ifelse(rep('metadata' %in% colnames(df.met), length(rna)), df.met[, 'metadata'], ''))
      id <- pre.id$id; id.url <- lis.url$par$ids 
      if (check_obj(id.url)) {
        if (all(id.url %in% rna)) id <- id.url
      }; req(all(id %in% rna))
      if ('metadata' %in% colnames(df.met)) sel <- paste0(id, ' ', df.met[id, 'metadata']) else sel <- paste0(id, ' ')
      sch.sgl$cho <- cho; sch.sgl$sel <- sel
      # updateSelectizeInput(session, 'sch.sgl', choices=cho, selected=sel, server=TRUE)
    })
    mul.box <- reactiveValues()
    observeEvent(list(datIn$v, pre.id$id, lis.url$par$ids), {
      se.scl <- datIn$v; if (!check_obj(list(se.scl))) return()
      rna <- rownames(se.scl)
      id <- pre.id$id; id.url <- lis.url$par$ids
      if (check_obj(id.url)) { 
        if (all(id.url %in% rna)) id <- id.url
      }; req(all(id %in% rna))
      mul.box$v <- id
    })
   observe({
      ids$but.mul <- input$sch.mul.but; ids$but.sgl <- input$sch.sgl.but 
   })
  output$sch.box <- renderUI({
    cat('Rendering search box ... \n')
    sch.mod <- input$sch.mode; init.v <- mul.box$v
    init.v <- paste0(init.v, collapse=',')
    # At the very beginning, two search buttons and two search boxes are all NULL, since they are not rendered. When toggling between search modes, one of two buttons is reverted to 0 not NULL, while the search content is reverted to NULL then strings.
    if (sch.mod=='Single') { 
      sbox <- selectizeInput(ns('sch.sgl'), p(style='margin-bottom:-10px;', lab.sgl, actionButton(ns("sch.sgl.but"), "Plot", style=conf.col)), choices='none', selected='none', multiple = TRUE, options=list(placeholder = 'Supports auto-completion.')) 
    } else if (sch.mod=='Multiple') { 
    sbox <- textInput(inputId=ns('sch.mul'), p(style='margin-bottom:-10px;', lab.mul, actionButton(ns("sch.mul.but"), "Plot", style=conf.col)), value=init.v, placeholder='Muliple IDs must ONLY be separated by space or comma.', width='100%')
   }; cat('Done! \n'); sbox

  })
  rna <- reactiveValues(val=NULL)
  observe({ # rna$val: accounts for filtered data.  
    se.scl <- datIn$v
    if (is.null(se.scl)) return()
    rna$val <- rownames(se.scl)
  })
  upd.sel <- reactiveValues()
  observeEvent(list(sch.sgl$cho, sch.sgl$sel, input$sch.mode, rna$val), {
    # if (is.null(input$sch.sgl)) return()
    sch.sgl.but <- input$sch.sgl.but
    pars <- list(sch.sgl$cho, sch.sgl$sel, input$sch.mode, rna$val)
    # if (is.null(sch.sgl.but) | 0 %in% sch.sgl.but) pars <- c(pars, list(sch.sgl.but))
    upd.sel$pars <- pars
  })
  observeEvent(upd.sel$pars, {
    updateSelectizeInput(session, 'sch.sgl', choices=sch.sgl$cho, selected=sch.sgl$sel, server=TRUE)
  }) 
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.sgl.but, rna$val), { 
    cat('Single search IDs ... \n')
    sch.sgl.but <- input$sch.sgl.but; req(sch.sgl.but)
    if (is.null(input$sch.sgl)) return()
    if (is.null(datIn$v)|'none' %in% input$sch.sgl) return()
    sel <- sub(' .*', '', input$sch.sgl)
    # validate: holds on til condition is met, can address issues in execution orders.
    validate(need(all(sel %in% rna$val), ''))
    lgc <- identical(sort(sel), sort(ids$sel)) 
    if (sch.sgl.but > 0) if (lgc) showModal(modal(msg = 'Please enter new gene ids to continue!')); req(!lgc) 
    ids$sel <- sel; ids$but.sgl <- sch.sgl.but
    cat('Done! \n')
  })

  observe({ # On-start IDs in single search mode.
    return() # Ensures a default gene is selected but SHMs are not plotted. 
    cat('Single search IDs: on-start ... \n')
    # print(list(input$sch.sgl, ids$sel, rna$val[1:3], input$sch.sgl.but))
    sch.sgl <- input$sch.sgl
    # Ensures on-start.
    if (sum(ids$sel %in% rna$val)>0|is.null(sch.sgl)) return()
    if (sch.sgl[1]=='none') return()
    ids$sel <- sub(' .*', '', sch.sgl); ids$but.sgl <- input$sch.sgl.but
    cat('Done! \n')
  })
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.mul.but, rna$val), {
    cat('Multiple search IDs ... \n')
    sch.mul.but <- input$sch.mul.but; sch.mul <- input$sch.mul
    req(sch.mul.but)
    if (is.null(rna$val)|!is.character(sch.mul)) return()
    lgc <- sch.mul==''
    if (lgc) { 
      # sel <- pre.id$id
      showModal(modal(msg = 'No genes ids entered!')); req(!lgc)
    } else {
      gens <- strsplit(gsub(' |,', '__', sch.mul), '__')[[1]]
      sel <- gens[tolower(gens) %in% tolower(rna$val)]
      # Invalid IDs.
      dif <- setdiff(gens, sel)
      msg <- paste0('ID(s) not detected: ', paste0(dif, collapse=', ')); cat(msg, '\n')
      lgc.dif <- length(dif)==0
      if (!lgc.dif) showModal(modal(msg = msg)); req(lgc.dif)
      lgc.sel <- length(sel)>0 
      if (!lgc.sel) { 
        # sel <- pre.id$id
        showModal(modal(msg = 'Genes ids are invalid!')); req(lgc.sel)
      }
      # Eleminate case difference.
      sel <- rna$val[tolower(rna$val) %in% tolower(sel)]
    }
    lgc.new <- identical(sort(sel), sort(ids$sel)) 
    if (sch.mul.but>0) if (lgc.new) showModal(modal(msg = 'Please enter new gene ids to continue!')); req(!lgc.new) 
    ids$sel <- sel; ids$but.mul <- sch.mul.but
    cat('Done! \n') 
  })
  # observeEvent: input$sch.mul and ids$sel affects each other and lead to infinite circles.  
  observeEvent(list(input$sch.mul, rna$val, pre.id$id), { # On-start IDs in multiple search mode.
    return() # Ensures a default gene is selected but SHMs are not plotted. 
    cat('Multiple search IDs: on-start ... \n')
    sch.mul <- input$sch.mul
    # Ensures on-start.
    if (sum(ids$sel %in% rna$val)>0|is.null(sch.mul)) return()
    if (sch.mul=='') return()
    # if (!all(ids$sel %in% sch.mul)) return()
    gens <- strsplit(gsub(' |,', '__', sch.mul), '__')[[1]]
    sel <- gens[tolower(gens) %in% tolower(rna$val)]
    # Invalid IDs.
    dif <- setdiff(gens, sel)
    msg <- paste0('ID(s) not detected: ', paste0(dif, collapse=', ')); cat(msg, '\n')
    lgc.dif <- length(dif)==0
    if (!lgc.dif) showModal(modal(msg = msg)); req(lgc.dif)
    lgc <- length(sel) >0
    if (!lgc) {
      # sel <- pre.id$id
      showModal(modal(msg = 'Genes ids are invalid!')); req(lgc)
    }
    # Eleminate case difference.
    sel <- rna$val[tolower(rna$val) %in% tolower(sel)]
    ids$sel <- sel; ids$but.mul <- input$sch.mul.but
    cat('Done! \n')
  })
  # observeEvent(ipt$fileIn, { ids$sel <- NULL })
  return(list(ids=ids))
  onBookmark(function(state) { state })
  })
}
