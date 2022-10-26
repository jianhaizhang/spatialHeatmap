# Module for searching gene ids.

search_server <- function(id, ids, lis.url, url.id, upl.mod.lis, dat.mod.lis, session) {
  moduleServer(id, function(input, output, session) {
    cat('ID search box ... \n'); ns <- session$ns 
    ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
    observeEvent(session$clientData$url_search, {
      id.sch.sgl <- url_val('sear-sch.sgl', lis.url)
      # lis <- reactiveValuesToList(lis.url)
      # In bookmarked url of single search box, if there is a comma in gene ID and annotation combination, the combined ID and annotation will be broken at the comma, resulting two strings.
      if (id.sch.sgl[1]!='null') {
        # Multiple ID-discs are concatenated and separated by comma. E.g [ID1 disc1,ID2 disc2].
        id.sch.sgl <- unlist(strsplit(gsub('^\\[|\\]$', '', id.sch.sgl), ','))
      }
      url.id$sch.sgl <- lis.url$par[['sear-sch.sgl']] <- id.sch.sgl
      url.id$sch.mul <- lis.url$par[['sear-sch.mul']] <- unlist(strsplit(gsub(' |,', '__', url_val('sear-sch.mul', lis.url)), '__'))
      cat('Parsing IDs from URL ... done! \n')
    })
    geneIn <- reactiveValues(lis=NULL)
    observe({
      if (is.null(dat.mod.lis)) return()
      gene.all <- dat.mod.lis$geneIn; geneIn$lis <- gene.all()
    })
    pre.id <- reactiveValues(id=NULL)
    observe({ # Pre-selected ids in config file.
      if (is.null(geneIn$lis)) return()
      df.aggr.tran <- geneIn$lis$df.aggr.tran 
      rna <- rownames(df.aggr.tran)
      id <- cfg$lis.par$data.matrix['selected.id', 'default']
      if (length(id)==0) { pre.id$id <- rna[1]; return() }
      id <- make.names(strsplit(id, ',')[[1]])
      id <- id[id %in% rna]
      if (length(id)==0) id <- rna[1]; pre.id$id <- id
    })
    sch.sgl <- reactiveValues(cho=NULL, sel=NULL)
    observe({ # Accounts for IDs in URL for single search box.
      if (is.null(geneIn$lis)) return()
      # geneIn <- dat.mod.lis$geneIn; gen.lis <- geneIn()
      df.aggr.tran <- geneIn$lis$df.aggr.tran
      rna <- rownames(df.aggr.tran); df.met <- geneIn$lis$df.met
      cho <- paste0(rna, ' ', ifelse(rep('metadata' %in% colnames(df.met), length(rna)), df.met[, 'metadata'], ''))
      id <- pre.id$id
      if (!is.null(url.id$sch.sgl)) if (url.id$sch.sgl[1]!='null') {
        # Only obtain the ID, no annotation.
        sel <- sub(' .*', '', url.id$sch.sgl)
        id <- sel[sel %in% rna]
      }
      sel <- paste0(id, ' ', df.met[id, 'metadata'])
      sch.sgl$cho <- cho; sch.sgl$sel <- sel
      # updateSelectizeInput(session, 'sch.sgl', choices=cho, selected=sel, server=TRUE)
    })
  output$sch.box <- renderUI({
    cat('Rendering search box ... \n')
    sch.mod <- input$sch.mode; mul.val <- pre.id$id
    if (length(url.id$sch.mul)!=0) if (url.id$sch.mul[1]!='null') mul.val <- url.id$sch.mul
    mul.val <- paste0(mul.val, collapse=',')
    # At the very beginning, two search buttons and two search boxes are all NULL, since they are not rendered. When toggling between search modes, one of two buttons is reverted to 0 not NULL, while the search content is reverted to NULL then strings.
    if (sch.mod=='Single') { 
      sbox <- selectizeInput(ns('sch.sgl'), p(lab.sgl, actionButton(ns("sch.sgl.but"), "Confirm selection")), choices='none', selected='none', multiple = TRUE, options=list(placeholder = 'Supports auto-completion.')) 
    } else if (sch.mod=='Multiple') { 
    sbox <- textInput(inputId=ns('sch.mul'), p(lab.mul, actionButton(ns("sch.mul.but"), "Confirm selection")), value=mul.val, placeholder='Muliple IDs must ONLY be separated by space or comma.', width='100%')
   }; cat('Done! \n'); sbox

  })
  rna <- reactiveValues(val=NULL)
  observe({ # rna$val: accounts for filtered data.  
    if (is.null(geneIn$lis)) return()
    rna$val <- rownames(geneIn$lis$df.aggr.tran)
  })
  upd.sel <- reactiveValues()
  observe({
    if (is.null(input$sch.sgl)) return()
    pars <- list(sch.sgl$cho, sch.sgl$sel, input$sch.mode, rna$val)
    upd.sel$pars <- pars
  })
  observeEvent(upd.sel$pars, {
    updateSelectizeInput(session, 'sch.sgl', choices=sch.sgl$cho, selected=sch.sgl$sel, server=TRUE)
  }) 
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.sgl.but, rna$val), { 
    cat('Single search IDs ... \n')
    if (is.null(input$sch.sgl)) return()
    if (is.null(geneIn$lis)|'none' %in% input$sch.sgl) return()
    sel <- sub(' .*', '', input$sch.sgl)
    # validate: holds on til condition is met, can address issues in execution orders.
    validate(need(all(sel %in% rna$val), ''))
    ids$sel <- sel; ids$but.sgl <- input$sch.sgl.but
    cat('Done! \n')
  })

  observe({ # On-start IDs in single search mode.
    cat('Single search IDs: on-start ... \n')
    # print(list(input$sch.sgl, ids$sel, rna$val[1:3], input$sch.sgl.but))
    sch.sgl <- input$sch.sgl
    if (sum(ids$sel %in% rna$val)>0|is.null(sch.sgl)) return()
    if (sch.sgl[1]=='none') return()
    ids$sel <- sub(' .*', '', sch.sgl); ids$but.sgl <- input$sch.sgl.but
    cat('Done! \n')
  })
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.mul.but, rna$val), {
    cat('Multiple search IDs ... \n')
    sch.mul <- input$sch.mul
    if (is.null(rna$val)|!is.character(sch.mul)) return()
    if (sch.mul=='') sel <- pre.id$id else {
      gens <- strsplit(gsub(' |,', '__', sch.mul), '__')[[1]]
      sel <- gens[tolower(gens) %in% tolower(rna$val)]
      # Invalid IDs.
      dif <- setdiff(gens, sel)
      msg <- paste0('ID(s) not detected: ', paste0(dif, collapse=', ')); cat(msg, '\n')
      if (length(dif)>0) showNotification(msg, duration=2)
      validate(need(length(dif)==0, ''))
      if (length(sel)==0) sel <- pre.id$id
      # Eleminate case difference.
      sel <- rna$val[tolower(rna$val) %in% tolower(sel)]
    }; ids$sel <- sel; ids$but.mul <- input$sch.mul.but
    cat('Done! \n') 
  })
  # observeEvent: input$sch.mul and ids$sel affects each other and lead to infinite circles.  
  observeEvent(list(input$sch.mul, rna$val, pre.id$id), { # On-start IDs in multiple search mode.
    cat('Multiple search IDs: on-start ... \n')
    sch.mul <- input$sch.mul
    if (sum(ids$sel %in% rna$val)>0|is.null(sch.mul)) return()
    if (sch.mul=='') return()
    # if (!all(ids$sel %in% sch.mul)) return()
    gens <- strsplit(gsub(' |,', '__', sch.mul), '__')[[1]]
    sel <- gens[tolower(gens) %in% tolower(rna$val)]
    # Invalid IDs.
    dif <- setdiff(gens, sel)
    msg <- paste0('ID(s) not detected: ', paste0(dif, collapse=', ')); cat(msg, '\n')
    if (length(dif)>0) showNotification(msg, duration=2)
    validate(need(length(dif)==0, ''))
    if (length(sel)==0) sel <- pre.id$id
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
