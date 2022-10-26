# Module for matrix heatmap and network analysis.
network_ui <- function(id) {
  ns <- NS(id)
  list(
  tabPanel(id='cluster', "Cluster", value='clus', icon=NULL,
      fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('1%', '15%', '1%', '15%'), '', 
      dropdownButton(inputId=ns('dpbMea'), label='Similarity/Dissimilarity', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      radioButtons(inputId=ns('measure'), label="Measure", choices=c('Correlation', 'Distance'), selected='Correlation', inline=TRUE, width='100%'), 
      div(title='Only applicable when "Correlation" is selected.',
      radioButtons(inputId=ns("cor.abs"), label="Absolute correlation", choices=c('No', 'Yes'), selected='No', inline=TRUE, width='100%')
      )
      ), '',
      dropdownButton(inputId=ns('dpbThr'), label='Subset', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=400,
        div(title='Subset the most similar neighbors at a cutoff by top proportion, specific number, or similarity/dissimilairty value.',
        radioButtons(inputId=ns("thr"), label="By", choices=c('Proportion'='p', 'Number'='n', 'Similarity/Dissimilarity'='v'), selected='p', inline=TRUE, width='100%'),
        numericInput(inputId=ns('mhm.v'), label='Cutoff', value=0.2, min=-Inf, max=Inf, step=0.1, width='100%')
        )
      )
      )),
      bsTooltip(id=ns('dpbMea'), title='Measure to subset the nearest neighbors for the target gene (s) in spatial heatmaps.', placement = "top", trigger = "hover"),
      navbarPage('', id=ns('clusNav'),
      tabPanel(strong('Matrix Heatmap (MHM)'), value='mhmPlot', 
      actionButton(ns("mhm.but"), "Click to show/update", icon=icon("sync"), style="color: #fff; background-color:#499fe9;border-color:#2e6da4"), plotlyOutput(ns("HMly"))
      ),
      tabPanel('Parameter (MHM)', value='mhmPar',
      fluidRow(splitLayout(cellWidths=c('1%', '20%'), '', 
      radioButtons(inputId=ns("mat.scale"), label="Scale by", choices=c("No", "Column", "Row"), selected='No', inline=TRUE)
      # radioButtons(inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected='No', inline=TRUE)
      ))
      ), # tabPanel('Plot', value='mhmPar'
    tabPanel(strong("Interactive Network (NET)"), value='netPlot', icon=NULL, 
      actionButton(ns("cpt.nw"), "Click to show/update", icon=icon("sync"), style="color: #fff; background-color:#499fe9;border-color:#2e6da4"),
      fluidRow(splitLayout(cellWidths=c("1%", "5%", "92%", "2%"), "", plotOutput(ns("bar.net")), visNetworkOutput(ns("vis")), ""))
    ),
    tabPanel("Parameter (NET)", value='netPar', icon=NULL, 
      #fluidRow( # If column widths are not integers, columns are vertically aligned.
      fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('0.5%', '9%', '0.5%', '10%', '0.5%', '14%', '0.5%', '15%', '0.5%', '11%'), '',
      dropdownButton(inputId=ns('dpwNetTar'), label='Target gene', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        selectInput(ns("gen.sel"), "", c("None"), selected='None', width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwNetType'), label='Network type', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      selectInput(inputId=ns("net.type"), label="", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected='signed', width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwModSize'), label='Module size/splitting', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
        numericInput(ns("min.size"), "Min module size", value=15, min=15, max=5000, width='100%'),
        selectInput(ns("ds"), HTML("Module splitting sensitivity level <br/> (Larger value results in more <br/> modules with smaller sizes)"), 3:2, selected=3, width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwNetAdj'), label='Adjacency/edges', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
        selectInput(ns("adj.in"), HTML("Adjacency threshold <br/> (the smaller, the more edges)"), sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width='100%'),
        numericInput(ns("max.edg"), HTML("Maximun edges <br/> (too many edges may crash the app)"), value=50, min=1, max=500, width='100%'), htmlOutput(ns("edge"))
      ), '',
      dropdownButton(inputId=ns('dpwColNet'), label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '',
      textInput(ns("color.net"), "Color scheme", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width='100%'),
      actionButton(ns("col.but.net"), "Confirm", icon=icon("sync"), style='margin-top:24px'))) # fluidRow(splitLayout
      )
      )) # fluidRow(splitLayout
     )
    ) # navbarPage('', id=ns('mhmNav')
  )
  ) # list
}
