log_exp_ui <- function(id) {
  ns <- NS(id)
  radioButtons(inputId=ns('log'), label='Log/exp-transform', choices=c("No", "log2", "exp2"), selected='No', inline=TRUE)
}

log_exp_server <- function(id, arg.lis) {
  moduleServer(id, function(input, output, session) {
  observe({
    arg.lis; cfg <- arg.lis$cfg
    updateRadioButtons(session, inputId='log', label='Log/exp transform', choices=c("No", "log2", "exp2"), selected=cfg$lis.par$data.matrix['log.exp', 'default'], inline=TRUE)
  })

  })
}


