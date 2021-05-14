# In server.R and ui.R, the relative directories are used such as 'config/config.yaml'. Though they are in the R directory, they are used by app.R, so these relative directories are valid. Since the relative directories are effective at the place where they are actualy called (i.e. app.R).
shinyApp(ui, server, enableBookmarking = "url")
