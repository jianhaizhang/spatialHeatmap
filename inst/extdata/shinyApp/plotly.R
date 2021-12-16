library(plotly)
#install.packages("dplyr")
library(dplyr)
#install.packages("shiny")
library(shiny)
#install.packages("htmlwidgets")
library(htmlwidgets)

## Read the fake csv data set about ice cream
ice.cream.df <- read.csv("~/Downloads/icecream-img.csv", stringsAsFactors = FALSE)

## Create the shiny UI layout
ui <- fluidPage(
  headerPanel("Price per Scoops"),
  sidebarPanel(
    sliderInput("priceRange", label = h3("Price Range"), min = 0, 
                max = 10, value = c(2, 8))
  ),
  mainPanel(
    plotlyOutput("icePlot"),
    h4("Click on the dots to learn more about the ice cream flavor."),
    h4("Use the lasso/box select to learn more about the ratings of each ice cream flavor."),
    plotlyOutput("ratingPlot"),
    uiOutput("imageLink")
  )
)

## Create the Shiny Server layout
server <- function(input, output) {
  ## Create the plotly plot that compares price vs scoops
  output$icePlot <- renderPlotly({
    range.ice.cream <- ice.cream.df %>% filter(prices >= input$priceRange[1]) %>% filter(prices <= input$priceRange[2])
    plot_ly(range.ice.cream, x = ~prices, y = ~scoops, type = "scatter", mode = "markers",
            text = ~paste("Flavor:", flavors), key = ~images, source = "imgLink") %>% 
            layout(title = paste("Ice Cream Price vs Scoops Given")) 
  })
  
  ## Create text paragraph of info on a selected point
  output$imageLink <- renderText({
    event.data <- event_data(event = "plotly_click", source = "imgLink")
    if (is.null(event.data)) {
      print("Click to see the link of the point.")
    } else { 
      ice.cream <- ice.cream.df %>% filter(images == event.data$key)
      HTML('<p>Flavor:',ice.cream$flavors, '</p>','<p>X Value:', event.data[["x"]], '</p>','<p>Y Value:', event.data[["y"]],'</p>',
           '<a href="', ice.cream$images,'">', ice.cream$images,'</a>','<p>','<img src="',ice.cream$images, '"/>','</p>')
    }
  })
  
  ## Create the plotly plot of price vs rating based on selection
  output$ratingPlot <- renderPlotly({
    event.data <- event_data(event = "plotly_selected", source = "imgLink")
    if (is.null(event.data)) {
      print("Click and drag events (i.e., select/lasso) to make the bar plot appear here")
      plot_ly(ice.cream.df, x = ~flavors, y = ~rating, type = "bar",
              text = ~paste("Flavor:", flavors)) %>%
              layout(title = paste("Ice Cream Ratings Given by Flavor"))
    } else {
      ice.cream <- ice.cream.df[ice.cream.df$images %in% event.data$key,]
      plot_ly(ice.cream, x = ~flavors, y = ~rating, type = "bar",
              text = ~paste("Flavor:", flavors), key = ~images, source = "imgLink") %>%
              layout(title = paste("Ice Cream Ratings Given by Flavor"))
    }
  })
}

shinyApp(ui = ui, server = server)
