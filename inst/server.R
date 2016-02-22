## server script for shiny GUI for scater
## Davis McCarthy
## 18 November 2015

library(shiny)

# Define server logic 
shinyServer(function(input, output) {
    
    # You can access the value of the widget with input$select, e.g.
    output$plottype <- renderPrint({input$plottype})
})

