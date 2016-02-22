## User-interface for shiny GUI for scater
## Davis McCarthy
## 18 November 2015

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    titlePanel("Welcome to scater!"),
    sidebarLayout(
        sidebarPanel(h2("Overview of scater"),
                     p("The scater package has a great deal of functionality for
    the pre-processing, quality control, normalisation and 
                        visualisation of single-cell RNA-seq datasets."),
                     br(),
                     p("This interactive app can be used to explore your data
                        and apply that functionality."),
                     p("Accessing the tabs from left to right will step you 
                       through the steps in pre-processing your data."),
                     selectInput("plottype", label = h3("Select plot type"), 
                                 choices = list(
                                     "cumulative expression plot" = 1, 
                                     "plotExpression" = 2,
                                     "plotQC" = 3,
                                     "plotPhenoData" = 4,
                                     "plotFeatureData" = 5,
                                     "plotPCA" = 6,
                                     "plotTSNE" = 7), 
                                 selected = 1),
                     helpText("Select which type of plot you would like to 
                              produce from your SCESet object."),
                     hr(),
                     p("Current value for the plot type:"),
                     verbatimTextOutput("plottype"),
                     p("Help for functions can be accessed by typing", 
                       code("?function_name"), "at the R prompt.")
        ),
        
        mainPanel(img(src = "scater_qc_workflow.png", width = 650)
        )
    )
))


## Ideas: 
### show code for plots below the plot so that users can learn the commands

