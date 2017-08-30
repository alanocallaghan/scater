#' scater GUI function
#'
#' scater shiny app GUI for workflow for less programmatically inclined users or
#' those who would like a quick and easy way to view multiple plots.
#'
#' @param object SinglCellExperiment object after running \code{\link{calculateQCMetrics}} 
#' on it
#'
#' @return Opens a browser window with an interactive shiny app and visualize
#' all possible plots included in the scater
#'
#' @import shiny shinydashboard
#' @author Davis McCarthy and Vladimir Kiselev
#' @export
#' 
#' @examples 
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#' example_sce <- calculateQCMetrics(example_sce, 
#' feature_controls = list(set1 = 1:40))
#' \dontrun{
#' scater_gui(example_sce)
#' }
scater_gui <- function(object) {
    
    pd <- colnames(colData(object))
    pd.plot <- pd[!grepl("filter_", pd) & !grepl("is_", pd)]
    featurenames <- rownames(object)
    
    exprs_values <- assayNames(object)
    exprs_values <- exprs_values[!grepl("is_exprs", exprs_values)]
    
    shinyApp(
        ui <- dashboardPage(
            dashboardHeader(title = "scater"),
            dashboardSidebar(
                sidebarMenu(
                    menuItem("plot", tabName = "plot"),
                    menuItem("plotQC", tabName = "plotQC"),
                    menuItem("plotPCA - QC metrics", tabName = "plotPCA_QC"),
                    menuItem("plotPCA - expression", tabName = "plotPCA"),
                    menuItem("plotTSNE", tabName = "plotTSNE"),
                    menuItem("plotDiffusionMap", tabName = "plotDiffusionMap"),
                    menuItem("plotExpression", tabName = "plotExpression")
                )
            ),
            dashboardBody(
                tabItems(
                    tabItem(tabName = "plot",
                            fluidRow(
                                box(HTML("<h4>Overview of expression for each cell</h4>
                                         Plots produced by this function are intended 
                                         to provide an overview of large-scale 
                                         differences between cells. For each cell, 
                                         the features are ordered from most-expressed 
                                         to least-expressed and the cumulative 
                                         proportion of the total expression for 
                                         the cell is computed across the top 
                                         nfeatures features. These plots can flag 
                                         cells with a very high proportion of the 
                                         library coming from a small number of features; 
                                         such cells are likely to be problematic 
                                         for analyses. Using the colour and blocking 
                                         arguments can flag overall differences in 
                                         cells under different experimental conditions 
                                         or affected by different batch and other variables."),
                                    width = 12,
                                    status = "success")
                            ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plot", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       selectInput("block1", "block1:",
                                                   pd.plot,
                                                   selected = pd.plot[2]),
                                       selectInput("block2", "block2:",
                                                   pd.plot,
                                                   selected = pd.plot[3]),
                                       selectInput("colour_by", "colour_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("exprs_values", "exprs_values:",
                                                   exprs_values,
                                                   selected = "logcounts")
                                )
                            )
                    ),
                    tabItem(tabName = "plotQC",
                            fluidRow(
                                box(HTML("<h4>General plots</h4>
                                         <b>highest-expression</b> shows features with 
                                         highest expression<br>
                                         <b>explanatory-variables</b> shows a set of 
                                         explanatory variables plotted against each other, 
                                         ordered by marginal variance explained<br>
                                         <b>exprs-mean-vs-freq</b> plots the mean expression 
                                         levels against the frequency of expression for a 
                                         set of features"),
                                    width = 12,
                                    status = "success")
                            ),
                            fluidRow(
                                box(plotOutput("plotQC", height = 600), width = 8),
                                box(
                                    radioButtons("QCtype",
                                                 label = "Choose a type of QC plot",
                                                 choices = c("highest-expression",
                                                             "explanatory-variables",
                                                             "exprs-freq-vs-mean"),
                                                 selected = "highest-expression"),
                                    width = 4
                                )
                            ),
                            fluidRow(
                                box(HTML("<h4>Find PCs</h4>
                                         This plot shows the most important principal 
                                         components for a given variable"),
                                    width = 12,
                                    status = "success"),
                                box(plotOutput("plotQCfindpc", height = 600), width = 8),
                                box(
                                    radioButtons("QCvar",
                                                 label = "Choose a variable of interest",
                                                 choices = pd.plot,
                                                 selected = "total_features"),
                                    width = 4
                                )
                            )
                    ),
                    tabItem(tabName = "plotPCA_QC",
                            fluidRow(
                                box(HTML("<h4>PCA using QC metrics and cell variables</h4>
                                         Principal component analysis plots using QC metrics and cell metadata variables rather than expression levels. PCA on QC metrics can be used to identify potentially problematic cells, distinct from biological effects captured when using feature expression levels. See output in R session for names of detected outlier cells.<br>
                                         <br>
                                         Points can be coloured either by cell metadata variables (see drop-down menus for shape_by and size_by for options) or feature expression levels, just enter a valid name into the text box below."),
                                    width = 12,
                                    status = "success")
                                ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plotPCA_QC", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       checkboxInput("pcaqc_detect_outliers", 
                                                     "detect outliers?", 
                                                     value = TRUE),
                                       selectInput("pcaqc_selected_vars", 
                                                   "variables to use for PCA:",
                                                   pd.plot,
                                                   selected = 
                                                       c("pct_counts_top_100_features",
                                                         "total_features",
                                                         "pct_counts_feature_controls",
                                                         "n_detected_feature_controls",
                                                         "log10_counts_endogenous_features",
                                                         "log10_counts_feature_controls"),
                                                   multiple = TRUE),
                                       textInput("pcaqc_colour_by",
                                                 "colour_by (either cell metadata or feature expression):",
                                                 pd.plot[3],
                                                 placeholder = "Gene_0082"),
                                       selectInput("pcaqc_shape_by", "shape_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("pcaqc_size_by", "size_by:",
                                                   pd.plot,
                                                   selected = pd.plot[7]),
                                       numericInput("pcaqc_ncomponents",
                                                    "number of components:",
                                                    2, min = 2, max = 15),
                                       checkboxInput("pcaqc_scale_features", 
                                                     "scale_features", 
                                                     value = TRUE)
                                )
                            )
                            ),
                    tabItem(tabName = "plotPCA",
                            fluidRow(
                                box(HTML("<h4>PCA</h4>
                                        Principal component analysis plots using feature expression levels. PCA is particularly good for QC purposes.<br>
                                        <br>
                                        Points can be coloured either by cell metadata variables (see drop-down menus for shape_by and size_by for options) or feature expression levels, just enter a valid name into the text box below."),
                                    width = 12,
                                    status = "success")
                            ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plotPCA", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       # selectInput("pca_colour_by", "colour_by (either cell metadata or feature expression):",
                                       #             pd.plot,
                                       #             selected = pd.plot[4]),
                                       textInput("pca_colour_by",
                                                 "colour_by (either cell metadata or feature expression):",
                                                 pd.plot[3],
                                                 placeholder = "Gene_0082"),
                                       selectInput("pca_shape_by", "shape_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("pca_size_by", "size_by:",
                                                   pd.plot,
                                                   selected = pd.plot[7]),
                                       selectInput("pca_exprs_values",
                                                   "exprs_values:",
                                                   exprs_values, 
                                                   selected = "logcounts"),
                                       numericInput("pca_ntop",
                                                    "number of most variable features to use:",
                                                    500,
                                                    min = 50, max = 10000,
                                                    step = 25),
                                       numericInput("pca_ncomponents",
                                                    "number of components:",
                                                    2, min = 2, max = 15),
                                       checkboxInput("pca_scale_features", 
                                                     "scale_features", 
                                                     value = TRUE)
                                )
                            )
                    ),
                    tabItem(tabName = "plotTSNE",
                            fluidRow(
                                box(HTML("<h4>t-SNE</h4>
                                         Show a t-distributed stochastic neighbour embedding plot of cells. t-SNE is particularly good for displaying multiple distinct cell types.<br>
                                         <br>
                                         Points can be coloured either by cell metadata variables (see drop-down menus for shape_by and size_by for options) or feature expression levels, just enter a valid name into the text box below."),
                                    width = 12,
                                    status = "success")
                                ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plotTSNE", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       textInput("tsne_colour_by",
                                                 "colour_by (either cell metadata or feature expression):",
                                                 pd.plot[3],
                                                 placeholder = "Gene_0082"),
                                       selectInput("tsne_shape_by", "shape_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("tsne_size_by", "size_by:",
                                                   pd.plot,
                                                   selected = pd.plot[7]),
                                       selectInput("tsne_exprs_values",
                                                   "exprs_values:",
                                                   exprs_values, 
                                                   selected = "logcounts"),
                                       numericInput("tsne_ntop",
                                                    "number of most variable features to use:",
                                                    500,
                                                    min = 50, max = 10000,
                                                    step = 25),
                                       numericInput("tsne_ncomponents",
                                                    "number of components:",
                                                    2, min = 2, max = 15),
                                       checkboxInput("tsne_scale_features", 
                                                     "scale_features", 
                                                     value = TRUE),
                                       numericInput("tsne_rand_seed",
                                                    "random seed to make plot reproducible:",
                                                    5000)
                                )
                            )
                    ),
                    tabItem(tabName = "plotDiffusionMap",
                            fluidRow(
                                box(HTML("<h4>Diffusion Map</h4>
                                         Show a diffusion map plot of cells. Diffusion maps are particularly good for displaying cells at various stages along a continuous differentiation process.<br>
                                         <br>
                                         Points can be coloured either by cell metadata variables (see drop-down menus for shape_by and size_by for options) or feature expression levels, just enter a valid name into the text box below."),
                                    width = 12,
                                    status = "success")
                                ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plotDiffusionMap", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       textInput("diffmap_colour_by",
                                                 "colour_by (either cell metadata or feature expression):",
                                                 pd.plot[3],
                                                 placeholder = "Gene_0082"),
                                       selectInput("diffmap_shape_by", "shape_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("diffmap_size_by", "size_by:",
                                                   pd.plot,
                                                   selected = pd.plot[7]),
                                       selectInput("diffmap_exprs_values",
                                                   "exprs_values:",
                                                   exprs_values, 
                                                   selected = "logcounts"),
                                       numericInput("diffmap_ntop",
                                                    "number of most variable features to use:",
                                                    500,
                                                    min = 50, max = 10000,
                                                    step = 25),
                                       numericInput("diffmap_ncomponents",
                                                    "number of components:",
                                                    2, min = 2, max = 15),
                                       checkboxInput("diffmap_scale_features", 
                                                     "scale_features", 
                                                     value = TRUE),
                                       numericInput("diffmap_rand_seed",
                                                    "random seed to make plot reproducible:",
                                                    5000),
                                       radioButtons("diffmap_distance",
                                                    label = "Choose a distance metric",
                                                    choices = c("euclidean",
                                                                "cosine",
                                                                "rankcor"),
                                                    selected = "euclidean")
                                )
                            )
                    ),
                    tabItem(tabName = "plotExpression",
                            fluidRow(
                                box(HTML("<h4>Feature-level expression</h4>
                                         Plot expression levels for a set of features.<br>
                                         
                                         The x-axis variable for the plot can either be a cell metadata variable or expression levels for another feature."
                                         ),
                                    width = 12,
                                    status = "success")
                                ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plotExpression", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       selectInput("exprs_features", 
                                                   "features:",
                                                   featurenames,
                                                   selected = featurenames[1:6],
                                                   multiple = TRUE),
                                       # selectInput("exprs_x", "x-axis variable:",
                                       #             pd.plot,
                                       #             selected = pd.plot[4]),
                                       textInput("exprs_x",
                                                 "x-axis variable (either cell metadata variable or feature name):",
                                                 pd.plot[4],
                                                 placeholder = "Gene_0082"),
                                       selectInput("exprs_colour_by", "colour_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("exprs_shape_by", "shape_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("exprs_size_by", "size_by:",
                                                   pd.plot,
                                                   selected = pd.plot[7]),
                                       selectInput("exprs_exprs_values",
                                                   "exprs_values:",
                                                   exprs_values, 
                                                   selected = "logcounts"),
                                       numericInput("exprs_ncols",
                                                    "number of columns:",
                                                    2, min = 1, max = 8),
                                       checkboxInput("exprs_show_median", 
                                                     "show median?", 
                                                     value = FALSE),
                                       checkboxInput("exprs_show_violin", 
                                                     "show violin?", 
                                                     value = TRUE),
                                       checkboxInput("exprs_show_smooth", 
                                                     "show smoothed fit?", 
                                                     value = FALSE),
                                       checkboxInput("exprs_log2", 
                                                     "transform expression values to log2 scale?", 
                                                     value = FALSE)
                                )
                            )
                    )
                )
            )
        ),
        server <- function(input, output, session) {
            output$plot <- renderPlot({
                plotScater(object, exprs_values = input$exprs_values,
                     block1 = input$block1,
                     block2 = input$block2,
                     colour_by = input$colour_by)
            })
            output$plotQC <- renderPlot({
                plotQC(object, type = input$QCtype)
            })
            output$plotQCfindpc <- renderPlot({
                plotQC(object, type = "find-pcs", variable = input$QCvar)
            })
            output$plotPCA_QC <- renderPlot({
                plotPCA(object, 
                        ncomponents = input$pcaqc_ncomponents,
                        pca_data_input = "pdata",
                        selected_variables = input$pcaqc_selected_vars,
                        detect_outliers = input$pcaqc_detect_outliers,
                        colour_by = input$pcaqc_colour_by,
                        size_by = input$pcaqc_size_by,
                        shape_by = input$pcaqc_shape_by,
                        scale_features = input$pcaqc_scale_features,
                        legend = "all") +
                    theme(legend.position = "bottom")
            })
            output$plotPCA <- renderPlot({
                plotPCA(object, 
                        ntop = input$pca_ntop,
                        ncomponents = input$pca_ncomponents,
                        exprs_values = input$pca_exprs_values,
                        colour_by = input$pca_colour_by,
                        size_by = input$pca_size_by,
                        shape_by = input$pca_shape_by,
                        scale_features = input$pca_scale_features) +
                    theme(legend.position = "bottom")
            })
            output$plotTSNE <- renderPlot({
                plotTSNE(object, 
                        ntop = input$tsne_ntop,
                        ncomponents = input$tsne_ncomponents,
                        exprs_values = input$tsne_exprs_values,
                        colour_by = input$tsne_colour_by,
                        size_by = input$tsne_size_by,
                        shape_by = input$tsne_shape_by,
                        scale_features = input$tsne_scale_features,
                        rand_seed = input$tsne_rand_seed) +
                    theme(legend.position = "bottom")
            })
            output$plotDiffusionMap <- renderPlot({
                plotDiffusionMap(object, 
                         ntop = input$diffmap_ntop,
                         ncomponents = input$diffmap_ncomponents,
                         exprs_values = input$diffmap_exprs_values,
                         colour_by = input$diffmap_colour_by,
                         size_by = input$diffmap_size_by,
                         shape_by = input$diffmap_shape_by,
                         scale_features = input$diffmap_scale_features,
                         rand_seed = input$diffmap_rand_seed,
                         distance = input$diffmap_distance) +
                    theme(legend.position = "bottom")
            })
            output$plotExpression <- renderPlot({
                plotExpression(object, 
                               features = input$exprs_features,
                               x = input$exprs_x,
                               exprs_values = input$exprs_exprs_values,
                               colour_by = input$exprs_colour_by,
                               size_by = input$exprs_size_by,
                               shape_by = input$exprs_shape_by,
                               ncol = input$exprs_ncols,
                               show_median = input$exprs_show_median,
                               show_violin = input$exprs_show_violin,
                               show_smooth = input$exprs_show_smooth,
                               log2_values = input$exprs_log2) +
                    theme(legend.position = "bottom")
            })
            session$onSessionEnded(function() {
                stopApp()
            })
        },
        options = list(launch.browser = TRUE)
    )
}
