.metadata_dispatcher <- function(object, mode, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
# Internal function to create the data frames, given an indication of 
# whether we are looking at row or column-level metadata.    
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define dataframe to pass to plotMetadata
    x_by_out <- .choose_vis_values(object, x, mode = mode, search = "metadata")
    x_lab <- x_by_out$name
    y_by_out <- .choose_vis_values(object, y, mode = mode, search = "metadata")
    y_lab <- y_by_out$name
    if (!is.null(x)) {
        df_to_plot <- data.frame(X=x_by_out$val, Y=y_by_out$val)
    } else {
        df_to_plot <- data.frame(Y=y_by_out$val, X=factor(character(nrow(object))))
    }

    ## checking visualization arguments
    colour_by_out <- .choose_vis_values(object, colour_by, mode = mode, search = "any", 
                                        exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val
    
    shape_by_out <- .choose_vis_values(object, shape_by, mode = mode, search = "any", 
                                       exprs_values = exprs_values, coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val
    
    size_by_out <- .choose_vis_values(object, size_by, mode = mode, search = "any", 
                                      exprs_values = exprs_values)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    df_to_plot$shape_by <- shape_by_vals
    df_to_plot$size_by <- size_by_vals
    df_to_plot$colour_by <- colour_by_vals

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
                     colour_by = colour_by, size_by = size_by, shape_by = shape_by, ...)
}

#' Plot cell phenotype data from an SingleCellExperiment object
#'
#' \code{plotPhenoData}, \code{plotColData} and \code{plotCellData} are
#' synonymous.
#'
#' @param object a \code{\link{SingleCellExperiment}} object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to plot
#' total_features against log10(total_counts).
#' @param ... arguments passed to \code{\link{plotPhenoData}} (if
#' \code{\link{plotColData}} or \code{\link{plotCellData}}) or to
#' \code{\link{plotMetadata}}, e.g.\code{theme_size}, \code{size},
#' \code{alpha}, \code{shape}.
#'
#' @details Plot phenotype data from a SingleCellExperiment object. If one variable is
#' supplied then a density plot will be returned. If both variables are
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid
#' will be returned. If both variables are discrete then a jitter plot will be
#' produced. The object returned is a ggplot object, so further layers and
#' plotting options (titles, facets, themes etc) can be added.
#'
#' @return a ggplot plot object
#'
#' @name plotColData
#' @rdname plotColData
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce)
#' example_sce <- normalize(example_sce)
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status",
#'    size_by = "Gene_0001", shape_by = "Treatment")
#'
#' plotColData(example_sce, y = "Treatment", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "Mutation_Status", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "Mutation_Status", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status",
#'    size_by = "Gene_0001", shape_by = "Treatment")
plotColData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
{
    .metadata_dispatcher(object, mode = "column", y = y, x = x,
                         colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                         exprs_values = exprs_values, ...)
}


#' @rdname plotColData
#' @export
plotPhenoData <- function(...) {
    plotColData(...)
}

#' @rdname plotColData 
#' @export
plotCellData <- function(...) {
    plotColData(...)
}

#' Plot feature (gene) data from a SingleCellExperiment object
#'
#' \code{plotFeatureData} and \code{plotRowData} are synonymous.
#'
#' @param object a \code{\link{SingleCellExperiment}} object containing
#' expression values and experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to produce
#' a density plot of number of cells expressing the feature (requires
#' \code{calculateQCMetrics} to have been run on the \code{SingleCellExperiment} object prior).
#' @param ... arguments passed to \code{\link{plotMetadata}}, e.g.
#' \code{theme_size}, \code{size}, \code{alpha}, \code{shape}.
#'
#' @details Plot feature (gene) data from an SingleCellExperiment object. If one variable is
#' supplied then a density plot will be returned. If both variables are
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid
#' will be returned. If both variables are discrete then a jitter plot will be
#' produced. The object returned is a ggplot object, so further layers and
#' plotting options (titles, facets, themes etc) can be added.
#'
#' @return a ggplot plot object
#'
#' @name plotRowData 
#' @rdname plotRowData
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce,
#'     feature_controls = list(ERCC=1:40))
#' example_sce <- normalize(example_sce)
#' 
#' plotRowData(example_sce, y="n_cells_counts", x="log10_total_counts")
#' plotRowData(example_sce, y="n_cells_counts", 
#'    size_by ="log10_total_counts",
#'    colour_by = "is_feature_control")
#'
plotRowData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
{
    .metadata_dispatcher(object, mode = "row", y = y, x = x,
                         colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                         exprs_values = exprs_values, ...)
}

#' @rdname plotRowData 
#' @export
plotFeatureData <- function(...) {
    plotRowData(...)
}



