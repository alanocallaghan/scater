#' Plot column metadata
#'
#' Plot column-level (i.e., cell) metadata in an SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object containing expression values and experimental information.
#' @param y Specification of the column-level metadata to show on the y-axis, see \code{?"\link{scater-vis-var}"} for possible values.
#' Note that only metadata fields will be searched, \code{assays} will not be used.
#' @param x Specification of the column-level metadata to show on the x-axis, see \code{?"\link{scater-vis-var}"} for possible values.
#' Again, only metadata fields will be searched, \code{assays} will not be used.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details 
#' If \code{y} is continuous and \code{x=NULL}, a violin plot is generated.
#' If \code{x} is categorical, a grouped violin plot will be generated, with one violin for each level of \code{x}.
#' If \code{x} is continuous, a scatter plot will be generated.
#'
#' If \code{y} is categorical and \code{x} is continuous, horizontal violin plots will be generated.
#' If \code{x} is missing or categorical, rectangule plots will be generated where the area of a rectangle is proportional to the number of points for a combination of factors.
#'
#' Note that \code{plotPhenoData} and \code{plotCellData} are synonyms for \code{plotColData}.
#' These are artifacts of the transition from the old SCESet class, and will be deprecated in future releases.
#'
#' @return A ggplot object.
#'
#' @name plotColData
#' @rdname plotColData
#' @export
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
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
    .Deprecated(new="plotColData")
    plotColData(...)
}

#' @rdname plotColData 
#' @export
plotCellData <- function(...) {
    .Deprecated(new="plotColData")
    plotColData(...)
}

#' Plot row metadata
#'
#' Plot row-level (i.e., gene) metadata from a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object containing expression values and experimental information.
#' @param y Specification of the row-level metadata to show on the y-axis, see \code{?"\link{scater-vis-var}"} for possible values.
#' Note that only metadata fields will be searched, \code{assays} will not be used.
#' @param x Specification of the row-level metadata to show on the x-axis, see \code{?"\link{scater-vis-var}"} for possible values.
#' Again, only metadata fields will be searched, \code{assays} will not be used.
#' @param colour_by Specification of a row metadata field or a cell to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a row metadata field or a cell to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a row metadata field or a cell to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details 
#' If \code{y} is continuous and \code{x=NULL}, a violin plot is generated.
#' If \code{x} is categorical, a grouped violin plot will be generated, with one violin for each level of \code{x}.
#' If \code{x} is continuous, a scatter plot will be generated.
#'
#' If \code{y} is categorical and \code{x} is continuous, horizontal violin plots will be generated.
#' If \code{x} is missing or categorical, rectangule plots will be generated where the area of a rectangle is proportional to the number of points for a combination of factors.
#'
#' Note that \code{plotFeatureData} is a synonym for \code{plotRowData}.
#' This is an artifact of the transition from the old SCESet class, and will be deprecated in future releases.
#'
#' @return A ggplot object.
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
#' plotRowData(example_sce, y="n_cells_by_counts", x="log10_total_counts")
#' plotRowData(example_sce, y="n_cells_by_counts", 
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
    .Deprecated(new="plotRowData")
    plotRowData(...)
}


.metadata_dispatcher <- function(object, mode, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", legend = "auto", ...)
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
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, mode = mode,
                                       colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                       by_exprs_values = exprs_values, legend = legend)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
    legend <- vis_out$legend

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
                     colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
                     legend=legend, ...)
}


