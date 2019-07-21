#' Plot column metadata
#'
#' Plot column-level (i.e., cell) metadata in an SingleCellExperiment object.
#'
#' @param object A \linkS4class{SingleCellExperiment} object containing expression values and experimental information.
#' @param y String specifying the column-level metadata field to show on the y-axis.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' @param x String specifying the column-level metadata to show on the x-axis.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' If \code{NULL}, nothing is shown on the x-axis.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see \code{?\link{retrieveCellInfo}} for details.
#' @param by_show_single Deprecated and ignored.
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
plotColData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        by_exprs_values = "logcounts", by_show_single = FALSE,
                        ...)
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define dataframe to pass to plotMetadata
    x_by_out <- retrieveCellInfo(object, x, search = "colData")
    x_lab <- x_by_out$name
    y_by_out <- retrieveCellInfo(object, y, search = "colData")
    y_lab <- y_by_out$name
    if (!is.null(x)) {
        df_to_plot <- data.frame(X=x_by_out$val, Y=y_by_out$val)
    } else {
        num <- ifelse(mode=="row", nrow(object), ncol(object))
        df_to_plot <- data.frame(Y=y_by_out$val, X=factor(character(num)))
    }

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis_col(df_to_plot, se = object, 
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        by_exprs_values = by_exprs_values)

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
        colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
        ..., point_FUN=NULL)
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
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see \code{?"\link{scater-vis-var}"} for details.
#' @param by_show_single Logical scalar specifying whether single-level factors should be used for point aesthetics, see \code{?"\link{scater-vis-var}"} for details.
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
                        by_exprs_values = "logcounts", by_show_single = FALSE,
                        ...)
{
    .metadata_dispatcher(object, mode = "row", y = y, x = x,
                         colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                         by_exprs_values = by_exprs_values, by_show_single = by_show_single,
                         ...)
}

.metadata_dispatcher <- function(object, mode, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        by_exprs_values = "logcounts", by_show_single = FALSE,
                        ...)
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
        num <- ifelse(mode=="row", nrow(object), ncol(object))
        df_to_plot <- data.frame(Y=y_by_out$val, X=factor(character(num)))
    }

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, 
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        by_exprs_values = by_exprs_values)

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
                     colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
                     ..., point_FUN=NULL)
}
