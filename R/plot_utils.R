#' Variable selection for visualization
#'
#' A number of \pkg{scater} functions accept a SingleCellExperiment object and extract (meta)data from it for use in a plot.
#' These values are then used on the x- or y-axes (e.g., \code{\link{plotColData}}) or for tuning visual parameters, e.g.,
#' \code{colour_by}, \code{shape_by}, \code{size_by}.
#' This page describes how the selection of these values can be controlled by the user,
#' by passing appropriate values to the arguments of the desired plotting function.
#' 
#' @section When plotting by cells:
#' Here, we assume that each visual feature of interest (e.g., point or line) corresponds to a cell in the SingleCellExperiment object \code{sce}.
#' We will also assume that the user wants to change the colour of each feature according to the cell (meta)data.
#' To do so, the user can pass as an argument:
#' \itemize{
#' \item An unnamed character vector of length 1, i.e., a string. 
#' This is initially assumed to be the name of a column-level metadata field.
#' The function will first search the column names of \code{colData(sce)}, and extract metadata for all cells if a matching field is found.
#' If no match is found, the function will assume that the string represents a gene name.
#' It will search \code{rownames(sce)} and extract gene expression values for any matching row across all cells.
#' Otherwise, an error is raised.
#' \item A named character vector of length 1, where the name is either \code{"exprs"} or \code{"metadata"}.
#' This forces the function to only search for the string in \code{rownames(sce)} or \code{colnames(colData(sce))}, respectively.
#' Adding an explicit name is useful when the same field exists in both the row names and column metadata names.
#' \item A character vector of length greater than 1.
#' This will search for nested fields in \code{colData(sce)}.
#' For example, supplying a character vector \code{c("A", "B", "C")} will retrieve \code{colData(sce)$A$B$C}, where both \code{A} and \code{B} contain nested DataFrames.
#' See \code{\link{calculateQCMetrics}} with \code{compact=TRUE} for an example of how these can be constructed. 
#' The concatenated name \code{"A:B:C"} will be used in the legend.
#' \item A character vector of length greater than 1 and the first element set to \code{NA}.
#' This will search for nested fields in the internal column data of a \linkS4class{SingleCellExperiment}, i.e., in \code{\link{int_colData}}.
#' For example, \code{c(NA, "size_factor")} would retrieve the values corresponding to \code{sizeFactors(object)}.
#' The concatenated name without the \code{NA} is used in the legend.
#' Note that internal fields are \emph{only} searched when \code{NA} is the first element.
#' \item A data frame with one column and number of rows equal to the number of cells.
#' This should contain values to use for visualization, e.g., for plotting on the x-/y-axis, or for colouring by.
#' In this manner, the user can use new information without manually adding it to the SingleCellExperiment object.
#' The column name of the data frame will be used in the legend.
#' }
#'
#' The same logic applies for other visualization parameters such as \code{shape_by} and \code{size_by}.
#' Other arguments may also use the same scheme, but this depends on the context; see the documentation for each function for details.
#' In particular, if an argument explicitly refers to a metadata field, any names for the character string will be ignored.
#' Similarly, a character vector of length > 1 is not allowed for an argument that explicitly refers to expression values.
#' 
#' @section When plotting by features:
#' Here, we assume that each visual feature of interest (e.g., point or line) corresponds to a feature in the SingleCellExperiment object \code{sce}.
#' The scheme is mostly the same as described above, with a few differences:
#' \itemize{
#' \item \code{rowData} is searched instead of \code{colData}, as we are extracting metadata for each feature.
#' \item When extracting expression values, the name of a single cell must be specified.
#' Visualization will then use the expression profile for all features in that cell.
#' (This tends to be a rather unusual choice for colouring.)
#' \item Character strings named with \code{"exprs"} will search for the string in \code{colnames(sce)}.
#' \item A data frame input should have number of rows equal to the number of features.
#' }
#'
#' @section Miscellaneous details:
#' Most functions will have a \code{by_exprs_values} parameter. 
#' This defines the assay of the SingleCellExperiment object from which expression values are extracted for use in colouring, shaping or sizing the points.
#' The setting of \code{by_exprs_values} will usually default to \code{"logcounts"}, or to the value of \code{exprs_values} in functions such as \code{\link{plotExpression}}.
#' However, it can be specified separately from \code{exprs_values}, which is useful for visualizing two different types of expression values on the same plot.
#'
#' Most functions will also have a \code{by_show_single} parameter.
#' If \code{FALSE}, variables with only one level are not used for visualization, i.e., the visual aspect (colour or shape or size) is set to the default for all points.
#' No guide is created for this aspect, avoiding clutter in the legend when that aspect provides no information.
#' If \code{TRUE}, all supplied variables are used for visualization, regardless of how many levels they have.
#'
#' @name scater-vis-var
#'
#' @seealso
#' \code{\link{plotColData}}, 
#' \code{\link{plotRowData}}, 
#' \code{\link{plotReducedDim}}, 
#' \code{\link{plotExpression}}, 
#' \code{\link{plotPlatePosition}},
#' and most other plotting functions.
NULL

#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom SingleCellExperiment int_elementMetadata int_colData
#' @importFrom BiocGenerics rownames colnames
.choose_vis_values <- function(x, by, mode=c("column", "row"), search=c("any", "metadata", "exprs"),
    exprs_values = "logcounts", coerce_factor = FALSE, level_limit = Inf, discard_solo = FALSE) 
# This function looks through the visualization data and returns the
# values to be visualized. Either 'by' itself, or a column of colData,
# or a column of rowData, or the expression values of a feature.
{
    vals <- NULL
    if (is.character(by) && length(by) > 0) { 
        mode <- match.arg(mode)
        search <- match.arg(search)
        internal_only <- FALSE

        # Determining what to check, based on input 'by'.
        if (length(by)>1) {
            if (search=="exprs") {
                stop("character vector of length > 1 not allowed for search='exprs'")
            }
            search <- "metadata"

            if (is.na(by[1])) { 
                internal_only <- TRUE
                by <- by[-1]
            }
        } else if (search=="any") { 
            cur_name <- names(by)
            if (!is.null(cur_name) && !is.na(cur_name)) { 
                if (cur_name=="metadata" || cur_name =="exprs") {
                    search <- cur_name
                } 
            }
        }
        names(by) <- NULL
           
        # Checking the metadata; note the loop to account for nesting.
        if (search=="any" || search=="metadata") {
            if (!internal_only) { 
                meta_data <- if (mode=="column") colData(x) else rowData(x)
            } else {
                meta_data <- if (mode=="column") int_colData(x) else int_elementMetadata(x)
            }
            for (field in by) {
                if (!field %in% colnames(meta_data)) {
                    break
                }
                vals <- meta_data[[field]]
                meta_data <- vals
            }
            by <- paste(by, collapse=":") # collapsing to a single string for output.
        }

        # Metadata takes priority, so we don't bother searching if 'vals' is non-NULL.
        if ((search=="any" || search=="exprs") && is.null(vals)) { 
            exprs <- assay(x, i = exprs_values, withDimNames=FALSE)
            if (mode=="column") {
                m <- match(by, rownames(x)) # coloring columns, so we take the row values.
                if (!is.na(m)) {
                    vals <- exprs[m,] 
                }
            } else if (mode=="row") {
                m <- match(by, colnames(x)) 
                if (!is.na(m)) {
                    vals <- exprs[,m] 
                }
            }
        }

        if (is.null(vals)) {
            stop(sprintf("cannot find '%s' in %s fields", by, search))
        }

    } else if (is.data.frame(by)) {
        if (ncol(by) != 1L) {
            stop("input data frame should only have one column")
        } else {
            if (mode=="column" && nrow(by) != ncol(x)) {
                stop("number of rows of input data frame should be equal to 'ncol(object)'")
            }
            if (mode=="row" && nrow(by) != nrow(x)) {
                stop("number of rows of input data frame should be equal to 'nrow(object)'")
            }
        }

        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)

    } else if (!is.null(by)) {
        # We have to allow by=NULL to pass through smoothly.
        stop("invalid value of 'by' supplied")
    }

    # Checking the level limit.
    if (coerce_factor && !is.null(vals)) {
        vals <- factor(vals)
        if (level_limit < nlevels(vals)) {
            stop(sprintf("number of unique levels exceeds %i", level_limit))
        }
    }

    # If only one level for the variable, set to NULL.
    if (length(unique(vals))<=1L && discard_solo) { 
        by <- NULL
        vals <- NULL
    }
    return(list(name = by, val = vals))
}

.incorporate_common_vis <- function(df, se, mode, colour_by, size_by, shape_by, by_exprs_values, by_show_single)
# A convenience wrapper to incorporate colour, size and shape arguments into the data.frame for plotting.
# Do NOT use the supplied names to name fields in 'df', as these may clash with internal names.
{
    ## check colour argument:
    colour_by_out <- .choose_vis_values(se, colour_by, mode = mode, search = "any", discard_solo = !by_show_single,
                                        exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    df$colour_by <- colour_by_out$val

    ## check shape argument (note the limiter):
    shape_by_out <- .choose_vis_values(se, shape_by, mode = mode, search = "any",  discard_solo = !by_show_single,
                                       exprs_values = by_exprs_values, coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    df$shape_by <- shape_by_out$val

    ## check size argument:
    size_by_out <- .choose_vis_values(se, size_by, mode = mode, search = "any", discard_solo = !by_show_single,
                                      exprs_values = by_exprs_values)
    size_by <- size_by_out$name
    df$size_by <- size_by_out$val

    return(list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by))
}

.incorporate_common_vis_col <- function(df, se, mode, colour_by, size_by, shape_by, by_exprs_values, by_show_single) {
    colour_by_out <- retrieveCellInfo(se, colour_by, assay.type = by_exprs_values)
    colour_by <- colour_by_out$name
    df$colour_by <- colour_by_out$val

    shape_by_out <- retrieveCellInfo(se, shape_by, assay.type = by_exprs_values)
    shape_by <- shape_by_out$name
    df$shape_by <- .coerce_to_factor(shape_by_out$value, 10, "shape_by")

    size_by_out <- retrieveCellInfo(se, size_by, assay.type = by_exprs_values)
    size_by <- size_by_out$name
    df$size_by <- size_by_out$val

    list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by)
}

################################################
## Creating pair plots.

.makePairs <- function(data_matrix) 
# with thanks to Gaston Sanchez, who posted this code online
# https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
{
    if ( is.null(names(data_matrix)) )
        names(data_matrix) <- paste0("row", 1:nrow(data_matrix))
    exp_grid <- expand.grid(x = 1:ncol(data_matrix), y = 1:ncol(data_matrix))
    exp_grid <- exp_grid[exp_grid$x != exp_grid$y,]
    all_panels <- do.call("rbind", lapply(1:nrow(exp_grid), function(i) {
        xcol <- exp_grid[i, "x"]
        ycol <- exp_grid[i, "y"]
        data.frame(xvar = names(data_matrix)[ycol], yvar = names(data_matrix)[xcol],
                   x = data_matrix[, xcol], y = data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels = names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels = names(data_matrix))
    densities <- do.call("rbind", lapply(1:ncol(data_matrix), function(i) {
        data.frame(xvar = names(data_matrix)[i], yvar = names(data_matrix)[i],
                   x = data_matrix[, i])
    }))
    list(all = all_panels, densities = densities)
}
