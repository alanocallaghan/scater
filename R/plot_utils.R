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
#' To do so, the user can pass to \code{colour_by}:
#' \itemize{
#' \item An unnamed character string.
#' This is initially assumed to be the name of a column-level metadata field.
#' The function will first search the column names of \code{colData(sce)}, and extract metadata for all cells if a matching field is found.
#' If no match is found, the function will assume that the string represents a gene name.
#' It will search \code{rownames(sce)} and extract gene expression values for any matching row across all cells.
#' Otherwise, an error is raised.
#' \item A named character string, where the name is either \code{"Feature"} or \code{"Metadata"}.
#' This forces the function to only search for the name in \code{rownames(sce)} or \code{colnames(colData(sce))}, respectively.
#' Adding an explicit name is useful when the same field exists in both the row names and column metadata names.
#' \item A character vector of length greater than 1.
#' This will search for nested fields in \code{colData(sce)}.
#' For example, supplying a character vector \code{c("A", "B", "C")} will retrieve \code{colData(sce)$A$B$C}, 
#' where both \code{A} and \code{B} contain nested DataFrames.
#' See \code{\link{calculateQCMetrics}} with \code{compact=TRUE} for an example of how these can be formed.
#' \item A data frame with one column and number of rows equal to the number of cells.
#' This should contain values to use for visualization (in this case, for colouring by).
#' In this manner, the user can use new information without manually adding it to the SingleCellExperiment object.
#' The column name of the data frame will be used in the legend.
#' }
#' Of course, the same logic applies for other visualization parameters such as \code{shape_by} and \code{size_by}.
#' Other arguments may also use the same scheme, but this depends on the context; see the documentation for each function for details.
#' 
#' @section When plotting by features:
#' Here, we assume that each visual feature of interest (e.g., point or line) corresponds to a feature in the SingleCellExperiment object \code{sce}.
#' The scheme is mostly the same as described above, with a few differences:
#' \itemize{
#' \item \code{rowData} is used instead of \code{colData}, as we are extracting metadata for each feature.
#' \item When extracting expression values, the name of a single cell must be specified.
#' Visualization will then use the expression profile for all features in that cell.
#' This tends to be a rather unusual choice for colouring, but we will not judge.
#' \item Named character strings should use \code{"Cell"} instead of \code{"Feature"}.
#' \item A data frame input should have number of rows equal to the number of features.
#' }
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

.choose_vis_values <- function(x, by, mode=c("column", "row"), search=c("any", "metadata", "feature"),
                               exprs_values = "logcounts", coerce_factor = FALSE, level_limit = NA,
                               discard_solo = FALSE) 
# This function looks through the visualization data and returns the
# values to be visualized. Either 'by' itself, or a column of colData,
# or a column of rowData, or the expression values of a feature.
{
    vals <- NULL
    if (is.character(by)) {
        mode <- match.arg(mode)
        search <- match.arg(search)

        # Determining what to check, based on input 'by'.
        if (search=="any") { 
            check_metadata <- check_features <- TRUE

            if (length(by)==0) {
                check_metadata <- FALSE 
                check_features <- FALSE
            } else if (length(by)>1) {
                check_metadata <- TRUE
                check_features <- FALSE
            } else if (length(by)==1L) {
                cur_name <- names(by)
                if (!is.null(cur_name) && !is.na(cur_name)) { 
                    if (cur_name=="Metadata") {
                        check_features <- FALSE
                        check_metadata <- TRUE
                    } else if ((mode=="column" && cur_name=="Feature")
                               || (mode=="row" && cur_name=="Cell")) {
                        check_features <- TRUE
                        check_metadata <- FALSE
                    } 
                }
            }
        } else {
            check_metadata <- (search=="metadata")
            check_features <- !check_metadata            
        }
           
        # Checking the metadata; note the loop to account for nesting.
        if (check_metadata) { 
            if (mode=="column") {
                meta_data <- colData(x)
            } else {
                meta_data <- rowData(x)
            }
            
            for (field in by) {
                if (!field %in% colnames(meta_data)) {
                    break
                }
                vals <- meta_data[[field]]
                metadata <- vals
            }
            by <- paste(by, collapse=":") # collapsing to a single string for output.
        }

        # Metadata takes priority, so we don't bother searching if 'vals' is non-NULL.
        if (check_features) {
            if (is.null(vals)){
                exprs <- assay(x, i = exprs_values)
                if (mode=="column") {
                    vals <- exprs[by,] # coloring columns, so we take the row values.
                } else {
                    vals <- exprs[,by]
                }
            }
        }

        if (is.null(vals) && (check_metadata || check_features)) {
            stop("cannot find the supplied '*_by' in features or metadata")
        }
    } else if (is.data.frame(by)) {
        if (ncol(by) != 1L) {
            stop("'*_by' should be a data frame with one column")
        } else {
            if (mode=="column" && nrow(by) != ncol(x)) {
                stop("'nrow(*_by)' should be equal to number of columns in 'x'")
            }
            if (mode=="row" && nrow(by) != nrow(x)) {
                stop("'nrow(*_by)' should be equal to number of rows in 'x'")
            }
        }

        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)

    } else if (!is.null(by)) {
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

.incorporate_common_vis <- function(df, se, mode, colour_by, size_by, shape_by, by_exprs_values, legend='auto') 
# A convenience wrapper to incorporate colour, size and shape arguments into the data.frame for plotting.
# Do NOT use the supplied names to name fields in 'df', as these may clash with internal names.
{
    legend <- match.arg(legend, c("auto", "none", "all"))
    discard_solo <- legend=="auto"

    ## check colour argument:
    colour_by_out <- .choose_vis_values(se, colour_by, mode = mode, search = "any", discard_solo = discard_solo,
                                        exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    df$colour_by <- colour_by_out$val

    ## check shape argument (note the limiter):
    shape_by_out <- .choose_vis_values(se, shape_by, mode = mode, search = "any",  discard_solo = discard_solo,
                                       exprs_values = by_exprs_values, coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    df$shape_by <- shape_by_out$val

    ## check size argument:
    size_by_out <- .choose_vis_values(se, size_by, mode = mode, search = "any", discard_solo = discard_solo,
                                      exprs_values = by_exprs_values)
    size_by <- size_by_out$name
    df$size_by <- size_by_out$val

    return(list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by, legend = legend))
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
