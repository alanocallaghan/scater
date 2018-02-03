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
        } else if (nrow(by) != ncol(x)) {
            stop("'nrow(*_by)' should be equal to number of columns in 'x'")
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
