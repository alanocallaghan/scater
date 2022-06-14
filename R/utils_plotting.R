.incorporate_common_vis_row <- function(df, se, mode, colour_by, size_by, shape_by, 
    by_exprs_values, by_show_single, other_fields, multiplier = NULL) 
{
    colour_by_out <- retrieveFeatureInfo(se, colour_by, exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    if (!is.null(multiplier)) {
        colour_by_out$value<- colour_by_out$value[multiplier]
    }
    df$colour_by <- colour_by_out$value

    shape_by_out <- retrieveFeatureInfo(se, shape_by, exprs_values = by_exprs_values)
    shape_by <- shape_by_out$name
    if (!is.null(multiplier)) {
        shape_by_out$value<- shape_by_out$value[multiplier]
    }
    df$shape_by <- .coerce_to_factor(shape_by_out$value, 10, "shape_by")

    size_by_out <- retrieveFeatureInfo(se, size_by, exprs_values = by_exprs_values)
    size_by <- size_by_out$name
    if (!is.null(multiplier)) {
        size_by_out$value<- size_by_out$value[multiplier]
    }
    df$size_by <- size_by_out$value

    for (o in other_fields) {
        other <- retrieveFeatureInfo(se, o, exprs_values=by_exprs_values)
        if (!is.null(multiplier)) {
            other$value<- other$value[multiplier]
        }
        df <- .add_other_or_warn(df, other)
    }

    list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by)
}

.incorporate_common_vis_col <- function(df, se, mode, colour_by, size_by, shape_by, 
    by_exprs_values, by_show_single, other_fields, multiplier = NULL, swap_rownames = NULL) 
{
    colour_by_out <- retrieveCellInfo(se, colour_by, exprs_values = by_exprs_values,
        swap_rownames=swap_rownames)
    colour_by <- colour_by_out$name
    if (!is.null(multiplier)) {
        colour_by_out$value <- colour_by_out$value[multiplier]
    }
    df$colour_by <- colour_by_out$value

    shape_by_out <- retrieveCellInfo(se, shape_by, exprs_values = by_exprs_values,
        swap_rownames=swap_rownames)
    shape_by <- shape_by_out$name
    if (!is.null(multiplier)) {
        shape_by_out$value <- shape_by_out$value[multiplier]
    }
    df$shape_by <- .coerce_to_factor(shape_by_out$value, 10, "shape_by")

    size_by_out <- retrieveCellInfo(se, size_by, exprs_values = by_exprs_values,
        swap_rownames=swap_rownames)
    size_by <- size_by_out$name
    if (!is.null(multiplier)) {
        size_by_out$value <- size_by_out$value[multiplier]
    }
    df$size_by <- size_by_out$value

    for (o in other_fields) {
        other <- retrieveCellInfo(se, o, exprs_values=by_exprs_values,
            swap_rownames=swap_rownames)
        if (!is.null(multiplier)) {
            other$value <- other$value[multiplier]
        }
        df <- .add_other_or_warn(df, other)
    }

    list(df = df, colour_by = colour_by, shape_by = shape_by, size_by = size_by)
}

.add_other_or_warn <- function(df, other) {
    if (!is.null(other$name)) {
        if (!other$name %in% colnames(df)) {
            df[[other$name]] <- other$value
        } else {
            warning(sprintf("not adding duplicated '%s' from 'other_fields'", other$name))
        }
    }
    df
}

################################################
## Creating pair plots.

.makePairs <- function(data_matrix) 
# with thanks to Gaston Sanchez, who posted this code online
# https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
{
    if ( is.null(names(data_matrix)) )
        names(data_matrix) <- paste0("row", seq_len(nrow(data_matrix)))
    exp_grid <- expand.grid(x = seq_len(ncol(data_matrix)), y = seq_len(ncol(data_matrix)))
    exp_grid <- exp_grid[exp_grid$x != exp_grid$y,]
    all_panels <- do.call("rbind", lapply(seq_len(nrow(exp_grid)), function(i) {
        xcol <- exp_grid[i, "x"]
        ycol <- exp_grid[i, "y"]
        data.frame(xvar = names(data_matrix)[ycol], yvar = names(data_matrix)[xcol],
                   x = data_matrix[, xcol], y = data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels = names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels = names(data_matrix))
    densities <- do.call("rbind", lapply(seq_len(ncol(data_matrix)), function(i) {
        data.frame(xvar = names(data_matrix)[i], yvar = names(data_matrix)[i],
                   x = data_matrix[, i])
    }))
    list(all = all_panels, densities = densities)
}

################################################
## Using non-standard gene IDs
.swap_rownames <- function(object, swap_rownames = NULL) {
    if (is.null(swap_rownames)) {
        return(object)
    }
    rownames(object) <- .get_rowData_column(object, swap_rownames)
    object
}

.handle_features <- function(features, object) {
    if (is.logical(features)) {
        if (length(features) != nrow(object)) {
            stop("logical features index should be the of length nrow(object).")
        }
        features <- rownames(object)[features]
    } else if (is.numeric(features)) {
        if (any(features > nrow(object)) | any(features <= 0) | any(features != round(features))) {
            stop("All features should be round numbers; > 0 and <= nrow(object).")
        }
        features <- rownames(object)[features]
    } else if (is.character(features) | is.factor(features)) {
        if (!all(as.character(features) %in% rownames(object))) {
            stop("Some features not in input object.")
        }
    }
    # if it's character by now, preserve the order
    if (is.character(features)) {
        features <- factor(features, levels = features)
    }
    # otherwise it's presumably a factor and should keep its levels
    features
}

.get_rowData_column <- function(object, column) {
    if (!column %in% colnames(rowData(object))) {
        stop("Cannot find column ", column, " in rowData")
    }
    rowData(object)[[column]]
}
