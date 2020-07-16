#' Plot heatmap of gene expression values
#'
#' Create a heatmap of expression values for each cell and specified features in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param features A character vector of row names, a logical vector of integer vector of indices specifying rows of \code{object} to show in the heatmap.
#' @param columns A vector specifying the subset of columns in \code{object} to show as columns in the heatmap. 
#' Also specifies the column order if \code{cluster_cols=FALSE} and \code{order_columns_by=NULL}.
#' By default, all columns are used.
#' @param exprs_values A string or integer scalar indicating which assay of \code{object} should be used as expression values for colouring in the heatmap.
#' @param center A logical scalar indicating whether each row should have its mean expression centered at zero prior to plotting. 
#' @param zlim A numeric vector of length 2, specifying the upper and lower bounds for the expression values. 
#' This winsorizes the expression matrix prior to plotting (but after centering, if \code{center=TRUE}). 
#' If \code{NULL}, it defaults to the range of the expression matrix.
#' @param symmetric A logical scalar specifying whether the default \code{zlim} should be symmetric around zero. 
#' If \code{TRUE}, the maximum absolute value of \code{zlim} will be computed and multiplied by \code{c(-1, 1)} to redefine \code{zlim}.
#' @param color A vector of colours specifying the palette to use for mapping expression values to colours. 
#' This defaults to the default setting in \code{\link[pheatmap]{pheatmap}}.
#' @param colour_columns_by A list of values specifying how the columns should be annotated with colours.
#' Each entry of the list can be any acceptable input to the \code{by} argument in \code{?\link{retrieveCellInfo}}.
#' A character vector can also be supplied and will be treated as a list of strings.
#' @param column_annotation_colors Passed to \code{\link[pheatmap]{pheatmap}} 
#' as the \code{annotation_colours} argument.
#' @param order_columns_by A list of values specifying how the columns should be ordered.
#' Each entry of the list can be any acceptable input to the \code{by} argument in \code{?\link{retrieveCellInfo}}.
#' A character vector can also be supplied and will be treated as a list of strings.
#' This argument is automatically appended to \code{colour_columns_by}.
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for colouring of column-level data - see the \code{exprs_values} argument in \code{?\link{retrieveCellInfo}}.
#' @param show_colnames,cluster_cols,... Additional arguments to pass to \code{\link[pheatmap]{pheatmap}}.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to 
#'  identify features instead of \code{rownames(object)} when labelling plot 
#'  elements.
#'
#' @details 
#' Setting \code{center=TRUE} is useful for examining log-fold changes of each cell's expression profile from the average across all cells.
#' This avoids issues with the entire row appearing a certain colour because the gene is highly/lowly expressed across all cells.
#'
#' Setting \code{zlim} preserves the dynamic range of colours in the presence  of outliers. 
#' Otherwise, the plot may be dominated by a few genes, which will \dQuote{flatten} the observed colours for the rest of the heatmap.
#'
#' Setting \code{order_columns_by} is useful for automatically ordering the heatmap by one or more factors of interest, e.g., cluster identity.
#' This the need to set \code{colour_columns_by}, \code{cluster_cols} and \code{columns} to achieve the same effect.
#' @return A heatmap is produced on the current graphics device. 
#' The output of \code{\link[pheatmap]{pheatmap}} is invisibly returned.
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @author Aaron Lun
#' 
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10])
#'
#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     center=TRUE, symmetric=TRUE)
#'
#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     colour_columns_by=c("Mutation_Status", "Cell_Cycle"))
#'
#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom viridis viridis
#' @importFrom SummarizedExperiment assay assayNames
plotHeatmap <- function(object, features, columns = NULL,
    exprs_values = "logcounts", center = FALSE, zlim = NULL, symmetric = FALSE,
    color = NULL, colour_columns_by = NULL, column_annotation_colors = list(),
    order_columns_by = NULL, by_exprs_values = exprs_values, 
    show_colnames = FALSE, cluster_cols = is.null(order_columns_by),
    swap_rownames = NULL, ...) 
{
    # Setting names, otherwise the downstream colouring fails.
    if (is.null(colnames(object))) { 
        colnames(object) <- seq_len(ncol(object)) 
    }

    # Pulling out the features. swap_rownames makes features index a rowdata col
    feats <- .swap_rownames(object, features, swap_rownames)
    heat.vals <- assay(object, exprs_values)[feats, , drop=FALSE]
    rownames(heat.vals) <- features
    if (!is.null(columns)) {
        columns <- .subset2index(columns, object, byrow=FALSE)
        heat.vals <- heat.vals[,columns,drop=FALSE]
    }
    if (center) {
        heat.vals <- heat.vals - rowMeans2(DelayedArray(heat.vals))
    }

    if (!is.null(order_columns_by)) {
        ordering <- list()
        for (i in seq_along(order_columns_by)) {
            vals <- retrieveCellInfo(object, order_columns_by[[i]], exprs_values = by_exprs_values)$val
            if (!is.null(columns)) {
                vals <- vals[columns]
            }
            ordering[[i]] <- vals
        }
        heat.vals <- heat.vals[,do.call(order, ordering),drop=FALSE]
        cluster_cols <- FALSE
        colour_columns_by <- c(colour_columns_by, order_columns_by)
    }

    # Winsorizing to preserve the dynamic range of colours.
    if (is.null(zlim)) {
        zlim <- range(heat.vals)
    }
    if (symmetric) {
        extreme <- max(abs(zlim))
        zlim <- c(-extreme, extreme)
    }
    heat.vals[heat.vals < zlim[1]] <- zlim[1]
    heat.vals[heat.vals > zlim[2]] <- zlim[2]

    if (is.null(color)) {
        color <- eval(formals(pheatmap::pheatmap)$color, envir=environment(pheatmap::pheatmap))
    }
    color.breaks <- seq(zlim[1], zlim[2], length.out=length(color)+1L) 

    # Collecting variables to colour_by.
    if (length(colour_columns_by)) {
        column_variables <- list()

        for (i in seq_along(colour_columns_by)) { 
            field <- colour_columns_by[[i]]
            colour_by_out <- retrieveCellInfo(object, field,
                exprs_values = by_exprs_values, swap_rownames = swap_rownames)

            if (is.null(colour_by_out$val)) { 
                next
            } else if (is.numeric(colour_by_out$val)) { 
                colour_fac <- colour_by_out$val
                col_scale <- viridis(25)
            } else {
                colour_fac <- as.factor(colour_by_out$val)

                nlevs_colour_by <- nlevels(colour_fac)
                if (nlevs_colour_by <= 10) {
                    col_scale <- .get_palette("tableau10medium")
                } else if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                    col_scale <- .get_palette("tableau20") 
                } else {
                    col_scale <- viridis(nlevs_colour_by)                    
                }

                col_scale <- col_scale[seq_len(nlevs_colour_by)]
                names(col_scale) <- levels(colour_fac)
            }

            col_name <- colour_by_out$name
            if (col_name=="") {
                col_name <- paste0("unnamed", i)
            }
            column_variables[[col_name]] <- colour_fac
            if (is.null(column_annotation_colors[[col_name]])) {
                column_annotation_colors[[col_name]] <- col_scale
            }
        }

        # No need to subset for 'columns' or 'order_columns_by',
        # as pheatmap::pheatmap uses the rownames to handle this for us.
        column_variables <- do.call(data.frame,
            c(column_variables, list(row.names=colnames(object))))
        column_annotation_colors <- column_annotation_colors[as.character(colour_columns_by)]
    } else {
        column_variables <- column_annotation_colors <- NULL
    }

    # Creating the heatmap as specified.
    pheatmap::pheatmap(heat.vals, color=color, breaks=color.breaks, 
        annotation_col=column_variables, annotation_colors=column_annotation_colors, 
        show_colnames=show_colnames, cluster_cols=cluster_cols, ...) 
}
