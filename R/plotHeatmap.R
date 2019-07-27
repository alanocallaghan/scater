#' Plot heatmap of gene expression values
#'
#' Create a heatmap of expression values for each cell and specified features in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param features A character vector of row names, a logical vector of integer vector of indices specifying rows of \code{object} to show in the heatmap.
#' @param columns A vector specifying the subset of columns in \code{object} to show as columns in the heatmp. 
#' By default, all columns are used in their original order.
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
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for colouring of column-level data - see the \code{exprs_values} argument in \code{?\link{retrieveCellInfo}}.
#' @param by_show_single Deprecated and ignored.
#' @param show_colnames Logical scalar specifying whether column names should be shown, if available in \code{object}.
#' @param ... Additional arguments to pass to \code{\link[pheatmap]{pheatmap}}.
#'
#' @details Setting \code{center=TRUE} is useful for examining log-fold changes of each cell's expression profile from the average across all cells.
#' This avoids issues with the entire row appearing a certain colour because the gene is highly/lowly expressed across all cells.
#'
#' Setting \code{zlim} preserves the dynamic range of colours in the presence  of outliers. 
#' Otherwise, the plot may be dominated by a few genes, which will \dQuote{flatten} the observed colours for the rest of the heatmap.
#'
#' @return A heatmap is produced on the current graphics device. 
#' The output of \code{\link[pheatmap]{pheatmap}} is invisibly returned.
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @author Aaron Lun
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
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
#' @importFrom SummarizedExperiment assay
plotHeatmap <- function(object, features, columns=NULL, exprs_values="logcounts",
    center=FALSE, zlim=NULL, symmetric=FALSE, color=NULL, 
    colour_columns_by=NULL, by_exprs_values = exprs_values, by_show_single=FALSE,
    show_colnames = TRUE, ...) 
{
    features_to_use <- .subset2index(features, object, byrow=TRUE)
    heat.vals <- assay(object, exprs_values, withDimnames=FALSE)[features_to_use,,drop=FALSE]

    rownames(heat.vals) <- rownames(object)[features_to_use]
    if (is.null(colnames(object))) { 
        colnames(heat.vals) <- seq_len(ncol(object)) # otherwise downstream colouring fails.
        show_colnames <- FALSE
    } else {
        colnames(heat.vals) <- colnames(object)
    }

    if (!is.null(columns)) { 
        heat.vals <- heat.vals[,columns,drop=FALSE]        
    }
    if (center) {
        heat.vals <- heat.vals - rowMeans2(DelayedArray(heat.vals))
    }

    # Winsorizing to preserve the dynamic range.
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
        column_variables <- column_colorings <- list()
        for (field in colour_columns_by) { 
            colour_by_out <- retrieveCellInfo(object, field, exprs_values = by_exprs_values)

            if (is.null(colour_by_out$val)) { 
                next
            } else if (is.numeric(colour_by_out$val)) { 
                colour_fac <- cut(colour_by_out$val, 25)
            } else {
                colour_fac <- as.factor(colour_by_out$val)
            } 

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

            column_variables[[colour_by_out$name]] <- colour_fac
            column_colorings[[colour_by_out$name]] <- col_scale
        }

        column_variables <- do.call(data.frame, c(column_variables, list(row.names=colnames(object))))
    } else {
        column_variables <- column_colorings <- NULL
    }

    # Creating the heatmap as specified.
    pheatmap::pheatmap(heat.vals, color=color, breaks=color.breaks, 
        annotation_col=column_variables, annotation_colors=column_colorings, 
        show_colnames=show_colnames, ...) 
}
