################################################################################
### Plot cells in plate positions

#' Plot cells in plate positions
#'
#' Plots cells in their position on a plate, coloured by phenotype data or
#' feature expression.
#'
#' @param object an \code{SingleCellExperiment} object. If \code{object$plate_position} is not
#' \code{NULL}, then this will be used to define each cell's position on the
#' plate, unless the \code{plate_position} argument is specified.
#' @param plate_position optional character vector providing a position on the
#' plate for each cell (e.g. A01, B12, etc, where letter indicates row and
#' number indicates column). Specifying this argument overrides any plate
#' position information extracted from the SingleCellExperiment object.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, a
#' data frame with one column containing values to map to colours for all cells.
#' @param x_position numeric vector providing x-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
#' @param y_position numeric vector providing y-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
#' @param exprs_values a string specifying the expression values to use for
#' colouring the points, if \code{colour_by} is set as a feature name.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#'
#' @details This function expects plate positions to be given in a charcter
#' format where a letter indicates the row on the plate and a numeric value
#' indicates the column. So each cell has a plate position such as "A01", "B12",
#' "K24" and so on. From these plate positions, the row is extracted as the
#' letter, and the column as the numeric part. If \code{object$plate_position}
#' or the \code{plate_position} argument are used to define plate positions,
#' then positions should be provided in this format. Alternatively, numeric
#' values to be used as x- and y-coordinates by supplying both the
#' \code{x_position} and \code{y_position} arguments to the function.
#'
#' @return
#' A ggplot object.
#'
#' @export
#'
#' @examples
#' ## prepare data
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## define plate positions
#' example_sce$plate_position <- paste0(
#' rep(LETTERS[1:5], each = 8), rep(formatC(1:8, width = 2, flag = "0"), 5))
#'
#' ## plot plate positions
#' plotPlatePosition(example_sce, colour_by = "Mutation_Status")
#'
#' ## Must have exprs slot defined in object
#' plotPlatePosition(example_sce, colour_by = "Gene_0004")
#'
plotPlatePosition <- function(object, plate_position = NULL,
                              colour_by = NULL,
                              x_position = NULL, y_position = NULL,
                              exprs_values = "logcounts", theme_size = 24, legend = "auto") {
    ## check object is SingleCellExperiment object
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be of class SingleCellExperiment")

    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ## Checking colour validity
    colour_by_out <- .choose_vis_values(object, colour_by, cell_control_default = TRUE,
                                        check_features = TRUE, exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    ## obtain well positions
    if ( !is.null(plate_position) ) {
        if ( length(plate_position) != ncol(object) )
            stop("Supplied plate_position argument must have same length as number of columns of SingleCellExperiment object.")
        plate_position_char <- plate_position

    } else
        plate_position_char <- object$plate_position

    if ( is.null(plate_position_char) ) {
        if ( is.null(x_position) || is.null(y_position) )
            stop("If plate_position is NULL then both x_position and y_position must be supplied.")
        plate_position_x <- x_position
        plate_position_y <- y_position
    } else {
        plate_position_y <- gsub("[0-9]*", "", plate_position_char)
        plate_position_y <- factor(plate_position_y,
                                   rev(sort(unique(plate_position_y))))
        plate_position_x <- gsub("[A-Z]*", "", plate_position_char)
        plate_position_x <- ordered(as.integer(plate_position_x))
    }

    ## Define data.frame for plotting
    df_to_plot <- data.frame(plate_position_x, plate_position_y)
    if ( !is.null(plate_position_char) )
        df_to_plot[["plate_position_char"]] <- plate_position_char
    df_to_plot$colour_by <- colour_by_vals

    ## make the plot
    aesth <- aes(x = plate_position_x, y = plate_position_y, fill = colour_by)
    if ( !is.null(plate_position_char) )
        aesth$label <- as.symbol("plate_position_char")

    plot_out <- ggplot(df_to_plot, aesth) +
        geom_point(shape = 21, size = theme_size, colour = "gray50")
    if ( !is.null(plate_position_char) )
        plot_out <- plot_out + geom_text(colour = "gray90")
    ## make sure colours are nice
    plot_out <- .resolve_plot_colours(plot_out,
                                      df_to_plot$colour_by,
                                      colour_by, fill = TRUE)

    ## Define plotting theme
    plot_out <- plot_out + theme_bw(theme_size) +
        theme(axis.title = element_blank(), axis.ticks = element_blank(),
              legend.text = element_text(size = theme_size / 2),
              legend.title = element_text(size = theme_size / 2)) +
        guides(fill = guide_legend(override.aes = list(size = theme_size / 2)))
    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")

    ## return plot
    plot_out
}

################################################################################
### Plot expression against transcript length

#' Plot expression against transcript length
#'
#' Plot expression values from a \code{\link{SingleCellExperiment}} object
#' against transcript length values defined in the SingleCellExperiment object
#' or supplied as an argument.
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param tx_length transcript lengths to plot on the x-axis. Can be one of: (1)
#' the name of a column of \code{rowData(object)} containing the transcript length
#' values, or (2) the name of an element of \code{assays(object)} containing
#' a matrix of transcript length values, or (3) a numeric vector of length equal
#' to the number of rows of \code{object} (number of features).
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assays"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable by which to colour
#' expression values on the plot. Alternatively, a data frame with
#' one column, containing a value for each feature to map to a colour.
#' @param shape_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot. Alternatively, a data frame
#' with one column containing values to map to shapes.
#' @param size_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot. Alternatively, a data frame
#' with one column containing values to map to sizes.
#' @param xlab label for x-axis; if \code{NULL} (default), then \code{x} will be
#' used as the x-axis label
#' @param show_exprs_sd logical, show the standard deviation of expression
#' values for each feature on the plot
#' @param show_smooth logical, show a smoothed fit through the expression values
#'  on the plot
#' @param alpha numeric value between 0 (completely transparent) and 1 (completely
#' solid) defining how transparent plotted points (cells) should be.
#' Points are jittered horizontally if the x-axis value is categorical rather
#' than numeric to avoid overplotting.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10)
#' @param log2_values should the expression values be transformed to the
#' log2-scale for plotting (with an offset of 1 to avoid logging zeroes)?
#' @param size numeric scalar optionally providing size for points if
#' \code{size_by} argument is not given. Default is \code{NULL}, in which case
#' \pkg{ggplot2} default is used.
#' @param se logical, should standard errors be shown (default \code{TRUE}) for
#' the smoothed fit through the cells. (Ignored if \code{show_smooth} is \code{FALSE}).
#'
#' @return a ggplot object
#' @export
#'
#' @importFrom DelayedMatrixStats rowMedians
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' rd <- DataFrame(gene_id = rownames(sc_example_counts),
#'         feature_id = paste("feature", rep(1:500, each = 4), sep = "_"),
#'      median_tx_length = rnorm(2000, mean = 5000, sd = 500))
#' rownames(rd) <- rownames(sc_example_counts)
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts),
#' colData = sc_example_cell_info, rowData = rd)
#' example_sce <- normalize(example_sce)
#'
#' plotExprsVsTxLength(example_sce, "median_tx_length")
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE)
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE,
#' show_exprs_sd = TRUE)
#'
#' ## using matrix of tx length values in assays(object)
#' mat <- matrix(rnorm(ncol(example_sce) * nrow(example_sce), mean = 5000,
#'  sd = 500), nrow = nrow(example_sce))
#' dimnames(mat) <- dimnames(example_sce)
#' assay(example_sce, "tx_len") <- mat
#'
#' plotExprsVsTxLength(example_sce, "tx_len", show_smooth = TRUE,
#' show_exprs_sd = TRUE)
#'
#' ## using a vector of tx length values
#' plotExprsVsTxLength(example_sce, rnorm(2000, mean = 5000, sd = 500))
#'
plotExprsVsTxLength <- function(object, tx_length = "median_feat_eff_len",
                                exprs_values = "logcounts",
                                colour_by = NULL, shape_by = NULL,
                                size_by = NULL, xlab = NULL,
                                show_exprs_sd = FALSE,
                                show_smooth = FALSE, alpha = 0.6,
                                theme_size = 10, log2_values = FALSE, size = NULL,
                                se = TRUE) {
    ## Check object is an SingleCellExperiment object
    if ( !is(object, "SingleCellExperiment") )
        stop("object must be an SingleCellExperiment")

    tx_length_values <- rep(NA, nrow(object))
    ## Check arguments are valid
    if ( length(tx_length) == 1 ) {
        if ( tx_length %in% colnames(rowData(object)) )
            tx_length_values <- rowData(object)[[tx_length]]
        else {
            if ( tx_length %in% SummarizedExperiment::assayNames(object) ) {
                tx_length_mat <- assay(object, tx_length)
                tx_length_values <- DelayedMatrixStats::rowMedians(DelayedArray(tx_length_mat))
            } else
                stop("the argument 'tx_length' should specify a column of rowData(object) or an element of assayNames(object) [see names(assayNames(object))")
        }
    } else {
        if ( length(tx_length) != nrow(object) )
            stop("If tx_length is a vector it must have length equal to nrow(object).")
        else {
            if ( !is.numeric(tx_length) )
                stop("If a vector, tx_length must contain numeric values.")
            tx_length_values <- tx_length
        }
    }

    exprs_mat <- assay(object, exprs_values)
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else
        ylab <- paste0("Expression (", exprs_values, ")")

    ## compute mean expression and sd of expression values
    exprs_mean <- rowMeans(exprs_mat)
    exprs_sd <- sqrt(.rowVars(exprs_mat))

    df_to_plot <- data.frame(tx_length_values, exprs_mean, exprs_sd,
                             ymin = exprs_mean - exprs_sd,
                             ymax = exprs_mean + exprs_sd)

    ## check colour, size, shape arguments
    colour_by_out <- .choose_vis_values(object, colour_by, check_coldata = FALSE)
    colour_by <- colour_by_out$name
    if (!is.null(colour_by)) df_to_plot[[colour_by]] <- colour_by_out$val

    shape_by_out <- .choose_vis_values(object, shape_by, check_coldata = FALSE,
                                       coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    if (!is.null(shape_by)) df_to_plot[[shape_by]] <- shape_by_out$val

    size_by_out <- .choose_vis_values(object, size_by, check_coldata = FALSE)
    size_by <- size_by_out$name
    if (!is.null(size_by)) df_to_plot[[size_by]] <- size_by_out$val

    ## Construct a ggplot2 aesthetic for the plot
    aesth <- aes()
    aesth$x <- as.symbol("tx_length_values")
    aesth$y <- as.symbol("exprs_mean")
    aesth$ymin <- as.symbol("ymin")
    aesth$ymax <- as.symbol("ymax")

    if ( !is.null(colour_by) )
        aesth$colour <- as.symbol(colour_by)
    if ( !is.null(shape_by) )
        aesth$shape <- as.symbol(shape_by)
    if ( !is.null(size_by) )
        aesth$size <- as.symbol(size_by)

    ## Define sensible x-axis label if NULL
    if ( is.null(xlab) )
        xlab <- "Median transcript length"

    ## Make the plot
    plot_out <- ggplot2::ggplot(df_to_plot, aesth) + xlab(xlab) + ylab(ylab)

    ## if colour aesthetic is defined, then choose sensible colour palette
    if ( !is.null(aesth$colour) )
        plot_out <- .resolve_plot_colours(plot_out,
                                          df_to_plot[[as.character(aesth$colour)]],
                                          as.character(aesth$colour))

    if ( is.null(aesth$size) & !is.null(size) ) {
        ## add SDs
        if ( show_exprs_sd )
            plot_out <- plot_out + geom_pointrange(size = size, alpha = 0.9 * alpha)
        ## add points to plot
        plot_out <- plot_out + geom_point(size = size, alpha = alpha)
    }  else {
        ## add SDs
        if ( show_exprs_sd )
            plot_out <- plot_out + geom_pointrange(alpha = 0.9 * alpha)
        ## add points to plot
        plot_out <- plot_out + geom_point(alpha = alpha)
    }

    ## show optional decorations on plot if desired
    if (show_smooth) {
        plot_out <- plot_out + stat_smooth(colour = "firebrick", linetype = 2,
                                           se = se)
    }

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    plot_out
}



