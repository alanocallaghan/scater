################################################################################
### Overview plot function for SingleCellExperiment

#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library accounted for by the most highly
#' expressed features for each cell for a \code{SingleCellExperiment} object. 
#'
#' @param x a \code{SingleCellExperiment} object
#' @param block1 character string defining the column of \code{colData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param block2 character string defining the column of \code{colData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param colour_by character string defining the column of \code{colData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively,
#' a data frame with one column containing a value for each cell, which will be
#' mapped to a corresponding colour.
#' @param nfeatures numeric scalar indicating the number of features to include
#' in the plot.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"counts"} (raw counts) [default], \code{"cpm"}
#' (counts per million), or \code{"fpkm"} (FPKM values).
#' @param linewidth numeric scalar giving the "size" parameter (in ggplot2
#' parlance) for the lines plotted. Default is 1.5.
#' @param ... arguments passed to \code{plotSCE}
#' @param ncol number of columns to use for \code{facet_wrap} if only one block is
#' defined.
#' @param theme_size numeric scalar giving font size to use for the plotting
#' theme
#'
#' @details Plots produced by this function are intended to provide an overview
#' of large-scale differences between cells. For each cell, the features are
#' ordered from most-expressed to least-expressed and the cumulative proportion
#' of the total expression for the cell is computed across the top
#' \code{nfeatures} features. These plots can flag cells with a very high
#' proportion of the library coming from a small number of features; such cells
#' are likely to be problematic for analyses. Using the colour and blocking
#' arguments can flag overall differences in cells under different experimental
#' conditions or affected by different batch and other variables.
#'
#' @return a ggplot plot object
#'
#' @importFrom dplyr mutate
#' @importFrom plyr aaply
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#'
#' plotScater(example_sce)
#' plotScater(example_sce, exprs_values = "counts", colour_by = "Cell_Cycle")
#' plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle")
#'
#' cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)
#' plotScater(example_sce, exprs_values = "cpm", block1 = "Treatment",
#' block2 = "Mutation_Status", colour_by = "Cell_Cycle")
#' # Error is thrown if chosen expression values are not available
#'
plotScater <- function(x, block1 = NULL, block2 = NULL, colour_by = NULL,
                    nfeatures = 500, exprs_values = "counts", ncol = 3,
                    linewidth = 1.5, theme_size = 10) {
    if (!is(x, "SingleCellExperiment"))
        stop("x must be of class SingleCellExperiment")
    if ( !is.null(block1) ) {
        if ( !(block1 %in% colnames(colData(x))) )
            stop("The block1 argument must either be NULL or a column of colData(x).")
    }
    if ( !is.null(block2) ) {
        if ( !(block2 %in% colnames(colData(x))) )
            stop("The block2 argument must either be NULL or a column of colData(x).")
    }

    ## Setting values to colour by.
    colour_by_out <- .choose_vis_values(x, colour_by)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    ## Define an expression matrix depending on which values we're using
    exprs_mat <- assay(x, i = exprs_values)
    nfeatures <- min(nfeatures, nrow(exprs_mat))

    ## Use C++ to get the sequencing real estate accounted for by features
    to_plot <- seq_len(nfeatures)
    ncells <- ncol(exprs_mat)
    seq_real_estate <- .Call(cxx_calc_top_features, exprs_mat, to_plot, NULL)
    seq_real_estate_long <- data.frame(Feature=rep(to_plot, each=ncells),
                                       Cell=rep(seq_len(ncells), nfeatures))
    seq_real_estate_long$Proportion_Library <- unlist(seq_real_estate) / 100

    ## Add block and colour_by information if provided
    if ( !is.null(block1) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block1 = as.factor(rep(x[[block1]],
                                                         each = nfeatures)))
    if ( !is.null(block2) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block2 = as.factor(rep(x[[block2]],
                                                         each = nfeatures)))
    if ( !is.null(colour_by) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, colour_by = rep(colour_by_vals,
                                                  each = nfeatures))

    ## Set up plot
    if ( is.null(colour_by) ) {
        plot_out <- ggplot(seq_real_estate_long,
                           aes_string(x = "Feature", y = "Proportion_Library",
                                      group = "Cell")) +
            geom_line(linetype = "solid", alpha = 0.3, size = linewidth)
    } else {
        plot_out <- ggplot(seq_real_estate_long,
                           aes_string(x = "Feature", y = "Proportion_Library",
                                      group = "Cell", colour = "colour_by")) +
            geom_line(linetype = "solid", alpha = 0.3, size = linewidth)
    }
    ## Deal with blocks for grid
    if ( !(is.null(block1) | is.null(block2)) )
        plot_out <- plot_out + facet_grid(block2 ~ block1)
    else {
        if ( !is.null(block1) && is.null(block2) ) {
            plot_out <- plot_out +
                facet_wrap(~block1, ncol = ncol)
        }
        if ( is.null(block1) && !is.null(block2) ) {
            plot_out <- plot_out +
                facet_wrap(~block2, ncol = ncol)
        }
    }
    ## Add extra plot theme and details
    if ( !is.null(seq_real_estate_long$colour_by) ) {
        plot_out <- .resolve_plot_colours(plot_out,
                                          seq_real_estate_long$colour_by,
                                          colour_by)
    }

    plot_out <- plot_out +
        xlab("Number of features") + ylab("Cumulative proportion of library")

    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot
    plot_out
}

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

#' Plot expression values for a set of features (e.g. genes or transcripts)
#'
#' @param object an \code{\link{SingleCellExperiment}} object containing expression values and
#' experimental information. Must have been appropriately prepared. For the
#' \code{plotExpressionDefault} method, the \code{object} argument is a
#' data.frame in 'long' format providing expression values for a set of features
#'  to plot, plus metadata used in the \code{aesth} argument, but this is not
#'  meant to be a user-level operation.
#' @param features a character vector of feature names or Boolean
#' vector or numeric vector of indices indicating which features should have
#' their expression values plotted
#' @param x character string providing a column name of \code{pData(object)} or
#' a feature name (i.e. gene or transcript) to plot on the x-axis in the
#' expression plot(s). If a feature name, then expression values for the feature
#' will be plotted on the x-axis for each subplot.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable by which to colour
#' expression values on the plot. Alternatively, a data frame with one column,
#' containing a value for each cell that will be mapped to a colour.
#' @param shape_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot. Alternatively, a data frame with
#' one column containing values to map to shapes.
#' @param size_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot. Alternatively, a data frame with
#' one column containing values to map to sizes.
#' @param ncol number of columns to be used for the panels of the plot
#' @param xlab label for x-axis; if \code{NULL} (default), then \code{x} will be
#' used as the x-axis label
#' @param show_median logical, show the median for each group on the plot
#' @param show_violin logical, show a violin plot for the distribution
#' for each group on the plot
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
#' @param scales character scalar, should scales be fixed ("fixed"),
#' free ("free"), or free in one dimension ("free_x"; "free_y", the default).
#' Passed to the \code{scales} argument in the \code{\link[ggplot2]{facet_wrap}} function
#' from the \code{ggplot2} package.
#' @param se logical, should standard errors be shown (default \code{TRUE}) for
#' the smoothed fit through the cells. (Ignored if \code{show_smooth} is \code{FALSE}).
#' @param jitter character scalar to define whether points are to be jittered
#' (\code{"jitter"}) or presented in a "beeswarm" style (if \code{"swarm"}; default).
#' "Beeswarm" style usually looks more attractive, but for datasets with a large
#' number of cells, or for dense plots, the jitter option may work better.
#' @param ... optional arguments (from those listed above) passed to
#' \code{plotExpressionDefault}
#'
#' @details Plot expression values (default log2(counts-per-million +
#' 1), if available) for a set of features.
#' @return a ggplot plot object
#'
#' @name plotExpression
#' @aliases plotExpression plotExpression,SingleCellExperiment-method plotExpression,data.frame-method
#' @import ggplot2
#' @importFrom Biobase varLabels
#' @importFrom Biobase featureNames
#' @importFrom ggbeeswarm geom_quasirandom
#' @export
#'
#' @examples
#' ## prepare data
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#  example_sce <- normalize(example_sce)
#' example_sce <- calculateQCMetrics(example_sce)
#' sizeFactors(example_sce) <- colSums(counts(example_sce))
#' example_sce <- normalize(example_sce)
#'
#' ## default plot
#' plotExpression(example_sce, 1:15)
#' plotExpression(example_sce, 1:15, jitter = "jitter")
#'
#' ## plot expression against an x-axis value
#' plotExpression(example_sce, 1:6, "Mutation_Status")
#'
#' ## explore options
#' plotExpression(example_sce, 1:6, x = "Mutation_Status", exprs_values = "logcounts",
#' colour_by = "Cell_Cycle", show_violin = TRUE, show_median = TRUE)
#' plotExpression(example_sce, 1:6, x = "Mutation_Status", exprs_values = "counts",
#' colour_by = "Cell_Cycle", show_violin = TRUE, show_median = TRUE)
#'
#' plotExpression(example_sce, "Gene_0001", x = "Mutation_Status")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Mutation_Status")
#' plotExpression(example_sce, "Gene_0001", x = "Gene_0002")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Gene_0002")
#' ## plot expression against expression values for Gene_0004
#' plotExpression(example_sce, 1:4, "Gene_0004")
#' plotExpression(example_sce, 1:4, "Gene_0004", show_smooth = TRUE)
#' plotExpression(example_sce, 1:4, "Gene_0004", show_smooth = TRUE, se = FALSE)
#'
plotExpression <- function(object, features, x = NULL,
                              exprs_values = "logcounts", log2_values = FALSE,
                              colour_by = NULL, shape_by = NULL, size_by = NULL,
                              ncol = 2, xlab = NULL, show_median = FALSE,
                           show_violin = TRUE, theme_size = 10, ...) {

    ## Define number of features to plot
    if (is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)

    if ( exprs_values == "exprs" && !(exprs_values %in% assayNames(object)) )
        exprs_values <- "logcounts"
    if ( !(exprs_values %in% assayNames(object)) )
        stop(paste(exprs_values, "not found in assays(object)."))
    exprs_mat <- assay(object, i = exprs_values)
    exprs_mat <- exprs_mat[features,,drop = FALSE]
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else
        ylab <- paste0("Expression (", exprs_values, ")")
    to_melt <- as.matrix(exprs_mat)

    ## check x-coordinates are valid
    xcoord <- NULL
    if ( !is.null(x) ) {
        if ( x %in% rownames(object) ) {
            xcoord <- assay(object, i = exprs_values)[x,]
            show_violin <- FALSE
            show_median <- FALSE
        } else if (x %in% colnames(colData(object))) {
            xcoord <- object[[x]]
        } else {
            stop("'x' should be a column of colData(object) or in rownames(object)")
        }
    }

    ## checking visualization arguments
    colour_by_out <- .choose_vis_values(object, colour_by)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    shape_by_out <- .choose_vis_values(object, shape_by, coerce_factor = TRUE,
                                       level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val

    size_by_out <- .choose_vis_values(object, size_by)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    ## Melt the expression data and metadata into a convenient form
    evals_long <- reshape2::melt(to_melt, value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")

    ## Prepare the samples information
    samps <- data.frame(row.names = colnames(object))
    if ( !is.null(xcoord) ) samps[[x]] <- xcoord

    ## Construct a ggplot2 aesthetic for the plot
    aesth <- aes()
    if ( is.null(x) ) {
        aesth$x <- as.symbol("Feature")
        one_facet <- TRUE
    } else {
        aesth$x <- as.symbol(x)
        one_facet <- FALSE
    }
    aesth$y <- as.symbol("evals")

    ## Define sensible x-axis label if NULL
    if ( is.null(xlab) )
        xlab <- x

    if ( !is.null(shape_by) ) {
        aesth$shape <- as.symbol(shape_by)
        samps[[shape_by]] <- shape_by_vals
    }
    if ( !is.null(size_by) ) {
        aesth$size <- as.symbol(size_by)
        samps[[size_by]] <- size_by_vals
    }
    if ( !is.null(colour_by) ) {
        aesth$colour <- as.symbol(colour_by)
        samps[[colour_by]] <- colour_by_vals
    }

    ## Extend the sample information, combine with the expression values
    samples_long <- samps[rep(seq_len(ncol(object)), each = nfeatures), , drop = FALSE]

    ## create plotting object
    object <- cbind(evals_long, samples_long)

    ## Make the plot
    plot_out <- plotExpressionDefault(object, aesth, ncol, xlab, ylab,
                                      show_median, show_violin, 
                                      one_facet = one_facet, ...)

    if ( is.null(x) ) { ## in this case, do not show x-axis ticks or labels
        plot_out <- plot_out + theme(
            axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(.03, .02, .05, .02), "npc"))
        if (is.null(colour_by))
            plot_out <- plot_out + guides(fill = "none", colour = "none")
    }
    plot_out
}


#' @param aesth an \code{aes} object to use in the call to \code{\link[ggplot2]{ggplot}}.
#' @param ylab character string defining a label for the y-axis (y-axes) of the
#' plot.
#' @param one_facet logical, should expression values for features be plotted in one facet
#' instead of mutiple facets, one per feature? Default if \code{x = NULL}.
#' @rdname plotExpression
#' @aliases plotExpression
#' @export
plotExpressionDefault <- function(object, aesth, ncol = 2, xlab = NULL,
                                  ylab = NULL, show_median = FALSE,
                                  show_violin = TRUE, show_smooth = FALSE,
                                  theme_size = 10,
                                  alpha = 0.6, size = NULL, scales = "fixed",
                                  one_facet = FALSE, se = TRUE, jitter = "swarm") {
    if ( !("Feature" %in% names(object)) )
        stop("object needs a column named 'Feature' to define the feature(s) by which to plot expression.")

    ## use x as group for violin plot if discrete
    group_by_x <- (show_violin &&
        (!is.numeric(object[[as.character(aesth$x)]]) ||
             nlevels(as.factor(object[[as.character(aesth$x)]])) <= 5))
    if (group_by_x)
        aesth$group <- aesth$x
    else
        aesth$group <- 1

    ## Define the plot
    if (one_facet) {
        if (is.null(aesth$colour))
            aesth$colour <- as.symbol("Feature")
        plot_out <- ggplot(object, aesth) + xlab(xlab) + ylab(ylab)
    } else {
        plot_out <- ggplot(object, aesth) +
            facet_wrap(~Feature, ncol = ncol, scales = scales) +
            xlab(xlab) +
            ylab(ylab)
    }

    ## if colour aesthetic is defined, then choose sensible colour palette
    if ( !is.null(aesth$colour) )
        plot_out <- .resolve_plot_colours(plot_out,
                                          object[[as.character(aesth$colour)]],
                                          as.character(aesth$colour))

    ## if x axis variable is not numeric, then jitter points horizontally
    if ( is.numeric(aesth$x) ) {
        if ( is.null(aesth$size) & !is.null(size) )
            plot_out <- plot_out + geom_point(size = size, alpha = alpha)
        else
            plot_out <- plot_out + geom_point(alpha = alpha)
    } else {
        if ( is.null(aesth$size) & !is.null(size) ) {
            if (jitter == "swarm")
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    alpha = alpha, size = size, groupOnX = TRUE)
            else
                plot_out <- plot_out + geom_jitter(
                    alpha = alpha, size = size, position = position_jitter(height = 0))
        }
        else {
            if (jitter == "swarm")
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    alpha = alpha, groupOnX = TRUE)
            else
                plot_out <- plot_out + geom_jitter(
                    alpha = alpha, position = position_jitter(height = 0))
        }
    }

    ## show optional decorations on plot if desired
    if (show_violin) {
        if (one_facet && (aesth$colour == as.symbol("Feature"))) {
            plot_out <- plot_out +
                geom_violin(aes_string(fill = "Feature"), colour = "gray60",
                            alpha = 0.2, scale = "width")
            plot_out <- .resolve_plot_colours(
                plot_out, object[[as.character(aesth$colour)]],
                as.character(aesth$colour), fill = TRUE)
        }
        else
            plot_out <- plot_out + geom_violin(colour = "gray60", alpha = 0.3,
                                               fill = "gray80", scale = "width")
    }
    if (show_median) {
        plot_out <- plot_out +
            stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                         geom = "crossbar", width = 0.3, alpha = 0.8)
    }
    if (show_smooth) {
        plot_out <- plot_out + stat_smooth(colour = "firebrick", linetype = 2, se = se)
    }

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)

    plot_out
}

################################################################################

#' Plot heatmap of gene expression values
#'
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param features A character vector containing row names of \code{object} 
#' to be shown as rows in the heatmap.
#' @param columns A vector specifying the subset of columns in \code{object}
#' to show as columns in the heatmp. All columns are used by default.
#' @param exprs_values A string indicating which assay of \code{object}
#' should be used as expression values for colouring in the heatmap.
#' @param center A logical scalar indicating whether each row should be 
#' centered at zero prior to plotting. Setting this to \code{TRUE} 
#' @param zlim A numeric vector of length 2, specifying the upper and lower
#' bounds for the expression values. This winsorizes the expression matrix
#' prior to plotting (but after centering, if \code{center=TRUE}. If 
#' \code{NULL}, it defaults to the range of the expression matrix.
#' @param symmetric A logical scalar specifying whether the default \code{zlim}
#' should be symmetric around zero. If \code{TRUE}, the maximum absolute
#' value of \code{zlim} will be computed and multiplied by \code{c(-1, 1)} to 
#' redefine \code{zlim}.
#' @param color A vector of colours specifying the palette to use for mapping 
#' expression values to colours. This defaults to the default setting in 
#' \code{\link[pheatmap]{pheatmap}}.
#' @param colour_columns_by A character vector of column metadata names 
#' in \code{colData(object)} or of gene names, to use to annotate the columns
#' by colour.
#' @param ... Additional arguments to pass to \code{\link[pheatmap]{pheatmap}}.
#'
#' @details Setting \code{center=TRUE} is useful for examining log-fold 
#' changes of each cell's expression profile from the average across all cells.
#' This avoids issues with the entire row appearing a certain colour because
#' the gene is highly/lowly expressed across all cells.
#'
#' Setting \code{zlim} preserves the dynamic range of colours in the presence 
#' of outliers. Otherwise, the plot may be dominated by a few genes,
#' which will \dQuote{flatten} the observed colours for the rest of the heatmap.
#'
#' @return A heatmap is produced on the current graphics device. The output
#' of \code{\link[pheatmap]{pheatmap}} is invisibly returned.
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @author Aaron Lun
#' 
#' @examples
#' example(normalizeSCE) # borrowing the example objects in here.
#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10])

#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     center=TRUE, symmetric=TRUE)
#'
#' plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     colour_columns_by=c("Mutation_Status", "Cell_Cycle"))
#'
#' @export
plotHeatmap <- function(object, features, columns=NULL, exprs_values="logcounts",
                        center=FALSE, zlim=NULL, symmetric=FALSE, color=NULL, 
                        colour_columns_by=NULL, ...) {

    heat.vals <- assay(object, exprs_values)[features,,drop=FALSE]
    if (!is.null(columns)) { 
        heat.vals <- heat.vals[,columns,drop=FALSE]        
    }
    if (center) {
        heat.vals <- heat.vals - rowMeans(heat.vals)
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
            colour_by_out <- .choose_vis_values(object, field, 
                cell_control_default = TRUE, check_features = TRUE,
                exprs_values = exprs_values)

            if (is.numeric(colour_by_out$val)) { 
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
                col_scale <- viridis::viridis(nlevs_colour_by)                    
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
        annotation_col=column_variables, annotation_colors=column_colorings, ...) 
}

################################################################################

#' Plot metadata for cells or features
#'
#' @param object a data.frame (or object that can be coerced to such) object
#' containing metadata in columns to plot.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to plot
#' total_features against log10(total_counts).
#' @param shape numeric scalar to define the plotting shape. Ignored if shape is
#' included in the \code{aesth} argument.
#' @param alpha numeric scalar (in the interval 0 to 1) to define the alpha
#' level (transparency) of plotted points. Ignored if alpha is included in the
#' \code{aesth} argument.
#' @param size numeric scalar to define the plotting size. Ignored if size is
#' included in the \code{aesth} argument.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10)
#'
#' @details Plot cell or feature metadata from an SingleCellExperiment object. If one variable
#'  is supplied then a density plot will be returned. If both variables are
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid
#' will be returned. If both variables are discrete then a jitter plot will be
#' produced. The object returned is a ggplot object, so further layers and
#' plotting options (titles, facets, themes etc) can be added.
#'
#' @return a ggplot plot object
#'
#' @import viridis
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce)
#' plotMetadata(colData(example_sce))
#'
plotMetadata <- function(object,
                         aesth = aes_string(x = "log10(total_counts)",
                                          y = "total_features"),
                         shape = NULL, alpha = NULL, size = NULL,
                         theme_size = 10) {
    object <- as.data.frame(object)
    ## Must have at least an x variable in the aesthetics
    if (is.null(aesth$x))
        stop("No x variable defined. Must have at least an x variable defined
             in the aesth argument.")

    ## Determine type of plot to make
    ### Define plot type if there is only an x variable but no y
    if (is.null(aesth$y)) {
        typeof_x <- typeof(aesth[[1]])
        plot_type <- "density"
        if (typeof_x == "symbol") {
            var_name <- as.character(aesth[[1]])
            x <- typeof(object[, var_name])
            if (is.character(x) | is.factor(x))
                plot_type <- "bar"
        }
    } else {
        ### Define plot type if both x and y variables are provided
        typeof_x <- typeof(aesth$x)
        typeof_y <- typeof(aesth$y)
        var_type_x <- var_type_y <- "continuous"
        if (typeof_x == "symbol") {
            var_name <- as.character(aesth$x)
            x <- object[, var_name]
            if (is.character(x) | is.factor(x))
                var_type_x <- "discrete"
        }
        if (typeof_y == "symbol") {
            var_name <- as.character(aesth$y)
            y <- object[, var_name]
            if (is.character(y) | is.factor(y))
                var_type_y <- "discrete"
        }
        if ( var_type_x == "continuous" && var_type_y == "continuous" )
            plot_type <- "scatter"
        else {
            if ( var_type_x == "discrete" && var_type_y == "discrete" )
                plot_type <- "jitter"
            else
                plot_type <- "violin"
        }
    }

    ## Setup plot
    show_size_guide <- show_alpha_guide <- show_shape_guide <- TRUE
    if ( is.null(aesth$size) ) {
        show_size_guide <- FALSE
        if ( is.null(size) )
            size <- 1
    }
    if ( is.null(aesth$alpha) ) {
        show_alpha_guide <- FALSE
        if ( is.null(alpha) )
            alpha <- 0.7
    }
    if ( is.null(aesth$shape) ) {
        show_shape_guide <- FALSE
        if ( is.null(shape) )
            shape <- 16
    }

    ## Set up basics of plot
    plot_out <- ggplot(object, aesth)

    if (plot_type == "bar") {
        plot_out <- plot_out + geom_bar(stat = "identity")
    } else if (plot_type == "density") {
        plot_out <- plot_out + geom_density(kernel = "rectangular", size = 2) +
            geom_rug(alpha = 0.5, size = 1)
    } else if (plot_type == "violin") {
        plot_out <- plot_out + geom_violin(size = 1, scale = "width")
    } else {
        plot_out <- plot_out + geom_rug(alpha = 0.5, size = 1)
        if (!show_shape_guide && !show_size_guide && !show_alpha_guide) {
            if (plot_type == "scatter") {
                plot_out <- plot_out + geom_point(
                    shape = shape, size = size, alpha = alpha)
            } else {
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    shape = shape, size = size, alpha = alpha)
            }
        } else {
            if (!show_shape_guide && !show_size_guide) {
                if (plot_type == "scatter") {
                    plot_out <- plot_out + geom_point(shape = shape, size = size)
                } else {
                    plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                        shape = shape, size = size)
                }
            } else if (!show_size_guide && !show_alpha_guide) {
                if (plot_type == "scatter") {
                    plot_out <- plot_out + geom_point(
                        alpha = alpha, size = size)
                } else {
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    size = size, alpha = alpha)
                }
            } else if (!show_shape_guide && !show_alpha_guide) {
                if (plot_type == "scatter") {
                    plot_out <- plot_out + geom_point(
                        shape = shape, alpha = alpha)
                } else {
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    shape = shape, alpha = alpha)
                }
            } else {
                if (!show_shape_guide) {
                    if (plot_type == "scatter") {
                        plot_out <- plot_out + geom_point(shape = shape)
                    } else {
                    plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                        shape = shape)
                    }
                } else if (!show_size_guide) {
                    if (plot_type == "scatter") {
                        plot_out <- plot_out + geom_point(size = size)
                    } else {
                        plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                            size = size)
                    }
                } else if (!show_alpha_guide) {
                    if (plot_type == "scatter") {
                        plot_out <- plot_out + geom_point(alpha = alpha)
                    } else {
                        plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                            alpha = alpha)
                    }
                } else {
                    if (plot_type == "scatter") {
                        plot_out <- plot_out + geom_point()
                    } else {
                        plot_out <- plot_out + ggbeeswarm::geom_quasirandom()
                    }
                }
            }
        }
    }
    
    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)

    ## Define plot colours
    if ( "colour" %in% names(aesth) || "color" %in% names(aesth) ) {
        if ( is.null(aesth$colour) )
            colvar <- as.character(aesth$color)
        else
            colvar <- as.character(aesth$colour)
        plot_out <- .resolve_plot_colours(plot_out, object[, colvar],
                                          as.character(colvar))
    }
    if ( !is.null(aesth$fill) ) {
        plot_out <- .resolve_plot_colours(plot_out, object[, aesth$fill],
                                          as.character(aesth$fill), fill = TRUE)
    }

    ## Tweak plot guides
    if ( !show_alpha_guide )
        plot_out <- plot_out + guides(alpha = FALSE)
    if ( !show_shape_guide )
        plot_out <- plot_out + guides(shape = FALSE)
    if ( !show_size_guide )
        plot_out <- plot_out + guides(size = FALSE)

    ## Return plot object
    plot_out
}


################################################################################

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
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce)
#' plotPhenoData(example_sce, aesth = aes_string(x = "log10(total_counts)",
#' y = "total_features", colour = "Mutation_Status"))
#'
#' plotColData(example_sce, aesth = aes_string(x = "log10(total_counts)",
#' y = "total_features", colour = "Mutation_Status"))
#'
#' plotCellData(example_sce, aesth = aes_string(x = "log10(total_counts)",
#' y = "total_features", colour = "Mutation_Status"))

#'
plotPhenoData <- function(object, aesth=aes_string(x = "log10(total_counts)",
                                                   y = "total_features"), ...) {
    ## We must have an SingleCellExperiment object
    if (!is(object, "SingleCellExperiment"))
        stop("object must be an SingleCellExperiment object.")

    ## Define dataframe to pass to plotMetadata
    df_to_plot <- colData(object)

    ## Check that aesthetics make sense for feature names if used
    for (item in unlist(aesth)) {
        item <- as.character(item)
        if ( !(item %in% colnames(colData(object))) &&
             (item %in% rownames(object)) ) {
            df_to_plot <- data.frame(df_to_plot, exprs(object)[item,])
            colnames(df_to_plot)[ncol(df_to_plot)] <- item
        }
    }

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(df_to_plot, aesth, ...)

    plot_out
}

#' @rdname plotPhenoData
#' @export
plotColData <- function(...) {
    plotPhenoData(...)
}

#' @rdname plotPhenoData
#' @export
plotCellData <- function(...) {
    plotPhenoData(...)
}


################################################################################

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
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce)
#' plotFeatureData(example_sce, aesth = aes(x = n_cells_counts, y = log10_total_counts))
#'
#' plotRowData(example_sce, aesth = aes(x = n_cells_counts, y = log10_total_counts))
#'
plotFeatureData <- function(object,
                            aesth = aes_string(x = "n_cells_counts",
                                               y = "log10_total_counts"), ...) {
    ## We must have an SingleCellExperiment object
    if (!is(object, "SingleCellExperiment"))
        stop("object must be an SingleCellExperiment object.")

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(as.data.frame(rowData(object)), aesth, ...)

    plot_out
}

#' @rdname plotFeatureData
#' @export
plotRowData <- function(...) {
    plotFeatureData(...)
}


################################################################################
### Multiplot function for ggplot2 plots

#' Multiple plot function for ggplot2 plots
#'
#' Place multiple \code{\link[ggplot2]{ggplot}} plots on one page.
#'
#' @param ...,plotlist ggplot objects can be passed in ..., or to plotlist (as
#' a list of ggplot objects)
#' @param cols numeric scalar giving the number of columns in the layout
#' @param layout a matrix specifying the layout. If present, \code{cols} is
#' ignored.
#'
#' @details If the layout is something like
#' \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the
#' upper left, 2 will go in the upper right, and 3 will go all the way across
#' the bottom. There is no way to tweak the relative heights or widths of the
#' plots with this simple function. It was adapted from
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'
#' @return a \code{ggplot} plot object
#'
#' @importFrom grid grid.newpage
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @importFrom grid grid.layout
#' @export
#' @examples
#' library(ggplot2)
#' ## This example uses the ChickWeight dataset, which comes with ggplot2
#' ## First plot
#' p1 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet, group = Chick)) +
#'    geom_line() +
#'    ggtitle("Growth curve for individual chicks")
#' ## Second plot
#' p2 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet)) +
#'    geom_point(alpha = .3) +
#'    geom_smooth(alpha = .2, size = 1) +
#'    ggtitle("Fitted growth curve per diet")
#' ## Third plot
#' p3 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, colour = Diet)) +
#'    geom_density() +
#'    ggtitle("Final weight, by diet")
#' ## Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, fill = Diet)) +
#'     geom_histogram(colour = "black", binwidth = 50) +
#'    facet_grid(Diet ~ .) +
#'    ggtitle("Final weight, by diet") +
#'    theme(legend.position = "none")        # No legend (redundant in this graph)
#' ## Combine plots and display
#' multiplot(p1, p2, p3, p4, cols = 2)
#'
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    num_plots <- length(plots)

    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
                         ncol = cols, nrow = ceiling(num_plots / cols))
    }

    if (num_plots == 1) {
        print(plots[[1]])
    } else {
        ## Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(
            layout = grid::grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:num_plots) {
            # Get i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
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



