## Suite of plotting functions

################################################################################
## define colour palettes
.get_palette <- function(palette_name) {
    switch(palette_name,
           tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                         "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                         "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                         "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
           tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                               "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                               "#CDCC5D", "#6DCCDA"),
           colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                            "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                            "#FFBC79", "#CFCFCF"),
           trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
                            "#FFC156", "#69B764", "#F26C64", "#FFDD71",
                            "#9FCD99"),
           purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
                            "#5F5A41", "#B4B19B", "#995688", "#D898BA",
                            "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
           bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
                         "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
                         "#BD0A36", "#F4737A"),
           greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                             "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                             "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
           cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
                      "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
                      "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
                      "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB")
    )
}

# Get nice plotting colour schemes for very general colour variables
.resolve_plot_colours <- function(plot_out, colour_by, colour_by_name,
                                  fill = FALSE) {
    ## if the colour_by object is NULL, return the plot_out object unchanged
    if ( is.null(colour_by) )
        return(plot_out)
    ## Otherwise, set a sensible colour scheme and return the plot_out object
    if ( fill ) {
        if ( is.numeric(colour_by) ) {
            plot_out <- plot_out +
                viridis::scale_fill_viridis(name = colour_by_name)
        } else {
            nlevs_colour_by <- nlevels(as.factor(colour_by))
            if (nlevs_colour_by <= 10) {
                plot_out <- plot_out + scale_fill_manual(
                    values = .get_palette("tableau10medium"),
                    name = colour_by_name)
            } else {
                if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                    plot_out <- plot_out + scale_fill_manual(
                        values = .get_palette("tableau20"),
                        name = colour_by_name)
                } else {
                    plot_out <- plot_out +
                        viridis::scale_fill_viridis(
                            name = colour_by_name, discrete = TRUE)
                }
            }
        }
    } else {
        if ( is.numeric(colour_by) ) {
            plot_out <- plot_out +
                viridis::scale_color_viridis(name = colour_by_name)
        } else {
            nlevs_colour_by <- nlevels(as.factor(colour_by))
            if (nlevs_colour_by <= 10) {
                plot_out <- plot_out + scale_colour_manual(
                    values = .get_palette("tableau10medium"),
                    name = colour_by_name)
            } else {
                if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                    plot_out <- plot_out + scale_colour_manual(
                        values = .get_palette("tableau20"),
                        name = colour_by_name)
                } else {
                    plot_out <- plot_out +
                        viridis::scale_color_viridis(
                            name = colour_by_name, discrete = TRUE)
                }
            }
        }
    }
    plot_out
}

################################################################################
### Generic plot function for SCESet

#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library accounted for by the most highly
#' expressed features for each cell for an \code{SCESet} dataset.
#'
#' @param x an \code{SCESet} object
#' @param block1 character string defining the column of \code{pData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param block2 character string defining the column of \code{pData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param nfeatures numeric scalar indicating the number of features to include
#' in the plot.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"cpm"}
#' (counts per million), \code{"fpkm"} (FPKM values),
#' \code{"counts"} (counts for each feature) or \code{"exprs"} (whatever is in
#' the \code{'exprs'} slot of the \code{SCESet} object; if already on the log2
#' scale, as indicated by the \code{logged} slot of the object, then exprs
#' values are set to the power of 2 (so they are back on the raw scale they were
#'  on) before making the plot).
#' @param linewidth numeric scalar giving the "size" parameter (in ggplot2
#' parlance) for the lines plotted. Default is 1.5.
#' @param y optional argument for generic \code{plot} functions, not used for
#' plotting an \code{SCESet} object
#' @param ... arguments passed to \code{plotSCESet}
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
#' @name plot
#' @aliases plot plot,SCESet-method plot,SCESet,ANY-method
#' @export
#' @method plot
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#'
#' plot(example_sceset, exprs_values = "exprs")
#' plot(example_sceset, exprs_values = "exprs", colour_by = "Cell_Cycle")
#' plot(example_sceset, exprs_values = "exprs", block1 = "Treatment",
#' colour_by = "Cell_Cycle")
#' plot(example_sceset, exprs_values = "exprs", block1 = "Treatment",
#' block2 = "Mutation_Status", colour_by = "Cell_Cycle")
#' # What happens if chosen expression values are not available?
#' plot(example_sceset, block1 = "Treatment", colour_by = "Cell_Cycle")
#'
setMethod("plot", signature("SCESet"),
          function(x, ...) {
              plotSCESet(x, ...)
          })

#' @rdname plot
#' @aliases plot
#' @importFrom dplyr mutate
#' @importFrom plyr aaply
#' @importFrom reshape2 melt
#' @export
plotSCESet <- function(x, block1 = NULL, block2 = NULL, colour_by = NULL,
                        nfeatures = 500, exprs_values = "tpm", ncol = 3,
                       linewidth = 1.5, theme_size = 10) {
    object <- x
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( !is.null(block1) ) {
        if ( !(block1 %in% colnames(pData(object))) )
            stop("The block1 argument must either be NULL or a column of pData(object).")
    }
    if ( !is.null(block2) ) {
        if ( !(block2 %in% colnames(pData(object))) )
            stop("The block2 argument must either be NULL or a column of pData(object).")
    }
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% colnames(pData(object))) )
            stop("The colour_by argument must either be NULL or a column of pData(object).")
    }

    ## Define an expression matrix depending on which values we're using
    exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", 
                                              "cpm", "counts"))
    exprs_mat <- switch(exprs_values,
                        exprs = exprs(object),
                        tpm = tpm(object),
                        cpm = cpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", exprs_values, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        exprs_values <- "exprs"
    }
    if ( exprs_values == "exprs" && object@logged )
        exprs_mat <- 2 ^ (exprs_mat) - object@logExprsOffset

    ## Use plyr to get the sequencing real estate accounted for by features
    nfeatures_total <- nrow(exprs_mat)
    seq_real_estate <- t(plyr::aaply(exprs_mat, 2, .fun = function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }))
    rownames(seq_real_estate) <- 1:nfeatures_total
    nfeatures_to_plot <- nfeatures
    seq_real_estate_long <- reshape2::melt(seq_real_estate[1:nfeatures_to_plot, ],
                                           value.name = "exprs")

    ## Get the proportion of the library accounted for by the top features
    prop_library <- reshape2::melt(t(t(seq_real_estate[1:nfeatures_to_plot, ]) /
                                         colSums(exprs_mat)),
                                   value.name = "prop_library")
    colnames(seq_real_estate_long) <- c("Feature", "Cell", "exprs")
    seq_real_estate_long$Proportion_Library <- prop_library$prop_library

    ## Add block and colour_by information if provided
    if ( !is.null(block1) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block1 = as.factor(rep(object[[block1]],
                                                       each = nfeatures_to_plot)))
    if ( !is.null(block2) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block2 = as.factor(rep(object[[block2]],
                                                       each = nfeatures_to_plot)))
    if ( !is.null(colour_by) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, colour_by = rep(object[[colour_by]],
                                             each = nfeatures_to_plot))
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

#' Plot PCA for an SCESet object
#'
#' Produce a principal components analysis (PCA) plot of two or more principal
#' components for an \code{SCESet} dataset.
#'
#' @param object an \code{SCESet} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the PCA. Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other named element of the \code{assayData} slot of the \code{SCESet}
#' object that can be accessed with the \code{get_exprs} function.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param return_SCESet logical, should the function return an \code{SCESet}
#' object with principal component values for cells in the
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param pca_data_input character argument defining which data should be used
#' as input for the PCA. Possible options are \code{"exprs"} (default), which
#' uses expression data to produce a PCA at the cell level; \code{"pdata"} which
#' uses numeric variables from \code{pData(object)} to do PCA at the cell level;
#' and \code{"fdata"} which uses numeric variables from \code{fData(object)} to
#' do PCA at the feature level.
#' @param selected_variables character vector indicating which variables in
#' \code{pData(object)} to use for the phenotype-data based PCA. Ignored if
#' the argument \code{pca_data_input} is anything other than \code{"pdata"}.
#' @param detect_outliers logical, should outliers be detected in the PC plot?
#' Only an option when \code{pca_data_input} argument is \code{"pdata"}. Default
#' is \code{FALSE}.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is 
#' \code{"auto"}, which hides legends that have only one level and shows others. 
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#'
#' @details The function \code{\link{prcomp}} is used internally to do the PCA.
#' The function checks whether the \code{object} has standardised
#' expression values (by looking at \code{stand_exprs(object)}). If yes, the
#' existing standardised expression values are used for the PCA. If not, then
#' standardised expression values are computed using \code{\link{scale}} (with
#' feature-wise unit variances or not according to the \code{scale_features}
#' argument), added to the object and PCA is done using these new standardised
#' expression values.
#' 
#' If the arguments \code{detect_outliers} and \code{return_SCESet} are both 
#' \code{TRUE}, then the element \code{$outlier} is added to the pData 
#' (phenotype data) slot of the \code{SCESet} object. This element contains
#' indicator values about whether or not each cell has been designated as an
#' outlier based on the PCA. These values can be accessed for filtering 
#' low quality cells with, foe example, \code{example_sceset$outlier}. 
#'
#' @return either a ggplot plot object or an SCESet object
#'
#' @name plotPCA
#' @aliases plotPCA plotPCA,SCESet-method
#' @importFrom BiocGenerics plotPCA
#' @importFrom matrixStats rowVars
#' @importFrom matrixStats colVars
#' @export
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#'
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sceset)
#' plotPCA(example_sceset, colour_by = "Cell_Cycle")
#' plotPCA(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment")
#' plotPCA(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' size_by = "Mutation_Status")
#' plotPCA(example_sceset, shape_by = "Treatment", size_by = "Mutation_Status")
#' plotPCA(example_sceset, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#' 
#' ## experiment with legend
#' example_subset <- example_sceset[, example_sceset$Treatment == "treat1"]
#' plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", legend = "all")
#'
#' plotPCA(example_sceset, shape_by = "Treatment", return_SCESet = TRUE)
#'
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sceset, ncomponents = 8)
#' plotPCA(example_sceset, ncomponents = 4, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
plotPCASCESet <- function(object, ntop=500, ncomponents=2,
                          exprs_values = "exprs", colour_by = NULL,
                          shape_by = NULL, size_by = NULL, feature_set = NULL,
                          return_SCESet = FALSE, scale_features = TRUE,
                          draw_plot = TRUE, pca_data_input = "exprs",
                          selected_variables = NULL, detect_outliers = FALSE,
                          theme_size = 10, legend = "auto") {
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ## Set up indicator for whether to use pData or features for size_by and
    ## colour_by
    colour_by_use_pdata <- TRUE
    size_by_use_pdata <- TRUE
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) &&
             !(colour_by %in% featureNames(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(colour_by %in% varLabels(object)) &&
             (colour_by %in% featureNames(object)) )
            colour_by_use_pdata <- FALSE
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            shape_by <- "is_cell_control"
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) &&
             !(size_by %in% featureNames(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(size_by %in% varLabels(object)) &&
             (size_by %in% featureNames(object)) )
            size_by_use_pdata <- FALSE
    }
    if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
    }

    ## Define an expression matrix depending on which values we're using
#     exprs_values <- match.arg(
#         exprs_values, c("exprs", "tpm", "fpkm", "counts", "cpm", "norm_exprs",
#                         "stand_exprs"))
    exprs_mat <- get_exprs(object, exprs_values)
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", exprs_values, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        exprs_values <- "exprs"
    }

    ## Define features to use: either ntop, or if a set of features is defined,
    ## then those
    if ( is.null(feature_set) ) {
        rv <- matrixStats::rowVars(exprs_mat)
        feature_set <-
            order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    }

    if ( pca_data_input == "pdata" ) {
        #use_variable <- sapply(pData(object), is.double)
        ## select pData features to use
        if ( is.null(selected_variables) ) {
            selected_variables <- c("pct_counts_top_100_features",
                                    "total_features",
                                    "pct_counts_feature_controls",
                                    "n_detected_feature_controls",
                                    "log10_counts_endogenous_features",
                                    "log10_counts_feature_controls")
        }
        use_variable <- varLabels(object) %in% selected_variables
        vars_not_found <- !(selected_variables %in% varLabels(object))
        if ( any(vars_not_found) )
            message(paste("The following selected_variables were not found in pData(object):", selected_variables[vars_not_found]))
        ## scale double variables
        exprs_to_plot <- scale(pData(object)[, use_variable], 
                               scale = scale_features)
    } else if ( pca_data_input == "fdata" ) {
        use_variable <- sapply(fData(object), is.double)
        ## scale double variables
        exprs_to_plot <- scale(fData(object)[, use_variable],
                               scale = scale_features)
    } else {
        # Subsetting to the desired features (do NOT move below 'scale()')
        exprs_to_plot <- exprs_mat[feature_set,,drop=FALSE]
        ## Standardise expression if scale_features argument is TRUE
        exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
    }

    ## Drop any features with zero variance
    keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[, keep_feature]

    ## Compute PCA
    pca <- prcomp(exprs_to_plot)
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    ## Define data.frame for plotting
    if ( pca_data_input == "fdata")
        df_to_plot <- data.frame(pca$x[, 1:ncomponents],
                                 row.names = featureNames(object))
    else
        df_to_plot <- data.frame(pca$x[, 1:ncomponents],
                                 row.names = sampleNames(object))

    if ( !is.null(colour_by) && pca_data_input != "fdata") {
        if ( colour_by_use_pdata )
            colour_by_vals <- pData(object)[[colour_by]]
        else
            colour_by_vals <- exprs(object)[colour_by,]
        if ( nlevels(as.factor(colour_by_vals)) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = colour_by_vals)
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(colour_by_vals))
    }
    if ( !is.null(shape_by) && pca_data_input != "fdata" )
        df_to_plot <- data.frame(
            df_to_plot,
            shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) && pca_data_input != "fdata" ) {
        if ( size_by_use_pdata )
            size_by_vals <- pData(object)[[size_by]]
        else
            size_by_vals <- exprs(object)[size_by,]
        df_to_plot <- data.frame(df_to_plot, size_by = size_by_vals)
    }

    ## conduct outlier detection
    if ( detect_outliers ) {
        if ( requireNamespace("mvoutlier", quietly = TRUE) ) {
            if ( !(pca_data_input == "pdata") ) {
                warning("outlier detection will only be done if pca_data_input
                        argument is 'pdata' (operating on QC metrics)")
            } else {
                outliers <- mvoutlier::pcout(exprs_to_plot, makeplot = FALSE,
                                             explvar = 0.5, crit.M1 = 0.9,
                                             crit.c1 = 5, crit.M2 = 0.9,
                                             crit.c2 = 0.99, cs = 0.25,
                                             outbound = 0.05)
                outlier <- !as.logical(outliers$wfinal01)
                object$outlier <- outlier
                df_to_plot$colour_by <- outlier
                colour_by <- "outlier"
                cat("The following cells/samples are detected as outliers:\n")
                cat(paste(sampleNames(object)[outlier], collapse = "\n"))
                cat("\nVariables with highest loadings for PC1 and PC2:")
                oo <- order(pca$rotation[, 1], decreasing = TRUE)
                if ( requireNamespace("knitr", quietly = TRUE))
                    print(knitr::kable(pca$rotation[oo, 1:2]))
                else
                    print(pca$rotation[oo, 1:2])
            }
        } else {
            warning("The package mvoutlier must be installed to do outlier
                    detection")
        }
    }

    ## Make reduced-dimension plot
    if ( pca_data_input == "fdata" )
        plot_out <- plotReducedDim.default(df_to_plot, ncomponents,
                                           percentVar = percentVar,
                                           legend = legend)
    else
        plot_out <- plotReducedDim.default(df_to_plot, ncomponents, colour_by,
                                           shape_by, size_by, percentVar,
                                           legend = legend)

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")
    
    ## Plot PCA and return appropriate object
    if (return_SCESet) {
        ncomp_out <- max(ncomponents, 10)
        ncomp_out <- min(ncomp_out, ncol(pca$x))
        df_out <- pca$x[, 1:ncomp_out]
        rownames(df_out) <- sampleNames(object)
        attr(df_out, "percentVar") <- percentVar[1:ncomp_out]
        reducedDimension(object) <- df_out
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        ## Return PCA plot
        return(plot_out)
    }
}


#' @rdname plotPCA
#' @param ... further arguments passed to \code{\link{plotPCASCESet}}
#' @aliases plotPCA
#' @export
setMethod("plotPCA", signature("SCESet"),
          function(object, ntop = 500, ncomponents = 2, exprs_values = "exprs",
                   colour_by = NULL, shape_by = NULL, size_by = NULL,
                   feature_set = NULL, return_SCESet = FALSE,
                   scale_features = TRUE, draw_plot = TRUE,
                   pca_data_input = "exprs", selected_variables = NULL,
                   detect_outliers = FALSE, theme_size = 10, legend = "auto") {
              plotPCASCESet(object, ntop, ncomponents, exprs_values, colour_by,
                            shape_by, size_by, feature_set, return_SCESet,
                            scale_features, draw_plot, pca_data_input,
                            selected_variables, detect_outliers, theme_size,
                            legend)
          })




.makePairs <- function(data_matrix) {
    ## with thanks to Gaston Sanchez, who posted this code online
    ## https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
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


################################################################################
### plotTSNE

#' Plot t-SNE for an SCESet object
#'
#' Produce a t-distributed stochastic neighbour embedding (t-SNE) plot of two
#' components for an \code{SCESet} dataset.
#'
#' @param object an \code{SCESet} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the t-SNE Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of t-SNE
#' components to plot, starting from the first t-SNE component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of component 1 vs component
#' 2 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#' top components is produced. NB: computing more than two components for t-SNE
#' can become very time consuming.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values), 
#' or any other named element of the \code{assayData} slot of the \code{SCESet}
#' object that can be accessed with the \code{get_exprs} function.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the t-SNE calculation. If character, entries must all be 
#' in \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param return_SCESet logical, should the function return an \code{SCESet}
#' object with principal component values for cells in the
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param perplexity numeric scalar value defining the "perplexity parameter"
#' for the t-SNE plot. Passed to \code{\link[Rtsne]{Rtsne}} - see documentation
#' for that package for more details.
#' @param legend character, specifying how the legend(s) be shown? Default is 
#' \code{"auto"}, which hides legends that have only one level and shows others. 
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to
#' compute the t-SNE.
#'
#' @return If \code{return_SCESet} is \code{TRUE}, then the function returns an
#' \code{SCESet} object, otherwise it returns a \code{ggplot} object.
#' @name plotTSNE
#' @aliases plotTSNE plotTSNE,SCESet-method
#'
#' @export
#' @seealso
#' \code{\link[Rtsne]{Rtsne}}
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International
#' Conference on Learning Representations, 2013.
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#'
#' ## Examples plotting PC1 and PC2
#' plotTSNE(example_sceset, perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' size_by = "Mutation_Status", perplexity = 10)
#' plotTSNE(example_sceset, shape_by = "Treatment", size_by = "Mutation_Status",
#' perplexity = 5)
#' plotTSNE(example_sceset, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status", perplexity = 5)
#'
#' plotTSNE(example_sceset, shape_by = "Treatment", return_SCESet = TRUE,
#' perplexity = 10)
#'
#'
setMethod("plotTSNE", signature("SCESet"),
          function(object, ntop = 500, ncomponents = 2, exprs_values = "exprs",
                   colour_by = NULL, shape_by = NULL, size_by = NULL,
                   feature_set = NULL, return_SCESet = FALSE,
                   scale_features = TRUE, draw_plot = TRUE, theme_size = 10,
                   rand_seed = NULL, perplexity = floor(ncol(object) / 5), 
                   legend = "auto", ...) {
              ##
              if ( !requireNamespace("Rtsne", quietly = TRUE) )
                  stop("This function requires the 'Rtsne' package.
                       Try: install.packages('Rtsne').")
              ## check legend argument
              legend <- match.arg(legend, c("auto", "none", "all"))
              ## Set up indicator for whether to use pData or features for size_by and
              ## colour_by
              colour_by_use_pdata <- TRUE
              size_by_use_pdata <- TRUE
              ## Check arguments are valid
              if ( !is.null(colour_by) ) {
                  if ( !(colour_by %in% varLabels(object)) &&
                       !(colour_by %in% featureNames(object)) )
                      stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
                  if ( !(colour_by %in% varLabels(object)) &&
                       (colour_by %in% featureNames(object)) )
                      colour_by_use_pdata <- FALSE
              } else {
                  if ( "is_cell_control" %in% varLabels(object) )
                      colour_by <- "is_cell_control"
              }
              if ( !is.null(shape_by) ) {
                  if ( !(shape_by %in% varLabels(object)) )
                      stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
                  if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
                      stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
              } else {
                  if ( "is_cell_control" %in% varLabels(object) )
                      shape_by <- "is_cell_control"
              }
              if ( !is.null(size_by) ) {
                  if ( !(size_by %in% varLabels(object)) &&
                       !(size_by %in% featureNames(object)) )
                      stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
                  if ( !(size_by %in% varLabels(object)) &&
                       (size_by %in% featureNames(object)) )
                      size_by_use_pdata <- FALSE
              }
              if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
                  if ( !(all(feature_set %in% featureNames(object))) )
                      stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
              }

              ## Define an expression matrix depending on which values we're
              ## using
              exprs_mat <- get_exprs(object, exprs_values)
              if ( is.null(exprs_mat) ) {
                  warning(paste0("The object does not contain ", exprs_values, " expression values. Using exprs(object) values instead."))
                  exprs_mat <- exprs(object)
                  exprs_values <- "exprs"
              }

              ## Define features to use: either ntop, or if a set of features is
              ## defined, then those
              if ( is.null(feature_set) ) {
                  rv <- matrixStats::rowVars(exprs_mat)
                  ntop <- min(ntop, length(rv))
                  feature_set <- order(rv, decreasing = TRUE)[seq_len(ntop)]
              }

              ## Drop any features with zero variance
              exprs_to_plot <- exprs_mat[feature_set,]
              keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
              keep_feature[is.na(keep_feature)] <- FALSE
              exprs_to_plot <- exprs_to_plot[keep_feature, ]
              
              ## Standardise expression if stand_exprs(object) is null
              exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))

              ## Compute t-SNE
              if ( !is.null(rand_seed) )
                  set.seed(rand_seed)
              tsne_out <- Rtsne::Rtsne(t(exprs_to_plot),
                                       initial_dims = max(50, ncol(object)),
                                       perplexity = perplexity, ...)


              ## Define data.frame for plotting
              df_to_plot <- data.frame(tsne_out$Y[, 1:ncomponents],
                                       row.names = sampleNames(object))
              if ( !is.null(colour_by) ) {
                  if ( colour_by_use_pdata )
                      colour_by_vals <- pData(object)[[colour_by]]
                  else
                      colour_by_vals <- exprs(object)[colour_by,]
                  if ( nlevels(as.factor(colour_by_vals)) > 10 )
                      df_to_plot <- data.frame(
                          df_to_plot, colour_by = colour_by_vals)
                  else
                      df_to_plot <- data.frame(
                          df_to_plot, colour_by = as.factor(colour_by_vals))
              }
              if ( !is.null(shape_by) )
                  df_to_plot <- data.frame(
                      df_to_plot,
                      shape_by = as.factor(pData(object)[[shape_by]]))
              if ( !is.null(size_by) ) {
                  if ( size_by_use_pdata )
                      size_by_vals <- pData(object)[[size_by]]
                  else
                      size_by_vals <- exprs(object)[size_by,]
                  df_to_plot <- data.frame(df_to_plot, size_by = size_by_vals)
              }

              ## Make reduced-dimension plot
              plot_out <- plotReducedDim.default(df_to_plot, ncomponents,
                                                 colour_by, shape_by, size_by,
                                                 legend = legend)

              ## Define plotting theme
              if ( requireNamespace("cowplot", quietly = TRUE) )
                  plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
              else
                  plot_out <- plot_out + theme_bw(theme_size)
              ## remove legend if so desired
              if ( legend == "none" )
                  plot_out <- plot_out + theme(legend.position = "none")
              
              ## Plot t-SNE and return appropriate object
              if (return_SCESet) {
                  df_out <- tsne_out$Y[, 1:ncomponents]
                  rownames(df_out) <- sampleNames(object)
                  reducedDimension(object) <- df_out
                  if ( draw_plot )
                      print(plot_out)
                  return(object)
              } else {
                  ## Return t-SNE plot
                  return(plot_out)
              }
          })

################################################################################
### plotDiffusionMap

#' Plot a diffusion map for an SCESet object
#'
#' Produce a diffusion map plot of two components for an \code{SCESet} dataset.
#'
#' @param object an \code{SCESet} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the diffusion map. Default is \code{500}, but any \code{ntop} 
#' argument is overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first diffusion map component. Default 
#' is 2. If \code{ncomponents} is 2, then a scatterplot of component 1 vs 
#' component 2 is produced. If \code{ncomponents} is greater than 2, a pairs 
#' plots for the top components is produced. NB: computing many components for 
#' the diffusion map can become time consuming.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other named element of the \code{assayData} slot of the \code{SCESet}
#' object that can be accessed with the \code{get_exprs} function.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the diffusion map. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param return_SCESet logical, should the function return an \code{SCESet}
#' object with principal component values for cells in the
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param sigma argument passed to \code{\link[destiny]{DiffusionMap}}
#' @param distance argument passed to \code{\link[destiny]{DiffusionMap}}
#' @param legend character, specifying how the legend(s) be shown? Default is 
#' \code{"auto"}, which hides legends that have only one level and shows others. 
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param ... further arguments passed to \code{\link[destiny]{DiffusionMap}}
#'
#' @details The function \code{\link[destiny]{DiffusionMap}} is used internally 
#' to compute the diffusion map.
#'
#' @return If \code{return_SCESet} is \code{TRUE}, then the function returns an
#' \code{SCESet} object, otherwise it returns a \code{ggplot} object.
#' @name plotDiffusionMap
#' @aliases plotDiffusionMap plotDiffusionMap,SCESet-method
#'
#' @export
#' @seealso
#' \code{\link[destiny]{destiny}}
#' @references
#' Haghverdi L, Buettner F, Theis FJ. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics. 2015; doi:10.1093/bioinformatics/btv325
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#'
#' ## Examples plotting diffusion maps
#' plotDiffusionMap(example_sceset)
#' plotDiffusionMap(example_sceset, colour_by = "Cell_Cycle")
#' plotDiffusionMap(example_sceset, colour_by = "Cell_Cycle",
#' shape_by = "Treatment")
#' plotDiffusionMap(example_sceset, colour_by = "Cell_Cycle",
#' shape_by = "Treatment", size_by = "Mutation_Status")
#' plotDiffusionMap(example_sceset, shape_by = "Treatment",
#' size_by = "Mutation_Status")
#' plotDiffusionMap(example_sceset, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
#' plotDiffusionMap(example_sceset, shape_by = "Treatment",
#' return_SCESet = TRUE)
#'
#'
plotDiffusionMapSCESet <- function(object, ntop = 500, ncomponents = 2,
                                   exprs_values = "exprs", colour_by = NULL,
                                   shape_by = NULL, size_by = NULL,
                                   feature_set = NULL, return_SCESet = FALSE,
                                   scale_features = TRUE, draw_plot = TRUE,
                                   theme_size = 10, rand_seed = NULL,
                                   sigma = NULL, distance = "euclidean",
                                   legend = "auto", ...) {
    ##
    if ( !requireNamespace("destiny", quietly = TRUE) )
        stop("This function requires the 'destiny' package.
                       Try from Bioconductor with:
                       source('https://bioconductor.org/biocLite.R')
                       biocLite('destiny').")
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ## Set up indicator for whether to use pData or features
    ## for size_by and
    ## colour_by
    colour_by_use_pdata <- TRUE
    size_by_use_pdata <- TRUE
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) &&
             !(colour_by %in% featureNames(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(colour_by %in% varLabels(object)) &&
             (colour_by %in% featureNames(object)) )
            colour_by_use_pdata <- FALSE
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            shape_by <- "is_cell_control"
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) &&
             !(size_by %in% featureNames(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(size_by %in% varLabels(object)) &&
             (size_by %in% featureNames(object)) )
            size_by_use_pdata <- FALSE
    }
    if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
    }

    ## Define an expression matrix depending on which values we're
    ## using
    exprs_mat <- get_exprs(object, exprs_values)
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", exprs_values,
                       " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        exprs_values <- "exprs"
    }

    ## Define features to use: either ntop, or if a set of features is
    ## defined, then those
    if ( is.null(feature_set) ) {
        rv <- matrixStats::rowVars(exprs_mat)
        feature_set <-
            order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    }

    ## Drop any features with zero variance
    exprs_to_plot <- exprs_mat
    exprs_to_plot <- exprs_to_plot[feature_set,]
    keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[keep_feature, ]
    
    ## Standardise expression if indicated by scale_features argument
    exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))

    ## Compute DiffusionMap
    if ( !is.null(rand_seed) )
        set.seed(rand_seed)
    difmap_out <- destiny::DiffusionMap(
        t(exprs_to_plot), sigma = sigma, distance = distance, ...)


    ## Define data.frame for plotting
    df_to_plot <- data.frame(difmap_out@eigenvectors[, 1:ncomponents],
                             row.names = sampleNames(object))
    if ( !is.null(colour_by) ) {
        if ( colour_by_use_pdata )
            colour_by_vals <- pData(object)[[colour_by]]
        else
            colour_by_vals <- exprs(object)[colour_by,]
        if ( nlevels(as.factor(colour_by_vals)) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = colour_by_vals)
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(colour_by_vals))
    }
    if ( !is.null(shape_by) )
        df_to_plot <- data.frame(
            df_to_plot,
            shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) ) {
        if ( size_by_use_pdata )
            size_by_vals <- pData(object)[[size_by]]
        else
            size_by_vals <- exprs(object)[size_by,]
        df_to_plot <- data.frame(df_to_plot, size_by = size_by_vals)
    }

    ## Make reduced-dimension plot
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents,
                                       colour_by, shape_by, size_by, 
                                       legend = legend)

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")
    
    ## Plot PCA and return appropriate object
    if (return_SCESet) {
        df_out <- difmap_out@eigenvectors[, 1:ncomponents]
        rownames(df_out) <- sampleNames(object)
        reducedDimension(object) <- df_out
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        ## Return PCA plot
        return(plot_out)
    }
}

#' @rdname plotDiffusionMap
#' @aliases plotDiffusionMap
#' @export
setMethod("plotDiffusionMap", signature("SCESet"),
          function(object, ntop = 500, ncomponents = 2, exprs_values = "exprs",
                   colour_by = NULL, shape_by = NULL, size_by = NULL,
                   feature_set = NULL, return_SCESet = FALSE,
                   scale_features = FALSE, draw_plot = TRUE, theme_size = 10,
                   rand_seed = NULL, sigma = NULL, distance = "euclidean",
                   legend = "auto", ...) {
              plotDiffusionMapSCESet(object, ntop, ncomponents, exprs_values,
                                     colour_by, shape_by, size_by,
                                     feature_set, return_SCESet,
                                     scale_features, draw_plot, theme_size,
                                     rand_seed, sigma, distance, legend, ...)
          })


################################################################################
### plotMDS

#' Produce a multidimensional scaling plot for an SCESet object
#'
#' #' Produce an MDS plot from the cell pairwise distance data in an
#' \code{SCESet} dataset.
#'
#' @param object an \code{SCESet} object
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced. NB: computing more than two components for t-SNE can become very
#' time consuming.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param return_SCESet logical, should the function return an \code{SCESet}
#' object with principal component values for cells in the
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned.
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is 
#' \code{"auto"}, which hides legends that have only one level and shows others. 
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param ... arguments passed to S4 plotMDS method
#'
#' @details The function \code{\link{cmdscale}} is used internally to
#' compute the multidimensional scaling components to plot.
#'
#' @return If \code{return_SCESet} is \code{TRUE}, then the function returns an
#' \code{SCESet} object, otherwise it returns a \code{ggplot} object.
#' @name plotMDS
#' @aliases plotMDS plotMDS,SCESet-method
#'
#' @export
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#' example_sceset <- calculateQCMetrics(example_sceset)
#'
#' ## define cell-cell distances
#' cellDist(example_sceset) <- as.matrix(dist(t(exprs(example_sceset))))
#'
#' ## Examples plotting
#' plotMDS(example_sceset)
#' plotMDS(example_sceset, colour_by = "Cell_Cycle")
#' plotMDS(example_sceset, colour_by = "Cell_Cycle",
#' shape_by = "Treatment")
#'
#' ## define cell-cell distances differently
#' cellDist(example_sceset) <- as.matrix(dist(t(counts(example_sceset)),
#' method = "canberra"))
#' plotMDS(example_sceset, colour_by = "Cell_Cycle",
#' shape_by = "Treatment", size_by = "Mutation_Status")
#'
plotMDSSCESet <- function(object, ncomponents = 2, colour_by = NULL,
                          shape_by = NULL, size_by = NULL,
                          return_SCESet = FALSE, draw_plot = TRUE,
                          theme_size = 10, legend = "auto") {
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ##
    cell_dist <- cellDist(object)
    ncells <- ncol(object)
    if ( !(dim(cellDist(object))[1] == ncells &&
           dim(cellDist(object))[2] == ncells) )
        stop("cellDist(object) is not of the correct dimensions. Please define cell pairwise distances and try again:
             e.g. cellDist(object) <- as.matrix(dist(t(exprs(object))))")
    ## Set up indicator for whether to use pData or features
    ## for size_by and
    ## colour_by
    colour_by_use_pdata <- TRUE
    size_by_use_pdata <- TRUE
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) &&
             !(colour_by %in% featureNames(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(colour_by %in% varLabels(object)) &&
             (colour_by %in% featureNames(object)) )
            colour_by_use_pdata <- FALSE
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            shape_by <- "is_cell_control"
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) &&
             !(size_by %in% featureNames(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(size_by %in% varLabels(object)) &&
             (size_by %in% featureNames(object)) )
            size_by_use_pdata <- FALSE
    }

    ## Compute multidimentional scaling
    mds_out <- cmdscale(cell_dist, k = ncomponents)

    ## Define data.frame for plotting
    df_to_plot <- data.frame(mds_out[, 1:ncomponents],
                             row.names = sampleNames(object))
    colnames(df_to_plot) <- paste0("Component_", 1:ncomponents)
    if ( !is.null(colour_by) ) {
        if ( colour_by_use_pdata )
            colour_by_vals <- pData(object)[[colour_by]]
        else
            colour_by_vals <- exprs(object)[colour_by,]
        if ( nlevels(as.factor(colour_by_vals)) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = colour_by_vals)
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(colour_by_vals))
    }
    if ( !is.null(shape_by) )
        df_to_plot <- data.frame(
            df_to_plot,
            shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) ) {
        if ( size_by_use_pdata )
            size_by_vals <- pData(object)[[size_by]]
        else
            size_by_vals <- exprs(object)[size_by,]
        df_to_plot <- data.frame(df_to_plot, size_by = size_by_vals)
    }

    ## Make reduced-dimension plot
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents,
                                       colour_by, shape_by, size_by, 
                                       legend = legend)

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")
    
    ## Plot PCA and return appropriate object
    if (return_SCESet) {
        df_out <- mds_out[, 1:ncomponents]
        rownames(df_out) <- sampleNames(object)
        reducedDimension(object) <- df_out
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        ## Return PCA plot
        return(plot_out)
    }
}

#' @rdname plotMDS
#' @aliases plotMDS
#' @export
setMethod("plotMDS", signature("SCESet"),
          function(object, ncomponents = 2, colour_by = NULL, shape_by = NULL,
                   size_by = NULL, return_SCESet = FALSE, draw_plot = TRUE,
                   theme_size = 10, legend = "auto") {
              plotMDSSCESet(object, ncomponents, colour_by, shape_by, size_by,
                            return_SCESet, draw_plot, theme_size, legend)
          })



################################################################################
### plotReducedDim

#' Plot reduced dimension representation of cells
#'
#' @param object an \code{SCESet} object with a non-NULL \code{reducedDimension}
#' slot.
#' @param df_to_plot data.frame containing a reduced dimension represenation of
#' cells and optional metadata for the plot.
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of Dimension 2 vs Dimension
#' 1 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#'  top dimensions is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param percentVar numeric vector giving the proportion of variance in
#' expression explained by each reduced dimension. Only expected to be used
#' internally in the \code{\link[scater]{plotPCA}} function.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is 
#' \code{"auto"}, which hides legends that have only one level and shows others. 
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param ... optional arguments (from those listed above) passed to
#' \code{plotReducedDim.SCESet} or \code{plotReducedDim.default}
#'
#' @details The function \code{plotReducedDim.default} assumes that the first
#' \code{ncomponents} columns of \code{df_to_plot} contain the reduced dimension
#'  components to plot, and that any subsequent columns define factors for
#'  \code{colour_by}, \code{shape_by} and \code{size_by} in the plot.
#'
#' @return a ggplot plot object
#'
#' @name plotReducedDim
#' @aliases plotReducedDim plotReducedDim,SCESet-method plotReducedDim,data.frame-method
#' @import viridis
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#'
#' reducedDimension(example_sceset) <- prcomp(t(exprs(example_sceset)), scale. = TRUE)$x
#' plotReducedDim(example_sceset)
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle")
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle", shape_by="Treatment")
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle", size_by="Treatment")
#' plotReducedDim(example_sceset, ncomponents=5)
#' plotReducedDim(example_sceset, ncomponents=5, colour_by="Cell_Cycle", shape_by="Treatment")
#'
#'
plotReducedDim.default <- function(df_to_plot, ncomponents=2, colour_by=NULL,
                           shape_by=NULL, size_by=NULL, percentVar=NULL,
                           theme_size = 10, legend = "auto") {
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"), several.ok = FALSE)
    ## Define plot
    if ( ncomponents > 2 ) {
        ## expanding numeric columns for pairs plot
        df_to_expand <- df_to_plot[, 1:ncomponents]
        if ( is.null(percentVar) ) {
            colnames(df_to_expand) <- colnames(df_to_plot)[1:ncomponents]
        } else {
            colnames(df_to_expand) <- paste0(
                colnames(df_to_plot)[1:ncomponents], ": ",
                round(percentVar[1:ncomponents] * 100), "% variance")
        }
        gg1 <- .makePairs(df_to_expand)
        ## new data frame
        df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -c(1:ncomponents)])
        colnames(df_to_plot_big)[-c(1:4)] <- colnames(df_to_plot)
        ## pairs plot
        plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) +
            facet_grid(xvar ~ yvar, scales = "free") +
            stat_density(aes_string(x = "x",
                                    y = "(..scaled.. * diff(range(x)) + min(x))"),
                         data = gg1$densities, position = "identity",
                         colour = "grey20", geom = "line") +
            xlab("") +
            ylab("") +
            theme_bw(theme_size)
    } else {
        comps <- colnames(df_to_plot)[1:2]
        if ( is.null(percentVar) ) {
            x_lab <- "Dimension 1"
            y_lab <- "Dimension 2"
        } else {
            x_lab <- paste0("Component 1: ", round(percentVar[1] * 100),
                            "% variance")
            y_lab <- paste0("Component 2: ", round(percentVar[2] * 100),
                            "% variance")
        }
        plot_out <- ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
            xlab(x_lab) +
            ylab(y_lab) +
            geom_rug(colour = "gray20", alpha = 0.65) +
            theme_bw(theme_size)
    }

    ## if only one level for the variable, set to NULL
    if ( legend == "auto" ) {
        if ( !is.null(colour_by) && length(unique(df_to_plot$colour_by)) == 1)
            colour_by <- NULL
        if ( !is.null(shape_by) && length(unique(df_to_plot$shape_by)) == 1)
            shape_by  <- NULL
        if ( !is.null(size_by) && length(unique(df_to_plot$size_by)) == 1)
            size_by <- NULL
    }
    ## Apply colour_by, shape_by and size_by variables if defined
    if ( !is.null(colour_by) && !is.null(shape_by) && !is.null(size_by) ) {
        plot_out <- plot_out +
            geom_point(aes_string(colour = "colour_by", shape = "shape_by",
                                  size = "size_by"), alpha = 0.65) +
            guides(size = guide_legend(title = size_by),
                   shape = guide_legend(title = shape_by))
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$colour_by,
                                          colour_by)
    } else {
        if  ( sum(is.null(colour_by) + is.null(shape_by) + is.null(size_by)) == 1 ) {
            if ( !is.null(colour_by) && !is.null(shape_by) ) {
                plot_out <- plot_out +
                    geom_point(aes_string(colour = "colour_by",
                                          shape = "shape_by"),
                               alpha = 0.65) +
                    guides(shape = guide_legend(title = shape_by))
                plot_out <- .resolve_plot_colours(plot_out,
                                                  df_to_plot$colour_by,
                                                  colour_by)
            }
            if ( !is.null(colour_by) && !is.null(size_by) ) {
                plot_out <- plot_out +
                    geom_point(aes_string(fill = "colour_by", size = "size_by"),
                               shape = 21, colour = "gray70", alpha = 0.65) +
                    guides(size = guide_legend(title = size_by))
                plot_out <- .resolve_plot_colours(plot_out,
                                                  df_to_plot$colour_by,
                                                  colour_by, fill = TRUE)
            }
            if ( !is.null(shape_by) && !is.null(size_by) ) {
                plot_out <- plot_out +
                    geom_point(aes_string(shape = "shape_by", size = "size_by"),
                               fill = "gray20", colour = "gray20",
                               alpha = 0.65) +
                    guides(size = guide_legend(title = size_by),
                           shape = guide_legend(title = shape_by))
            }
        } else {
            if ( sum(is.null(colour_by) + is.null(shape_by) +
                     is.null(size_by)) == 2 ) {
                if ( !is.null(colour_by) ) {
                    plot_out <- plot_out +
                        geom_point(aes_string(fill = "colour_by"),
                                   shape = 21, colour = "gray70", alpha = 0.65)
                    plot_out <- .resolve_plot_colours(plot_out,
                                                      df_to_plot$colour_by,
                                                      colour_by, fill = TRUE)
                }
                if ( !is.null(shape_by) ) {
                    plot_out <- plot_out +
                        geom_point(aes_string(shape = "shape_by"),
                                   colour = "gray20", alpha = 0.65) +
                        guides(shape = guide_legend(title = shape_by))
                }
                if ( !is.null(size_by) ) {
                    plot_out <- plot_out +
                        geom_point(aes_string(size = "size_by"),
                                   fill = "gray20", shape = 21,
                                   colour = "gray70", alpha = 0.65) +
                        guides(size = guide_legend(title = size_by))
                }
            } else {
                plot_out <- plot_out +
                    geom_point(fill = "gray20", shape = 21,
                               colour = "gray70", alpha = 0.65)
            }
        }
    }

    ## Return plot
    plot_out
}

#' @rdname plotReducedDim
#' @aliases plotReducedDim
#' @export
plotReducedDim.SCESet <- function(object, ncomponents=2, colour_by=NULL,
                                  shape_by=NULL, size_by=NULL, theme_size = 10,
                                  legend = "auto") {
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ## Set up indicator for whether to use pData or features for size_by and
    ## colour_by
    colour_by_use_pdata <- TRUE
    size_by_use_pdata <- TRUE
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) &&
             !(colour_by %in% featureNames(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(colour_by %in% varLabels(object)) &&
             (colour_by %in% featureNames(object)) )
            colour_by_use_pdata <- FALSE
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            shape_by <- "is_cell_control"
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) &&
             !(size_by %in% featureNames(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(size_by %in% varLabels(object)) &&
             (size_by %in% featureNames(object)) )
            size_by_use_pdata <- FALSE
    }

    ## Extract reduced dimension representation of cells
    if ( is.null(reducedDimension(object)) )
        stop("reducedDimension slot of object is NULL. Need non null reducedDimension to plot.")
    red_dim <- redDim(object)
    if ( ncomponents > ncol(red_dim) )
        stop("ncomponents to plot is larger than number of columns of reducedDimension(object)")

    ## Define data.frame for plotting
    df_to_plot <- data.frame(red_dim[, 1:ncomponents])
    if ( !is.null(colour_by) ) {
        if ( colour_by_use_pdata )
            colour_by_vals <- pData(object)[[colour_by]]
        else
            colour_by_vals <- exprs(object)[colour_by,]
        if ( nlevels(as.factor(colour_by_vals)) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = colour_by_vals)
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(colour_by_vals))
    }
    if ( !is.null(shape_by) )
        df_to_plot <- data.frame(df_to_plot,
                                 shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) ) {
        if ( size_by_use_pdata )
            size_by_vals <- pData(object)[[size_by]]
        else
            size_by_vals <- exprs(object)[size_by,]
        df_to_plot <- data.frame(df_to_plot, size_by = size_by_vals)
    }

    ## Call default method to make the plot
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents, colour_by,
                                       shape_by, size_by, percentVar = NULL,
                                       legend = legend)

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")

    ## return plot
    plot_out
    
}

#' @rdname plotReducedDim
#' @aliases plotReducedDIm
#' @export
setMethod("plotReducedDim", signature("SCESet"),
          function(object, ncomponents=2, colour_by=NULL, shape_by=NULL,
                   size_by=NULL, theme_size = 10, legend="auto") {
              plotReducedDim.SCESet(object, ncomponents, colour_by, shape_by,
                                    size_by, theme_size, legend)
          })

#' @rdname plotReducedDim
#' @aliases plotReducedDim
#' @export
setMethod("plotReducedDim", signature("data.frame"),
          function(object, ncomponents=2, colour_by=NULL, shape_by=NULL,
                   size_by=NULL, percentVar=NULL, legend="auto") {
              plotReducedDim.default(object, ncomponents, colour_by, shape_by,
                                     size_by, percentVar, legend)
          })




################################################################################
### Plot cells in plate positions

#' Plot cells in plate positions
#' 
#' Plots cells in their position on a plate, coloured by phenotype data or 
#' feature expression.
#' 
#' @param object an \code{SCESet} object. If \code{object$plate_position} is not
#' \code{NULL}, then this will be used to define each cell's position on the 
#' plate, unless the \code{plate_position} argument is specified.
#' @param plate_position optional character vector providing a position on the
#' plate for each cell (e.g. A01, B12, etc, where letter indicates row and 
#' number indicates column). Specifying this argument overrides any plate 
#' position information extracted from the SCESet object.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param x_position numeric vector providing x-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
#' @param y_position numeric vector providing y-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
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
#' @export
#' 
#' @examples 
#' ## prepare data
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#'
#' ## define plate positions
#' example_sceset$plate_position <- paste0(
#' rep(LETTERS[1:5], each = 8), rep(formatC(1:8, width = 2, flag = "0"), 5))
#'
#' ## plot plate positions
#' plotPlatePosition(example_sceset, colour_by = "Mutation_Status")
#' 
#' plotPlatePosition(example_sceset, colour_by = "Gene_0004")
#' 
plotPlatePosition <- function(object, plate_position = NULL, 
                              colour_by = NULL,
                              x_position = NULL, y_position = NULL,
                              theme_size = 24, legend = "auto") {
    ## check object is SCESet object
    if ( !is(object, "SCESet") )
        stop("Object must be of class SCESet.")
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    ## Set up indicator for whether to use pData or features
    ## for colour_by
    colour_by_use_pdata <- TRUE
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) &&
             !(colour_by %in% featureNames(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)]")
        if ( !(colour_by %in% varLabels(object)) &&
             (colour_by %in% featureNames(object)) )
            colour_by_use_pdata <- FALSE
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    
    ## obtain well positions
    if ( !is.null(plate_position) ) {
        if ( length(plate_position) != ncol(object) )
            stop("Supplied plate_position argument must have same length as number of columns of SCESet object.")
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
        
    if ( !is.null(colour_by) ) {
        if ( colour_by_use_pdata )
            colour_by_vals <- pData(object)[[colour_by]]
        else
            colour_by_vals <- exprs(object)[colour_by,]
        df_to_plot <- data.frame(df_to_plot, colour_by = colour_by_vals)
    }
    
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
#' @param object an SCESet object containing expression values and
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
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable by which to colour
#' expression values on the plot.
#' @param shape_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot.
#' @param size_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot.
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
#' \code{plotExpressionSCESet} or \code{plotExpressionDefault}
#'
#' @details Plot expression values (default log2(transcripts-per-million +
#' 1), if available) for a set of features.
#' @return a ggplot plot object
#'
#' @name plotExpression
#' @aliases plotExpression plotExpression,SCESet-method plotExpression,data.frame-method
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#'
#' ## default plot
#' plotExpression(example_sceset, 1:15)
#' plotExpression(example_sceset, 1:15, jitter = "jitter")
#' 
#' ## plot expression against an x-axis value
#' plotExpression(example_sceset, 1:6, "Mutation_Status")
#'
#' ## explore options
#' plotExpression(example_sceset, 1:6, x="Mutation_Status", exprs_values="exprs",
#' colour_by="Cell_Cycle", show_violin=TRUE, show_median=TRUE)
#' plotExpression(example_sceset, 1:6, x="Mutation_Status", exprs_values="counts",
#' colour_by="Cell_Cycle", show_violin=TRUE, show_median=TRUE)
#'
#' ## plot expression against expression values for Gene_0004
#' plotExpression(example_sceset, 1:4, "Gene_0004")
#' plotExpression(example_sceset, 1:4, "Gene_0004", show_smooth = TRUE)
#' plotExpression(example_sceset, 1:4, "Gene_0004", show_smooth = TRUE, se = FALSE)
#'
plotExpressionSCESet <- function(object, features, x = NULL, exprs_values = "exprs",
                                 colour_by = NULL, shape_by = NULL,
                                 size_by = NULL, ncol = 2, xlab = NULL,
                                 show_median = FALSE, show_violin = TRUE,
                                 show_smooth = FALSE, alpha = 0.6, 
                                 theme_size = 10, log2_values = FALSE, size = NULL, 
                                 scales = "fixed", se = TRUE, jitter = "swarm") {
    ## Check object is an SCESet object
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet")

    ## Define number of features to plot
    if (is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)

    ## Check arguments are valid
    if ( is.null(x) ) {
        x_is_feature <- FALSE
    } else {
        if ( !(x %in% varLabels(object)) && !(x %in% featureNames(object)))
            stop("the argument 'x' should specify a column of pData(object) [see varLabels(object)] or a feature [see featureNames(object)")
        if ( x %in% featureNames(object) ) {
            x_is_feature <- TRUE
            show_violin <- FALSE
            show_median <- FALSE
        } else
            x_is_feature <- FALSE
    }
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if ( typeof(features) == "character" ) {
        if ( !(all(features %in% featureNames(object))) )
            stop("when the argument 'features' is of type character, all features must be in featureNames(object)")
    }

    exprs_mat <- get_exprs(object, exprs_values)
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else
        ylab <- paste0("Expression (", exprs_values, ")")
    to_melt <- as.matrix(exprs_mat[features, , drop = FALSE])

    ## Melt the expression data and metadata into a convenient form
    evals_long <- reshape2::melt(to_melt, value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")
    ## Extend the samples information
    samps <-  pData(object)
    if ( x_is_feature )
        samps[[x]] <- exprs_mat[x, ]
    samples_long <- samps[rep(seq_len(ncol(object)), each = nfeatures),]

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
    if ( !is.null(colour_by) )
        aesth$colour <- as.symbol(colour_by)
    else {
        if ( !is.null(is_exprs(object)) ) {
            ## Colour by is_exprs if we can (i.e. is_exprs is not NULL)
            isexpr_long <- reshape2::melt(is_exprs(object)[features,],
                                          value.name = "is_exprs")
            evals_long <- dplyr::mutate(
                evals_long, Is_Expressed = as.vector(isexpr_long$is_exprs))
            aesth$colour <- as.symbol("Is_Expressed")
        }
    }
    if ( !is.null(shape_by) )
        aesth$shape <- as.symbol(shape_by)
    ## Define sensible x-axis label if NULL
    if ( is.null(xlab) )
        xlab <- x

    ## Combine the expression values and sample information
    object <- cbind(evals_long, samples_long)

    ## Make the plot
    plot_out <- plotExpressionDefault(object, aesth, ncol, xlab, ylab,
                                      show_median, show_violin, show_smooth,
                                      alpha, size, scales, one_facet, se, jitter)

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
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
plotExpressionDefault <- function(object, aesth, ncol=2, xlab = NULL,
                                  ylab = NULL, show_median = FALSE,
                                  show_violin = TRUE, show_smooth = FALSE, 
                                  alpha = 0.6, size = NULL, scales = "fixed", 
                                  one_facet = FALSE, se = TRUE, jitter = "swarm") {
    if ( !("Feature" %in% names(object)) )
        stop("object needs a column named 'Feature' to define the feature(s) by which to plot expression.")

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
                    alpha = alpha, size = size)
            else
                plot_out <- plot_out + geom_jitter(
                    alpha = alpha, size = size, position = position_jitter(height = 0))
        }
        else {
            if (jitter == "swarm")
                plot_out <- plot_out + ggbeeswarm::geom_quasirandom(
                    alpha = alpha)
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
    plot_out
}


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature(object = "SCESet"),
          function(object, ...) {
              plotExpressionSCESet(object, ...)
          })


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature("data.frame"),
          function(object, ...) {
              plotExpressionDefault(object, ...)
          })

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
#' @details Plot cell or feature metadata from an SCESet object. If one variable
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotMetadata(pData(example_sceset))
#'
plotMetadata <- function(object,
                         aesth=aes_string(x = "log10(total_counts)",
                                          y = "total_features"),
                         shape = NULL, alpha = NULL, size = NULL,
                         theme_size = 10) {
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
                var_type_x <- "discrete"
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
            aesth$size <- 4
        else
            aesth$size <- size
    }
    if ( is.null(aesth$alpha) ) {
        show_alpha_guide <- FALSE
        if ( is.null(alpha) )
            aesth$alpha <- 0.7
        else
            aesth$alpha <- alpha
    }
    if ( is.null(aesth$shape) ) {
        show_shape_guide <- FALSE
        if ( is.null(shape) )
            shape <- 16
    }

    ## Set up basics of plot
    plot_out <- ggplot(object, aesth)

    ## Density plot
    if (plot_type == "bar") {
        plot_out <- plot_out + geom_bar(stat = "identity")
    }
    if (plot_type == "density") {
        plot_out <- plot_out + geom_density(kernel = "rectangular", size = 2) +
            geom_rug(alpha = 0.5, size = 1)
    }
    if (plot_type == "jitter") {
        if ( !show_shape_guide )
            plot_out <- plot_out +  ggbeeswarm::geom_quasirandom(shape = shape)
        else
            plot_out <- plot_out +  ggbeeswarm::geom_quasirandom()
    }
    if (plot_type == "scatter") {
        if ( !show_shape_guide )
            plot_out <- plot_out +  ggbeeswarm::geom_quasirandom(shape = shape)
        else
            plot_out <- plot_out +  ggbeeswarm::geom_quasirandom()
        plot_out <- plot_out + geom_rug(alpha = 0.5, size = 1)
    }
    if (plot_type == "violin") {
        if ( !show_shape_guide )
            plot_out  <- plot_out +  ggbeeswarm::geom_quasirandom(shape = shape)
        else
            plot_out  <- plot_out +  ggbeeswarm::geom_quasirandom()
        plot_out <- plot_out + geom_violin(size = 1, scale = "width")

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

#' Plot phenotype data from an SCESet object
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to plot
#' total_features against log10(total_counts).
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param ... arguments passed to \code{\link{plotMetadata}}.
#'
#' @details Plot phenotype data from an SCESet object. If one variable is
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotPhenoData(example_sceset, aesth = aes_string(x = "log10(total_counts)",
#' y = "total_features", colour = "Mutation_Status"))
#'
plotPhenoData <- function(object, aesth=aes_string(x = "log10(total_counts)",
                                                   y = "total_features"),
                          theme_size = 10, ...) {
    ## We must have an SCESet object
    if (!is(object, "SCESet"))
        stop("object must be an SCESet object.")

    ## Define dataframe to pass to plotMetadata
    df_to_plot <- pData(object)

    ## Check that aesthetics make sense for feature names if used
    for (item in unlist(aesth)) {
        item <- as.character(item)
        if ( !(item %in% varLabels(object)) &&
             (item %in% featureNames(object)) ) {
            df_to_plot <- data.frame(df_to_plot, exprs(object)[item,])
            colnames(df_to_plot)[ncol(df_to_plot)] <- item
        }
    }

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(df_to_plot, aesth, ...)

    ## Define plotting theme
#     if ( requireNamespace("cowplot", quietly = TRUE) )
#         plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
#     else
#         plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot object
    plot_out
}


################################################################################

#' Plot feature (gene) data from an SCESet object
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to produce
#' a density plot of number of cells expressing the feature (requires
#' \code{calculateQCMetrics} to have been run on the SCESet object prior).
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param ... arguments passed to \code{\link{plotMetadata}}.
#'
#' @details Plot feature (gene) data from an SCESet object. If one variable is
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotFeatureData(example_sceset, aesth=aes(x=n_cells_exprs, y=pct_total_counts))
#'
plotFeatureData <- function(object,
                            aesth = aes_string(x = "n_cells_exprs",
                                               y = "prop_total_counts"),
                            theme_size = 10, ...) {
    ## We must have an SCESet object
    if (!is(object, "SCESet"))
        stop("object must be an SCESet object.")

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(fData(object), aesth, ...)

    ## Define plotting theme
#     if ( requireNamespace("cowplot", quietly = TRUE) )
#         plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
#     else
#         plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot object
    plot_out
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
#' Plot expression values from an SCESet object against transcript length values
#' defined in the SCESet object or supplied as an argument.
#'
#' @param object an \code{\link{SCESet}} object
#' @param tx_length transcript lengths to plot on the x-axis. Can be one of: (1) 
#' the name of a column of \code{fData(object)} containing the transcript length 
#' values, or (2) the name of an element of \code{assayData(object)} containing 
#' a matrix of transcript length values, or (3) a numeric vector of length equal
#' to the number of rows of \code{object} (number of features).
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{fData(object)} which will be used as a variable by which to colour
#' expression values on the plot.
#' @param shape_by optional character string supplying name of a column of
#' \code{fData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot.
#' @param size_by optional character string supplying name of a column of
#' \code{fData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot.
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
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' fd <- new("AnnotatedDataFrame", data = 
#' data.frame(gene_id = rownames(sc_example_counts), 
#'         feature_id = paste("feature", rep(1:500, each = 4), sep = "_"),
#'      median_tx_length = rnorm(2000, mean = 5000, sd = 500)))
#' rownames(fd) <- rownames(sc_example_counts)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd,
#' featureData = fd)
#' 
#' plotExprsVsTxLength(example_sceset, "median_tx_length")
#' plotExprsVsTxLength(example_sceset, "median_tx_length", show_smooth = TRUE)
#' plotExprsVsTxLength(example_sceset, "median_tx_length", show_smooth = TRUE,
#' show_exprs_sd = TRUE)
#' 
#' ## using matrix of tx length values in assayData(object)
#' mat <- matrix(rnorm(ncol(example_sceset) * nrow(example_sceset), mean = 5000,
#'  sd = 500), nrow = nrow(example_sceset))
#' dimnames(mat) <- dimnames(example_sceset)
#' set_exprs(example_sceset, "tx_len") <- mat
#' 
#' plotExprsVsTxLength(example_sceset, "tx_len", show_smooth = TRUE,
#' show_exprs_sd = TRUE)
#' 
#' ## using a vector of tx length values
#' plotExprsVsTxLength(example_sceset, rnorm(2000, mean = 5000, sd = 500))
#' 
plotExprsVsTxLength <- function(object, tx_length = "median_feat_eff_len", 
                                exprs_values = "exprs",
                                colour_by = NULL, shape_by = NULL,
                                size_by = NULL, xlab = NULL, 
                                show_exprs_sd = FALSE,
                                show_smooth = FALSE, alpha = 0.6, 
                                theme_size = 10, log2_values = FALSE, size = NULL, 
                                se = TRUE) {
    ## Check object is an SCESet object
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet")
    
    tx_length_values <- rep(NA, nrow(object))
    ## Check arguments are valid
    if ( length(tx_length) == 1 ) {
        if ( tx_length %in% Biobase::fvarLabels(object) )
            tx_length_values <- fData(object)[[tx_length]]
        else {
            if ( tx_length %in% names(Biobase::assayData(object)) ) {
                tx_length_mat <- Biobase::assayData(object)[[tx_length]]
                tx_length_values <- matrixStats::rowMedians(tx_length_mat)
            } else
                stop("the argument 'tx_length' should specify a column of fData(object) [see Biobase::fvarLabels(object)] or an element of assayData(object) [see names(assayData(object))")
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
    
    exprs_mat <- get_exprs(object, exprs_values)
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else
        ylab <- paste0("Expression (", exprs_values, ")")
    
    ## compute mean expression and sd of expression values
    exprs_mean <- rowMeans(exprs_mat)
    exprs_sd <- matrixStats::rowSds(exprs_mat)
    
    df_to_plot <- data.frame(tx_length_values, exprs_mean, exprs_sd,
                             ymin = exprs_mean - exprs_sd, 
                             ymax = exprs_mean + exprs_sd)
    
    ## check colour, size, shape arguments
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% Biobase::fvarLabels(object)) )
            stop("the argument 'colour_by' should specify a column of fData(object) [see fvarLabels(object)]")
        df_to_plot[[colour_by]] <- fData(object)[[colour_by]]
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% Biobase::fvarLabels(object)) )
            stop("the argument 'shape_by' should specify a column of fData(object) [see fvarLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")
        df_to_plot[[shape_by]] <- fData(object)[[shape_by]]
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% Biobase::fvarLabels(object)) )
            stop("the argument 'size_by' should specify a column of fData(object) [see fvarLabels(object)]")
        df_to_plot[[size_by]] <- fData(object)[[size_by]]
    }
    
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



