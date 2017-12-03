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

.choose_vis_values <- function(x, by, check_coldata = TRUE,
                               cell_control_default = FALSE,
                               check_features = FALSE,
                               exprs_values = "logcounts",
                               coerce_factor = FALSE, level_limit = NA) {
    ## This function looks through the visualization data and returns the
    ## values to be visualized. Either 'by' itself, or a column of colData,
    ## or a column of rowData, or the expression values of a feature.

    if (is.character(by)) {
        if (length(by) != 1L) {
            stop("'by' should be a character vector of length 1")
        }

        ## checking if it refers to a field in colData/rowData.
        vals <- NULL
        if (check_coldata) {
            if (by %in% colnames(colData(x))) {
                vals <- colData(x)[[by]]
            }
            unfound <- "colnames(colData(x))"
        } else {
            if (by %in% colnames(rowData(x))) {
                vals <- rowData(x)[[by]]
            }
            unfound <- "colnames(rowData(x))"
        }

        ## checking if it refers to a feature
        if (check_features) {
            if (is.null(vals) && by %in% rownames(x)) {
                vals <- assay(x, i = exprs_values)[by,]
            }
            unfound <- append(unfound, "rownames(x)")
        }

        ## throwing an error if still unfound
        if (is.null(vals)) {
            stop(sprintf("'%s' not found %s", by,
                         paste(sprintf("in '%s'", unfound), collapse = " or ")))
        }

    } else if (is.data.frame(by)) {
        if (ncol(by) != 1L) {
            stop("'by' should be a data frame with one column")
        } else if (nrow(by) != ncol(x)) {
            stop("'nrow(by)' should be equal to number of columns in 'x'")
        }

        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)

    } else {
        if (!is.null(by)) stop("invalid value of 'by' supplied")
        vals <- NULL

        ##  Switching to cell controls if desired.
        if (cell_control_default) {
            if ( "is_cell_control" %in% colnames(colData(x)) ) {
                by <- "is_cell_control"
                vals <- colData(x)[[by]]
            }
        }
    }

    # Checking the level limit.
    if (coerce_factor && !is.null(vals)) {
        vals <- factor(vals)
        if (level_limit < nlevels(vals))
            stop(sprintf("number of unique levels exceeds %i", level_limit))
    }

    return(list(name = by, val = vals))
}


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
#' cpm(example_sce) <- calculateCPM(example_sce, use.size.factors = FALSE)
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

    ## Use plyr to get the sequencing real estate accounted for by features
    nfeatures_total <- nrow(exprs_mat)
    seq_real_estate <- t(plyr::aaply(exprs_mat, 2, .fun = function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }))
    rownames(seq_real_estate) <- seq_len(nfeatures_total)
    nfeatures_to_plot <- nfeatures
    to_plot <- seq_len(nfeatures_to_plot)
    seq_real_estate_long <- reshape2::melt(seq_real_estate[to_plot, ],
                                           value.name = exprs_values)

    ## Get the proportion of the library accounted for by the top features
    prop_library <- reshape2::melt(t(t(seq_real_estate[to_plot, ]) /
                                         .general_colSums(exprs_mat)),
                                   value.name = "prop_library")
    colnames(seq_real_estate_long) <- c("Feature", "Cell", exprs_values)
    seq_real_estate_long$Proportion_Library <- prop_library$prop_library

    ## Add block and colour_by information if provided
    if ( !is.null(block1) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block1 = as.factor(rep(x[[block1]],
                                                         each = nfeatures_to_plot)))
    if ( !is.null(block2) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block2 = as.factor(rep(x[[block2]],
                                                         each = nfeatures_to_plot)))
    if ( !is.null(colour_by) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, colour_by = rep(colour_by_vals,
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

#' Run PCA for a SingleCellExperiment object
#'
#' Perform a principal components analysis (PCA) based on the data stored in 
#' a \code{\link{SingleCellExperiment}} object. 
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the PCA. Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal
#' components to obtain from \code{\link{prcomp}}.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other named element of the \code{assays} slot of the \code{SingleCellExperiment}
#' object that can be accessed with the \code{assay} function.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param pca_data_input character argument defining which data should be used
#' as input for the PCA. Possible options are \code{"logcounts"} (default), which
#' uses log-count data to produce a PCA at the cell level; \code{"coldata"} or
#' \code{"pdata"} (for backwards compatibility) which uses numeric variables
#' from \code{colData(object)} to do PCA at the cell level; and
#' \code{"rowdata"} which uses numeric variables from \code{rowData(object)} to
#' do PCA at the feature level.
#' @param selected_variables character vector indicating which variables in
#' \code{colData(object)} to use for the phenotype-data based PCA. Ignored if
#' the argument \code{pca_data_input} is anything other than \code{"pdata"}
#' or \code{"coldata"}.
#' @param detect_outliers logical, should outliers be detected based on PCA
#' coordinates generated from column-level metadata? Only an option when 
#' \code{pca_data_input} argument is \code{"pdata"} or \code{"coldata"}. 
#' Default is \code{FALSE}.
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
#' If the arguments \code{detect_outliers} and \code{return_SCE} are both
#' \code{TRUE}, then the element \code{$outlier} is added to the pData
#' (phenotype data) slot of the \code{SingleCellExperiment} object. This element contains
#' indicator values about whether or not each cell has been designated as an
#' outlier based on the PCA. These values can be accessed for filtering
#' low quality cells with, for example, \code{example_sce$outlier}.
#'
#' When \code{pca_data_input="pdata"} or \code{"coldata"}, the selected variables 
#' default to a vector containing:
#' \itemize{
#' \item \code{"pct_counts_top_100_features"}
#' \item \code{"total_features"}
#' \item \code{"pct_counts_feature_control"}
#' \item \code{"total_features_feature_control"}
#' \item \code{"log10_total_counts_endogenous"}
#' \item \code{"log10_total_counts_feature_control"}
#' }
#' These metrics were chosen due to their utility in distinguishing low-quality
#' libraries. However, they can be overriden by setting \code{selected_variables}
#' manually. In particular, \code{"log10_total_counts"} is more useful than 
#' the \code{_endogenous} and \code{_control} metrics when spike-ins are not
#' available.
#'
#' @rdname runPCA
#' @seealso \code{\link[scater]{plotPCA}}
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runPCA(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runPCA <- function(object, ntop=500, ncomponents=2, exprs_values = "logcounts",
       feature_set = NULL, scale_features = TRUE, pca_data_input = "logcounts",
       selected_variables = NULL, detect_outliers = FALSE) {

    if ( pca_data_input == "pdata" || pca_data_input == "coldata" ) {
        #use_variable <- sapply(pData(object), is.double)
        ## select pData features to use
        if ( is.null(selected_variables) ) {
            selected_variables <- c("pct_counts_top_100_features",
                                    "total_features",
                                    "pct_counts_feature_control",
                                    "total_features_feature_control",
                                    "log10_total_counts_endogenous",
                                    "log10_total_counts_feature_control")
        }

        col_data_names <- colnames(colData(object)) 
        use_variable <- col_data_names %in% selected_variables
        vars_not_found <- !(selected_variables %in% col_data_names)
        if ( any(vars_not_found) ) {
            for (missing_var in selected_variables[vars_not_found]) {
                warning(sprintf("selected variable '%s' not found in 'colData(object)'", 
                                missing_var))
            }
        }
        ## scale double variables
        exprs_to_plot <- scale(colData(object)[, use_variable],
                               scale = scale_features)
    } else {
        exprs_mat <- assay(object, i = exprs_values)
        
        # Choosing a set of features, if null.
        if (is.null(feature_set)) {
            rv <- .general_rowVars(exprs_mat)
            o <- order(rv, decreasing = TRUE)
            feature_set <- o[seq_len(min(ntop, length(rv)))]
        }

        # Subsetting to the desired features (do NOT move below 'scale()')
        exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]
        ## Standardise expression if scale_features argument is TRUE
        exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
    }

    ## Drop any features with zero variance
    keep_feature <- .general_colVars(exprs_to_plot) > 0.001
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[, keep_feature]

    ## conduct outlier detection
    if ( detect_outliers ) {
        if ( !(pca_data_input == "pdata" || pca_data_input == "coldata") ) {
            warning("outlier detection requires 'pca_data_input=\"coldata\"'")
        } else {
            outliers <- mvoutlier::pcout(exprs_to_plot, makeplot = FALSE,
                                         explvar = 0.5, crit.M1 = 0.9,
                                         crit.c1 = 5, crit.M2 = 0.9,
                                         crit.c2 = 0.99, cs = 0.25,
                                         outbound = 0.05)
             outlier <- !as.logical(outliers$wfinal01)
             object$outlier <- outlier
        }
    }

    ## Compute PCA
    pca <- prcomp(exprs_to_plot, rank. = ncomponents)
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar

    # Saving the stuff in the object.
    reducedDim(object, "PCA") <- pcs
    return(object)
}


#' Plot PCA for a SingleCellExperiment object
#'
#' Produce a principal components analysis (PCA) plot of two or more principal
#' components for a \code{\link{SingleCellExperiment}} object. 
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced.
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}. 
#' @param rerun logical, should PCA be recomputed even if \code{object} contains a
#' \code{"PCA"} element in the \code{reducedDims} slot?
#' @param return_SCE logical, should the function return a \code{SingleCellExperiment}
#' object with principal component values for cells in the
#' \code{reducedDim} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned. This will be deprecated in the next 
#' development cycle in favour of directly calling \code{\link{runPCA}}.
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param run_args Arguments to pass to \code{\link{runPCA}} when \code{rerun=TRUE}
#' or if there is no existing \code{"PCA"} element in the \code{reducedDims} slot.
#'
#' @details For back-compatabibility purposes, users can specify arguments to 
#' \code{\link{runPCA}} in \code{...}. This will trigger a warning as it will be
#' deprecated in the next development cycle.
#'
#' @return either a ggplot plot object or an SingleCellExperiment object
#'
#' @name plotPCA
#' @rdname plotPCA
#' @aliases plotPCA plotPCA,SingleCellExperiment-method
#' @importFrom BiocGenerics plotPCA
#' @seealso \code{\link{runPCA}}
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sce)
#' plotPCA(example_sce, colour_by = "Cell_Cycle")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' size_by = "Mutation_Status")
#' plotPCA(example_sce, shape_by = "Treatment", size_by = "Mutation_Status")
#' plotPCA(example_sce, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
#' ## experiment with legend
#' example_subset <- example_sce[, example_sce$Treatment == "treat1"]
#' plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", legend = "all")
#'
#' plotPCA(example_sce, shape_by = "Treatment", return_SCE = TRUE)
#'
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sce, ncomponents = 8)
#' plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
plotPCASCE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, 
                       ncomponents = 2, run_args=list()) {
    
    new_args <- .disambiguate_args(...)    
    run_args <- c(run_args, new_args$run)
    plot_args <- new_args$plot

    ## Running PCA if necessary.
    if (!("PCA" %in% names(reducedDims(object))) || rerun) {
        object <- do.call(runPCA, c(list(object = object, ncomponents = ncomponents),
                                    run_args))
    }

    plot_out <- do.call(plotReducedDim, 
                        c(list(object = object, ncomponents = ncomponents, 
                               use_dimred = "PCA"), plot_args))

    ## Plot PCA and return appropriate object
    if (return_SCE) {
        .Deprecated(msg="'return_SCE=TRUE' is deprecated, use 'runPCA' instead")
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

.disambiguate_args <- function(...) 
# This function is only necessary to provide some protection in the transition 
# from having running arguments in "..." to plotting arguments in "...". It can
# be removed in the next development cycle. 
{
    plot_arg_names <- union(names(formals(plotReducedDim)), 
                            names(formals(plotReducedDimDefault)))
    extra_args <- list(...)
    for_plotting <- !is.na(pmatch(names(extra_args), plot_arg_names))
    if (!all(for_plotting)) { 
        warning(sprintf("non-plotting arguments like '%s' should go in 'run_args'", 
                        names(extra_args)[!for_plotting][1]))
    }
    return(list(plot=extra_args[for_plotting],
                run=extra_args[!for_plotting]))
}

#' @rdname plotPCA
#' @aliases plotPCA
#' @export
setMethod("plotPCA", "SingleCellExperiment", plotPCASCE)

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

#' Run t-SNE for a SingleCellExperiment object
#'
#' Perform t-stochastic neighbour embedding (t-SNE) based on the data stored in 
#' a \code{\link{SingleCellExperiment}} object.
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the t-SNE Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of t-SNE
#' components to obtain.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values),
#' or any other named element of the \code{assayData} slot of the \code{SingleCellExperiment}
#' object that can be accessed with the \code{assay} function.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the t-SNE calculation. If character, entries must all be
#' in \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param perplexity numeric scalar value defining the "perplexity parameter"
#' for the t-SNE plot. Passed to \code{\link[Rtsne]{Rtsne}} - see documentation
#' for that package for more details.
#' @param ... Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#'
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to
#' compute the t-SNE. Note that the algorithm is not deterministic, so different
#' runs of the function will produce differing plots (see \code{\link{set.seed}}
#' to set a random seed for replicable results). The value of the
#' \code{perplexity} parameter can have a large effect on the resulting plot, so
#' it can often be worthwhile to try multiple values to find the most appealing
#' visualisation and to ensure that the conclusions are robust.
#'
#' @rdname runTSNE
#' @seealso 
#' \code{\link[Rtsne]{Rtsne}},
#' \code{\link[scater]{plotTSNE}}
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runTSNE(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runTSNE <- function(object, ntop = 500, ncomponents = 2, exprs_values = "logcounts",
        feature_set = NULL, use_dimred = NULL, n_dimred = NULL, scale_features = TRUE,
        rand_seed = NULL, perplexity = floor(ncol(object) / 5), ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr
        do_pca <- FALSE
        pca_dims <- ncol(vals)

    } else {
        ## Define an expression matrix depending on which values we're
        ## using
        exprs_mat <- assay(object, i = exprs_values)

        ## Define features to use: either ntop, or if a set of features is
        ## defined, then those
        if ( is.null(feature_set) ) {
            rv <- .general_rowVars(exprs_mat)
            ntop <- min(ntop, length(rv))
            feature_set <- order(rv, decreasing = TRUE)[seq_len(ntop)]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat[feature_set,,drop = FALSE]
        keep_feature <- .general_rowVars(vals) > 0.001
        keep_feature[is.na(keep_feature)] <- FALSE
        vals <- vals[keep_feature,,drop = FALSE]

        ## Standardise expression if stand_exprs(object) is null
        vals <- t(vals)
        if (scale_features) {
            vals <- scale(vals, scale = TRUE)
        }
        do_pca <- TRUE
        pca_dims <- max(50, ncol(object))
    }

    # Actually running the Rtsne step.
    if ( !is.null(rand_seed) )
        set.seed(rand_seed)
    tsne_out <- Rtsne::Rtsne(vals, initial_dims = pca_dims, pca = do_pca,
                             perplexity = perplexity, dims = ncomponents,...)
    reducedDim(object, "TSNE") <- tsne_out$Y
    return(object)
}


#' Plot t-SNE for an SingleCellExperiment object
#'
#' Produce a t-distributed stochastic neighbour embedding (t-SNE) plot of two
#' components for a \code{SingleCellExperiment} object.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}.
#' @param ncomponents numeric scalar indicating the number of t-SNE
#' components to plot, starting from the first t-SNE component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of component 1 vs component
#' 2 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#' top components is produced. NB: computing more than two components for t-SNE
#' can become very time consuming.
#' @param return_SCE logical, should the function return a \code{SingleCellExperiment}
#' object with principal component values for cells in the
#' \code{reducedDims} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned. This will be deprecated in the next 
#' development cycle in favour of directly calling \code{\link{runTSNE}}.
#' @param rerun logical, should PCA be recomputed even if \code{object} contains a
#' \code{"TSNE"} element in the \code{reducedDims} slot?
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param run_args Arguments to pass to \code{\link{runTSNE}} when \code{rerun=TRUE}
#' or if there is no existing \code{"TSNE"} element in the \code{reducedDims} slot.
#'
#' @details For back-compatabibility purposes, users can specify arguments to 
#' \code{\link{runTSNE}} in \code{...}. This will trigger a warning as it will be
#' deprecated in the next development cycle.
#'
#' @return If \code{return_SCE} is \code{TRUE}, then the function returns a
#' \code{SingleCellExperiment} object, otherwise it returns a \code{ggplot} object.
#' @name plotTSNE
#' @rdname plotTSNE
#' @aliases plotTSNE plotTSNE,SingleCellExperiment-method
#'
#' @export
#' @seealso
#' \code{\link{runTSNE}}
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International
#' Conference on Learning Representations, 2013.
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' ## Examples plotting t-SNE
#' plotTSNE(example_sce, perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' size_by = "Mutation_Status", perplexity = 10)
#' plotTSNE(example_sce, shape_by = "Treatment", size_by = "Mutation_Status",
#' perplexity = 5)
#' plotTSNE(example_sce, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status", perplexity = 5)
#'
#' plotTSNE(example_sce, shape_by = "Treatment", return_SCE = TRUE,
#' perplexity = 10)
#'
#'
plotTSNE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE,
                   rerun = FALSE, ncomponents = 2, run_args=list()) {

    new_args <- .disambiguate_args(...)    
    run_args <- c(run_args, new_args$run)
    plot_args <- new_args$plot

    # Re-running t-SNE if necessary.
    if ( !("TSNE" %in% names(reducedDims(object))) || rerun) {
        object <- do.call(runTSNE, c(list(object=object, ncomponents = ncomponents),
                                     run_args))
    }

    plot_out <- do.call(plotReducedDim, 
                        c(list(object = object, ncomponents = ncomponents, 
                               use_dimred = "TSNE"), plot_args))

    if (return_SCE) {
        .Deprecated(msg="'return_SCE=TRUE' is deprecated, use 'runTSNE' instead")
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

################################################################################
### plotDiffusionMap

#' Create a diffusion map for an SingleCellExperiment object
#'
#' Produce a diffusion map plot using data stored in a \code{SingleCellExperiment} 
#' object.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the diffusion map. Default is \code{500}, but any \code{ntop}
#' argument is overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of diffusion
#' components to obtain.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other named element of the \code{assayData} slot of the \code{SingleCellExperiment}
#' object that can be accessed with the \code{assay} function.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the diffusion map. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param ... Additional arguments to pass to \code{\link[destiny]{DiffusionMap}}.
#'
#' @details The function \code{\link[destiny]{DiffusionMap}} is used internally
#' to compute the diffusion map.
#'
#' @export
#' @rdname runDiffusionMap
#' @seealso 
#' \code{\link[destiny]{destiny}},
#' \code{\link[scater]{plotDiffusionMap}}
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runDiffusionMap(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runDiffusionMap <- function(object, ntop = 500, ncomponents = 2, feature_set = NULL,
        exprs_values = "logcounts", scale_features = TRUE, use_dimred=NULL, n_dimred=NULL,
        rand_seed = NULL, ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results.
        vals <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            vals <- vals[,seq_len(n_dimred),drop = FALSE]
        }
    } else {
        ## Define an expression matrix depending on which values we're
        ## using
        exprs_mat <- assay(object, i = exprs_values)

        ## Define features to use: either ntop, or if a set of features is
        ## defined, then those
        if ( is.null(feature_set) ) {
            rv <- .general_rowVars(exprs_mat)
            feature_set <-
                order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat
        vals <- vals[feature_set,,drop = FALSE]
        keep_feature <- .general_rowVars(vals) > 0.001
        keep_feature[is.na(keep_feature)] <- FALSE
        vals <- vals[keep_feature,,drop = FALSE]

        ## Standardise expression if indicated by scale_features argument
        vals <- t(vals)
        if (scale_features) {
            vals <- scale(vals, scale = TRUE)
        }
    }

    ## Compute DiffusionMap
    if ( !is.null(rand_seed) )
        set.seed(rand_seed)
    difmap_out <- destiny::DiffusionMap(vals, ...)

    reducedDim(object, "DiffusionMap") <- difmap_out@eigenvectors[, seq_len(ncomponents), drop = FALSE]
    return(object)
}

#' Plot a diffusion map for a SingleCellExperiment object
#'
#' Produce a diffusion map plot of two components for a \code{SingleCellExperiment} object. 
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}.
#' @param ncomponents numeric scalar indicating the number of diffusion 
#' components to plot, starting from the first diffusion component. Default
#' is 2. If \code{ncomponents} is 2, then a scatterplot of component 1 vs
#' component 2 is produced. If \code{ncomponents} is greater than 2, a pairs
#' plots for the top components is produced. NB: computing many components for
#' the diffusion map can become time consuming.
#' @param return_SCE logical, should the function return a \code{SingleCellExperiment}
#' object with principal component values for cells in the
#' \code{reducedDims} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned. This will be deprecated in the next 
#' development cycle in favour of directly calling \code{\link{runDiffusionMap}}.
#' @param rerun logical, should PCA be recomputed even if \code{object} contains a
#' \code{"DiffusionMap"} element in the \code{reducedDims} slot?
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param run_args Arguments to pass to \code{\link{runDiffusionMap}} when \code{rerun=TRUE}
#' or if there is no existing \code{"DiffusionMap"} element in the \code{reducedDims} slot.
#'
#' @details For back-compatabibility purposes, users can specify arguments to 
#' \code{\link{runDiffusionMap}} in \code{...}. This will trigger a warning as it will be
#' deprecated in the next development cycle.
#'
#' @return If \code{return_SCE} is \code{TRUE}, then the function returns an
#' \code{SingleCellExperiment} object, otherwise it returns a \code{ggplot} object.
#' @name plotDiffusionMap
#' @rdname plotDiffusionMap
#' @aliases plotDiffusionMap plotDiffusionMap,SingleCellExperiment-method
#'
#' @export
#' @seealso
#' \code{\link{runDiffusionMap}}
#' @references
#' Haghverdi L, Buettner F, Theis FJ. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics. 2015; doi:10.1093/bioinformatics/btv325
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' \dontrun{
#' ## Examples plotting diffusion maps
#' plotDiffusionMap(example_sce)
#' plotDiffusionMap(example_sce, colour_by = "Cell_Cycle")
#' plotDiffusionMap(example_sce, colour_by = "Cell_Cycle",
#' shape_by = "Treatment")
#' plotDiffusionMap(example_sce, colour_by = "Cell_Cycle",
#' shape_by = "Treatment", size_by = "Mutation_Status")
#' plotDiffusionMap(example_sce, shape_by = "Treatment",
#' size_by = "Mutation_Status")
#' plotDiffusionMap(example_sce, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
#' plotDiffusionMap(example_sce, shape_by = "Treatment",
#' return_SCE = TRUE)
#' }
#'
plotDiffusionMap <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, 
      rerun = FALSE, ncomponents = 2, run_args=list()) {

    new_args <- .disambiguate_args(...)    
    run_args <- c(run_args, new_args$run)
    plot_args <- new_args$plot

    # Re-running the diffusion map if necessary.
    if ( !("DiffusionMap" %in% names(reducedDims(object))) || rerun) {
        object <- do.call(runDiffusionMap, c(list(object, ncomponents = ncomponents),
                                             run_args))
    }

    plot_out <- do.call(plotReducedDim, 
                        c(list(object = object, ncomponents = ncomponents, 
                               use_dimred = "DiffusionMap"), plot_args))

    if (return_SCE) {
        .Deprecated(msg="'return_SCE=TRUE' is deprecated, use 'runDiffusionMap' instead")
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

################################################################################
### plotMDS

#' Run MDS for a SingleCellExperiment object
#'
#' Perform multi-dimensional scaling using data stored in a
#' \code{SingleCellExperiment} object. 
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the diffusion map. Default is \code{500}, but any \code{ntop}
#' argument is overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of MDS dimensions
#' to obtain.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other named element of the \code{assayData} slot of the \code{SingleCellExperiment}
#' object that can be accessed with the \code{assay} function.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the diffusion map. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param method string specifying the type of distance to be computed between cells.
#'
#' @details The function \code{\link{cmdscale}} is used internally to
#' compute the multidimensional scaling components to plot.
#'
#' @export
#' @rdname runMDS
#' @seealso \code{\link[scater]{plotMDS}}
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runMDS(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runMDS <- function(object, ntop = 500, ncomponents = 2, feature_set = NULL,
        exprs_values = "logcounts", scale_features = TRUE, use_dimred=NULL, n_dimred=NULL,
        method = "euclidean") {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results.
        vals <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            vals <- vals[,seq_len(n_dimred),drop = FALSE]
        }
    } else {
        ## Define an expression matrix depending on which values we're
        ## using
        exprs_mat <- assay(object, i = exprs_values)

        ## Define features to use: either ntop, or if a set of features is
        ## defined, then those
        if ( is.null(feature_set) ) {
            rv <- .general_rowVars(exprs_mat)
            feature_set <-
                order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat
        vals <- vals[feature_set,,drop = FALSE]
        keep_feature <- .general_rowVars(vals) > 0.001
        keep_feature[is.na(keep_feature)] <- FALSE
        vals <- vals[keep_feature,,drop = FALSE]

        ## Standardise expression if indicated by scale_features argument
        vals <- t(vals)
        if (scale_features) {
            vals <- scale(vals, scale = TRUE)
        }
    }

    cell_dist <- stats::dist(vals, method = method)
    mds_out <- cmdscale(cell_dist, k = ncomponents)
    reducedDim(object, "MDS") <- mds_out
    object
}

#' Produce a multidimensional scaling plot for a SingleCellExperiment object
#'
#' Produce an MDS plot from the cell pairwise distance data in a
#' \code{SingleCellExperiment} object. 
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}.
#' @param ncomponents numeric scalar indicating the number of MDS dimensions
#' to plot, starting from the first dimension. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced. NB: computing more than two components for t-SNE can become very
#' time consuming.
#' @param return_SCE logical, should the function return a \code{SingleCellExperiment}
#' object with principal component values for cells in the
#' \code{reducedDims} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned. This will be deprecated in the next 
#' development cycle in favour of directly calling \code{\link{runMDS}}.
#' @param rerun logical, should PCA be recomputed even if \code{object} contains a
#' \code{"MDS"} element in the \code{reducedDims} slot?
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param run_args Arguments to pass to \code{\link{runMDS}} when \code{rerun=TRUE}
#' or if there is no existing \code{"MDS"} element in the \code{reducedDims} slot.
#'
#' @details For back-compatabibility purposes, users can specify arguments to 
#' \code{\link{runMDS}} in \code{...}. This will trigger a warning as it will be
#' deprecated in the next development cycle.
#'
#' @return If \code{return_SCE} is \code{TRUE}, then the function returns an
#' \code{SingleCellExperiment} object, otherwise it returns a \code{ggplot} object.
#' @name plotMDS
#' @aliases plotMDS plotMDS,SingleCellExperiment-method
#' @seealso \code{\link{runMDS}}
#'
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' ## Examples plotting
#' plotMDS(example_sce)
#' plotMDS(example_sce, colour_by = "Cell_Cycle")
#' plotMDS(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#'
#' ## define cell-cell distances differently
#' plotMDS(example_sce, colour_by = "Cell_Cycle",
#' shape_by = "Treatment", size_by = "Mutation_Status", method = "canberra")
#'
plotMDS <- function(object, ..., ncomponents = 2, return_SCE = FALSE,
                    rerun = FALSE, draw_plot = TRUE, run_args=list()) {

    new_args <- .disambiguate_args(...)    
    run_args <- c(run_args, new_args$run)
    plot_args <- new_args$plot

    # Re-running the MDS if necessary.
    if ( !("MDS" %in% names(reducedDims(object))) || rerun) {
        object <- do.call(runMDS, c(list(object = object, ncomponents = ncomponents),
                                    run_args))
    }

    plot_out <- do.call(plotReducedDim,
                        c(list(object = object, ncomponents = ncomponents, 
                               use_dimred = "MDS"), plot_args))

    if (return_SCE) {
        .Deprecated(msg="'return_SCE=TRUE' is deprecated, use 'runMDS' instead")
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

################################################################################
### plotReducedDim

#' Plot reduced dimension representation of cells
#'
#' @param object an \code{SingleCellExperiment} object with a non-NULL \code{reducedDimension}
#' slot.
#' @param use_dimred character, name of reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object to plot (e.g. "PCA", "TSNE", etc).
#' @param df_to_plot data.frame containing a reduced dimension represenation of
#' cells and optional metadata for the plot.
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of Dimension 2 vs Dimension
#' 1 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#'  top dimensions is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, a
#' data frame with one column containing values to map to colours for all cells.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param percentVar numeric vector giving the proportion of variance in
#' expression explained by each reduced dimension. Only expected to be used
#' internally in the \code{\link[scater]{plotPCA}} function.
#' @param exprs_values a string specifying the expression values to use for
#' colouring the points, if \code{colour_by} or \code{size_by} are set as feature names.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param add_ticks logical scalar indicating whether ticks should be drawn
#' on the axes corresponding to the location of each point.
#' @param ... optional arguments (from those listed above) passed to
#' \code{plotReducedDimDefault}
#'
#' @details The function \code{plotReducedDim.default} assumes that the first
#' \code{ncomponents} columns of \code{df_to_plot} contain the reduced dimension
#'  components to plot, and that any subsequent columns define factors for
#'  \code{colour_by}, \code{shape_by} and \code{size_by} in the plot.
#'
#' @return a ggplot plot object
#'
#' @name plotReducedDim
#' @aliases plotReducedDim plotReducedDim,SingleCellExperiment-method plotReducedDim,data.frame-method
#' @import viridis
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' reducedDim(example_sce, "PCA") <- prcomp(t(exprs(example_sce)), scale. = TRUE)$x
#' plotReducedDim(example_sce, "PCA")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", shape_by="Treatment")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", size_by="Treatment")
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle", shape_by="Treatment")
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#'
plotReducedDimDefault <- function(df_to_plot, ncomponents=2, percentVar=NULL,
                           colour_by=NULL, shape_by=NULL, size_by=NULL,
                           theme_size = 10, legend = "auto", add_ticks=TRUE) {

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
            theme_bw(theme_size)
        if (add_ticks) {
            plot_out <- plot_out + geom_rug(colour = "gray20", alpha = 0.65)
        }
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

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)

    ## remove legend if so desired
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")

    ## Return plot
    plot_out
}

#' @rdname plotReducedDim
#' @aliases plotReducedDim
#' @export
plotReducedDim <- function(object, use_dimred, ncomponents = 2,
                              colour_by = NULL, shape_by = NULL, size_by = NULL,
                              exprs_values = "logcounts", percentVar = NULL, ...) {

    ## Check arguments are valid
    colour_by_out <- .choose_vis_values(
        object, colour_by, cell_control_default = TRUE, check_features = TRUE,
        exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    shape_by_out <- .choose_vis_values(
        object, shape_by, cell_control_default = TRUE, coerce_factor = TRUE,
        level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val

    size_by_out <- .choose_vis_values(
        object, size_by, check_features = TRUE, exprs_values = exprs_values)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    ## Extract reduced dimension representation of cells
    red_dim <- reducedDim(object, use_dimred)
    if ( ncomponents > ncol(red_dim) )
        stop("ncomponents to plot is larger than number of columns of reducedDimension(object)")
    if (is.null(percentVar)) {
        percentVar <- attr(red_dim, "percentVar")
    }

    ## Define data.frame for plotting
    df_to_plot <- data.frame(red_dim[, seq_len(ncomponents),drop=FALSE])
    df_to_plot$colour_by <- colour_by_vals
    df_to_plot$shape_by <- shape_by_vals
    df_to_plot$size_by <- size_by_vals

    ## Call default method to make the plot
    plotReducedDimDefault(df_to_plot, ncomponents = ncomponents, percentVar = percentVar,
                colour_by = colour_by, shape_by = shape_by, size_by = size_by, ...)
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
                tx_length_values <- matrixStats::rowMedians(tx_length_mat)
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
    exprs_sd <- sqrt(.general_rowVars(exprs_mat))

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



