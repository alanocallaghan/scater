## Convenience function for computing QC metrics and adding to colData & rowData
### This file contains definitions for the following functions:
### * findImportantPCs
### * plotExplanatoryVariables
###
### * .calculateSilhouetteWidth
### * .getRSquared
### * .getTypeOfVariable

################################################################################

#' Find most important principal components for a given variable
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param variable character scalar providing a variable name (column from
#' \code{colData(object)}) for which to determine the most important PCs.
#' @param plot_type character string, indicating which type of plot to produce.
#' Default, \code{"pairs-pcs"} produces a pairs plot for the top 5 PCs based on
#' their R-squared with the variable of interest. A value of
#' \code{"pcs-vs-vars"} produces plots of the top PCs against the variable of
#' interest.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "counts",
#' "tpm", "fpkm" and "logcounts" (default), or anything else in the object added manually by 
#' the user.
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the PCA. Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in
#' \code{rownames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param theme_size numeric scalar providing base font size for ggplot theme.
#'
#' @details Plot the top 5 or 6 most important PCs (depending on the 
#' \code{plot_type} argument for a given variable. Importance here is defined as
#' the R-squared value from a linear model regressing each PC onto the variable 
#' of interest. 
#'
#' @return a \code{\link{ggplot}} plot object
#'
#' @import viridis
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#' example_sce <- calculateQCMetrics(example_sce)
#' findImportantPCs(example_sce, variable="total_features")
#'
findImportantPCs <- function(object, variable="total_features",
                             plot_type = "pcs-vs-vars", exprs_values = "logcounts",
                             ntop = 500, feature_set = NULL, 
                             scale_features = TRUE, theme_size = 10) {
    if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% rownames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in rownames(object)")
    }
    df_for_pca <- assay(object, exprs_values)
    if ( is.null(df_for_pca) )
        stop("The supplied 'exprs_values' argument not found in assayData(object). Try 'exprs' or similar.")
    if ( is.null(feature_set) ) {
        rv <- .rowVars(df_for_pca)
        feature_set <-
            order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    }
    df_for_pca <- df_for_pca[feature_set,]
    df_for_pca <- t(df_for_pca)

    ## Drop any features with zero variance
    keep_feature <- .colVars(df_for_pca) > 0.001
    keep_feature[is.na(keep_feature)] <- FALSE
    df_for_pca <- df_for_pca[, keep_feature]
    ## compute PCA
    pca <- prcomp(df_for_pca, retx = TRUE, center = TRUE, 
                  scale. = scale_features)
    colnames(pca$x) <- paste("component", 1:ncol(pca$x))
    if (!(variable %in% colnames(colData(object))))
        stop("variable not found in colData(object).
             Please make sure colData(object)[, variable] exists.")
    x <- colData(object)[, variable]
    x_na <- is.na(x)
    x <- x[!x_na]
    if (length(unique(x)) <= 1)
        stop("variable only has one unique value, so cannot determine important
             principal components.")
    ## Determine type of variable
    typeof_x <- .getTypeOfVariable(object, variable)
    if ( typeof_x == "discrete" ) {
        ## If x is a discrete variable
        x_int <- as.factor(x)
        ## Compute R-squared for each PC
        design <- model.matrix(~x_int)
    } else {
        ## If x is a continuous variable - use as a continuous variable
        design <- model.matrix(~x)
    }
    ## Get R-squared for each PC for the variable of interest
    pca_r_squared <- .getRSquared(t(pca$x[!x_na,]), design)
    ## Tidy up names and choose top 5 most important PCs for the variable
    # names(ave_sil_width) <- colnames(pca$x)
    names(pca_r_squared) <- colnames(pca$x)
    colnames(pca$x) <- paste0(colnames(pca$x), "\n(R-squared ",
                              formatC(signif(pca_r_squared, digits = 2),
                                      digits = 2, format = "fg", flag = "#"), ")")
    top5 <- order(pca_r_squared, decreasing = TRUE)[1:5]
    if ( plot_type == "pairs-pcs" ) {
        ## Define colours for points
        colour_by <- colData(object)[, variable]
        ## Generate a larger data.frame for pairs plot
        df_to_expand <- pca$x[, top5]
#         colnames(df_to_expand) <- colnames(pca$x)[, top5]
#         rownames(df_to_expand) <- sampleNames(object)
        names(df_to_expand) <- colnames(df_to_expand)
        gg1 <- .makePairs(df_to_expand)
        ## new data frame
        df_to_plot_big <- data.frame(gg1$all, colour_by)
        # colnames(df_to_plot_big)[-c(1:4)] <- get("variable")
        ## pairs plot
        plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) +
            geom_point(aes_string(fill = "colour_by"), colour = "gray40",
                       shape = 21, alpha = 0.65) +
            facet_grid(xvar ~ yvar, scales = "free") +
            stat_density(aes_string(x = "x",
                                    y = "(..scaled.. * diff(range(x)) + min(x))"),
                         data = gg1$densities, position = "identity",
                         colour = "grey20", geom = "line") +
            xlab("") +
            ylab("") +
            theme_bw(theme_size)
        plot_out <- .resolve_plot_colours(plot_out, colour_by, get("variable"),
                                          fill = TRUE)
        return(plot_out)
    } else {
        top6 <- order(pca_r_squared, decreasing = TRUE)[1:6]
        df_to_plot <- reshape2::melt(pca$x[, top6])
        xvar <- colData(object)[, variable]
        df_to_plot$xvar <- rep(xvar, 6)
        pcs_vars_plot <- ggplot(df_to_plot, aes_string(x = "xvar", y = "value"),
                                colour = "black") +
            facet_wrap(~ Var2, nrow = 3, scales = "free_y") +
            xlab(variable) +
            ylab("Principal component value") +
            theme_bw(theme_size)
        if ( typeof_x == "discrete") {
            pcs_vars_plot <- pcs_vars_plot +
                geom_violin(fill = "aliceblue", colour = "gray60",
                            alpha = 0.6, scale = "width") +
                geom_boxplot(width = 0.25, outlier.size = 0)
            if ( ncol(object) <= 150 ) {
                pcs_vars_plot <- pcs_vars_plot +
                    geom_dotplot(fill = "gray10", alpha = 0.6, binaxis = 'y',
                                 stackdir = 'center', dotsize = 1)
            }
        } else {
            pcs_vars_plot <- pcs_vars_plot +
                geom_point(fill = "gray10", alpha = 0.6, shape = 21) +
                stat_smooth(aes(group = 1), method = "lm", alpha = 0.3)
        }
        return(pcs_vars_plot)
    }
}

.getTypeOfVariable <- function(object, variable) {
    ## Extract variable
    x <- colData(object)[, variable]
    ## Get type
    if (is.character(x) || is.factor(x) || is.logical(x)) {
        typeof_x <- "discrete"
    } else {
        if (is.integer(x)) {
            if (length(unique(x)) > 10)
                typeof_x <- "continuous"
            else
                typeof_x <- "discrete"
        } else {
            if (is.numeric(x))
                typeof_x <- "continuous"
            else {
                x <- as.character(x)
                typeof_x <- "discrete"
                warning(paste0("Unrecognised variable type for ", variable,
". Variable being coerced to discrete. Please make sure colData(object)[, variable] is a proper discrete or continuous variable"))
            }
        }
    }
    typeof_x
}
