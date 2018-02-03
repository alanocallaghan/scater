## Convenience function for computing QC metrics and adding to colData & rowData
### This file contains definitions for the following functions:
### * findImportantPCs
### * plotExplanatoryVariables
### * plotHighestExprs
### * plotQC
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



#' @importFrom limma lmFit
.getRSquared <- function(y, design, chunk=NULL) { 
    QR <- qr(design)
    if (!is.null(chunk)) {
        ngenes <- nrow(y)
        by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
        sst <- ssr <- numeric(nrow(y))
        for (element in levels(by.chunk)) {
            current <- by.chunk==element
            out <- .getRSquared_internal(QR, y[current,, drop = FALSE])
            sst[current] <- out$sst
            ssr[current] <- out$ssr
        }
    } else {
        out <- .getRSquared_internal(QR, y)
        sst <- out$sst
        ssr <- out$ssr
    }
        
    # Return proportion of variance explained    
    (ssr/sst)
}

.getRSquared_internal <- function(QR, y) {
    ## Compute total sum of squares
    sst <- .rowVars(y) * (ncol(y)-1)    
    ## Compute residual sum of squares
    effects <- qr.qty(QR, as.matrix(t(y)))
    ssr <- sst - colSums(effects[-seq_len(QR$rank),, drop = FALSE] ^ 2) # no need for .colSums, as this is always dense.
    return(list(sst = sst, ssr = ssr))
}

################################################################################

#' Plot the features with the highest expression values
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param col_by_variable variable name (must be a column name of colData(object))
#' to be used to assign colours to cell-level values.
#' @param n numeric scalar giving the number of the most expressed features to
#' show. Default value is 50.
#' @param drop_features a character, logical or numeric vector indicating which
#' features (e.g. genes, transcripts) to drop when producing the plot. For
#' example, control genes might be dropped to focus attention on contribution
#' from endogenous rather than synthetic genes.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "counts" (default),
#' "tpm", "fpkm" and "logcounts".
#' @param feature_names_to_plot character scalar indicating which column of the 
#' rowData slot in the \code{object} is to be used for the feature names 
#' displayed on the plot. Default is \code{NULL}, in which case 
#' \code{rownames(object)} is used.
#' @param as_percentage logical scalar indicating whether percentages should be
#' plotted. If \code{FALSE}, the raw \code{exprs_values} are shown instead.
#'
#' @details Plot the percentage of counts accounted for by the top n most highly
#' expressed features across the dataset.
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce, 
#' feature_controls = list(set1 = 1:500))
#' plotHighestExprs(example_sce, col_by_variable="total_features")
#' plotHighestExprs(example_sce, col_by_variable="Mutation_Status")
#' plotQC(example_sce, type = "highest-express")
#'
plotHighestExprs <- function(object, col_by_variable = "total_features", n = 50,
                             drop_features = NULL, exprs_values = "counts",
                             feature_names_to_plot = NULL, as_percentage = TRUE) {
    ## Figure out which features to drop
    if ( !(is.null(drop_features) || length(drop_features) == 0) ) {
        if (is.character(drop_features)) {
            drop_features <- which(rownames(object) %in% drop_features)
        }
        if (is.logical(drop_features)) {
            object <- object[!drop_features,]
        } else {
            object <- object[-drop_features,]
        }
    }

    ## Define expression values to be used
    ## Find the most highly expressed features in this dataset
    exprs_mat <- assay(object, exprs_values, withDimnames=FALSE)
    ave_exprs <- .rowSums(exprs_mat)
    oo <- order(ave_exprs, decreasing=TRUE)
    chosen <- oo[seq_len(n)]

    ## define feature names for plot
    rdata <- rowData(object)
    if (is.null(feature_names_to_plot) || 
        is.null(rowData(object)[[feature_names_to_plot]])) {
        feature_names <- rownames(object)
    } else {
        feature_names <- rowData(object)[[feature_names_to_plot]]
        feature_names <- as.character(feature_names)
    }
    rownames(exprs_mat) <- feature_names 
    rdata$Feature <- factor(feature_names, levels=feature_names[rev(oo)])

    ## Check if is_feature_control is defined
    if ( is.null(rdata$is_feature_control) ) { 
        rdata$is_feature_control <- rep(FALSE, nrow(rdata))
    }

    ## Determine percentage expression accounted for by top features across all
    ## cells, and determine percentage of counts for top features by cell
    df_exprs_by_cell <- t(exprs_mat[chosen,])
    if (as_percentage) { 
        total_exprs <- sum(ave_exprs)
        top50_pctage <- 100 * sum(ave_exprs[chosen]) / total_exprs
        df_exprs_by_cell <- 100 * df_exprs_by_cell / .colSums(exprs_mat)
        pct_total <- 100 * ave_exprs / total_exprs
        rdata[["pct_total"]] <- pct_total
    } else {
        rdata[[exprs_values]] <- ave_exprs
    }
    df_exprs_by_cell <- as.matrix(df_exprs_by_cell)

    ## Melt dataframe so it is conducive to ggplot
    if ( is.null(rownames(rdata)) ) { 
        rownames(rdata) <- as.character(rdata$Feature)
    }
    df_exprs_by_cell_long <- reshape2::melt(df_exprs_by_cell)
    colnames(df_exprs_by_cell_long) <- c("Cell", "Tags", "value")
    df_exprs_by_cell_long$Feature <- factor(
        rdata[as.character(df_exprs_by_cell_long$Tags), "Feature"],
        levels = as.character(rdata$Feature[rev(chosen)]))
    
    ## Check that variable to colour points exists
    if (!(col_by_variable %in% colnames(colData(object)))) {
        warning(sprintf("'%s' not found in colData(object)", col_by_variable))
        plot_x <- FALSE
        aes_to_use <- aes_string(y="Feature", x="value")
    } else {
        plot_x <- TRUE
        x <- colData(object)[, col_by_variable]
        typeof_x <- .getTypeOfVariable(object, col_by_variable)

        ## Add colour variable information
        if (typeof_x == "discrete") {
            df_exprs_by_cell_long$colour_by <- factor(x)
        } else {
            df_exprs_by_cell_long$colour_by <- x
        }

        aes_to_use <- aes_string(y="Feature", x="value", colour="colour_by")
    }

    ## Make plot
    plot_most_expressed <- ggplot(df_exprs_by_cell_long, aes_to_use) +
        geom_point(alpha = 0.6, shape = 124)

    if (as_percentage) { 
        plot_most_expressed <- plot_most_expressed + ggtitle(paste0("Top ", n, 
            " account for ", format(top50_pctage, digits = 3), "% of total")) +
            xlab(paste0("% of total ", exprs_values))
        legend_val <- "as.numeric(pct_total)"
    } else {
        plot_most_expressed <- plot_most_expressed + xlab(exprs_values)
        legend_val <- sprintf("as.numeric(%s)", exprs_values)
    }

    plot_most_expressed <- plot_most_expressed + ylab("Feature") + theme_bw(8) +
        theme(legend.position = c(1, 0), legend.justification = c(1, 0),
              axis.text.x = element_text(colour = "gray35"),
              axis.text.y = element_text(colour = "gray35"),
              axis.title.x = element_text(colour = "gray35"),
              axis.title.y = element_text(colour = "gray35"),
              title = element_text(colour = "gray35"))

    ## Sort of colouring of points
    if (plot_x) { 
        if (typeof_x == "discrete") {
            plot_most_expressed <- .resolve_plot_colours(
                plot_most_expressed, df_exprs_by_cell_long$colour_by,
                col_by_variable)
#         plot_most_expressed <- plot_most_expressed +
#             ggthemes::scale_colour_tableau(name = col_by_variable)
        } else {
            plot_most_expressed <- plot_most_expressed +
                scale_colour_gradient(name = col_by_variable, low = "lightgoldenrod",
                                      high = "firebrick4", space = "Lab")
        }
    }

    plot_most_expressed + geom_point(
        aes_string(x = legend_val, y = "Feature", fill = "is_feature_control"),
        data = as.data.frame(rdata[chosen,]), colour = "gray30", shape = 21) +
        scale_fill_manual(values = c("aliceblue", "wheat")) +
        guides(fill = guide_legend(title = "Feature control?"))
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

################################################################################


#' Plot explanatory variables ordered by percentage of phenotypic variance explained
#'
#' @param object an SingleCellExperiment object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param method character scalar indicating the type of plot to produce. If
#' "density", the function produces a density plot of R-squared values for each
#' variable when fitted as the only explanatory variable in a linear model. If
#' "pairs", then the function produces a pairs plot of the explanatory variables
#' ordered by the percentage of feature expression variance (as measured by
#' R-squared in a marginal linear model) explained.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "logcounts" (default),
#' "tpm", "fpkm", "cpm", and "counts".
#' @param nvars_to_plot integer, the number of variables to plot in the pairs
#' plot. Default value is 10.
#' @param min_marginal_r2 numeric scalar giving the minimal value required for
#' median marginal R-squared for a variable to be plotted. Only variables with a
#' median marginal R-squared strictly larger than this value will be plotted.
#' @param variables optional character vector giving the variables to be plotted.
#' Default is \code{NULL}, in which case all variables in \code{colData(object)}
#' are considered and the \code{nvars_to_plot} variables with the highest median
#' marginal R-squared are plotted.
#' @param return_object logical, should an \code{SingleCellExperiment} object with median
#' marginal R-squared values added to \code{varMetadata(object)} be returned?
#' @param theme_size numeric scalar giving font size to use for the plotting
#' theme
#' @param ... parameters to be passed to \code{\link{pairs}}.
#'
#' @details If the \code{method} argument is "pairs", then the function produces
#' a pairs plot of the explanatory variables ordered by the percentage of
#' feature expression variance (as measured by R-squared in a marginal linear
#' model) explained by variable. Median percentage R-squared is reported on the
#' plot for each variable. Discrete variables are coerced to a factor and
#' plotted as integers with jittering. Variables with only one unique value are
#' quietly ignored.
#'
#' @return A ggplot object
#' @importFrom Biobase varMetadata<-
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
#' vars <- names(colData(example_sce))[c(2:3, 5:14)]
#' plotExplanatoryVariables(example_sce, variables=vars)
#'
plotExplanatoryVariables <- function(object, method = "density",
                                     exprs_values = "logcounts", nvars_to_plot = 10,
                                     min_marginal_r2 = 0, variables = NULL,
                                     return_object = FALSE, theme_size = 10,
                                     ...) {
    ## Check method argument
    method <- match.arg(method, c("density", "pairs"))
    ## Checking arguments for expression values
    # exprs_values <- match.arg(
    #     exprs_values,
    #     choices = c("logcounts", "norm_exprs", "stand_exprs", "norm_exprs",
    #                 "counts", "norm_counts", "tpm", "norm_tpm", "fpkm",
    #                 "norm_fpkm", "cpm", "norm_cpm"))
    exprs_mat <- assay(object, exprs_values)
    if ( is.null(exprs_mat) )
        stop("The supplied 'exprs_values' argument not found in assayData(object). Try 'exprs' or similar.")

    ## Check that variables are defined
    if ( is.null(variables) ) {
        variables_to_plot <- colnames(colData(object))
    } else {
        variables_to_plot <- NULL
        for (var in variables) {
            if ( !(var %in% colnames(colData(object))) ) {
                warning(paste("variable", var, "not found in colData(object).
                     Please make sure colData(object)[, variable] exists. This variable will not be plotted."))
            } else {
                variables_to_plot <- c(variables_to_plot, var)
            }
        }
    }

    ## Initialise matrix to store R^2 values for each feature for each variable
    rsquared_mat <- matrix(NA_real_, nrow = nrow(object),
                           ncol = length(variables_to_plot))
    val_to_plot_mat <- matrix(NA_real_, nrow = ncol(object),
                              ncol = length(variables_to_plot))
    colnames(rsquared_mat) <- colnames(val_to_plot_mat) <- variables_to_plot
    rownames(rsquared_mat) <- rownames(object)
    rownames(val_to_plot_mat) <- colnames(object)

    ## Get R^2 values for each feature and each variable
    for (var in variables_to_plot) {
        if ( var %in% variables_to_plot ) {
            if (length(unique(colData(object)[, var])) <= 1) {
                message(paste("The variable", var, "only has one unique value, so R^2 is not meaningful.
This variable will not be plotted."))
                rsquared_mat[, var] <- NA
            } else {
                x <- colData(object)[, var]
                #     x_na <- is.na(x)
                #     x <- x[!x_na]
                ## Determine type of variable
                typeof_x <- .getTypeOfVariable(object, var)
                if ( typeof_x == "discrete" ) {
                    x <- factor(x)
                    val_to_plot_mat[, var] <- jitter(as.integer(x))
                } else {
                    val_to_plot_mat[, var] <- x
                }
                design <- model.matrix(~x)
                rsquared_mat[, var] <- .getRSquared(exprs_mat, design, chunk=500)
#                 rsq_base <- apply(exprs_mat, 1, function(y) {
#                     lm.first <- lm(y ~ -1 + design); summary(lm.first)$r.squared})
#                 all(abs(rsq_base - rsquared_mat[, var]) < 0.000000000001)
            }
        }
    }

    ## Get median R^2 for each variable, add to labels and order by median R^2
    median_rsquared <- apply(rsquared_mat, 2, median, na.rm=TRUE)
    oo_median <- order(median_rsquared, decreasing = TRUE)
    nvars_to_plot <- min(sum(median_rsquared > min_marginal_r2, na.rm = TRUE),
                         nvars_to_plot)

    if ( method == "pairs" ) {
        if (nvars_to_plot == 1) 
            stop("Only one variable to plot, which does not make sense for a pairs plot.")
        ## Generate a larger data.frame for pairs plot
        df_to_expand <- val_to_plot_mat[, oo_median[1:nvars_to_plot], drop = FALSE]
        names(df_to_expand) <- colnames(df_to_expand)
        gg1 <- .makePairs(df_to_expand)
        diag_labs <-  paste0("Median R-sq = \n",
                             formatC(signif(100*median_rsquared, digits = 3),
                                     digits = 3, format = "fg", flag = "#"),
                             "%")[oo_median[1:nvars_to_plot]]
        centres <- apply(df_to_expand, 2,
                         function(x) {diff(range(x))/2 + min(x)})
        gg1$diags <- data.frame(xvar = colnames(df_to_expand),
                                yvar = colnames(df_to_expand),
                                x = centres, y = centres,
                                xmax = apply(df_to_expand, 2, max),
                                xmin = apply(df_to_expand, 2, min),
                                label = diag_labs)
        ## Plot these bad boys
        plot_out <- ggplot(gg1$all, aes_string(x = "x", y = "y")) +
            geom_point(fill = "gray60", colour = "gray40",
                       shape = 21, alpha = 0.65) +
            facet_grid(xvar ~ yvar, scales = "free") +
            geom_rect(aes_string(xmin = "xmin", ymin = "xmin", xmax = "xmax",
                                 ymax = "xmax"), colour = "white",
                      fill = "white", data = gg1$diags) +
            geom_text(aes_string(x = "x", y = "y", label = "label"),
                      size = theme_size / 3, data = gg1$diags) +
            xlab("") +
            ylab("") +
            theme_bw(theme_size) +
            theme(legend.position = "none")
    } else {
        df_to_plot <- suppressMessages(reshape2::melt(
            rsquared_mat[, oo_median[1:nvars_to_plot], drop = FALSE]))
        colnames(df_to_plot) <- c("Feature", "Expl_Var", "R_squared")
        df_to_plot$Pct_Var_Explained <- 100 * df_to_plot$R_squared
        df_to_plot$Expl_Var <- factor(
            df_to_plot$Expl_Var,
            levels = colnames(rsquared_mat)[oo_median[1:nvars_to_plot]])
        plot_out <- ggplot(df_to_plot, aes_string(x = "Pct_Var_Explained",
                                                  colour = "Expl_Var")) +
            geom_line(stat = "density", alpha = 0.7, size = 2, trim = TRUE) +
            geom_vline(xintercept = 1, linetype = 2) +
            scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
            xlab(paste0("% variance explained (log10-scale)")) +
            ylab("Density") +
            coord_cartesian(xlim = c(10 ^ (-3), 100))
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$Expl_Var, "")
        if ( requireNamespace("cowplot", quietly = TRUE) )
            plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
        else
            plot_out <- plot_out + theme_bw(theme_size)
    }

    if ( return_object ) {
        ## Return object so that marginal R^2 are added to varMetadata
        varMetadata(object) <- data.frame(
            labelDescription = paste("Median marginal R-squared =",
                                     median_rsquared))
        rdata <- rowData(object)
        rsq_out <- rsquared_mat[, oo_median[1:nvars_to_plot], drop = FALSE]
        colnames(rsq_out) <- paste0("Rsq_", colnames(rsq_out))
        rdata_new <- DataFrame(cbind(rdata, rsq_out))
        rowData(object) <- rdata_new
        print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

################################################################################

#' Produce QC diagnostic plots
#'
#' @param object an SingleCellExperiment object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param type character scalar providing type of QC plot to compute:
#' "highest-expression" (showing features with highest expression), "find-pcs" (showing
#' the most important principal components for a given variable),
#' "explanatory-variables" (showing a set of explanatory variables plotted
#' against each other, ordered by marginal variance explained), or
#' "exprs-mean-vs-freq" (plotting the mean expression levels against the
#' frequency of expression for a set of features).
#' @param ... arguments passed to \code{\link{plotHighestExprs}},
#' \code{\link{findImportantPCs}}, \code{\link{plotExplanatoryVariables}} and
#' \code{{plotExprsMeanVsFreq}} as appropriate.
#'
#' @details Display useful quality control plots to help with pre-processing
#' of data and identification of potentially problematic features and cells.
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#' example_sce <- calculateQCMetrics(example_sce)
#' plotQC(example_sce, type="high", col_by_variable="Mutation_Status")
#' plotQC(example_sce, type="find", variable="total_features")
#' vars <- names(colData(example_sce))[c(2:3, 5:14)]
#' plotQC(example_sce, type="expl", variables=vars)
#'
plotQC <- function(object, type = "highest-expression", ...) {
    type <- match.arg(type, c("highest-expression", "find-pcs",
                              "explanatory-variables", "exprs-freq-vs-mean"))
    if (type == "highest-expression") {
        plot_out <- plotHighestExprs(object, ...)
        return(plot_out)
    }
    if (type == "find-pcs") {
        plot_out <- findImportantPCs(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }
    if (type == "explanatory-variables") {
        plot_out <- plotExplanatoryVariables(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }
    if (type == "exprs-freq-vs-mean") {
        plot_out <- plotExprsFreqVsMean(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }

}


################################################################################

#' Plot a relative log expression (RLE) plot
#'
#' Produce a relative log expression (RLE) plot of one or more transformations of 
#' cell expression values.
#'
#' @param object an \code{SingleCellExperiment} object
#' @param exprs_mats named list of expression matrices. Entries can either be a 
#' character string, in which case the corresponding expression matrix will be 
#' extracted from the SingleCellExperiment \code{object}, or a matrix of expression values.
#' @param exprs_logged logical vector of same length as \code{exprs_mats} indicating
#' whether the corresponding entry in \code{exprs_mats} contains logged expression
#' values (\code{TRUE}) or not (\code{FALSE}).
#' @param colour_by character string defining the column of \code{colData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, 
#' a data frame with one column, containing values to map to colours for all cells.
#' @param style character(1), either \code{"minimal"} (default) or \code{"full"},
#' defining the boxplot style to use. \code{"minimal"} uses Tufte-style boxplots and
#' is fast for large numbers of cells. \code{"full"} uses the usual 
#' \code{\link{ggplot2}} and is more detailed and flexible, but can take a long 
#' time to plot for large datasets.
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternative is "none" (hide all legends).
#' @param order_by_colour logical, should cells be ordered (grouped) by the 
#' \code{colour_by} variable? Default is \code{TRUE}. Useful for visualising 
#' differences between batches or experimental conditions.
#' @param ncol integer, number of columns for the facetting of the plot. 
#' Default is 1.
#' @param ... further arguments passed to \code{\link[ggplot2]{geom_boxplot}}.
#'
#' @return a ggplot plot object
#'
#' @details 
#' Unwanted variation can be highly problematic and so its detection is often crucial.
#' Relative log expression (RLE) plots are a powerful tool for visualising such 
#' variation in high dimensional data. RLE plots are particularly useful for
#' assessing whether a procedure aimed at removing unwanted variation, i.e. a 
#' normalisation procedure, has been successful. These plots, while originally 
#' devised for gene expression data from microarrays, can also be used to reveal 
#' unwanted variation in single-cell expression data, where such variation can be 
#' problematic.
#' 
#' If style is "full", as usual with boxplots, the box shows the inter-quartile 
#' range and whiskers extend no more than 1.5 * IQR from the hinge (the 25th or 
#' 75th percentile). Data beyond the whiskers are called outliers and are plotted 
#' individually. The median (50th percentile) is shown with a white bar.
#' 
#' If style is "minimal", then median is shown with a circle, the IQR in a grey
#' line, and "whiskers" (as defined above) for the plots are shown with coloured 
#' lines. No outliers are shown for this plot style.
#'
#' @references 
#' Gandolfo LC, Speed TP. RLE Plots: Visualising Unwanted Variation in High Dimensional Data. 
#' arXiv [stat.ME]. 2017. Available: http://arxiv.org/abs/1704.03590
#'
#' @author 
#' Davis McCarthy
#'
#' @name plotRLE
#' @aliases plotRLE plotRLE,SingleCellExperiment-method
#' @export
#'
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colQuantiles rowMedians
#' 
#' @examples 
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#'  example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(logcounts(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' plotRLE(example_sce, list(logcounts= "logcounts", counts = "counts"), c(TRUE, FALSE), 
#'        colour_by = "Mutation_Status", style = "minimal")
#'
#' plotRLE(example_sce, list(logcounts = "logcounts", counts = "counts"), c(TRUE, FALSE), 
#'        colour_by = "Mutation_Status", style = "full",
#'        outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)
#' 
plotRLE <- function(object, exprs_mats = list(logcounts = "logcounts"), exprs_logged = c(TRUE),
                    colour_by = NULL, style = "minimal", legend = "auto", 
                    order_by_colour = TRUE, ncol = 1,  ...) {
    for (i in seq_len(length(exprs_mats))) {
        
        if (is.character(exprs_mats[[i]]) && exprs_mats[[i]] == "exprs") 
            exprs_mats[[i]] <- "logcounts"
    }
    .plotRLE(object, exprs_mats = exprs_mats, exprs_logged = exprs_logged,
             colour_by = colour_by, legend = legend, 
             order_by_colour = order_by_colour, ncol = ncol, style = style, ...)
}

.plotRLE <- function(object, exprs_mats = list(logcounts = "logcounts"), exprs_logged = c(TRUE),
                     colour_by = NULL, legend = "auto", order_by_colour = TRUE, ncol = 1,
                     style = "minimal", ...) {
    if (any(is.null(names(exprs_mats))) || any(names(exprs_mats) == ""))
        stop("exprs_mats must be a named list, with all names non-NULL and non-empty.")
    ## check legend argument
    legend <- match.arg(legend, c("auto", "none", "all"))
    style <- match.arg(style, c("full", "minimal"))
    ## Check arguments are valid
    colour_by_out <- .choose_vis_values(object, colour_by, mode = "column", search = "any", exprs_values = "logcounts")
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val
    ncells <- NULL
    ## calculate RLE
    rle_mats <- list()
    for (i in seq_along(exprs_mats)) {
        rle_mats[[i]] <- .calc_RLE(.get_exprs_for_RLE(object, exprs_mats[[i]]), 
                                   exprs_logged[i])
        names(rle_mats)[i] <- names(exprs_mats)[i]
        if (is.null(ncells))
            ncells <- ncol(rle_mats[[i]])
        else {
            if (ncol(rle_mats[[i]]) != ncells)
                stop(paste("Number of cells for", names(rle_mats)[i], "does not match other exprs matrices."))
        }
    }
    ## get into df for ggplot
    df_to_plot <- NULL
    if (order_by_colour) {
        oo <- order(colour_by_vals)
        colour_by_vals <- colour_by_vals[oo]
    }
    for (i in seq_along(rle_mats)) {
        tmp_df <- dplyr::as_data_frame(rle_mats[[i]])
        if (order_by_colour)
            tmp_df <- tmp_df[, oo]
        if (style == "full") {
            tmp_df[["source"]] <- names(rle_mats)[i]
            tmp_df <- reshape2::melt(tmp_df, id.vars = c("source"), value.name = "rle")
            tmp_df[[colour_by]] <- rep(colour_by_vals, each = nrow(rle_mats[[i]]))
            tmp_df[["x"]] <- rep(seq_len(ncells), each = nrow(rle_mats[[i]]))
        } else if (style == "minimal") {
            boxstats <- .rle_boxplot_stats(as.matrix(tmp_df))
            boxstats[[colour_by]] <- colour_by_vals
            boxstats[["x"]] <- seq_len(ncells)
            tmp_df <- boxstats
        } else
            stop("style argument must be either 'full' or 'minimal'.")
        tmp_df[["source"]] <- names(rle_mats)[i]
        if (is.null(df_to_plot)) {
            df_to_plot <- tmp_df
        } else {
            df_to_plot <- dplyr::bind_rows(df_to_plot, tmp_df)
        }
    }
    if (style == "full") {
        aesth <- aes_string(x = "x", group = "x", y = "rle", 
                            colour = colour_by, fill = colour_by)
        plot_out <- .plotRLE_full(df_to_plot, aesth, ncol, ...)
    } else if (style == "minimal") {
        plot_out <- .plotRLE_minimal(df_to_plot, colour_by, ncol)
    } 
    plot_out <- .resolve_plot_colours(plot_out, colour_by_vals, colour_by,
                                      fill = FALSE)
    plot_out <- .resolve_plot_colours(plot_out, colour_by_vals, colour_by,
                                      fill = TRUE)
    if ( legend == "none" )
        plot_out <- plot_out + theme(legend.position = "none")
    plot_out
}

.rle_boxplot_stats <- function(mat) {
    boxstats <- colQuantiles(DelayedArray(mat))
    colnames(boxstats) <- c("q0", "q25", "q50", "q75", "q100")
    boxdf <- dplyr::as_data_frame(boxstats)
    interqr <- boxstats[, 4] - boxstats[, 2]
    boxdf[["whiskMin"]] <- pmax(boxdf[["q0"]], 
                                boxdf[["q25"]] - 1.5 * interqr)
    boxdf[["whiskMax"]] <- pmin(boxdf[["q100"]], 
                                boxdf[["q75"]] + 1.5 * interqr)
    boxdf[["variable"]] <- colnames(mat)
    boxdf
}

.plotRLE_minimal <- function(df, colour_by, ncol, ...) {
    plot_out <- ggplot(df, aes_string(x = "x", fill = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "q75"), 
                     colour = "gray60") +
        geom_segment(aes_string(xend = "x", y = "q75", yend = "whiskMax", 
                                colour = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "whiskMin",
                                colour = colour_by)) +
        geom_point(aes_string(y = "q50"), shape = 21) +
        geom_hline(yintercept = 0, colour = "gray40", alpha = 0.5) +
        facet_wrap(~source, ncol = ncol) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.line.x = element_blank())
    plot_out
}


.plotRLE_full <- function(df, aesth, ncol, ...) {
    plot_out <- ggplot(df, aesth) +
        geom_boxplot(...) + # geom_boxplot(notch=T) to compare groups
        stat_summary(geom = "crossbar", width = 0.65, fatten = 0, color = "white", 
                     fun.data = function(x){ 
                         return(c(y = median(x), ymin = median(x), ymax = median(x))) }) +
        facet_wrap(~source, ncol = ncol) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.line.x = element_blank())
    plot_out
}

.get_exprs_for_RLE <- function(object, exprs_mat) {
    if (is.matrix(exprs_mat)) {
        return(exprs_mat)
    } else {
        if (is.character(exprs_mat))
            return(assay(object, exprs_mat))
        else
            stop("exprs_mat must be either a matrix of expression values or a character string giving the name of an expression data element of the SCESet object.")
    } 
}

.calc_RLE <- function(exprs_mat, logged = TRUE) {
    if (!logged)
        exprs_mat <- log2(exprs_mat + 1)
    features_meds <- rowMedians(DelayedArray(exprs_mat))
    med_devs <- as.matrix(exprs_mat - features_meds)
    med_devs
}

