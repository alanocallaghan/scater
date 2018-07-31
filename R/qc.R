## Convenience function for computing QC metrics and adding to colData & rowData
### This file contains definitions for the following functions:
### * findImportantPCs
### * plotExplanatoryVariables
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
            xlab("% variance explained") +
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
