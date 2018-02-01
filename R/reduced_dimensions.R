#' Run PCA for a SingleCellExperiment object
#'
#' Perform a principal components analysis (PCA) based on the data stored in 
#' a \code{\link{SingleCellExperiment}} object. 
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param ntop numeric scalar indicating the number of most variable features to use for the PCA. 
#' Default is 500, but any \code{ntop} argument is overridden if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal components to obtain from \code{\link{prcomp}}.
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
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
#' @return A \code{SingleCellExperiment} object containing the first 
#' \code{ncomponent} principal coordinates for each cell in the \code{"PCA"}
#' entry of the \code{reducedDims} slot.
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
            rv <- .rowVars(exprs_mat)
            o <- order(rv, decreasing = TRUE)
            feature_set <- o[seq_len(min(ntop, length(rv)))]
        }

        # Subsetting to the desired features (do NOT move below 'scale()')
        exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]
        ## Standardise expression if scale_features argument is TRUE
        exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
    }

    ## Drop any features with zero variance
    keep_feature <- .colVars(exprs_to_plot) > 0.001
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
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
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
#' @return A \code{SingleCellExperiment} object containing the coordinates of
#' the first \code{ncomponent} t-SNE dimensions for each cell in the \code{"TSNE"}
#' entry of the \code{reducedDims} slot.
#'
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to
#' compute the t-SNE. Note that the algorithm is not deterministic, so different
#' runs of the function will produce differing plots (see \code{\link{set.seed}}
#' to set a random seed for replicable results). The value of the
#' \code{perplexity} parameter can have a large effect on the resulting plot, so
#' it can often be worthwhile to try multiple values to find the most appealing
#' visualisation and to ensure that the conclusions are robust.
#'
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International Conference on Learning Representations, 2013.
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
            rv <- .rowVars(exprs_mat)
            ntop <- min(ntop, length(rv))
            feature_set <- order(rv, decreasing = TRUE)[seq_len(ntop)]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat[feature_set,,drop = FALSE]
        keep_feature <- .rowVars(vals) > 0.001
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
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
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
#' @return A \code{SingleCellExperiment} object containing the coordinates of the first 
#' \code{ncomponent} diffusion map components for each cell in the \code{"DiffusionMap"} 
#' entry of the \code{reducedDims} slot.
#'
#' @export
#' @rdname runDiffusionMap
#' @seealso 
#' \code{\link[destiny]{destiny}},
#' \code{\link[scater]{plotDiffusionMap}}
#'
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
            rv <- .rowVars(exprs_mat)
            feature_set <-
                order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat
        vals <- vals[feature_set,,drop = FALSE]
        keep_feature <- .rowVars(vals) > 0.001
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
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
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
#' @return A \code{SingleCellExperiment} object containing the coordinates of the 
#' first \code{ncomponent} MDS dimensions for each cell in the #' \code{"MDS"} 
#' entry of the \code{reducedDims} slot.
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
            rv <- .rowVars(exprs_mat)
            feature_set <-
                order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat
        vals <- vals[feature_set,,drop = FALSE]
        keep_feature <- .rowVars(vals) > 0.001
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


