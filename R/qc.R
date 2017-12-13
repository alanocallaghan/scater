## Convenience function for computing QC metrics and adding to pData & rowData
### This file contains definitions for the following functions:
### * calculateQCMetrics
### * findImportantPCs
### * plotExplanatoryVariables
### * plotHighestExprs
### * plotQC
###
### * .calculateSilhouetteWidth
### * .getRSquared
### * .getTypeOfVariable


################################################################################
#' Calculate QC metrics
#'
#' @param object an SingleCellExperiment object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param exprs_values character(1), indicating slot of the \code{assays} of the \code{object}
#' should be used to define expression? Valid options are "counts" [default; recommended],
#' "tpm", "fpkm" and "logcounts", or anything else in the object added manually by 
#' the user.
#' @param feature_controls a named list containing one or more vectors 
#' (character vector of feature names, logical vector, or a numeric vector of
#' indices are all acceptable) used to identify feature controls 
#' (for example, ERCC spike-in genes, mitochondrial genes, etc). 
#' @param cell_controls a character vector of cell (sample) names, or a logical
#' vector, or a numeric vector of indices used to identify cell controls (for
#' example, blank wells or bulk controls).
#' @param nmads numeric scalar giving the number of median absolute deviations 
#' to be used to flag potentially problematic cells based on total_counts (total
#' number of counts for the cell, or library size) and total_features (number of
#' features with non-zero expression). For total_features, cells are flagged for
#' filtering only if total_features is \code{nmads} below the median. Default 
#' value is 5.
#' @param pct_feature_controls_threshold numeric scalar giving a threshold for
#' percentage of expression values accounted for by feature controls. Used as to
#' flag cells that may be filtered based on high percentage of expression from
#' feature controls.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic features and cells. 
#' 
#' The following QC metrics are computed:
#' \describe{
#'  \item{total_counts:}{Total number of counts for the cell (aka ``library 
#'  size'')}
#'  \item{log10_total_counts:}{Total counts on the log10-scale}
#'  \item{total_features:}{The number of endogenous features (i.e. not control
#'  features) for the cell that have expression above the detection limit
#'   (default detection limit is zero)}
#'  \item{filter_on_depth:}{Would this cell be filtered out based on its
#'  log10-depth being (by default) more than 5 median absolute deviations from
#'  the median log10-depth for the dataset?}
#'  \item{filter_on_coverage:}{Would this cell be filtered out based on its
#'  coverage being (by default) more than 5 median absolute deviations from the
#'  median coverage for the dataset?}
#'  \item{filter_on_pct_counts_feature_controls:}{Should the cell be filtered
#'  out on the basis of having a high percentage of counts assigned to control
#'  features? Default threshold is 80 percent (i.e. cells with more than 80
#'  percent of counts assigned to feature controls are flagged).}
#'  \item{counts_feature_controls:}{Total number of counts for the cell
#'  that come from (one or more sets of user-defined) control features. Defaults
#'   to zero if no control features are indicated. If more than one set of
#'   feature controls are defined (for example, ERCC and MT genes are defined
#'   as controls), then this metric is produced for all sets, plus the union of
#'   all sets (so here, we get columns
#'   \code{counts_feature_controls_ERCC},
#'   \code{counts_feature_controls_MT} and
#'   \code{counts_feature_controls}).}
#'  \item{log10_counts_feature_controls:}{Just as above, the total
#'   number of counts from feature controls, but on the log10-scale. Defaults
#'   to zero (i.e.~log10(0 + 1), offset to avoid negative infinite values) if
#'   no feature control are indicated.}
#'  \item{pct_counts_feature_controls:}{Just as for the counts
#'   described above, but expressed as a percentage of the total counts.
#'   Defined for all control sets and their union, just like the raw counts.
#'   Defaults to zero if no feature controls are defined.}
#'  \item{filter_on_pct_counts_feature_controls:}{Would this cell be
#'   filtered out on the basis that the percentage of counts from feature
#'   controls is higher than a defined threhold (default is 80\%)? Just as with
#'   \code{counts_feature_controls}, this is defined for all control sets
#'   and their union.}
#'  \item{pct_counts_top_50_features:}{What percentage of the total counts is accounted for by the 50 highest-count features? Also computed for the top 100 and top 200 features, with the obvious changes to the column names. Note that the top ``X'' percentage will not be computed if the total number of genes is less than ``X''.}
#'  \item{pct_dropout:}{Percentage of features that are not ``detectably 
#'  expressed'', i.e. have expression below the \code{lowerDetectionLimit} 
#'  threshold.}
#'  \item{counts_endogenous_features:}{Total number of counts for the cell
#'   that come from endogenous features (i.e. not control features). Defaults
#'   to `depth` if no control features are indicated.}
#'  \item{log10_counts_endogenous_features:}{Total number of counts from
#'   endogenous features on the log10-scale. Defaults to all counts if no 
#'   control features are indicated.}
#'  \item{n_detected_feature_controls:}{Number of defined feature controls
#'    that have expression greater than the threshold defined in the object
#'    (that is, they are ``detectably expressed''; see
#'    \code{object@lowerDetectionLimit} to check the threshold). As with other
#'    metrics for feature controls, defined for all sets of feature controls
#'    (set names appended as above) and their union. So we might commonly get
#'    columns \code{n_detected_feature_controls_ERCC},
#'    \code{n_detected_feature_controls_MT} and
#'    \code{n_detected_feature_controls} (ERCC and MT genes detected).}
#'  \item{is_cell_control:}{Has the cell been defined as a cell control? If
#'    more than one set of cell controls are defined (for example, blanks and
#'    bulk libraries are defined as cell controls), then this metric is produced
#'     for all sets, plus the union of all sets (so we could typically get
#'     columns \code{is_cell_control_Blank},
#'     \code{is_cell_control_Bulk}, and \code{is_cell_control}, the latter
#'     including both blanks and bulks as cell controls).}
#' }
#' These cell-level QC metrics are added as columns to the ``phenotypeData''
#' slot of the \code{\link{SingleCellExperiment}} object so that they can be inspected and are
#' readily available for other functions to use. Furthermore, wherever
#' ``counts'' appear in the above metrics, the same metrics will also be
#' computed for ``exprs'', ``tpm'' and ``fpkm'' values (if TPM and FPKM values
#' are present in the \code{SingleCellExperiment} object), with the appropriate term
#' replacing ``counts'' in the name. The following feature-level QC metrics are
#' also computed:
#' \describe{
#' \item{mean_exprs:}{The mean expression level of the  gene/feature.}
#' \item{exprs_rank:}{The rank of the feature's mean expression level in the
#' cell.}
#' \item{n_cells_counts:}{The number of cells for which the expression level of
#' the feature is above the detection limit (default detection limit is zero).}
#' \item{total_feature_counts:}{The total number of counts assigned to that
#' feature across all cells.}
#' \item{log10_total_feature_counts:}{Total feature counts on the log10-scale.}
#' \item{pct_total_counts:}{The percentage of all counts that are accounted for
#' by the counts assigned to the feature.}
#' \item{pct_dropout:}{The percentage of all cells that have no detectable 
#' expression (i.e. \code{is_exprs(object)} is \code{FALSE}) for the feature.}
#' \item{is_feature_control:}{Is the feature a control feature? Default is
#' `FALSE` unless control features are defined by the user. If more than one
#' feature control set is defined (as above), then a column of this type is
#' produced for each control set (e.g. here, \code{is_feature_control_ERCC} and
#' \code{is_feature_control_MT}) as well as the column named
#' \code{is_feature_control}, which indicates if the feature belongs to any of
#' the control sets.}
#' }
#' These feature-level QC metrics are added as columns to the ``rowData''
#' slot of the \code{SingleCellExperiment} object so that they can be inspected and are
#' readily available for other functions to use. As with the cell-level metrics,
#'  wherever ``counts'' appear in the above, the same metrics will also be
#'  computed for ``exprs'', ``tpm'' and ``fpkm'' values (if TPM and FPKM values
#'  are present in the \code{SingleCellExperiment} object), with the appropriate term
#'  replacing ``counts'' in the name.
#'
#' @return an SingleCellExperiment object
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase sampleNames<- sampleNames assayDataElement assayDataElement<-
#' @importFrom stats cmdscale coef mad median model.matrix nls prcomp quantile var dist
#' @importFrom methods is new
#' @importFrom utils read.table
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames rowData rowData<- colData colData<-
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## with a set of feature controls defined
#' example_sce <- calculateQCMetrics(example_sce, 
#' feature_controls = list(set1 = 1:40))
#' 
#' ## with a named set of feature controls defined
#' example_sce <- calculateQCMetrics(example_sce, 
#'                                      feature_controls = list(ERCC = 1:40))
#' 
calculateQCMetrics <- function(object, exprs_values="counts", 
                               feature_controls = NULL, cell_controls = NULL, 
                               nmads = 5, pct_feature_controls_threshold = 80) {
    if ( !methods::is(object, "SingleCellExperiment"))
        stop("object must be a SingleCellExperiment")
    exprs_mat <- assay(object, i = exprs_values)
    if (exprs_values == "counts" || exprs_values == "cpm") { 
        linear <- TRUE
    } else { 
        linear <- FALSE
    }

    ## Adding general metrics for each cell.
    cd <- .get_qc_metrics_per_cell(exprs_mat, exprs_type = exprs_values,
            subset_row = NULL, subset_type = NULL, linear = TRUE)
    rd <- DataFrame(is_feature_control = logical(nrow(exprs_mat)), 
            row.names = rownames(exprs_mat))

    ## Adding metrics for the technical controls.
    n_feature_sets <- length(feature_controls)
    if (n_feature_sets) { 
        if (is.null(names(feature_controls))) {
            stop("feature_controls should be named")
        }
        
        # Converting to integer indices for all applications.
        reindexed <- lapply(feature_controls, FUN = .subset2index, 
                            target = exprs_mat)
        is_fcon <- Reduce(union, reindexed)
        rd$is_feature_control[is_fcon] <- TRUE

        # Adding feature controls.
        for (f in seq_len(n_feature_sets)) {
            cur.index <- logical(nrow(exprs_mat))
            cur.index[reindexed[[f]]] <- TRUE
            rd[[paste0("is_feature_control_", names(reindexed)[f])]] <- cur.index
        }

        # Running through all endogenous genes.
        is_endog <- which(!rd$is_feature_control)
        cd_endog <- .get_qc_metrics_per_cell(exprs_mat, exprs_type = exprs_values,
                subset_row = is_endog, subset_type = "endogenous", linear = linear)

        # Running through all feature controls.
        cd_fcon <- .get_qc_metrics_per_cell(exprs_mat, exprs_type = exprs_values,
                subset_row = is_fcon, subset_type = "feature_control", 
                linear = linear)

        # Running through each of the feature controls.
        cd_per_fcon <- vector("list", n_feature_sets)
        for (f in seq_len(n_feature_sets)) {
            cd_per_fcon[[f]] <- .get_qc_metrics_per_cell(
                exprs_mat, exprs_type = exprs_values,
                subset_row = reindexed[[f]], subset_type = names(reindexed)[f], 
                linear = linear)
        }

        cd <- do.call(cbind, c(list(cd, cd_endog, cd_fcon), cd_per_fcon))
    }
    
    ## Define cell controls
    ### Determine if vector or list
    rd_all <- .get_qc_metrics_per_gene(exprs_mat, exprs_type = exprs_values,
            subset_col = NULL, subset_type = NULL, linear = linear)
    rd <- cbind(rd, rd_all)
    cd$is_cell_control <- logical(ncol(exprs_mat))

    n_cell_sets <- length(cell_controls)
    if (n_cell_sets) {
        # Converting indices to integer.
        reindexed <- lapply(cell_controls, FUN = .subset2index, 
                            target = exprs_mat, byrow = FALSE)
        is_ccon <- Reduce(union, reindexed)
        cd$is_cell_control[is_ccon] <- TRUE
        
        # Adding sets to the colData.
        for (cx in seq_len(n_cell_sets)) {
            current_control <- logical(ncol(exprs_mat))
            current_control[reindexed[[cx]]] <- TRUE
            cd[[paste0("is_cell_control_", names(reindexed)[cx])]] <- current_control
        }

        # Adding statistics for non-control cells.
        is_noncon <- which(!cd$is_cell_control)
        rd_noncon <- .get_qc_metrics_per_gene(exprs_mat, exprs_type = exprs_values,
                subset_col = is_noncon, subset_type = "non_control", linear = linear)

        # Adding statistics for all control cells.
        rd_con <- .get_qc_metrics_per_gene(exprs_mat, exprs_type = exprs_values,
                subset_col = is_ccon, subset_type = "cell_control", linear = linear)

        # Adding statistics for each set of control cells.
        rd_collected <- vector("list", n_cell_sets)
        for (cx in seq_len(n_cell_sets)) {
            rd_current <- .get_qc_metrics_per_gene(exprs_mat, exprs_type = exprs_values,
                    subset_col = reindexed[[cx]], subset_type = names(reindexed)[cx], linear = linear)
            rd_collected[[cx]] <- rd_current
        }

        rd <- do.call(cbind, c(list(rd, rd_con, rd_noncon), rd_collected))
    }

    ### Remove columns to be replaced
    old_rd <- rowData(object)
    old_rd <- old_rd[, !(colnames(old_rd) %in% colnames(rd)), drop = FALSE]
    rd <- cbind(old_rd, rd)
    rowData(object) <- rd

    old_cd <- colData(object)
    old_cd <- old_cd[, !(colnames(old_cd) %in% colnames(cd)), drop = FALSE]
    cd <- cbind(old_cd, cd)
    colData(object) <- cd

    return(object)
}


.get_qc_metrics_per_cell <- function(exprs_mat, exprs_type = "counts",
        subset_row = NULL, subset_type = NULL, linear = TRUE) {
    ## Many thanks to Aaron Lun for suggesting efficiency improvements
    ## for this function.
    ## Get total expression from feature controls
    if (is.null(subset_type)) {
        subset_type <- ""
    } else {
        subset_type <- paste0("_", subset_type)
    }

    subset_row <- .subset2index(subset_row, target = exprs_mat, byrow = TRUE)
    nfeatures <- nexprs(exprs_mat, rows = subset_row, byrow = FALSE)
    rd <- DataFrame(nfeatures, log10(nfeatures + 1), row.names = colnames(exprs_mat))
    colnames(rd) <- paste0(c("", "log10_"), "total_features", subset_type)

    if (linear) {
        ## Adding the total sum.
        libsize <- .colSums(exprs_mat, rows = subset_row)
        rd[[paste0("total_", exprs_type, subset_type)]] <- libsize
        rd[[paste0("log10_total_", exprs_type, subset_type)]] <- log10(libsize + 1)

        if (!is.null(subset_row)) {
            ## Computing percentages of actual total.
            rd[[paste0("pct_", exprs_type, subset_type)]] <- 
                (100 * libsize / .colSums(exprs_mat))
        }
        
        ## Computing total percentages.
        pct_top <- .calc_top_prop(exprs_mat, subset_row = subset_row, 
                subset_type = subset_type, exprs_type = exprs_type)
        rd <- cbind(rd, pct_top)
    }

    return(rd)
}


.calc_top_prop <- function(exprs_mat, top.number = c(50L, 100L, 200L, 500L),
        subset_row = NULL, subset_type="", exprs_type="features") {
    ## Calculate the proportion of expression belonging to the top set of genes.
    ## Produces a matrix of proportions for each top number.
    
    if (is.null(subset_row)) { 
        total_nrows <- nrow(exprs_mat)
    } else {
        subset_row <- .subset2index(subset_row, exprs_mat) 
        subset_row <- subset_row - 1L # zero indexing needed for this C++ code.
        total_nrows <- length(subset_row)
    }

    can.calculate <- top.number <= total_nrows
    if (any(can.calculate)) { 
        top.number <- top.number[can.calculate]
        pct_exprs_top_out <- .Call(cxx_calc_top_features, exprs_mat, top.number, subset_row)
        names(pct_exprs_top_out) <- paste0("pct_", exprs_type, "_top_", top.number, "_features", subset_type)
        return(do.call(data.frame, pct_exprs_top_out))
    }
    return(data.frame(row.names = seq_len(ncol(exprs_mat))))
}


.get_qc_metrics_per_gene <- function(exprs_mat, exprs_type="counts", 
        subset_col=NULL, subset_type=NULL, linear=TRUE) {

    if (is.null(subset_type)) {
        subset_type <- ""
    } else {
        subset_type <- paste0("_", subset_type)
    }
    if (is.null(subset_col)) {
        total.cells <- ncol(exprs_mat)
    } else {
        subset_col <- .subset2index(subset_col, target = exprs_mat, byrow = FALSE) 
        total.cells <- length(subset_col)
    }

    sum_exprs <- .rowSums(exprs_mat, col = subset_col)
    ave <- sum_exprs/total.cells
    fd <- DataFrame(ave, log10(ave + 1), rank(ave), row.names = rownames(exprs_mat))
    colnames(fd) <- paste0(c("mean", "log10_mean", "rank"), "_", exprs_type, subset_type)

    ncells.exprs <- nexprs(exprs_mat, col = subset_col, byrow = TRUE)
    fd[[paste0("n_cells_", exprs_type, subset_type)]] <- ncells.exprs
    fd[[paste0("pct_dropout_", exprs_type, subset_type)]] <- 100 * (1 - ncells.exprs/total.cells)

    if (linear) {
        fd[[paste0("total_", exprs_type, subset_type)]] <- sum_exprs
        fd[[paste0("log10_total_", exprs_type, subset_type)]] <- log10(sum_exprs + 1)

        if (!is.null(subset_col)) { 
            total_exprs <- .rowSums(exprs_mat)
            fd[[paste0("pct_", exprs_type, subset_type)]] <- sum_exprs/total_exprs * 100
        }
    }

    return(fd) 
}

################################################################################

#' Identify if a cell is an outlier based on a metric
#' 
#' @param metric numeric or integer vector of values for a metric
#' @param nmads scalar, number of median-absolute-deviations away from median
#' required for a value to be called an outlier
#' @param type character scalar, choice indicate whether outliers should be 
#' looked for at both tails (default: "both") or only at the lower end ("lower") 
#' or the higher end ("higher")
#' @param log logical, should the values of the metric be transformed to the 
#' log10 scale before computing median-absolute-deviation for outlier detection?
#' @param subset logical or integer vector, which subset of values should be 
#' used to calculate the median/MAD? If \code{NULL}, all values are used.
#' Missing values will trigger a warning and will be automatically ignored. 
#' @param batch factor of length equal to \code{metric}, specifying the batch
#' to which each observation belongs. A median/MAD is calculated for each batch,
#' and outliers are then identified within each batch.
#' @param min.diff numeric scalar indicating the minimum difference from the 
#' median to consider as an outlier. The outlier threshold is defined from the 
#' larger of \code{nmads} MADs and \code{min.diff}, to avoid calling many 
#' outliers when the MAD is very small. If \code{NA}, it is ignored.
#' 
#' @description Convenience function to determine which values for a metric are
#' outliers based on median-absolute-deviation (MAD).
#' 
#' @return a logical vector of the same length as the \code{metric} argument
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## with a set of feature controls defined
#' example_sce <- calculateQCMetrics(example_sce, 
#' feature_controls = list(set1 = 1:40))
#' isOutlier(example_sce$total_counts, nmads = 3)
#' 
isOutlier <- function(metric, nmads = 5, type = c("both", "lower", "higher"), 
                      log = FALSE, subset = NULL, batch = NULL, min.diff = NA) {
    if (log) {
        metric <- log10(metric)
    }
    if (any(is.na(metric))) { 
        warning("missing values ignored during outlier detection")
    }

    if (!is.null(batch)) {
        N <- length(metric)
        if (length(batch) != N) { 
            stop("length of 'batch' must equal length of 'metric'")
        }

        # Coercing non-NULL subset into a logical vector.
        if (!is.null(subset)) { 
            new.subset <- logical(N)
            names(new.subset) <- names(metric)
            new.subset[subset] <- TRUE
            subset <- new.subset
        }
   
        # Computing QC metrics for each batch. 
        by.batch <- split(seq_len(N), batch)
        collected <- logical(N)
        for (b in by.batch) {
            collected[b] <- Recall(metric[b], nmads = nmads, type = type,
                                   log = FALSE, subset = subset[b], 
                                   batch = NULL, min.diff = min.diff)
        }
        return(collected)
    }

    # Computing median/MAD (possibly based on subset of the data).
    if (!is.null(subset)) {
        submetric <- metric[subset]
        if (length(submetric) == 0L) {
            warning("no observations remaining after subsetting")
        }
    } else {
        submetric <- metric
    }
    cur.med <- median(submetric, na.rm = TRUE)
    cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min.diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val 
    lower.limit <- cur.med - diff.val 
    
    type <- match.arg(type)
    if (type == "lower") {
        upper.limit <- Inf
    } else if (type == "higher") {
        lower.limit <- -Inf
    }

    return(metric < lower.limit | upper.limit < metric)
}

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
        rv <- matrixStats::rowVars(df_for_pca)
        feature_set <-
            order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    }
    df_for_pca <- df_for_pca[feature_set,]
    df_for_pca <- t(df_for_pca)

    ## Drop any features with zero variance
    keep_feature <- (matrixStats::colVars(df_for_pca) > 0.001)
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
    effects <- qr.qty(QR, t(y))
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
### Plot expression frequency vs mean for feature controls

#' Plot frequency of expression against mean expression level
#'
#' @param object an \code{SingleCellExperiment} object.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to plot. If character, entries must all be in
#' \code{rownames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}. If \code{NULL}, then the function checks if
#' feature controls are defined. If so, then only feature controls are plotted,
#' if not, then all features are plotted.
#' @param feature_controls character, numeric or logical vector indicating a set of
#' features to be used as feature controls for computing technical dropout
#' effects. If character, entries must all be in \code{rownames(object)}. If
#' numeric, values are taken to be indices for features. If logical, vector is
#' used to index features and should have length equal to \code{nrow(object)}.
#' If \code{NULL}, then the function checks if feature controls are defined. If
#' so, then these feature controls are used.
#' @param shape (optional) numeric scalar to define the plotting shape.
#' @param alpha (optional) numeric scalar (in the interval 0 to 1) to define the
#'  alpha level (transparency) of plotted points.
#' @param show_smooth logical, should a smoothed fit through feature controls
#' (if available; all features if not) be shown on the plot? Lowess used if a 
#' small number of feature controls. For details see 
#' \code{\link[ggplot2]{geom_smooth}}.
#' @param se logical, should standard error (confidence interval) be shown for 
#' smoothed fit?
#' @param ... further arguments passed to \code{\link{plotMetadata}} (should
#' only be \code{size}, if anythin).
#'
#' @details This function plots gene expression frequency versus mean
#' expression level, which can be useful to assess the effects of technical
#' dropout in the dataset. We fit a non-linear least squares curve for the
#' relationship between expression frequency and mean expression and use this to
#' define the number of genes above high technical dropout and the numbers of
#' genes that are expressed in at least 50% and at least 25% of cells. A subset
#' of genes to be treated as feature controls can be specified, otherwise any
#' feature controls previously defined are used.
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
#' example_sce <- calculateQCMetrics(example_sce, 
#'                                  feature_controls = list(set1 = 1:500))
#' plotExprsFreqVsMean(example_sce)
#'
#' example_sce <- calculateQCMetrics(
#' example_sce, feature_controls = list(controls1 = 1:20, 
#'                                       controls2 = 500:1000),
#'                                       cell_controls = list(set_1 = 1:5, 
#'                                       set_2 = 31:40))
#' plotExprsFreqVsMean(example_sce)
#'
plotExprsFreqVsMean <- function(object, feature_set = NULL,
                                feature_controls = NULL, shape = 1, alpha = 0.7,
                                show_smooth = TRUE, se = TRUE, ...) {
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be an SingleCellExperiment")
    if ( is.null(rowData(object)$n_cells_counts) ) {
        stop("rowData(object) does not have a 'n_cells_counts' column. Try running 'calculateQCMetrics' on this object (ensuring that the object contains counts data) and then rerun this command.")
    }
    if ( !is.null(feature_set) && feature_set != "feature_controls" &&
         typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% rownames(object))) )
            stop("when the argument 'feature_set' is of type character and not 'feature_controls', all features must be in rownames(object)")
    }
    if ( is.null(feature_set) ) {
        feature_set_logical <- rep(TRUE, nrow(object))
        x_lab <- "Mean normalised counts (all features; log2-scale)"
    } else {
        if ( length(feature_set) == 1 && feature_set == "feature_controls" ) {
            feature_set_logical <- rowData(object)$is_feature_control
            x_lab <- "Mean normalised counts (feature controls; log2-scale)"
        } else {
            if ( is.character(feature_set) )
                feature_set_logical <- rownames(object) %in% feature_set
            else {
                feature_set_logical <- rep(FALSE, nrow(object))
                if ( is.numeric(feature_set) )
                    feature_set_logical[feature_set] <- TRUE
                else
                    feature_set_logical <- feature_set
            }
               x_lab <- "Mean normalised counts (supplied feature set; log2-scale)"
        }
    }
    ## check that feature controls, if defined, are sensible
    if (!is.null(feature_controls) && typeof(feature_controls) == "character") {
        if ( !(all(feature_controls %in% rownames(object))) )
            stop("when the argument 'feature_controls' is of type character all features must be in rownames(object)")
    }
    if ( is.null(feature_controls) )
        feature_controls_logical <- rowData(object)$is_feature_control
    else {
        if ( is.character(feature_controls) )
            feature_controls_logical <- (rownames(object) %in%
                                             feature_controls)
        else {
            feature_controls_logical <- rep(FALSE, nrow(object))
            if ( is.numeric(feature_controls) )
                feature_controls_logical[feature_controls] <- TRUE
            else
                feature_controls_logical <- feature_controls
        }
        rowData(object)$is_feature_control <- feature_controls_logical
    }

    ## define percentage of cells expressing a gene
    rowData(object)$pct_cells_exprs <- (100 * rowData(object)$n_cells_counts /
                                          ncol(object))
    y_lab <- paste0("Frequency of expression (% of ", ncol(object), " cells)")

    rowData(object)$mean_exprs <- log2(calcAverage(object) + 1)
    ## Plot this
    if ( any(rowData(object)$is_feature_control[feature_set_logical]) &&
         !all(rowData(object)$is_feature_control[feature_set_logical]) ) {
        plot_out <- plotMetadata(
            as.data.frame(rowData(object)[feature_set_logical,]),
            aesth = aes_string(x = "mean_exprs",
                               y = "pct_cells_exprs",
                               colour = "is_feature_control",
                               shape = "is_feature_control"),
            alpha = alpha, ...) +
            scale_shape_manual(values = c(1, 17)) +
            ylab(y_lab) +
            xlab(x_lab)
    } else {
        plot_out <- plotMetadata(
            as.data.frame(rowData(object)[feature_set_logical,]),
            aesth = aes_string(x = "mean_exprs",
                               y = "pct_cells_exprs"),
            alpha = alpha, shape = shape, ...) +
            ylab(y_lab) +
            xlab(x_lab)
    }
    
    ## data frame with expression mean and frequency for feature controls
    mn_vs_fq <- data.frame(
        mn = rowData(object)$mean_exprs,
        fq = rowData(object)$pct_cells_exprs / 100,
        is_feature_control = feature_controls_logical)
    text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))
    
    if ( show_smooth ) {
        if ( any(feature_controls_logical) )
            plot_out <- plot_out +
                geom_smooth(aes_string(x = "mn", y = "100 * fq"), 
                            data = mn_vs_fq[feature_controls_logical,], 
                            colour = "firebrick", size = 1, se = se)
        else
            plot_out <- plot_out +
                geom_smooth(aes_string(x = "mn", y = "100 * fq"), data = mn_vs_fq, 
                            colour = "firebrick", size = 1, se = se)
    }        
    
    ## estimate 50% spike-in dropout
    if ( any(feature_controls_logical) ) {
        dropout <- nls(fq ~ (1 / (1 + exp(-(-i + 1 * mn)))),
                       start = list(i = 5),
                       data = mn_vs_fq[feature_controls_logical,])
        # nl_fit_df <- data.frame(
        #     x = quantile(mn_vs_fq$mn, probs = seq(0.01, 0.999, by = 0.001)))
        # nl_fit_df$y <- 100 * (1/(1 + exp(-(-coef(dropout) + 1 * nl_fit_df$x))))
        # nl_fit_df$is_feature_control <- FALSE
        ## annotate plot
        plot_out <- plot_out + 
            geom_vline(xintercept = coef(dropout), linetype = 2) +
            annotate("text", x = text_x_loc, y = 60, label = paste(
                sum(mn_vs_fq[!feature_controls_logical, "mn"] > coef(dropout)),
                " genes with mean expression\nhigher than value for 50% dropout of feature controls", sep = ""))
    }
    
    ## add annotations to existing plot
    plot_out <- plot_out +
        geom_hline(yintercept = 50, linetype = 2) + # 50% dropout
        annotate("text", x = text_x_loc, y = 40, label =  paste(
            sum(mn_vs_fq$fq >= 0.5),
            " genes are expressed\nin at least 50% of cells", sep = "" )) +
        annotate("text", x = text_x_loc, y = 20, label =  paste(
            sum(mn_vs_fq$fq >= 0.25),
            " genes are expressed\nin at least 25% of cells", sep = "" ))
    
    ## return the plot object
    plot_out
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
    colour_by_out <- .choose_vis_values(object, colour_by, cell_control_default = TRUE,
                                        check_features = TRUE, exprs_values = "logcounts")
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
    boxstats <- matrixStats::colQuantiles(mat)
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
    features_meds <- matrixStats::rowMedians(exprs_mat)
    med_devs <- exprs_mat - features_meds
    med_devs
}

