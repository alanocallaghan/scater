## Convenience function for computing QC metrics and adding to pData & fData
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
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param feature_controls a character vector of feature names, or a logical 
#' vector, or a numeric vector of indices used to identify feature controls 
#' (for example, ERCC spike-in genes, mitochondrial genes, etc).
#' @param cell_controls a character vector of cell (sample) names, or a logical
#' vector, or a numeric vector of indices used to identify cell controls (for
#' example, blank wells or bulk controls).
#' @param nmads numeric scalar giving the number of median absolute deviations 
#' to be used to flag potentially problematic cells based on total_counts (total
#' number of counts for the cell, or library size) and total_features (number of
#' features with non-zero expression). Default value is 5.
#' @param pct_feature_controls_threshold numeric scalar giving a threshold for
#' percentage of expression values accounted for by feature controls. Used as to
#' flag cells that may be filtered based on high percentage of expression from
#' feature controls.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic features and cells. The
#'  following QC metrics are computed:
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
#'  \item{pct_counts_top_50_features:}{What percentage of the total counts is accounted for by the 50 highest-count features? Also computed for the top 100 and top 200 features, with the obvious changes to the column names.}
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
#' slot of the \code{SCESet} object so that they can be inspected and are
#' readily available for other functions to use. Furthermore, wherever
#' ``counts'' appear in the above metrics, the same metrics will also be
#' computed for ``exprs'', ``tpm'' and ``fpkm'' values (if TPM and FPKM values
#' are present in the \code{SCESet} object), with the appropriate term
#' replacing ``counts'' in the name. The following feature-level QC metrics are
#' also computed:
#' \describe{
#' \item{mean_exprs:}{The mean expression level of the  gene/feature.}
#' \item{exprs_rank:}{The rank of the feature's mean expression level in the
#' cell.}
#' \item{n_cells_exprs:}{The number of cells for which the expression level of
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
#' These feature-level QC metrics are added as columns to the ``featureData''
#' slot of the \code{SCESet} object so that they can be inspected and are
#' readily available for other functions to use. As with the cell-level metrics,
#'  wherever ``counts'' appear in the above, the same metrics will also be
#'  computed for ``exprs'', ``tpm'' and ``fpkm'' values (if TPM and FPKM values
#'  are present in the \code{SCESet} object), with the appropriate term
#'  replacing ``counts'' in the name.
#'
#' @return an SCESet object
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase sampleNames<-
#' @importFrom matrixStats colCumsums
#' @importFrom stats cmdscale coef mad median model.matrix nls prcomp quantile var
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#'
#' ## with a set of feature controls define
#' example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:40)
calculateQCMetrics <- function(object, feature_controls = NULL,
                               cell_controls = NULL, nmads = 5,
                               pct_feature_controls_threshold = 80) {
    ## We must have an SCESet object
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet object.")
    ## the object must have some samples
    if ( ncol(object) < 1 )
        stop("object must have at least one sample (column)")
    if ( nrow(object) < 1 )
        stop("object must have at least one feature (row)")
    ## Compute cell-level metrics
    if ( is.null(is_exprs(object)) ) {
        if (is.null(counts(object))) {
            stop("Please define is_exprs(object).
                 E.g. use is_exprs(object) <- exprs(object) > 0.1")
        } else {
            warning("is_exprs(object) is null.
Defining 'is_exprs' using count data and
a lower count threshold of 0.")
            isexprs <- calcIsExprs(object, lowerDetectionLimit = 0)
            rownames(isexprs) <- rownames(counts(object))
            colnames(isexprs) <- colnames(counts(object))
            is_exprs(object) <- isexprs
        }
    }

    ## See what versions of the expression data are available in the object
    exprs_mat <- exprs(object)
    if ( object@logged )
        exprs_mat <- 2 ^ exprs_mat - object@logExprsOffset
    counts_mat <- counts(object)
    tpm_mat <- tpm(object)
    fpkm_mat <- fpkm(object)

    ## Contributions from control features
    ### Determine if vector or list
    if ( is.null(feature_controls) | length(feature_controls) == 0 ) {
        exprs_feature_controls <- pct_exprs_feature_controls <-
            rep(0, ncol(object))
        feature_controls_pdata <- data.frame(exprs_feature_controls,
                                          pct_exprs_feature_controls)
        if ( !is.null(counts_mat) ) {
            counts_feature_controls <- pct_counts_feature_controls <-
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata,
                                         data.frame(counts_feature_controls,
                                              pct_counts_feature_controls))
        }
        if ( !is.null(tpm_mat) ) {
            tpm_feature_controls <- pct_tpm_feature_controls <-
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata,
                                         data.frame(tpm_feature_controls,
                                                    pct_tpm_feature_controls))
        }
        if ( !is.null(fpkm_mat) ) {
            fpkm_feature_controls <- pct_fpkm_feature_controls <-
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata,
                                         data.frame(fpkm_feature_controls,
                                                    pct_fpkm_feature_controls))
        }
        is_feature_control <- rep(FALSE, nrow(object))
        feature_controls_fdata <- data.frame(is_feature_control)
        n_sets_feature_controls <- 1
    } else {
        if ( is.list(feature_controls) ) {
            feature_controls_list <- feature_controls
            n_sets_feature_controls <- length(feature_controls)
        }
        else {
            feature_controls_list <- list(feature_controls)
            n_sets_feature_controls <- 1
        }
        ## Cycle through the feature_controls list and add QC info
        for (i in seq_len(length(feature_controls_list)) ) {
            gc_set <- feature_controls_list[[i]]
            set_name <- names(feature_controls_list)[i]
            if ( is.logical(gc_set) ) {
                is_feature_control <- gc_set
                gc_set <- which(gc_set)
            } else {
                is_feature_control <- rep(FALSE, nrow(object))
            }
            if (is.character(gc_set))
                gc_set <- which(rownames(object) %in% gc_set)
            df_pdata_this <- .get_qc_metrics_exprs_mat(
                exprs_mat, gc_set, pct_feature_controls_threshold,
                calc_top_features = (n_sets_feature_controls == 1),
                exprs_type = "exprs")
            if ( !is.null(counts_mat) ) {
                df_pdata_counts <- .get_qc_metrics_exprs_mat(
                    counts_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1),
                    exprs_type = "counts")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_counts)
            }
            if ( !is.null(tpm_mat) ) {
                df_pdata_tpm <- .get_qc_metrics_exprs_mat(
                    tpm_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1),
                    exprs_type = "tpm")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_tpm)
            }
            if ( !is.null(fpkm_mat) ) {
                df_pdata_fpkm <- .get_qc_metrics_exprs_mat(
                    fpkm_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1),
                    exprs_type = "fpkm")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_fpkm)
            }
            is_feature_control[gc_set] <- TRUE
            ## Define number of feature controls expressed
            n_detected_feature_controls <- colSums(is_exprs(object)[gc_set,])
            df_pdata_this$n_detected_feature_controls <-
                n_detected_feature_controls
            if ( n_sets_feature_controls > 1 )
                colnames(df_pdata_this) <- paste(colnames(df_pdata_this),
                                                 set_name, sep = "_")
            if ( exists("feature_controls_pdata") )
                feature_controls_pdata <- cbind(feature_controls_pdata,
                                                df_pdata_this)
            else
                feature_controls_pdata <- df_pdata_this
            ## Construct data.frame for fData from this feature control set
            df_fdata_this <- data.frame(is_feature_control)
            colnames(df_fdata_this) <- paste(colnames(df_fdata_this), set_name,
                                          sep = "_")
            if ( exists("feature_controls_fdata") )
                feature_controls_fdata <- cbind(feature_controls_fdata,
                                                df_fdata_this)
            else
                feature_controls_fdata <- df_fdata_this
        }
    }

    ## Fix column names and define feature controls across all control sets
    if ( n_sets_feature_controls == 1 ) {
        colnames(feature_controls_pdata) <- gsub("_$", "",
                                              colnames(feature_controls_pdata))
        colnames(feature_controls_fdata) <- gsub("_$", "",
                                              colnames(feature_controls_fdata))
    } else {
        ## Combine all feature controls
        is_feature_control <- apply(feature_controls_fdata, 1, any)
        feature_controls_fdata <- cbind(feature_controls_fdata,
                                        is_feature_control)
        ## Compute metrics using all feature controls
        df_pdata_this <- .get_qc_metrics_exprs_mat(
            exprs_mat, is_feature_control, pct_feature_controls_threshold,
            calc_top_features = TRUE, exprs_type = "exprs")
#         exprs_feature_controls <- colSums(exprs_mat[is_feature_control,])
#         pct_exprs_feature_controls <- (100 * exprs_feature_controls /
#                                              colSums(exprs_mat))
#         filter_on_pct_exprs_feature_controls <-
#             (pct_exprs_feature_controls > pct_feature_controls_threshold)
#         df_pdata_this <- data.frame(exprs_feature_controls,
#                                     pct_exprs_feature_controls,
#                                     filter_on_pct_exprs_feature_controls)
        if ( !is.null(counts_mat) ) {
            df_pdata_counts <- .get_qc_metrics_exprs_mat(
                counts_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "counts")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_counts)
        }
        if ( !is.null(tpm_mat) ) {
            df_pdata_tpm <- .get_qc_metrics_exprs_mat(
                tpm_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "tpm")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_tpm)
        }
        if ( !is.null(fpkm_mat) ) {
            df_pdata_fpkm <- .get_qc_metrics_exprs_mat(
                fpkm_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "fpkm")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_fpkm)
        }
        n_detected_feature_controls <- colSums(
            is_exprs(object)[is_feature_control,])
        df_pdata_this$n_detected_feature_controls <- n_detected_feature_controls
        feature_controls_pdata <- cbind(feature_controls_pdata, df_pdata_this)
    }

    ## Compute total_features and find outliers
    total_features <- colSums(is_exprs(object)[!is_feature_control,])
    mad_total_features <- mad(total_features)
    med_total_features <- median(total_features)
    keep_total_features <- c(med_total_features - nmads * mad_total_features,
                             med_total_features + nmads * mad_total_features)
    filter_on_total_features <- (findInterval(total_features,
                                              keep_total_features) != 1)
    ## Compute total_counts if counts are present
    if ( !is.null(counts_mat) )
        total_counts <- colSums(counts_mat)
    else
        total_counts <- colSums(exprs_mat)
    mad_total_counts <- mad(log10(total_counts))
    med_total_counts <- median(log10(total_counts))
    keep_total_counts <- c(med_total_counts - nmads * mad_total_counts,
                           med_total_counts + nmads * mad_total_counts)
    filter_on_total_counts <- (findInterval(log10(total_counts),
                                            keep_total_counts) != 1)

    ## Define counts from endogenous features
    qc_pdata <- feature_controls_pdata
    qc_pdata$exprs_endogenous_features <- colSums(exprs_mat) -
        feature_controls_pdata$exprs_feature_controls
    if ( !is.null(counts_mat) ) {
        qc_pdata$counts_endogenous_features <- total_counts -
            feature_controls_pdata$counts_feature_controls
    }
    if ( !is.null(tpm_mat) ) {
        qc_pdata$tpm_endogenous_features <- colSums(tpm_mat) -
            feature_controls_pdata$tpm_feature_controls
    }
    if ( !is.null(fpkm_mat) ) {
        qc_pdata$fpkm_endogenous_features <- colSums(fpkm_mat) -
            feature_controls_pdata$fpkm_feature_controls
    }
    ## Define log10 read counts from feature controls
    cols_to_log <- grep("^counts_|^exprs_|^tpm_|^fpkm_",
                        colnames(qc_pdata))
    log10_cols <- log10(qc_pdata[, cols_to_log, drop = FALSE] + 1)
    colnames(log10_cols) <- paste0("log10_", colnames(qc_pdata)[cols_to_log])
    ## Combine into a big pdata object
    qc_pdata <- cbind(qc_pdata, log10_cols)

    ## Define cell controls
    ### Determine if vector or list
    if ( is.null(cell_controls) | length(cell_controls) == 0 ) {
        is_cell_control <- rep(FALSE, ncol(object))
        cell_controls_pdata <- data.frame(is_cell_control)
        n_sets_cell_controls <- 1
    } else {
        if ( is.list(cell_controls) ) {
            cell_controls_list <- cell_controls
            n_sets_cell_controls <- length(cell_controls)
        }
        else {
            cell_controls_list <- list(cell_controls)
            n_sets_cell_controls <- 1
        }
        for (i in seq_len(n_sets_cell_controls) ) {
            cc_set <- cell_controls_list[[i]]
            set_name <- names(cell_controls_list)[i]
            if ( is.logical(cc_set) ) {
                is_cell_control <- cc_set
                cc_set <- which(cc_set)
            } else {
                is_cell_control <- rep(FALSE, ncol(object))
            }
            if (is.character(cc_set))
                cc_set <- which(cellNames(object) %in% cc_set)
            is_cell_control[cc_set] <- TRUE
            ## Construct data.frame for pData from this feature control set
            is_cell_control <- as.data.frame(is_cell_control)
            colnames(is_cell_control) <- paste0("is_cell_control_", set_name)
            if ( exists("cell_controls_pdata") ) {
                cell_controls_pdata <- data.frame(cell_controls_pdata,
                                                  is_cell_control)
            } else
                cell_controls_pdata <- is_cell_control
        }
    }

    ## Check column names and get cell controls across all sets
    if ( n_sets_cell_controls == 1 ) {
        colnames(cell_controls_pdata) <- "is_cell_control"
    } else {
        ## Combine all cell controls
        is_cell_control <- apply(cell_controls_pdata, 1, any)
        cell_controls_pdata <- cbind(cell_controls_pdata, is_cell_control)
    }

    ## Add cell-level QC metrics to pData
    new_pdata <- as.data.frame(pData(object))
    ### Remove columns to be replaced
    to_replace <- colnames(new_pdata) %in%
        c(colnames(qc_pdata), colnames(cell_controls_pdata))
    new_pdata <- new_pdata[, !to_replace, drop = FALSE]
    ### Add new QC metrics
    new_pdata$total_counts <- total_counts
    new_pdata$log10_total_counts <- log10(total_counts)
    new_pdata$filter_on_total_counts <- filter_on_total_counts
    new_pdata$total_features <- total_features
    new_pdata$log10_total_features <- log10(total_features)
    new_pdata$filter_on_total_features <- filter_on_total_features
    new_pdata$pct_dropout <- 100 * colSums(!is_exprs(object)) / nrow(object)
    new_pdata <- cbind(new_pdata, qc_pdata, cell_controls_pdata)
    pData(object) <-  new("AnnotatedDataFrame", new_pdata)

    ## Add feature-level QC metrics to fData
    new_fdata <- as.data.frame(fData(object))
    ### Remove columns that are to be replaced
    to_replace <- colnames(new_fdata) %in% colnames(feature_controls_fdata)
    new_fdata <- new_fdata[, !to_replace, drop = FALSE]
    ### Add new QC information
    new_fdata$mean_exprs <- rowMeans(exprs(object))
    new_fdata$exprs_rank <- rank(rowMeans(exprs(object)))
    new_fdata$n_cells_exprs <- rowSums(is_exprs(object))
    total_exprs <- sum(exprs_mat)
    new_fdata$total_feature_exprs <- rowSums(exprs_mat)
    new_fdata$log10_total_feature_exprs <-
        log10(new_fdata$total_feature_exprs + 1)
    new_fdata$pct_total_exprs <- 100 * rowSums(exprs_mat) / total_exprs
    new_fdata$pct_dropout <- 100 * rowSums(!is_exprs(object)) / ncol(object)
    ## indicate if feature is feature control across any set
    if ( is.list(feature_controls) ) {
        feat_controls_cols <- grep("^is_feature_control",
                                   colnames(feature_controls_fdata))
        feature_controls_fdata$is_feature_control <- (
            rowSums(feature_controls_fdata[, feat_controls_cols, drop = FALSE])
            > 0)
    }

    if ( !is.null(counts_mat) ) {
        total_counts <- sum(counts_mat)
        new_fdata$total_feature_counts <- rowSums(counts_mat)
        new_fdata$log10_total_feature_counts <-
            log10(new_fdata$total_feature_counts + 1)
        new_fdata$pct_total_counts <- 100 * rowSums(counts_mat) / total_counts
    }
    if ( !is.null(tpm_mat) ) {
        total_tpm <- sum(tpm_mat)
        new_fdata$total_feature_tpm <- rowSums(tpm_mat)
        new_fdata$log10_total_feature_tpm <-
            log10(new_fdata$total_feature_tpm + 1)
        new_fdata$pct_total_tpm <- 100 * rowSums(tpm_mat) / total_tpm
    }
    if ( !is.null(fpkm_mat) ) {
        total_fpkm <- sum(fpkm_mat)
        new_fdata$total_feature_fpkm <- rowSums(fpkm_mat)
        new_fdata$log10_total_feature_fpkm <-
            log10(new_fdata$total_feature_fpkm + 1)
        new_fdata$pct_total_fpkm <- 100 * rowSums(fpkm_mat) / total_fpkm
    }
    ## Add new fdata to object
    new_fdata <- cbind(new_fdata, feature_controls_fdata)
    fData(object) <- new("AnnotatedDataFrame", new_fdata)

    ## Ensure sample names are correct and return object
    sampleNames(object) <- colnames(exprs(object))
    object
}

.get_qc_metrics_exprs_mat <- function(exprs_mat, is_feature_control,
                                           pct_feature_controls_threshold,
                                           calc_top_features = FALSE,
                                           exprs_type = "exprs") {
    ## Many thanks to Aaron Lun for suggesting efficiency improvements
    ## for this function.
    ## Get total expression from feature controls
    exprs_feature_controls <- colSums(exprs_mat[is_feature_control,])
    ## Get % expression from feature controls
    pct_exprs_feature_controls <- (100 * exprs_feature_controls /
                                         colSums(exprs_mat))
    ## Indicate whether or not to filter on percentage from controls
    filter_on_pct_exprs_feature_controls <-
        (pct_exprs_feature_controls > pct_feature_controls_threshold)
    ## Make a data frame
    df_pdata_this <- data.frame(exprs_feature_controls,
                                pct_exprs_feature_controls,
                                filter_on_pct_exprs_feature_controls)
    if (calc_top_features) { ## Do we want to calculate exprs accounted for by
        ## top features?
        ## Determine percentage of counts for top features by cell
        pct_exprs_top_out <- rep(list(list()), ncol(exprs_mat))
        for (cell in seq_along(pct_exprs_top_out)) {
            exprs_by_cell_sorted <- sort(exprs_mat[,cell], decreasing = TRUE)
            pct_exprs_top <- cumsum(exprs_by_cell_sorted[1:500]) / 
                sum(exprs_by_cell_sorted) * 100
            pct_exprs_top_out[[cell]] <- pct_exprs_top[c(50, 100, 200)]
        }
        pct_exprs_top_out <- do.call(rbind, pct_exprs_top_out)
        colnames(pct_exprs_top_out) <- paste0("pct_exprs_top_",
                                              c(50, 100, 200), "_features")
        df_pdata_this <- cbind(df_pdata_this, pct_exprs_top_out)
    }
    colnames(df_pdata_this) <- gsub("exprs", exprs_type, colnames(df_pdata_this))
    df_pdata_this
}



################################################################################

#' Find most important principal components for a given variable
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param variable character scalar providing a variable name (column from
#' \code{pData(object)}) for which to determine the most important PCs.
#' @param plot_type character string, indicating which type of plot to produce.
#' Default, \code{"pairs-pcs"} produces a pairs plot for the top 5 PCs based on
#' their R-squared with the variable of interest. A value of
#' \code{"pcs-vs-vars"} produces plots of the top PCs against the variable of
#' interest.
#' @param theme_size numeric scalar providing base font size for ggplot theme.
#'
#' @details Plot the top 5 most important PCs for a given variable. Importance
#' here is defined as the R-squared value from a linear model regressing each PC
#' onto the variable of interest. If the variable is discrete, unique values of
#' the variable define  the clusters. If the variable is continuous, it is coerced into three
#' discrete values, "low" (bottom quartile), "medium" (middle two quartiles)
#' and "high" (upper quartile), and these levels are used as a factor in the
#' linear model. Cells are coloured by their cluster in the plot.
#'
#' @return a \code{\link{ggplot}} plot object
#'
#' @import viridis
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#' example_sceset <- calculateQCMetrics(example_sceset)
#' findImportantPCs(example_sceset, variable="total_features")
#'
findImportantPCs <- function(object, variable="total_features",
                             plot_type = "pcs-vs-vars", theme_size = 10) {
    df_for_pca <- t(exprs(object))
    col_zero_var <- (apply(df_for_pca, 2, var) == 0)
    pca <- prcomp(df_for_pca[, !col_zero_var], retx = TRUE, center = TRUE,
                  scale. = TRUE)
    colnames(pca$x) <- paste("component", 1:ncol(pca$x))
    if (!(variable %in% colnames(pData(object))))
        stop("variable not found in pData(object).
             Please make sure pData(object)[, variable] exists.")
    x <- pData(object)[, variable]
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
        ## If x is a continuous variable - turn into an ordinal variable with three
        ## values - low, medium and high - or five values based on quintiles
        x_int <- cut(x, c(min(x) - 1, quantile(x, probs = c(0.2, 0.4, 0.6, 0.8)),
                                    max(x) + 1), right = TRUE)
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
        colour_by <- pData(object)[, variable]
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
        xvar <- pData(object)[, variable]
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
.getRSquared <- function(y, design) {
    ## Mean-centre rows to get correct R-squared values with the limma formula below
    y0 <- t(scale(t(y), center = TRUE, scale = FALSE))
    ## Get linear model fit
    fit <- limma::lmFit(y0, design = design)
    ## Compute total sum of squares
    sst <- rowSums(y0 ^ 2)
    ## Compute residual sum of squares
    ssr <- sst - fit$df.residual * fit$sigma ^ 2
    (ssr/sst)
}


################################################################################

#' Plot the features with the highest expression values
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param col_by_variable variable name (must be a column name of pData(object))
#' to be used to assign colours to cell-level values.
#' @param n numeric scalar giving the number of the most expressed features to
#' show. Default value is 50.
#' @param drop_features a character, logical or numeric vector indicating which
#' features (e.g. genes, transcripts) to drop when producing the plot. For
#' example, control genes might be dropped to focus attention on contribution
#' from endogenous rather than synthetic genes.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "counts" (default),
#' "tpm", "fpkm" and "exprs".
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:500)
#' plotHighestExprs(example_sceset, col_by_variable="total_features")
#' plotHighestExprs(example_sceset, col_by_variable="Mutation_Status")
#'
plotHighestExprs <- function(object, col_by_variable = "total_features", n = 50,
                              drop_features = NULL, exprs_values = "counts") {
    ## Check that variable to colour points exists
    if (!(col_by_variable %in% colnames(pData(object)))) {
        warning("col_by_variable not found in pData(object).
             Please make sure pData(object)[, variable] exists. Colours will not be plotted.")
        plot_cols <- FALSE
    } else
        plot_cols <- TRUE
    x <- pData(object)[, col_by_variable]
    #     x_na <- is.na(x)
    #     x <- x[!x_na]
    ## Determine type of variable
    typeof_x <- .getTypeOfVariable(object, col_by_variable)
    ## Figure out which features to drop
    if ( !(is.null(drop_features) | length(drop_features) == 0) ) {
        if (is.character(drop_features))
            drop_features <- which(rownames(object) %in% drop_features)
        if (is.logical(drop_features))
            object <- object[!drop_features,]
        else
            object <- object[-drop_features,]
    }
    ## Compute QC metrics on the (possibly) modified SCESet object to make sure
    ## we have the relevant values for this set of features
    if ( !is.null(fData(object)$is_feature_control) )
        object <- calculateQCMetrics(
            object, feature_controls = fData(object)$is_feature_control)
    else
        object <- calculateQCMetrics(object)

    ## Define expression values to be used
    exprs_values <- match.arg(exprs_values,
                              c("exprs", "tpm", "cpm", "fpkm", "counts"))
    exprs_mat <- get_exprs(object, exprs_values)
    if ( is.null(exprs_mat) && !is.null(counts(object)) ) {
        exprs_mat <- counts(object)
        message("Using counts as expression values.")
        exprs_values <- "counts"
    } else if ( is.null(exprs_mat) ) {
        exprs_mat <- exprs(object)
        message("Using exprs(object) values as expression values.")
        exprs_values <- "exprs"
    }
    if ( exprs_values == "exprs" && object@logged )
        exprs_mat <- 2 ^ (exprs_mat) - object@logExprsOffset

    ## Find the most highly expressed features in this dataset
    ### Order by total feature counts across whole dataset
    fdata <- fData(object)
    if ( paste0("total_feature_", exprs_values) %in% colnames(fdata) )
        oo <- order(fdata[[paste0("total_feature_", exprs_values)]],
                    decreasing = TRUE)
    else {
        if ( "total_feature_counts" %in% colnames(fdata) ) {
            oo <- order(fdata[["total_feature_counts"]], decreasing = TRUE)
            exprs_values <- "counts"
            message("Using counts to order total expression of features.")
        }
        else {
            exprs_values <- "exprs"
            oo <- order(fdata[["total_feature_exprs"]], decreasing = TRUE)
            message("Using 'exprs' to order total expression of features.")
        }
    }
    fdata$feature <- factor(featureNames(object),
                            levels = rownames(object)[rev(oo)])
    ## Check if is_feature_control is defined
    if ( is.null(fdata$is_feature_control) )
        fdata$is_feature_control <- rep(FALSE, nrow(fdata))
    if ( is.null(fdata$Feature) )
        fdata$Feature <- featureNames(object)

    ## Determine percentage expression accounted for by top features across all
    ## cells
    total_exprs <- sum(exprs_mat)
    total_feature_exprs <- fdata[[paste0("total_feature_", exprs_values)]]
    top50_pctage <- 100 * sum(total_feature_exprs[oo[1:n]]) / total_exprs
    ## Determine percentage of counts for top features by cell
    df_pct_exprs_by_cell <- (100 * t(exprs_mat[oo[1:n],]) / colSums(exprs_mat))

    ## Melt dataframe so it is conducive to ggplot
    df_pct_exprs_by_cell_long <- reshape2::melt(df_pct_exprs_by_cell)
    df_pct_exprs_by_cell_long$Var2 <- factor(
        df_pct_exprs_by_cell_long$Var2, levels = rownames(object)[rev(oo[1:n])])
    ## Add colour variable information
    if (typeof_x == "discrete")
        df_pct_exprs_by_cell_long$colour_by <- factor(x)
    else
        df_pct_exprs_by_cell_long$colour_by <- x
    ## Make plot
    plot_most_expressed <- ggplot(df_pct_exprs_by_cell_long,
                                  aes_string(y = "Var2", x = "value",
                                             colour = "colour_by")) +
        geom_point(alpha = 0.6, shape = 124) +
        ggtitle(paste0("Top ", n, " account for ",
                       format(top50_pctage, digits = 3), "% of total")) +
        ylab("Feature") +
        xlab(paste0("% of total ", exprs_values)) +
        theme_bw(8) +
        theme(legend.position = c(1, 0), legend.justification = c(1, 0),
              axis.text.x = element_text(colour = "gray35"),
              axis.text.y = element_text(colour = "gray35"),
              axis.title.x = element_text(colour = "gray35"),
              axis.title.y = element_text(colour = "gray35"),
              title = element_text(colour = "gray35"))
    ## Sort of colouring of points
    if (typeof_x == "discrete") {
        plot_most_expressed <- .resolve_plot_colours(
            plot_most_expressed, df_pct_exprs_by_cell_long$colour_by,
            col_by_variable)
#         plot_most_expressed <- plot_most_expressed +
#             ggthemes::scale_colour_tableau(name = col_by_variable)
    } else {
        plot_most_expressed <- plot_most_expressed +
            scale_colour_gradient(name = col_by_variable, low = "lightgoldenrod",
                                  high = "firebrick4", space = "Lab")
    }
    plot_most_expressed + geom_point(
        aes_string(x = paste0("as.numeric(pct_total_", exprs_values, ")"),
                   y = "Feature", fill = "is_feature_control"),
        data = fdata[oo[1:n],], colour = "gray30", shape = 21) +
        scale_fill_manual(values = c("aliceblue", "wheat")) +
        guides(fill = guide_legend(title = "Feature control?"))
}


.getTypeOfVariable <- function(object, variable) {
    ## Extract variable
    x <- pData(object)[, variable]
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
". Variable being coerced to discrete. Please make sure pData(object)[, variable] is a proper discrete or continuous variable"))
            }
        }
    }
    typeof_x
}

################################################################################


#' Plot explanatory variables ordered by percentage of phenotypic variance explained
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param method character scalar indicating the type of plot to produce. If
#' "density", the function produces a density plot of R-squared values for each
#' variable when fitted as the only explanatory variable in a linear model. If
#' "pairs", then the function produces a pairs plot of the explanatory variables
#' ordered by the percentage of feature expression variance (as measured by
#' R-squared in a marginal linear model) explained.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "exprs" (default),
#' "tpm", "fpkm", "cpm", and "counts".
#' @param nvars_to_plot integer, the number of variables to plot in the pairs
#' plot. Default value is 10.
#' @param min_marginal_r2 numeric scalar giving the minimal value required for
#' median marginal R-squared for a variable to be plotted. Only variables with a
#' median marginal R-squared strictly larger than this value will be plotted.
#' @param variables optional character vector giving the variables to be plotted.
#' Default is \code{NULL}, in which case all variables in \code{pData(object)}
#' are considered and the \code{nvars_to_plot} variables with the highest median
#' marginal R-squared are plotted.
#' @param return_object logical, should an \code{SCESet} object with median
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
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#' example_sceset <- calculateQCMetrics(example_sceset)
#' vars <- names(pData(example_sceset))[c(2:3, 5:14)]
#' plotExplanatoryVariables(example_sceset, variables=vars)
#'
plotExplanatoryVariables <- function(object, method = "density",
                                     exprs_values = "exprs", nvars_to_plot = 10,
                                     min_marginal_r2 = 0, variables = NULL,
                                     return_object = FALSE, theme_size = 10,
                                     ...) {
    ## Check method argument
    method <- match.arg(method, c("density", "pairs"))
    ## Checking arguments for expression values
    exprs_values <- match.arg(
        exprs_values,
        choices = c("exprs", "norm_exprs", "stand_exprs", "norm_exprs",
                    "counts", "norm_counts", "tpm", "norm_tpm", "fpkm",
                    "norm_fpkm", "cpm", "norm_cpm"))
    exprs_mat <- get_exprs(object, exprs_values)

    ## Check that variables are defined
    if ( is.null(variables) ) {
        variables_to_plot <- varLabels(object)
    } else {
        variables_to_plot <- NULL
        for (var in variables) {
            if ( !(var %in% colnames(pData(object))) ) {
                warning(paste("variable", var, "not found in pData(object).
                     Please make sure pData(object)[, variable] exists. This variable will not be plotted."))
            } else {
                variables_to_plot <- c(variables_to_plot, var)
            }
        }
    }
    variables_all <- varLabels(object)

    ## Initialise matrix to store R^2 values for each feature for each variable
    rsquared_mat <- matrix(NA, nrow = nrow(object),
                           ncol = length(variables_all))
    val_to_plot_mat <- matrix(NA, nrow = ncol(object),
                              ncol = length(variables_all))
    colnames(rsquared_mat) <- colnames(val_to_plot_mat) <- variables_all
    rownames(rsquared_mat) <- rownames(object)
    rownames(val_to_plot_mat) <- colnames(object)

    ## Get R^2 values for each feature and each variable
    for (var in variables_all) {
        if ( var %in% variables_to_plot ) {
            if (length(unique(pData(object)[, var])) <= 1) {
                message(paste("The variable", var, "only has one unique value, so R^2 is not meaningful.
This variable will not be plotted."))
                rsquared_mat[, var] <- NA
            } else {
                x <- pData(object)[, var]
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
                rsquared_mat[, var] <- .getRSquared(exprs_mat, design)
#                 rsq_base <- apply(exprs_mat, 1, function(y) {
#                     lm.first <- lm(y ~ -1 + design); summary(lm.first)$r.squared})
#                 all(abs(rsq_base - rsquared_mat[, var]) < 0.000000000001)
            }
        } else {
            rsquared_mat[, var] <- NA
        }
    }

    ## Get median R^2 for each variable, add to labels and order by median R^2
    median_rsquared <- apply(rsquared_mat, 2, median)
    oo_median <- order(median_rsquared, decreasing = TRUE)
    nvars_to_plot <- min(sum(median_rsquared > min_marginal_r2, na.rm = TRUE),
                         nvars_to_plot)

    if ( method == "pairs" ) {
        ## Generate a larger data.frame for pairs plot
        df_to_expand <- val_to_plot_mat[, oo_median[1:nvars_to_plot]]
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
            rsquared_mat[, oo_median[1:nvars_to_plot]]))
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
            ylab("") +
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
        fdata <- fData(object)
        rsq_out <- rsquared_mat[, oo_median[1:nvars_to_plot]]
        colnames(rsq_out) <- paste0("Rsq_", colnames(rsq_out))
        fdata_new <- new("AnnotatedDataFrame", cbind(fdata, rsq_out))
        fData(object) <- fdata_new
        print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}


# ### Emergency test code
# rsq <- fData(sce_simmons_filt_strict)[, grep("^Rsq", names(fData(sce_simmons_filt_strict)))]
# head(rsq)
# colSums(rsq > 0.99)
#
# rsq[rsq$Rsq_batch > 0.99,]
#
# plotExpression(sce_simmons_filt_strict, rownames(rsq)[rsq$Rsq_batch > 0.99], x = "batch",
#                ncol = 6, colour_by = "tissue_type", exprs_values = "norm_exprs")
#
# design <- model.matrix(~sce_simmons_filt_strict$batch)
# lm.first <- lm(norm_exprs(sce_simmons_filt_strict)["ERCC-00002",] ~ -1 + design)
# summary(lm.first)$r.squared
# summary(lm.first)$residuals
# plot(lm.first)
# lm2 <- lm(norm_exprs(sce_simmons_filt_strict)["ERCC-00002",] ~ sce_simmons_filt_strict$batch)
# summary(lm2)$r.squared
#
# rsq_limma <- .getRSquared(norm_exprs(sce_simmons_filt_strict), design)
#
# rsq_base2 <- apply(norm_exprs(sce_simmons_filt_strict), 1, function(y) {
#     lm.first <- lm(y ~ sce_simmons_filt_strict$batch); summary(lm.first)$r.squared})
#
# all(abs(rsq_base2 - rsq_limma) < 0.00000000000001)
#
# rsq_base <- apply(norm_exprs(sce_simmons_filt_strict), 1, function(y) {
#     lm.first <- lm(y ~ -1 + design); summary(lm.first)$r.squared})
# all(abs(rsq_base - rsq$Rsq_batch) < 0.000000000001)


################################################################################
### Plot expression frequency vs mean for feature controls

#' Plot frequency of expression against mean expression level
#'
#' @param object an \code{SCESet} object.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to plot. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}. If \code{NULL}, then the function checks if
#' feature controls are defined. If so, then only feature controls are plotted,
#' if not, then all features are plotted.
#' @param feature_controls character, numeric or logical vector indicating a set of
#' features to be used as feature controls for computing technical dropout
#' effects. If character, entries must all be in \code{featureNames(object)}. If
#' numeric, values are taken to be indices for features. If logical, vector is
#' used to index features and should have length equal to \code{nrow(object)}.
#' If \code{NULL}, then the function checks if feature controls are defined. If
#' so, then these feature controls are used.
#' @param shape (optional) numeric scalar to define the plotting shape.
#' @param alpha (optional) numeric scalar (in the interval 0 to 1) to define the
#'  alpha level (transparency) of plotted points.
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
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' ex_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' ex_sceset <- calculateQCMetrics(ex_sceset)
#' plotExprsFreqVsMean(ex_sceset)
#'
plotExprsFreqVsMean <- function(object, feature_set = NULL,
                                feature_controls = NULL, shape = 1, alpha = 0.7,
                                ...) {
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( is.null(fData(object)$n_cells_exprs) ||
         is.null(fData(object)$mean_exprs)) {
        stop("fData(object) does not have both 'n_cells_exprs' and 'mean_exprs' columns. Try running 'calculateQCMetrics' on this object, and then rerun this command.")
    }
    if ( !is.null(feature_set) && feature_set != "feature_controls" &&
         typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character and not 'feature_controls', all features must be in featureNames(object)")
    }
    if ( is.null(feature_set) ) {
        feature_set_logical <- rep(TRUE, nrow(object))
        x_lab <- "Mean expression level (all features)"
    } else {
        if ( length(feature_set) == 1 && feature_set == "feature_controls" ) {
            feature_set_logical <- fData(object)$is_feature_control
            x_lab <- "Mean expression level (feature controls)"
        } else {
            if ( is.character(feature_set) )
                feature_set_logical <- featureNames(object) %in% feature_set
            else {
                feature_set_logical <- rep(FALSE, nrow(object))
                if ( is.numeric(feature_set) )
                    feature_set_logical[feature_set] <- TRUE
                else
                    feature_set_logical <- feature_set
            }
               x_lab <- "Mean expression level (supplied feature set)"
        }
    }
    ## check that feature controls, if defined, are sensible
    if (!is.null(feature_controls) && typeof(feature_controls) == "character") {
        if ( !(all(feature_controls %in% featureNames(object))) )
            stop("when the argument 'feature_controls' is of type character all features must be in featureNames(object)")
    }
    if ( is.null(feature_controls) )
        feature_controls_logical <- fData(object)$is_feature_control
    else {
        if ( is.character(feature_controls) )
            feature_controls_logical <- (featureNames(object) %in%
                                             feature_controls)
        else {
            feature_controls_logical <- rep(FALSE, nrow(object))
            if ( is.numeric(feature_controls) )
                feature_controls_logical[feature_controls] <- TRUE
            else
                feature_controls_logical <- feature_controls
        }
        fData(object)$is_feature_control <- feature_controls_logical
    }

    ## define percentage of cells expressing a gene
    fData(object)$pct_cells_exprs <- (100 * fData(object)$n_cells_exprs /
                                          ncol(object))
    y_lab <- paste0("Frequency of expression (% of ", ncol(object), " cells)")

    ## Plot this
    if ( any(fData(object)$is_feature_control[feature_set_logical]) &&
         !all(fData(object)$is_feature_control[feature_set_logical]) ) {
        plot_out <- plotFeatureData(object[feature_set_logical,],
                                    aes_string(x = "mean_exprs",
                                               y = "pct_cells_exprs",
                                               colour = "is_feature_control",
                                               shape = "is_feature_control"),
                                    alpha = alpha, ...) +
            scale_shape_manual(values = c(1, 17)) +
            ylab(y_lab) +
            xlab(x_lab)
    } else {
        plot_out <- plotFeatureData(object[feature_set_logical,],
                                    aes_string(x = "mean_exprs",
                                               y = "pct_cells_exprs"),
                                    alpha = alpha, shape = shape, ...) +
            ylab(y_lab) +
            xlab(x_lab)
    }

    ## add dropout information if feature controls are present
    if ( any(feature_controls_logical) ) {
        ## data frame with expression mean and frequency for feature controls
        mn_vs_fq <- data.frame(
            mn = fData(object)$mean_exprs,
            fq = fData(object)$pct_cells_exprs / 100)

        ## estimate 50% spike-in dropout
        dropout <- nls(fq ~ (1 / (1 + exp(-(-i + 1 * mn)))),
                       start = list(i = 5),
                       data = mn_vs_fq[feature_controls_logical,])
        text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))

        ## add annotations to existing plot
        plot_out <- plot_out +
            geom_hline(yintercept = 50, linetype = 2) + # 50% dropout
            geom_vline(xintercept = coef(dropout), linetype = 2) +
            # expression at 50% dropout
            annotate("text", x = text_x_loc, y = 60, label = paste(
                sum(mn_vs_fq[!feature_controls_logical, "mn"] > coef(dropout)),
                " genes above\nhigh technical dropout", sep = "")) +
            annotate("text", x = text_x_loc, y = 40, label =  paste(
                sum(mn_vs_fq$fq >= 0.5),
                " genes are expressed\nin at least 50% of cells", sep = "" )) +
            annotate("text", x = text_x_loc, y = 20, label =  paste(
                sum(mn_vs_fq$fq >= 0.25),
                " genes are expressed\nin at least 25% of cells", sep = "" ))
    }

    ## return the plot object
    plot_out
}




################################################################################

#' Produce QC diagnostic plots
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param type character scalar providing type of QC plot to compute:
#' "highest-expression" (showing features with highest expression), "find-pcs" (showing
#' the most important principal components for a given variable),
#' "explanatory-variables" (showing a set of explanatory variables plotted
#' against each other, ordered by marginal variance explained), or
#' "exprs-mean-vs-freq" (plotting the mean expression levels against the
#' frequency of expression for a set of features).
#' @param ... arguments passed to \code{plotHighestExprs},
#' \code{plotImportantPCs}, \code{plotExplanatoryVariables} and
#' \code{plotExprsMeanVsFreq} as appropriate.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic features and cells.
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotQC(example_sceset, type="high", col_by_variable="Mutation_Status")
#' plotQC(example_sceset, type="find", variable="total_features")
#' vars <- names(pData(example_sceset))[c(2:3, 5:14)]
#' plotQC(example_sceset, type="expl", variables=vars)
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

