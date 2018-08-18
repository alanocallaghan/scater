#' Calculate QC metrics
#'
#' Compute quality control (QC) metrics for each feature and cell in a SingleCellExperiment object, accounting for specified control sets.
#'
#' @param object A SingleCellExperiment object containing expression values, usually counts.
#' @param exprs_values A string indicating which \code{assays} in the \code{object} should be used to define expression.
#' @param feature_controls A named list containing one or more vectors (a character vector of feature names, a logical vector, or a numeric vector of indices),
#' used to identify feature controls such as ERCC spike-in sets or mitochondrial genes. 
#' @param cell_controls A named list containing one or more vectors (a character vector of cell (sample) names, a logical vector, or a numeric vector of indices),
#' used to identify cell controls, e.g., blank wells or bulk controls.
#' @param percent_top An integer vector. 
#' Each element is treated as a number of top genes to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' See \code{pct_X_top_Y_features} below for more details.
#' @param detection_limit A numeric scalar to be passed to \code{\link{nexprs}}, specifying the lower detection limit for expression.
#' @param use_spikes A logical scalar indicating whether existing spike-in sets in \code{object} should be automatically added to \code{feature_controls}, see \code{?\link{isSpike}}.
#' @param compact A logical scalar indicating whether the metrics should be returned in a compact format as a nested DataFrame.
#' @param BPPARAM A BiocParallelParam object specifying whether the QC calculations should be parallelized. 
#'
#' @details 
#' This function calculates useful quality control metrics to help with pre-processing of data and identification of potentially problematic features and cells. 
#' 
#' Underscores in \code{assayNames(object)} and in \code{feature_controls} or \code{cell_controls} can cause theoretically cause ambiguities in the names of the output metrics.
#' While problems are highly unlikely, users are advised to avoid underscores when naming their controls/assays.
#'
#' @section Cell-level QC metrics:
#' Denote the value of \code{exprs_values} as \code{X}. 
#' Cell-level metrics are:
#' \describe{
#' \item{\code{total_X}:}{Sum of expression values for each cell (i.e., the library size, when counts are the expression values).}
#' \item{\code{log10_total_X}:}{Log10-transformed \code{total_X} after adding a pseudo-count of 1.}
#' \item{\code{total_features_by_X}:}{The number of features that have expression values above the detection limit.}
#' \item{\code{log10_total_features_by_X}:}{Log10-transformed \code{total_features_by_X} after adding a pseudo-count of 1.}
#' \item{\code{pct_X_in_top_Y_features}:}{The percentage of the total that is contained within the top \code{Y} most highly expressed features in each cell.
#' This is only reported when there are more than \code{Y} features. 
#' The top numbers are specified via \code{percent_top}.}
#' }
#' 
#' If any controls are specified in \code{feature_controls}, the above metrics will be recomputed using only the features in each control set. 
#' The name of the set is appended to the name of the recomputed metric, e.g., \code{total_X_F}.
#' A \code{pct_X_F} metric is also calculated for each set, representing the percentage of expression values assigned to features in \code{F}.
#' 
#' In addition to the user-specified control sets, two other sets are automatically generated when \code{feature_controls} is non-empty.
#' The first is the \code{"feature_control"} set, containing a union of all feature control sets;
#' and the second is an \code{"endogenous"} set, containing all genes not in any control set. 
#' Metrics are also computed for these sets in the same manner described above, suffixed with \code{_feature_control} and \code{_endogenous} instead of \code{_F}.
#'
#' Finally, there is the \code{is_cell_control} field, which indicates whether each cell has been defined as a cell control by \code{cell_controls}. 
#' If multiple sets of cell controls are defined (e.g., blanks or bulk libraries), a metric \code{is_cell_control_C} is produced for each cell control set \code{C}.
#' The union of all sets is stored in \code{is_cell_control}.
#'
#' All of these cell-level QC metrics are added as columns to the \code{colData} slot of the SingleCellExperiment object.
#' This allows them to be inspected by the user and makes them readily available for other functions to use. 
#'
#' @section Feature-level QC metrics:
#' Denote the value of \code{exprs_values} as \code{X}. 
#' Feature-level metrics are:
#' \describe{
#' \item{\code{mean_X}:}{Mean expression value for each gene across all cells.}
#' \item{\code{log10_mean_X}:}{Log10-mean expression value for each gene across all cells.}
#' \item{\code{n_cells_by_X}:}{Number of cells with expression values above the detection limit for each gene.}
#' \item{\code{pct_dropout_by_X}:}{Percentage of cells with expression values  below the detection limit for each gene.}
#' \item{\code{total_X}:}{Sum of expression values for each gene across all cells.}
#' \item{\code{log10_total_X}:}{Log10-sum of expression values for each gene across all cells.}
#' }
#'
#' If any controls are specified in \code{cell_controls}, the above metrics will be recomputed using only the cells in each control set. 
#' The name of the set is appended to the name of the recomputed metric, e.g., \code{total_X_C}.
#' A \code{pct_X_C} metric is also calculated for each set, representing the percentage of expression values assigned to cells in \code{C}.
#' 
#' In addition to the user-specified control sets, two other sets are automatically generated when \code{cell_controls} is non-empty. 
#' The first is the \code{"cell_control"} set, containing a union of all cell control sets;
#' and the second is an \code{"non_control"} set, containing all genes not in any control set. 
#' Metrics are computed for these sets in the same manner described above, suffixed with \code{_cell_control} and \code{_non_control} instead of\code{_C}.
#'
#' Finally, there is the \code{is_feature_control} field, which indicates whether each feature has been defined as a control by \code{feature_controls}. 
#' If multiple sets of feature controls are defined (e.g., ERCCs, mitochondrial genes), a metric \code{is_feature_control_F} is produced for each feature control set \code{F}.
#' The union of all sets is stored in \code{is_feature_control}.
#'
#' These feature-level QC metrics are added as columns to the \code{rowData} slot of the SingleCellExperiment object.
#' They can be inspected by the user and are readily available for other functions to use. 
#'
#' @section Compacted output:
#' If \code{compact=TRUE}, the QC metrics are stored in the \code{"scater_qc"} field of the \code{colData} and \code{rowData} as a nested DataFrame. 
#' This avoids cluttering the metadata with QC metrics, especially if many results are to be stored in a single SingleCellExperiment object. 
#'
#' Assume we have a feature control set \code{F} and a cell control set \code{C}.
#' The nesting structure in \code{scater_qc} in the \code{colData} is:
#' \preformatted{  scater_qc
#'   |-- is_cell_control
#'   |-- is_cell_control_C
#'   |-- all
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'   |   \-- ...
#'   +-- endogenous
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'       |-- pct_counts
#'   |   \-- ...
#'   +-- feature_control
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'       |-- pct_counts
#'   |   \-- ...
#'   \-- feature_control_F
#'       |-- total_counts
#'       |-- total_features_by_counts
#'       |-- pct_counts
#'       \-- ...
#' }
#' The nesting in \code{scater_qc} in the \code{rowData} is:
#' \preformatted{  scater_qc
#'   |-- is_feature_control
#'   |-- is_feature_control_F
#'   |-- all
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'   |   \-- ...
#'   +-- non_control 
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'       |-- pct_counts
#'   |   \-- ...
#'   +-- cell_control
#'   |   |-- total_counts
#'   |   |-- total_features_by_counts
#'       |-- pct_counts
#'   |   \-- ...
#'   \-- cell_control_C
#'       |-- total_counts
#'       |-- total_features_by_counts
#'       |-- pct_counts
#'       \-- ...
#' }
#'
#' No suffixing of the metric names by the control names is performed here. 
#' This is not necessary when each control set has its own nested DataFrame.
#'
#' @section Renamed metrics:
#' Several metric names have been changed in \pkg{scater} 1.7.5:
#' \itemize{
#'   \item \code{total_features} was changed to \code{total_features_by_X} 
#'   where \code{X} is the \code{exprs_values}. This avoids ambiguities if 
#'   \code{calculateQCMetrics} is called multiple times with different \code{exprs_values}.
#'   \item \code{n_cells_X} was changed to \code{n_cells_by_X}, to provide
#'   a more sensible name for the metric.
#'   \item \code{pct_dropout_X} was changed to \code{pct_dropout_by_X}.
#'   \item \code{pct_X_top_Y_features} was changed to \code{pct_X_in_top_Y_features}.
#' }
#' The old metric names have been removed in version 1.9.10.
#' 
#' @return A SingleCellExperiment object containing QC metrics in the row and column metadata.
#'
#' @author Davis McCarthy, with (many!) modifications by Aaron Lun
#'
#' @export
#' @importFrom methods is 
#' @importFrom SummarizedExperiment assay rowData rowData<- colData colData<-
#' @importFrom SingleCellExperiment isSpike spikeNames
#' @importFrom BiocParallel SerialParam bpmapply bpnworkers
#' @importFrom utils head tail
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## with a set of feature controls defined
#' example_sce <- calculateQCMetrics(example_sce, 
#' feature_controls = list(set1 = 1:40))
#' 
#' ## with a named set of feature controls defined
#' example_sce <- calculateQCMetrics(example_sce, 
#'      feature_controls = list(ERCC = 1:40))
#' 
calculateQCMetrics <- function(object, exprs_values="counts", feature_controls = NULL, cell_controls = NULL,
        percent_top = c(50, 100, 200, 500), detection_limit = 0, use_spikes = TRUE, compact = FALSE, 
        BPPARAM=SerialParam()) {

    if ( !is(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment")
    }
    exprs_mat <- assay(object, i = exprs_values)
    percent_top <- sort(as.integer(percent_top))
    
    # Add any existing spike-ins to the feature control set.
    existing_spikes <- spikeNames(object)
    if (use_spikes && length(existing_spikes)) {
        existing <- vector("list", length(existing_spikes))
        names(existing) <- existing_spikes
        for (spset in existing_spikes) {
            existing[[spset]] <- isSpike(object, type=spset)
        }

        already_there <- names(existing) %in% names(feature_controls)
        if (any(already_there)) {
            warning(sprintf("spike-in set '%s' overwritten by feature_controls set of the same name",
                    names(existing)[already_there][1]))
        }
        feature_controls <- c(feature_controls, existing[!already_there])
    }

    # Assemble controls.
    feat_out <- .organize_controls(feature_controls, "feature_control", "endogenous", exprs_mat, byrow=TRUE)
    all_feature_sets <- feat_out$sets
    feature_set_rdata <- feat_out$metadata

    cell_out <- .organize_controls(cell_controls, "cell_control", "non_control", exprs_mat, byrow=FALSE)
    all_cell_sets <- cell_out$sets
    cell_set_cdata <- cell_out$metadata

    # Computing all QC metrics, with cells split across workers. 
    n.cores <- bpnworkers(BPPARAM)
    boundaries <- as.integer(seq(from = 0L, to = ncol(exprs_mat), length.out = n.cores + 1L)) # zero indexed.
    work.ends <- tail(boundaries, -1L)
    work.starts <- head(boundaries, -1L)

    bp.out <- bpmapply(.compute_qc_metrics, start=work.starts, end=work.ends,
            MoreArgs=list(exprs_mat=exprs_mat, 
                all_feature_sets=all_feature_sets, 
                all_cell_sets=all_cell_sets, 
                percent_top=percent_top, 
                detection_limit=detection_limit),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    # Aggregating across cores (concatenating pre-cell statistics, summing per-feature statistics).
    cell_stats_by_feature_set <- bp.out[[1]][[1]]
    feature_stats_by_cell_set <- bp.out[[1]][[2]]
    if (n.cores > 1L) {
        for (i in seq_along(cell_stats_by_feature_set)) {
            current <- lapply(bp.out, FUN=function(sublist) { sublist[[1]][[i]] })
            cell_stats_by_feature_set[[i]] <- list(
                unlist(lapply(current, "[[", i=1L)),  # total count
                unlist(lapply(current, "[[", i=2L)),  # total features 
                do.call(cbind, lapply(current, "[[", i=3L)) # percentage in top X.
            )
        }
        
        for (i in seq_along(feature_stats_by_cell_set)) {
            current <- lapply(bp.out, FUN=function(sublist) { sublist[[2]][[i]] })
            feature_stats_by_cell_set[[i]] <- list(
                Reduce("+", lapply(current, "[[", i=1)), # total count
                Reduce("+", lapply(current, "[[", i=2))  # total non-zero cells
            )
        }
    }

    # Formatting the output of the loop.
    names(cell_stats_by_feature_set) <- c("all", names(all_feature_sets))
    total_libsize <- cell_stats_by_feature_set$all[[1]]

    for (fset in names(cell_stats_by_feature_set)) {
        if (fset!="all"){ 
            overall_total <- total_libsize
        } else {
            overall_total <- NULL
        }
        cell_stats_by_feature_set[[fset]] <- .format_per_cell(
            cell_stats_by_feature_set[[fset]], 
            exprs_type=exprs_values, percent_top=percent_top, 
            overall_total=overall_total)
    }

    names(feature_stats_by_cell_set) <- c("all", names(all_cell_sets))
    total_per_feature <- feature_stats_by_cell_set$all[[1]]

    for (cset in names(feature_stats_by_cell_set)) {
        if (cset!="all"){ 
            overall_total <- total_per_feature 
            total_ncells <- length(all_cell_sets[[cset]])
        } else {
            overall_total <- NULL
            total_ncells <- ncol(exprs_mat)
        }
        feature_stats_by_cell_set[[cset]] <- .format_per_feature(
            feature_stats_by_cell_set[[cset]], 
            exprs_type=exprs_values, total_ncells=total_ncells,
            overall_total=overall_total)
    }

    # Formatting output depending on whether we're compacting or not.
    if (compact) {
        scater_cd <- .convert_to_nested_DataFrame(colData(object)$scater_qc, 
            cell_set_cdata, cell_stats_by_feature_set)
        scater_rd <- .convert_to_nested_DataFrame(rowData(object)$scater_qc,
            feature_set_rdata, feature_stats_by_cell_set)
        colData(object)$scater_qc <- scater_cd
        rowData(object)$scater_qc <- scater_rd
    } else {
        scater_cd <- .convert_to_full_DataFrame(colData(object), cell_set_cdata, 
            cell_stats_by_feature_set, trim.fun=function(x) sub("^feature_control_", "", x))
        scater_rd <- .convert_to_full_DataFrame(rowData(object), feature_set_rdata, 
            feature_stats_by_cell_set, trim.fun=function(x) sub("^cell_control_", "", x))
        colData(object) <- scater_cd
        rowData(object) <- scater_rd
    }
    return(object)
}

##################################################

.organize_controls <- function(control_list, control_name, other_name, mat, byrow) { 
    dimlen <- if (byrow) nrow(mat) else ncol(mat)

    if (length(control_list)) { 
        if (is.null(names(control_list))) {
            stop(sprintf("%ss should be named", control_name))
        }
        
        # Converting to integer indices for all applications.
        reindexed <- lapply(control_list, FUN = .subset2index, target = mat, byrow = byrow)
        names(reindexed) <- sprintf("%s_%s", control_name, names(reindexed))

        is_con <- Reduce(union, reindexed)
        not_con <- seq_len(dimlen)[-is_con]
        extras <- list(not_con, is_con)
        names(extras) <- c(other_name, control_name)
        all_sets <- c(extras, reindexed)

        # But storing logical vectors in the metadata.
        as_logical <- vector("list", length(reindexed)+1)
        names(as_logical) <- c(control_name, names(reindexed))
        for (set in names(as_logical)) { 
            new_set <- logical(dimlen)
            new_set[all_sets[[set]]] <- TRUE
            as_logical[[set]] <- new_set
        }

    } else {
        all_sets <- list()
        as_logical <- list(logical(dimlen))
        names(as_logical) <- control_name
    }

    return(list(sets=all_sets, metadata=as_logical))
}

.compute_qc_metrics <- function(exprs_mat, start, end, all_feature_sets, all_cell_sets, percent_top, detection_limit) 
# A helper function defined in the scater namespace.
# This avoids the need to reattach scater in bpmapply for SnowParam().
{
    .Call(cxx_combined_qc, exprs_mat, start, end, all_feature_sets, all_cell_sets, percent_top, detection_limit)
}

##################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics colnames<-
.format_per_cell <- function(values, exprs_type, percent_top, overall_total=NULL) {
    total <- values[[1]]
    nfeatures <- values[[2]]
    pct_top <- t(values[[3]])

    cell_data <- DataFrame(nfeatures, log10(nfeatures + 1))
    colnames(cell_data) <- paste0(c("", "log10_"), "total_features_by_", exprs_type)

    cell_data[[paste0("total_", exprs_type)]] <- total
    cell_data[[paste0("log10_total_", exprs_type)]] <- log10(total + 1)
    
    if (!is.null(overall_total)) { # Computing percentages of actual total.
        cell_data[[paste0("pct_", exprs_type)]] <- 100 * total / overall_total
    }

    colnames(pct_top) <- sprintf("pct_%s_in_top_%i_features", exprs_type, percent_top)
    cbind(cell_data, pct_top)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics colnames<-
.format_per_feature <- function(values, exprs_type, total_ncells, overall_total=NULL) {
    sum_exprs <- values[[1]]
    ncells <- values[[2]]

    ave <- sum_exprs/total_ncells
    feature_data <- DataFrame(ave, log10(ave + 1))
    colnames(feature_data) <- paste0(c("mean", "log10_mean"), "_", exprs_type)

    feature_data[[paste0("n_cells_by_", exprs_type)]] <- ncells
    feature_data[[paste0("pct_dropout_by_", exprs_type)]] <- 100 * (1 - ncells/total_ncells)

    feature_data[[paste0("total_", exprs_type)]] <- sum_exprs
    feature_data[[paste0("log10_total_", exprs_type)]] <- log10(sum_exprs + 1)

    if (!is.null(overall_total)) { 
        feature_data[[paste0("pct_", exprs_type)]] <- sum_exprs/overall_total * 100
    }

    feature_data
}

##################################################

#' @importFrom BiocGenerics cbind
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.convert_to_nested_DataFrame <- function(existing, set_list, stat_list, exprs_values = "counts") {
    n_values <- length(stat_list[[1]][[1]]) # There should be at least one statistic.
    output <- .create_outer_DataFrame(set_list, n_values)

    sub_output <- new("DataFrame", nrows=n_values)
    for (x in names(stat_list)) { 
        current <- stat_list[[x]] # need to do it via "[[<-" to store DataFrames as columns.
        current <- .cbind_overwrite_DataFrames(existing[[x]], current)
        sub_output[[x]] <- current
    }
    
    output <- cbind(output, sub_output)
    .cbind_overwrite_DataFrames(existing, output)
}

#' @importFrom BiocGenerics colnames<- cbind
.convert_to_full_DataFrame <- function(existing, set_list, stat_list, trim.fun=identity) {
    n_values <- length(stat_list[[1]][[1]]) # There should be at least one statistic.
    output <- .create_outer_DataFrame(set_list, n_values)

    collected <- stat_list
    for (x in names(stat_list)) { 
        current <- stat_list[[x]]

        # For consistency with old output.
        if (x!="all") { 
            colnames(current) <- sprintf("%s_%s", colnames(current), trim.fun(x))
        }

        collected[[x]] <- current
    }

    combined <- do.call(cbind, c(list(output), unname(collected)))
    .cbind_overwrite_DataFrames(existing, combined)
}

#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.create_outer_DataFrame <- function(set_list, n_values) {
    if (length(set_list)) { 
        output <- do.call(DataFrame, set_list)
        colnames(output) <- sprintf("is_%s", colnames(output)) 
    } else {
        output <- new("DataFrame", nrows=n_values)
    }
    return(output)
}

#' @importFrom BiocGenerics colnames cbind
.cbind_overwrite_DataFrames <- function(existing, updated) {
    if (is.null(existing)) { 
        return(updated)
    }
    existing <- existing[, !(colnames(existing) %in% colnames(updated)), drop = FALSE]
    cbind(existing, updated)   
} 
