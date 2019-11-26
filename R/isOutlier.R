#' Identify outlier values 
#' 
#' Convenience function to determine which values in a numeric vector are outliers based on the median absolute deviation (MAD).
#'
#' @param metric Numeric vector of values.
#' @param nmads A numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier.
#' @param type String indicating whether outliers should be looked for at both tails (\code{"both"}), only at the lower tail (\code{"lower"}) or the upper tail (\code{"higher"}).
#' @param log Logical scalar, should the values of the metric be transformed to the log10 scale before computing MADs?
#' @param subset Logical or integer vector, which subset of values should be used to calculate the median/MAD? 
#' If \code{NULL}, all values are used.
#' @param batch Factor of length equal to \code{metric}, specifying the batch to which each observation belongs. 
#' A median/MAD is calculated for each batch, and outliers are then identified within each batch.
#' @param share_medians Logical scalar indicating whether the median calculation should be shared across batches.
#' Only used if \code{batch} is specified.
#' @param share_mads Logical scalar indicating whether the MAD calculation should be shared across batches.
#' Only used if \code{batch} is specified.
#' @param min_diff A numeric scalar indicating the minimum difference from the median to consider as an outlier. 
#' Ignored if \code{NA}.
#' 
#' @return A logical vector of the same length as the \code{metric} argument, specifying the observations that are considered as outliers.
#'
#' @details
#' Lower and upper thresholds are stored in the \code{"threshold"} attribute of the returned vector.
#' By default, this is a numeric vector of length 2 for the threshold on each side.
#' If \code{type="lower"}, the higher limit is \code{Inf}, while if \code{type="higher"}, the lower limit is \code{-Inf}.
#' 
#' If \code{min_diff} is not \code{NA}, the minimum distance from the median required to define an outlier is set as the larger of \code{nmads} MADs and \code{min_diff}.
#' This aims to avoid calling many outliers when the MAD is very small, e.g., due to discreteness of the metric.
#' 
#' If \code{subset} is specified, the median and MAD are computed from a subset of cells and the values are used to define the outlier threshold that is applied to all cells.
#' In a quality control context, this can be handy for excluding groups of cells that are known to be low quality (e.g., failed plates) so that they do not distort the outlier definitions for the rest of the dataset.
#' 
#' Missing values trigger a warning and are automatically ignored during estimation of the median and MAD.
#' The corresponding entries of the output vector are also set to \code{NA} values.
#'
#' @section Handling batches:
#' If \code{batch} is specified, outliers are defined within each batch separately using batch-specific median and MAD values.
#' This gives the same results as if the input metrics were subsetted by batch and \code{isOutlier} was run on each subset,
#' and is often useful when batches are known \emph{a priori} to have technical differences (e.g., in sequencing depth).
#' 
#' If \code{share_medians=TRUE}, the median is computed across all cells.
#' If \code{shared_mads=TRUE}, the MAD is computed using all cells (from either a batch-specific or shared median).
#' These settings are useful to enforce a common location or spread across batches, 
#' e.g., we might set \code{shared_mads=TRUE} for log-library sizes if coverage varies across batches but the variance across cells is expected to be consistent across batches.
#' 
#' If \code{batch} is specified, the \code{"threshold"} attribute in the returned vector is a matrix with one named column per level of \code{batch} and two rows (one per threshold).
#' 
#' If a batch does not have sufficient cells to compute the median or MAD (e.g., after applying \code{subset}),
#' all cells in that batch will have \code{NA} in the output.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' stats <- perCellQCMetrics(example_sce)
#'
#' str(isOutlier(stats$sum))
#' str(isOutlier(stats$sum, type="lower"))
#' str(isOutlier(stats$sum, type="higher"))
#' 
#' str(isOutlier(stats$sum, log=TRUE))
#'
#' b <- sample(LETTERS[1:3], ncol(example_sce), replace=TRUE)
#' str(isOutlier(stats$sum, log=TRUE, batch=b))
#' 
#' @seealso
#' \code{\link{quickPerCellQC}}, a convenience wrapper to perform outlier-based quality control.
#'
#' \code{\link{perCellQCMetrics}}, to compute potential QC metrics.
#' @export
#' @importFrom stats mad median
isOutlier <- function(metric, nmads = 3, type = c("both", "lower", "higher"), 
    log = FALSE, subset = NULL, batch = NULL, share_medians=FALSE, 
    share_mads=FALSE, min_diff = NA) 
{
    if (log) {
        metric <- log10(metric)
    }

    N <- length(metric)
    if (nobatch <- is.null(batch)) {
        batch <- rep("1", N)
    } else {
        if (length(batch) != N) { 
            stop("length of 'batch' must equal length of 'metric'")
        }
        batch <- as.character(batch)
    }

    if (!is.null(subset)) { 
        M <- metric[subset]
        B <- batch[subset]
    } else {
        M <- metric
        B <- batch
    }
    if (any(na.drop <- is.na(M))) { 
        M <- M[!na.drop]
        B <- B[!na.drop]
        warning("missing values ignored during outlier detection")
    }
    
    # Defining the mad or median, possibly by sharing information across batches.
    by.batch <- split(M, B)

    if (share_medians) {
        cur.med <- rep(median(M), length(by.batch))
        names(cur.med) <- names(by.batch)
    } else {
        cur.med <- vapply(by.batch, median, 0)
    }

    if (share_mads) {
        cur.mad <- median(abs(M - cur.med[B])) * formals(mad)$constant
        cur.mad <- rep(cur.mad, length(by.batch))
        names(cur.mad) <- names(by.batch)
    } else {
        cur.mad <- mapply(mad, x=by.batch, center=cur.med)
    }

    # Defining the thresholds.
    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val 
    lower.limit <- cur.med - diff.val 
    
    type <- match.arg(type)
    if (type == "lower") {
        upper.limit[] <- Inf
    } else if (type == "higher") {
        lower.limit[] <- -Inf
    }

    collected <- (metric < lower.limit[batch] | upper.limit[batch] < metric)
    names(collected) <- names(metric)

    # Formatting output.
    all.threshold <- rbind(lower=lower.limit, higher=upper.limit)
    if (nobatch) {
        all.threshold <- all.threshold[,1]
    }
    .store_thresholds(collected, all.threshold, logged=log)
}

.store_thresholds <- function(x, val, logged=FALSE) {
    if (logged) val <- 10^val
    attr(x, "thresholds") <- val
    x
}

.get_thresholds <- function(x) {
    attr(x, "thresholds")
}
