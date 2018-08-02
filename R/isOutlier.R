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
#' Missing values will trigger a warning and will be automatically ignored. 
#' @param batch Factor of length equal to \code{metric}, specifying the batch to which each observation belongs. 
#' A median/MAD is calculated for each batch, and outliers are then identified within each batch.
#' @param min_diff A numeric scalar indicating the minimum difference from the median to consider as an outlier. 
#' The outlier threshold is defined from the larger of \code{nmads} MADs and \code{min_diff}, to avoid calling many outliers when the MAD is very small. 
#' If \code{NA}, it is ignored.
#' 
#' @return A logical vector of the same length as the \code{metric} argument, specifying the observations that are considered as outliers.
#' 
#' @author Aaron Lun
#'
#' @export
#' @importFrom stats mad median
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
#' isOutlier(example_sce$total_counts, nmads = 3)
#' 
isOutlier <- function(metric, nmads = 5, type = c("both", "lower", "higher"), 
                      log = FALSE, subset = NULL, batch = NULL, min_diff = NA) {
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
                                   batch = NULL, min_diff = min_diff)
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

    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
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
