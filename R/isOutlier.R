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
#' @param min_diff A numeric scalar indicating the minimum difference from the median to consider as an outlier. 
#' The outlier threshold is defined from the larger of \code{nmads} MADs and \code{min_diff}, to avoid calling many outliers when the MAD is very small. 
#' Ignored if \code{NA}.
#' 
#' @return A logical vector of the same length as the \code{metric} argument, specifying the observations that are considered as outliers.
#'
#' @details
#' Lower and upper thresholds are stored in the \code{"threshold"} attribute of the returned vector.
#' This is a numeric vector of length 2 when \code{batch=NULL} for the threshold on each side.
#' Otherwise, it is a matrix with one named column per level of \code{batch} and two rows (one per threshold).
#' 
#' Missing values trigger a warning and are automatically ignored during estimation of the median and MAD.
#' The corresponding entries of the output vector are also set to \code{NA} values.
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
#' @export
#' @importFrom stats mad median
isOutlier <- function(metric, nmads = 3, type = c("both", "lower", "higher"), 
    log = FALSE, subset = NULL, batch = NULL, min_diff = NA) 
{
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
        all.threshold <- vector("list", length(by.batch))
        for (b in seq_along(by.batch)) {
            bdx <- by.batch[[b]]
            current <- Recall(metric[bdx], nmads = nmads, type = type, log = FALSE, 
                subset = subset[bdx], batch = NULL, min_diff = min_diff)
            all.threshold[[b]] <- .get_thresholds(current)
            collected[bdx] <- current
        }

        all.threshold <- do.call(cbind, all.threshold)
        colnames(all.threshold) <- names(by.batch)
        return(.store_thresholds(collected, all.threshold, logged=log))
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

    .store_thresholds( 
        (metric < lower.limit | upper.limit < metric),
        c(lower=lower.limit, higher=upper.limit),
        logged = log
    )
}

.store_thresholds <- function(x, val, logged=FALSE) {
    if (logged) val <- 10^val
    attr(x, "thresholds") <- val
    x
}

.get_thresholds <- function(x) {
    attr(x, "thresholds")
}
