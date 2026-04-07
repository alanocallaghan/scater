#' @export
scuttle::addPerCellQC

#' @export
scuttle::addPerFeatureQC

#' @export
scuttle::aggregateAcrossCells

#' @export
scuttle::aggregateAcrossFeatures

#' @export
scuttle::calculateAverage

#' @export
scuttle::calculateCPM

#' @export
scuttle::calculateFPKM

#' @export
scuttle::calculateTPM

#' @export
scuttle::computeLibraryFactors

#' @export
scuttle::computeMedianFactors

#' @export
scuttle::isOutlier

#' @export
scuttle::librarySizeFactors

#' @export
scuttle::logNormCounts

#' @export
scuttle::makePerCellDF

#' @export
scuttle::makePerFeatureDF

#' @export
scuttle::medianSizeFactors

#' @export
scuttle::mockSCE

#' @export
scuttle::normalizeCounts

#' @export
scuttle::numDetectedAcrossCells

#' @export
scuttle::numDetectedAcrossFeatures

#' @export
scuttle::perCellQCMetrics

#' @export
scuttle::perFeatureQCMetrics

#' @export
scuttle::quickPerCellQC

#' @export
scuttle::readSparseCounts

#' @export
scuttle::sumCountsAcrossCells

#' @export
scuttle::sumCountsAcrossFeatures

#' @export
scuttle::uniquifyFeatureNames

# Vendored from scuttle given that it's otherwise deprecated and we don't want to add a dependency on scrapper.
quick_means_by_group <- function(x, num.groups, group.id) {
    multiplier <- Matrix::sparseMatrix(i = seq_along(group.id), j = group.id, x = rep(1, length(group.id)))
    sums <- x %*% multiplier
    group.sizes <- tabulate(group.id, nbins=num.groups)
    t(t(as.matrix(sums)) / group.sizes)
}
