################################################################################
### updating an old SCESet object

#' @rdname toSingleCellExperiment
#' @export
#' @examples
#' \dontrun{
#' updateSCESet(example_sceset)
#' }
updateSCESet <- function(object) {
    toSingleCellExperiment(object)
}

#' Convert an SCESet object to a SingleCellExperiment object
#'
#' Convert an SCESet object produced with an older version of the
#' package to a SingleCellExperiment object compatible with the current version.
#'
#' @param object an \code{\link{SCESet}} object to be updated
#'
#' @return a \code{\link{SingleCellExperiment}} object
#' 
#' @name toSingleCellExperiment
#' @rdname toSingleCellExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @import SingleCellExperiment
#' @export
#' @examples
#' \dontrun{
#' toSingleCellExperiment(example_sceset)
#' }
toSingleCellExperiment <- function(object) {
    stopifnot(methods::is(object, "SCESet"))
    new.assay <- list()
    ass.names <- Biobase::assayDataElementNames(object)
    for (x in ass.names) {
        if (x == "exprs" || x == "logcounts") { 
            new.name <- "logcounts"
            if (new.name %in% names(new.assay)) {
                warning("'exprs' renamed to 'logcounts' will ", 
                        ifelse(x == "exprs", "overwrite existing value",
                               "be overwritten"))
            }
        } else {
            new.name <- x
        }
        new.assay[[new.name]] <- Biobase::assayDataElement(object, x)
    }
    
    new.coldata <- DataFrame(row.names = Biobase::sampleNames(object))
    pdat <- Biobase::phenoData(object)
    for (x in Biobase::varLabels(pdat)) { 
        new.coldata[[x]] <- pdat[[x]]
    }
    
    new.rowdata <- DataFrame(row.names = Biobase::featureNames(object)) 
    fdat <- Biobase::featureData(object)
    for (x in Biobase::varLabels(fdat)) { 
        new.rowdata[[x]] <- fdat[[x]]
    }

    sce <- SingleCellExperiment(new.assay, colData = new.coldata, 
                         rowData = new.rowdata)
        
    if (nrow(object@reducedDimension) > 0 && 
        ncol(object@reducedDimension) > 0) {
        new.reddim <- SimpleList(redDim = object@reducedDimension)
        reducedDims(sce) <- new.reddim
    }

    sce
}
