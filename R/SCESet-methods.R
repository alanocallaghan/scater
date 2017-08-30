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
    for (x in varLabels(pdat)) { 
        new.coldata[[x]] <- pdat[[x]]
    }
    
    new.rowdata <- DataFrame(row.names = Biobase::featureNames(object)) 
    fdat <- Biobase::featureData(object)
    for (x in varLabels(fdat)) { 
        new.rowdata[[x]] <- fdat[[x]]
    }

    sce <- SingleCellExperiment(new.assay, colData = new.coldata, 
                         rowData = new.rowdata)
        
    if (nrow(object@reducedDimension) > 0 && 
        ncol(object@reducedDimension) > 0) {
        new.reddim <- SimpleList(redDim = reducedDim(object))
        reducedDim(sce) <- new.reddim
    }

    sce
}


################################################################################

#' Create a new SCESet object
#'
#' Deprecated from scater version 1.3.29; the package now uses the 
#' \code{\link{SingleCellExperiment}} class. To convert an SCESet object to 
#' SingleCellExperiment see the \code{\link{toSingleCellExperiment}} function. 
#' This function is retained for backwards compatibility.
#'
#' @param exprsData expression data matrix for an experiment (features x cells)
#' @param countData data matrix containing raw count expression values
#' @param tpmData matrix of class \code{"numeric"} containing
#' transcripts-per-million (TPM) expression values
#' @param fpkmData matrix of class \code{"numeric"} containing fragments per
#' kilobase of exon per million reads mapped (FPKM) expression values
#' @param cpmData matrix of class \code{"numeric"} containing counts per
#' million (CPM) expression values (optional)
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param experimentData MIAME class object containing metadata data and details
#' about the experiment and dataset.
#' @param is_exprsData matrix of class \code{"logical"}, indicating whether
#'    or not each observation is above the \code{lowerDetectionLimit}.
#' @param cellPairwiseDistances object of class \code{"dist"} (or a class that
#' extends "dist") containing cell-cell distance or dissimilarity values.
#' @param featurePairwiseDistances object of class \code{"dist"} (or a class that
#' extends "dist") containing feature-feature distance or dissimilarity values.
#' @param lowerDetectionLimit the minimum expression level that constitutes true
#'  expression (defaults to zero and uses count data to determine if an
#'  observation is expressed or not).
#' @param logExprsOffset numeric scalar, providing the offset used when doing
#' log2-transformations of expression data to avoid trying to take logs of zero.
#' Default offset value is \code{1}.
#'
#' @details
#' This function now returns a \code{\link{SingleCellExperiment}} object, 
#' whereas earlier versions produced an \code{SCESet} object. The \code{scater}
#' package now uses \code{SingleCellExperiment} as its data structure instead
#' of \code{SCESet}.
#' 
#' @importFrom Biobase pData
#' @importFrom methods as
#' @return a \code{\link{SingleCellExperiment}} object
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' \dontrun{
#' example_sce <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' }
newSCESet <- function(exprsData = NULL, countData = NULL, tpmData = NULL, 
                      fpkmData = NULL, cpmData = NULL, phenoData = NULL, featureData = NULL, 
                      experimentData = NULL, is_exprsData = NULL, cellPairwiseDistances = dist(vector()), 
                      featurePairwiseDistances = dist(vector()), lowerDetectionLimit = NULL, 
                      logExprsOffset = NULL) {
    .Deprecated(old = "newSCESet", "SingleCellExperiment")

    ass.data <- list()
    if (!is.null(exprsData)) { 
        ass.data$logcounts <- exprsData
    }
    if (!is.null(countData)) { 
        ass.data$counts <- countData
    } 
    if (!is.null(tpmData)) {
        ass.data$tpm <- tpmData
    }
    if (!is.null(fpkmData)) {
        ass.data$fpkm <- fpkmData
    }
    if (!is.null(cpmData)) {
        ass.data$cpm <- cpmData
    }
    if (length(ass.data) == 0) {
        stop("one set of expression values should be supplied")
    }

    sce <- SingleCellExperiment(ass.data, rowData = featureData, 
                                metadata = experimentData)

    if (!is.null(phenoData)) { 
        colData(sce) <- as(pData(phenoData), "DataFrame")
    }
    if (!is.null(logExprsOffset)) {
        metadata(sce)$log.exprs.offset <- logExprsOffset
    }
    return(sce)
}



################################################################################

#' Retrieve a representation of gene expression
#'
#' Deprecated from scater version 1.3.29.
#'
#' @param object An object of type \code{SCESet}
#' @return A matrix representation of expression values.
#'
getExprs <- function(object) {
    stop("Deprecated from scater 1.3.29")
}


################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}. 
#' 
#' Deprecated from scater version 1.5.2.
#'
#' @param sce An \code{SCESet} object
#' @param exprs_values What should \code{exprs(cds)} be mapped from in the
#' \code{SCESet}? Should be one of "exprs", "tpm", "fpkm", "counts"
#'
#' @export
#' @rdname toCellDataSet
#' @name toCellDataSet
#' @return An object of class \code{CellDataSet}
#' @examples
#' \dontrun{"Deprecated"}
toCellDataSet <- function(sce, exprs_values = "exprs") {
    stop("Deprecated from scater 1.5.2.")
}

#' Convert a \code{CellDataSet} to an \code{SCESet}
#' 
#' Deprecated from scater version 1.5.2.
#'
#' @param cds A \code{CellDataSet} from the \code{monocle} package
#' @param exprs_values What should \code{exprs(cds)} be mapped to in the 
#' \code{SCESet}? Should be one of "exprs", "tpm", "fpkm", "counts"
#' @param logged logical, if \code{exprs_values="exprs"}, are
#'  the expression values already on the log2 scale, or not?
#' @param logExprsOffset numeric, value to add prior to log-transformation.
#'
#' @export
#' @importFrom Biobase featureData
#' @importFrom Biobase phenoData
#' @rdname fromCellDataSet
#' @name fromCellDataSet
#' @return An object of class \code{SCESet}
#' @examples
#' \dontrun{"Deprecated"}
fromCellDataSet <- function(cds, exprs_values = "tpm", logged = FALSE, 
                            logExprsOffset = 1) {
    stop("Deprecated from scater 1.5.2.")
}

