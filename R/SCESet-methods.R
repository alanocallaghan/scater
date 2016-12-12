### Methods for the SCESet class

################################################################################
### constructor function for SCESet class

#' Create a new SCESet object.
#'
#' Create a new SCESet object (the basic data container class in scater) from a
#' supplied matrix of expression values, plus cell and feature metadata. The
#' expression matrix have rows representing features (usually genes) and columns
#' representing cells.
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
#' @return a new SCESet object
#'
#' @details
#' Scater requires that all data be housed in SCESet objects. SCESet extends
#' Bioconductor's ExpressionSet class, and the same basic interface is
#' supported. newSCESet() expects a single matrix of expression values of a
#' nominated type to be provided, for example a matrix of counts or a matrix of
#' transcripts-per-million values. There is a hierarchy applied to the
#' expression data: counts > transcripts-per-million (tpm) > counts-per-million
#' (cpm) > fragments-per-kilobase-per-million-mapped (fpkm) > generic expression
#' values on the log2 scale (exprs). Data types higher in the higher are
#' preferred. Data types lower in the hierarchy will be computed from values
#' higher in the hierarchy - e.g. counts-per-million and expression values
#' (as log2(cpm + offset)) will be computed from counts. Data types higher in
#' the hierarchy will never be computed from types lower in the hierarchy (e.g.
#' counts will never be computed from exprs values). At a
#' minimum, an SCESet object will contain exprs values; these will be computed
#' as log2(*pm + offset) values if a data type higher in the hierarchy is
#' supplied as the expression matrix.
#'
#' Per-feature and per-cell metadata can be supplied with the featureData and
#' phenoData arguments, respectively. Use of these optional arguments is strongly encouraged.
#'
#' Many methods are provided in the package that operate on SCESet objects.
#'
#' Aside from the hierarchy of data types described above, scater is relatively
#' agnostic with respect to data the nature of the expression values. Most
#' frequently used values are feature counts or transcripts-per-million (tpm),
#' but any valid output from a program that calculates expression values from
#' RNA-Seq reads is supported. For example, expression values could also be
#' values from a single cell qPCR run or some other type of assay.
#'
#' In some cases it may be desirable to have both tpm and counts in an SCESet
#' object. In such cases, expression matrices can be added to an SCESet object
#' after it has been produced by using the \code{\link{set_exprs}} function to
#' add the expression matrix to the SCESet object.
#'
#' In many downstream functions it is most convenient if the
#' \code{'exprs'} values are on the log2-scale, so this is done by default.
#'
#' @importFrom Biobase annotatedDataFrameFrom
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase assayDataNew
#' @importFrom Biobase assayDataElement
#' @importFrom edgeR cpm.default
#' @import methods
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset
newSCESet <- function(exprsData = NULL,
                      countData = NULL,
                      tpmData = NULL,
                      fpkmData = NULL,
                      cpmData = NULL,
                      phenoData = NULL,
                      featureData = NULL,
                      experimentData = NULL,
                      is_exprsData = NULL,
                      cellPairwiseDistances = dist(vector()),
                      featurePairwiseDistances = dist(vector()),
                      lowerDetectionLimit = NULL,
                      logExprsOffset = NULL)
{
    ## Checking which value to use, in the hierarchy specified.
    have.data <- NULL
    for (dataname in c("countData", "tpmData", "cpmData", "fpkmData",
                       "exprsData")) {
        eData <- get(dataname)
        if (!is.null(eData)) {
            if (!is.null(have.data)) {
                warning(sprintf("'%s' provided, '%s' will be ignored",
                                have.data, dataname))
                assign(dataname, NULL)
            } else {
                assign(dataname, as.matrix(eData))
                have.data <- dataname
            }
        }
    }

    if (is.null(have.data)) {
        stop("one set of expression values should be supplied")
    }

    if (!is.null(is_exprsData)) {
        if (have.data != "exprsData") {
            warning(sprintf("'%s' provided, 'is_exprsData' will be ignored",
                            have.data))
            is_exprsData <- NULL
        } else {
            is_exprsData <- as.matrix(is_exprsData)
        }
    }

    ## Setting logExprsOffset and lowerDetectionLimit.
    if (is.null(logExprsOffset)) {
        logExprsOffset <- 1
        if (have.data != "countData") {
            warning("'logExprsOffset' should be set manually for non-count data")
        }
    }

    if (is.null(lowerDetectionLimit)) {
        lowerDetectionLimit <- 0
        if (have.data == "exprsData") {
            warning("'lowerDetectionLimit' should be set manually for log-expression values")
        }
    }

    ## If no exprsData provided, define it from counts or T/C/FPKMs
    if (have.data == "countData") {
        exprsData <- .compute_exprs(countData, size_factors = colSums(countData),
                                    log = TRUE, sum = FALSE,
                                    logExprsOffset = logExprsOffset)
        dimnames(exprsData) <- dimnames(countData)
    } else if (have.data != "exprsData") {
        exprsData <- log2(get(have.data) + logExprsOffset)
    }

    ## Generate valid phenoData and featureData if not provided
    if ( is.null(phenoData) )
        phenoData <- annotatedDataFrameFrom(exprsData, byrow = FALSE)
    if ( is.null(featureData) )
        featureData <- annotatedDataFrameFrom(exprsData, byrow = TRUE)

    ## Check experimentData
    expData_null <- new("MIAME",
                        name = "<your name here>",
                        lab = "<your lab here>",
                        contact = "<email address>",
                        title = "<title for this dataset>",
                        abstract = "An SCESet",
                        url = "<your website here>",
                        other = list(
                            notes = "This dataset created from ...",
                            coauthors = c("")
                        ))
    if ( !is.null( experimentData ) ) {
        if ( is(experimentData, "MIAME") )
            expData <- experimentData
        else {
            expData <- expData_null
            warning("'experimentData' is not an 'MIAME' object, setting to an empty object")
        }
    } else {
        expData <- expData_null
    }

    ## Generate new SCESet object
    assaydata <- assayDataNew("lockedEnvironment", exprs = exprsData)
    sceset <- new( "SCESet",
                   assayData = assaydata,
                   phenoData = phenoData,
                   featureData = featureData,
                   experimentData = expData,
                   cellPairwiseDistances = cellPairwiseDistances,
                   featurePairwiseDistances = featurePairwiseDistances,
                   lowerDetectionLimit = lowerDetectionLimit,
                   logExprsOffset = logExprsOffset,
                   logged = TRUE,
                   featureControlInfo = AnnotatedDataFrame(),
                   useForExprs = "exprs")

    ## Add non-null slots to assayData for SCESet object, omitting null slots
    if ( !is.null(is_exprsData) )
        is_exprs(sceset) <- is_exprsData
    if ( !is.null(tpmData) )
        tpm(sceset) <- tpmData
    if ( !is.null(fpkmData) )
        fpkm(sceset) <- fpkmData
    if ( !is.null(countData) )
        counts(sceset) <- countData
    if ( !is.null(cpmData) )
        cpm(sceset) <- cpmData

    ## Check validity of object
    validObject(sceset)
    sceset
}


################################################################################
### Define validity check for SCESet class object

setValidity("SCESet", function(object) {
    msg <- NULL
    valid <- TRUE

    ## Check that the dimensions and names of the bootstraps slot are sensible
    if ( (length(object@bootstraps) != 0) && (nrow(object@bootstraps)
                                              != nrow(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Number of boostrapped genes doesn't match number of genes in SCESet")
    }
    if ( (length(object@bootstraps) != 0) && (ncol(object@bootstraps)
                                              != ncol(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Number of boostrapped samples doesn't match number of samples in SCESet")
    }
    if (  (length(object@bootstraps) != 0) &&
          !identical(rownames(object@bootstraps), featureNames(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Boostrap row names must match SCESet featureNames")
    }
    if (  (length(object@bootstraps) != 0) &&
          !identical(colnames(object@bootstraps), sampleNames(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Boostrap column names must match SCESet sampleNames")
    }
    ## Check that the dimensions of the reducedDimension slot are sensible
    if ( (nrow(object@reducedDimension) != 0) &&
         (nrow(object@reducedDimension) != ncol(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Number of samples in reducedDimension doesn't match number of samples in SCESet")
    }
    if ( (nrow(object@reducedDimension) != 0) &&
         !identical(rownames(object@reducedDimension), sampleNames(object)) ) {
        valid <- FALSE
        msg <- c(msg, "Row names of reducedDimension don't match sampleNames of SCESet")
    }
    ## Check that the dimensions of the cellPairwiseDistances slot are sensible
    if ( (length(object@cellPairwiseDistances) != 0) &&
         (length(object@cellPairwiseDistances) != choose(ncol(object), 2)) ) {
        valid <- FALSE
        msg <- c(msg, "cellPairwiseDistances must be of length compatible with ncol(SCESet)")
    }
    if ( (length(object@cellPairwiseDistances) != 0) &&
         !identical(labels(object@cellPairwiseDistances), sampleNames(object)) ) {
        valid <- FALSE
        msg <- c(msg,
                 "Labels of cellPairwiseDistances must be identical to sampleNames(SCESet)")
    }
    ## Check that the dimensions of the featurePairwiseDistances slot are sensible
    if ( (length(object@featurePairwiseDistances) != 0) &&
         (length(object@featurePairwiseDistances) != choose(ncol(object), 2)) ) {
        valid <- FALSE
        msg <- c(msg, "featurePairwiseDistances must be of length compatible with nrow(SCESet)")
    }
    if ( (length(object@featurePairwiseDistances) != 0) &&
          !identical(labels(object@featurePairwiseDistances), featureNames(object)) ) {
        valid <- FALSE
        msg <- c(msg,
                 "Label names of featurePairwiseDistances must be identical to featureNames(SCESet)")
    }
    ## Check that we have sensible values for the counts
    if( .checkedCall(cxx_missing_exprs, exprs(object)) ) {
        warning( "The exprs data contain NA values." )
    }
    if ( (!is.null(counts(object))) && .checkedCall(cxx_negative_counts, counts(object)) )
        warning( "The count data contain negative values." )
    if ( !(object@useForExprs %in% c("exprs", "tpm", "fpkm", "counts")) ) {
        valid <- FALSE
        msg <- c(msg, "object@useForExprs must be one of 'exprs', 'tpm', 'fpkm', 'counts'")
    }

    if (valid) TRUE else msg
})


################################################################################
### updating an old SCESet object

#' Update an SCESet object to the current version
#'
#' It can be necessary to update an SCESet produced with an older version of the
#' package to be compatible with the current version of the package.
#'
#' @param object an \code{\link{SCESet}} object to be updated
#'
#' @return an updated \code{\link{SCESet}} object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' updateSCESet(example_sceset)
updateSCESet <- function(object) {
    if (!is(object, "SCESet"))
        stop("Object must be an SCESet.")
    sceset <- new( "SCESet",
                   assayData = object@assayData,
                   phenoData = object@phenoData,
                   featureData = object@featureData,
                   experimentData = object@experimentData,
                   cellPairwiseDistances = as.dist(object@cellPairwiseDistances),
                   featurePairwiseDistances = as.dist(object@featurePairwiseDistances),
                   consensus = ifelse(.hasSlot(object, "consensus"),
                                      object@consensus, list()),
                   lowerDetectionLimit = object@lowerDetectionLimit,
                   logExprsOffset = object@logExprsOffset,
                   logged = object@logged,
                   featureControlInfo = object@featureControlInfo,
                   useForExprs = object@useForExprs)

    ## Check validity of object
    validObject(sceset)
    sceset
}


################################################################################
### subsetting an SCESet object

#' Subsetting SCESet Objects
#'
#' Subset method for SCESet objects, which subsets both the expression data,
#' phenotype data, feature data and other slots in the object.
#'
#' @return an SCESet object
#' @rdname SCESet-subset
#' @name SCESet-subset
NULL

#' @inheritParams base::Extract
#' @param i,j,... indices specifying elements to extract or replace. Indices
#' are numeric or character vectors or empty (missing) or \code{NULL}. Numeric
#' values are coerced to integer as by \code{\link[base]{as.integer}} (and hence
#' truncated towards zero). Character vectors will be matched to the names of
#' the object (or for matrices/arrays, the dimnames): see
#' \code{\link[base]{Extract}} for further details.
#'
#' For \code{[}-indexing only: \code{i, j, ...} can be logical vectors, indicating
#' elements/slices to select. Such vectors are recycled if necessary to match
#' the corresponding extent. \code{i, j, ...} can also be negative integers,
#' indicating elements/slices to leave out of the selection. When indexing
#' arrays by \code{[} a single argument i can be a matrix with as many columns
#' as there are dimensions of \code{x}; the result is then a vector with
#' elements corresponding to the sets of indices in each row of \code{i}. An
#' index value of \code{NULL} is treated as if it were \code{integer(0)}.
#'
#' @aliases [,SCESet,ANY-method [,SCESet,ANY,ANY-method [,SCESet,ANY,ANY,ANY-method
#' @rdname SCESet-subset
#' @export
#' @seealso \code{\link[base]{Extract}}
#'
setMethod('[', 'SCESet', function(x, i, j, ..., drop=FALSE) {
    if ( !missing(i) && missing(j) ) {
        ## Subsetting features only
        x <- selectMethod('[', 'ExpressionSet')(x, i, , drop = drop)
        if ( length(x@featurePairwiseDistances) != 0 )
            x@featurePairwiseDistances <-
                as.dist(as.matrix(x@featurePairwiseDistances)[i, i, drop = drop])
        if ( !missing(...) ) {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if (  nrow(x@bootstraps) != 0 && length(dim(x@bootstraps)) > 2 )
                    x@bootstraps <- x@bootstraps[i, , ..., drop = drop]
            }
        } else {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if ( nrow(x@bootstraps) != 0 )
                    x@bootstraps <- x@bootstraps[i, , , drop = drop]
            }
        }
    } else if ( missing(i) && !missing(j) ) {
        ## Subsetting cells only
        x <- selectMethod('[', 'ExpressionSet')(x, , j, drop = drop)
        if ( length(x@cellPairwiseDistances) != 0 )
            x@cellPairwiseDistances <-
                as.dist(as.matrix(x@cellPairwiseDistances)[j, j, drop = drop])
        if ( nrow(x@reducedDimension) != 0 )
            x@reducedDimension <- x@reducedDimension[j, , drop = drop]
        if ( !missing(...) ) {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if ( ncol(x@bootstraps) != 0 && length(dim(x@bootstraps)) > 2 )
                    x@bootstraps <- x@bootstraps[, j, ..., drop = drop]
            }
        } else {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if ( ncol(x@bootstraps) != 0 )
                    x@bootstraps <- x@bootstraps[, j, , drop = drop]
            }
        }
    } else if ( !missing(i) && !missing(j) ) {
        ## Subsetting features (i) and cells (j)
        x <- selectMethod('[', 'ExpressionSet')(x, i, j, drop = drop)
        if ( length(x@featurePairwiseDistances) != 0 )
            x@featurePairwiseDistances <-
                as.dist(as.matrix(x@featurePairwiseDistances)[i, i, drop = drop])
        if ( length(x@cellPairwiseDistances) != 0 )
            x@cellPairwiseDistances <-
                as.dist(as.matrix(x@cellPairwiseDistances)[j, j, drop = drop])
        if ( nrow(x@reducedDimension) != 0 )
            x@reducedDimension <- x@reducedDimension[j, , drop = drop]
        if ( !missing(...) ) {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if ( nrow(x@bootstraps) != 0 && ncol(x@bootstraps) != 0 && length(dim(x@bootstraps)) > 2 )
                    x@bootstraps <- x@bootstraps[i, j, ..., drop = drop]
            }
        } else {
            if ( !is.na(ncol(x@bootstraps)) ) {
                if ( nrow(x@bootstraps) != 0 && ncol(x@bootstraps) != 0 )
                    x@bootstraps <- x@bootstraps[i, j, , drop = drop]
            }
        }
    } else{
        ## All missing: possibly not missing ... for subsetting bootstrap samples
        if ( !is.na(ncol(x@bootstraps)) ) {
            if ( nrow(x@bootstraps) != 0 && ncol(x@bootstraps) != 0 && length(dim(x@bootstraps)) > 2 )
                x@bootstraps <- x@bootstraps[, , ..., drop = drop]
        }
    }
    ## Check validity of object
    validObject(x)
    x
})


################################################################################
## cellNames

#' Get or set cell names from an SCESet object
#'
#' @param object An \code{\link{SCESet}} object.
#'
#' @return A vector of cell names.
#'
#' @details Simply a wrapper to \code{\link[Biobase]{sampleNames}}.
#' @importFrom Biobase sampleNames
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' cellNames(example_sceset)
#'
cellNames <- function(object) {
    sampleNames(object)
}


#' @usage
#' \S4method{cellNames}{SCESet,vector}(object)<-value
#'
#' @docType methods
#' @name cellNames
#' @rdname cellNames
#' @aliases cellNames cellNames<-,SCESet,vector-method
#'
#' @param value a vector of cell names to apply to the \code{SCESet} object.
#' @author Davis McCarthy
#'
#' @exportMethod "cellNames<-"
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' cellNames(example_sceset) <- 1:ncol(example_sceset)
#'
setReplaceMethod("cellNames", signature(object = "SCESet", value = "vector"),
                 function(object, value) {
                     sampleNames(object) <- value
                     validObject(object)
                     object
                 })


################################################################################
### Replacer methods for slots in an SCESet object

#' Replaces featureData in an SCESet object
#'
#' SCESet objects contain feature information (inherited from the ExpressionSet
#' class). This function conveniently replaces the feature data with the
#' value supplied, which must be an AnnotatedDataFrame.
#' @param object An SCESet object.
#' @param value an AnnotatedDataFrame with updated featureData to replace
#' existing
#' @return A matrix of expression count data, where rows correspond to features
#' (e.g. genes) and columns correspond to cells.
#'
#' @importFrom Biobase fData<-
#' @export
#' @rdname fData
#' @aliases fData fData,SCESet-method fData<-,SCESet,AnnotatedDataFrame-method fData<-,SCESet,data.frame-method
#'
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' fData(example_sceset)
#' }
setReplaceMethod("fData", signature(object = "SCESet", value = "AnnotatedDataFrame"),
                 function(object, value) {
                     object@featureData <- value
                     object
                 } )

#' @name fData
#' @rdname fData
#' @exportMethod "fData<-"
setReplaceMethod("fData", signature(object = "SCESet", value = "data.frame"),
                 function(object, value) {
                     object@featureData <- new("AnnotatedDataFrame", value)
                     object
                 } )



#' Replaces phenoData in an SCESet object
#'
#' SCESet objects contain phenotype information (inherited from the
#' ExpressionSet class). This function conveniently replaces the phenotype data
#' with the value supplied, which must be an AnnotatedDataFrame.
#' @param object An SCESet object.
#' @param value an AnnotatedDataFrame with updated phenoData to replace
#' existing
#' @return A matrix of expression count data, where rows correspond to features
#' (e.g. genes) and columns correspond to cells.
#'
#' @importFrom Biobase pData<-
#' @exportMethod "pData<-"
#' @rdname pData
#' @aliases pData pData,SCESet-method pData<-,SCESet,AnnotatedDataFrame-method pData<-,SCESet,data.frame-method
#'
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' pData(example_sceset)
#' }
setReplaceMethod("pData", signature(object = "SCESet", value = "AnnotatedDataFrame"),
                 function(object, value) {
                     object@phenoData <- value
                     object
                 } )


#' @name pData
#' @rdname pData
#' @exportMethod "pData<-"
setReplaceMethod("pData", signature(object = "SCESet", value = "data.frame"),
                 function(object, value) {
                     object@phenoData <- new("AnnotatedDataFrame", value)
                     object
                 } )


################################################################################
### get_exprs

#' Generic accessor for expression data from an SCESet object.
#'
#' Access by name a matrix of expression values, one row for each feature (gene,
#' exon, region, etc), and one column for each cell stored an element of the
#' assayData slot of the SCESet object.
#'
#' @usage
#' \S4method{get_exprs}{SCESet}(object, exprs_values, warning = TRUE)
#'
#' @docType methods
#' @name get_exprs
#' @rdname get_exprs
#' @aliases get_exprs get_exprs,SCESet-method
#'
#' @param object a \code{SCESet} object.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (default; transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param warning a logical scalar specifying whether a warning should be
#' raised, and \code{NULL} returned, if the requested expression values are
#' not present in \code{object}. Otherwise, an error will be thrown.
#' @param ... further arguments passed to \code{get_exprs.SCESet}
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' get_exprs(example_sceset, "counts")
#'
#' ## new slots can be defined and accessed
#' set_exprs(example_sceset, "scaled_counts") <- t(t(counts(example_sceset)) /
#' colSums(counts(example_sceset)))
#' get_exprs(example_sceset, "scaled_counts")[1:6, 1:6]
#'
get_exprs.SCESet <- function(object, exprs_values = "exprs", warning = TRUE) {
    exprs_mat <- object@assayData[[exprs_values]]

    if ( is.null(exprs_mat) ) {
        msg <- sprintf("'object' does not contain '%s' values", exprs_values)
        if (warning) {
            warning(paste0(msg, ", returning NULL"))
        } else {
            stop(msg)
        }
    }

    exprs_mat
}

#' @rdname get_exprs
#' @export
setMethod("get_exprs", signature(object = "SCESet"),
          function(object, exprs_values = "exprs", warning = TRUE) {
              get_exprs.SCESet(object, exprs_values, warning = warning)
          })


#' Assignment method for the new elements of an SCESet object.
#'
#' The assayData slot of an SCESet object holds the expression data matrices.
#' This functions makes it convenient to add new transformations of the
#' expression data to the assayData slot.
#'
#' @usage
#' \S4method{set_exprs}{SCESet,ANY,matrix}(object,name)<-value
#'
#' @docType methods
#' @name set_exprs
#' @rdname set_exprs
#' @aliases set_exprs set_exprs<-,SCESet,ANY,matrix-method set_exprs<-,SCESet,ANY,NULL-method
#'
#' @param object a \code{SCESet} object.
#' @param name character string giving the name of the slot to which the data
#' matrix is to be assigned (can already exist or not).
#' @param value a numeric or integer matrix matching the dimensions of the other
#' elements of the \code{assayData} slot of the \code{SCESet} object.
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' set_exprs(example_sceset, "scaled_counts") <- t(t(counts(example_sceset)) /
#' colSums(counts(example_sceset)))
#' get_exprs(example_sceset, "scaled_counts")[1:6, 1:6]
#'
#' ## get rid of scaled counts
#' set_exprs(example_sceset, "scaled_counts") <- NULL
#'
#' @name set_exprs
#' @rdname set_exprs
#' @exportMethod "set_exprs<-"
setReplaceMethod("set_exprs", signature(object = "SCESet", value = "matrix"),
                 function(object, name, value) {
                     if (!identical(dimnames(object), dimnames(value)))
                         stop("dimnames for new expression matrix must be
                              identical to dimnames for object")
                     Biobase::assayDataElement(object, name) <- value
                     validObject(object)
                     object
                 })

#' @name set_exprs
#' @rdname set_exprs
#' @exportMethod "set_exprs<-"
setReplaceMethod("set_exprs", signature(object = "SCESet", value = "NULL"),
                 function(object, name, value) {
                     Biobase::assayDataElement(object, name) <- value
                     validObject(object)
                     object
                 })


################################################################################
### counts

#' Accessors for the 'counts' element of an SCESet object.
#'
#' The counts element holds the count data as a matrix of non-negative integer
#' count values, one row for each feature (gene, exon, region, etc), and one
#' column for each cell. It is an element of the assayData slot of the SCESet
#' object.
#'
#' @usage
#' \S4method{counts}{SCESet}(object)
#'
#' \S4method{counts}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @importFrom BiocGenerics counts
#' @aliases counts counts,SCESet-method counts<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' counts(example_sceset)
#'
counts.SCESet <- function(object) {
    object@assayData$counts
}

#' @rdname counts
#' @export
setMethod("counts", signature(object = "SCESet"), counts.SCESet)

#' @name counts
#' @rdname counts
#' @importFrom BiocGenerics counts<-
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "counts") <- value
                     validObject(object)
                     object
                 })

################################################################################
### norm_counts

#' Accessors for the 'norm_counts' element of an SCESet object.
#'
#' The norm_counts element holds normalised count data as a matrix of
#' non-negative values, one row for each feature (gene, exon, region, etc), and one
#' column for each cell. It is an element of the assayData slot of the SCESet
#' object.
#'
#' @usage
#' \S4method{norm_counts}{SCESet}(object)
#'
#' \S4method{norm_counts}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_counts
#' @rdname norm_counts
#' @aliases norm_counts norm_counts,SCESet-method norm_counts<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' norm_counts(example_sceset)
#'
norm_counts.SCESet <- function(object) {
    object@assayData$norm_counts
}

#' @rdname norm_counts
#' @export
setMethod("norm_counts", signature(object = "SCESet"), norm_counts.SCESet)

#' @name norm_counts
#' @rdname norm_counts
#' @exportMethod "norm_counts<-"
setReplaceMethod("norm_counts", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "norm_counts") <- value
                     validObject(object)
                     object
                 })

################################################################################
### is_exprs

#' Accessors for the 'is_exprs' element of an SCESet object.
#'
#' The is_exprs element holds a logical matrix indicating whether or not each
#' observation is above the defined lowerDetectionLimit in the SCESet object. It
#' has the same dimensions as the 'exprs' and 'counts' elements, which hold the
#' transformed expression data and count data, respectively.
#'
#' @usage
#' \S4method{is_exprs}{SCESet}(object)
#'
#' \S4method{is_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name is_exprs
#' @rdname is_exprs
#' @aliases is_exprs is_exprs,SCESet-method is_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' is_exprs(example_sceset)
#'
is_exprs.SCESet <- function(object) {
    object@assayData$is_exprs
}

#' @rdname is_exprs
#' @export
setMethod("is_exprs", signature(object = "SCESet"), is_exprs.SCESet)

#' @name is_exprs<-
#' @rdname is_exprs
#' @exportMethod "is_exprs<-"
setReplaceMethod("is_exprs", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "is_exprs") <- value
                     validObject(object)
                     object
                 })

################################################################################
### norm_exprs

#' Accessors for the 'norm_exprs' (normalised expression) element of an SCESet object.
#'
#' The \code{norm_exprs} element of the arrayData slot in an SCESet object holds
#' a matrix containing normalised expression values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{norm_exprs}{SCESet}(object)
#'
#' \S4method{norm_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_exprs
#' @rdname norm_exprs
#' @aliases norm_exprs norm_exprs,SCESet-method norm_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#'
#' @details The default for normalised expression values is mean-centred and
#' variance-standardised expression data from the \code{exprs} slot of the
#' \code{SCESet} object. The function \code{normaliseExprs} (or
#' \code{normalizeExprs}) provides more options and functionality for
#' normalising expression data.
#'
#' @author Davis McCarthy
#' @export
#' @aliases norm_exprs norm_exprs,SCESet-method, norm_exprs<-,SCESet,matrix-method
#'
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' norm_exprs(example_sceset)
#'
norm_exprs.SCESet <- function(object) {
    object@assayData$norm_exprs
}

#' @rdname norm_exprs
#' @export
setMethod("norm_exprs", signature(object = "SCESet"), norm_exprs.SCESet)

#' @name norm_exprs<-
#' @rdname norm_exprs
#' @exportMethod "norm_exprs<-"
setReplaceMethod("norm_exprs", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "norm_exprs") <- value
                     validObject(object)
                     object
                 })


################################################################################
### stand_exprs

#' Accessors for the 'stand_exprs' (standardised expression) element of an
#' SCESet object.
#'
#' The \code{stand_exprs} element of the arrayData slot in an SCESet object holds
#' a matrix containing standardised (mean-centred, variance standardised, by
#' feature) expression values. It has the same dimensions as the 'exprs' and
#' 'counts' elements, which hold the transformed expression data and count data,
#'  respectively.
#'
#' @usage
#' \S4method{stand_exprs}{SCESet}(object)
#'
#' \S4method{stand_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name stand_exprs
#' @rdname stand_exprs
#' @aliases stand_exprs stand_exprs,SCESet-method stand_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#'
#' @details The default for normalised expression values is mean-centred and
#' variance-standardised expression data from the \code{exprs} slot of the
#' \code{SCESet} object. The function \code{normaliseExprs} (or
#' \code{normalizeExprs}) provides more options and functionality for
#' normalising expression data.
#'
#' @author Davis McCarthy
#' @export
#' @aliases stand_exprs stand_exprs,SCESet-method, stand_exprs<-,SCESet,matrix-method
#'
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' stand_exprs(example_sceset)
#'
stand_exprs.SCESet <- function(object) {
    object@assayData$stand_exprs
}

#' @rdname stand_exprs
#' @export
setMethod("stand_exprs", signature(object = "SCESet"), stand_exprs.SCESet)

#' @name stand_exprs<-
#' @rdname stand_exprs
#' @exportMethod "stand_exprs<-"
setReplaceMethod("stand_exprs", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "stand_exprs") <- value
                     validObject(object)
                     object
                 })



################################################################################
### tpm

#' Accessors for the 'tpm' (transcripts per million) element of an SCESet object.
#'
#' The \code{tpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing transcripts-per-million values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{tpm}{SCESet}(object)
#'
#' \S4method{tpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name tpm
#' @rdname tpm
#' @aliases tpm tpm,SCESet-method tpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#' @aliases tpm tpm,SCESet-method tpm<-,SCESet,matrix-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' tpm(example_sceset)
#'
tpm.SCESet <- function(object) {
    object@assayData$tpm
}

#' @name tpm
#' @rdname tpm
#' @export
#' @aliases tpm,SCESet-method
setMethod("tpm", signature(object = "SCESet"), tpm.SCESet)

#' @name tpm<-
#' @rdname tpm
#' @exportMethod "tpm<-"
#' @aliases tpm<-,SCESet,matrix-method
setReplaceMethod("tpm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "tpm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### norm_tpm

#' Accessors for the 'norm_tpm' (transcripts per million) element of an SCESet object.
#'
#' The \code{norm_tpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing normalised transcripts-per-million values. It has the
#' same dimensions as the 'exprs' and 'counts' elements, which hold the
#' transformed expression data and count data, respectively.
#'
#' @usage
#' \S4method{norm_tpm}{SCESet}(object)
#'
#' \S4method{norm_tpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_tpm
#' @rdname norm_tpm
#' @aliases norm_tpm norm_tpm,SCESet-method norm_tpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#' @aliases norm_tpm norm_tpm,SCESet-method norm_tpm<-,SCESet,matrix-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' norm_tpm(example_sceset)
#'
norm_tpm.SCESet <- function(object) {
    object@assayData$norm_tpm
}

#' @name norm_tpm
#' @rdname norm_tpm
#' @export
#' @aliases norm_tpm,SCESet-method
setMethod("norm_tpm", signature(object = "SCESet"), norm_tpm.SCESet)

#' @name norm_tpm<-
#' @rdname norm_tpm
#' @exportMethod "norm_tpm<-"
#' @aliases norm_tpm<-,SCESet,matrix-method
setReplaceMethod("norm_tpm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "norm_tpm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### cpm

#' Accessors for the 'cpm' (counts per million) element of an SCESet object.
#'
#' The \code{cpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing counts-per-million values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{cpm}{SCESet}(object)
#'
#' \S4method{cpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name cpm
#' @rdname cpm
#' @aliases cpm cpm,SCESet-method cpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#' @aliases cpm cpm,SCESet-method cpm<-,SCESet,matrix-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' cpm(example_sceset)[1:10, 1:6]
cpmSCESet <- function(object) {
    object@assayData$cpm
}

#' @name cpm
#' @rdname cpm
#' @export
#' @aliases cpm,SCESet-method
setMethod("cpm", signature(object = "SCESet"), cpmSCESet)

#' @name cpm<-
#' @rdname cpm
#' @exportMethod "cpm<-"
#' @aliases cpm<-,SCESet,matrix-method
setReplaceMethod("cpm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "cpm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### norm_cpm

#' Accessors for the 'norm_cpm' (normalised counts per million) element of an SCESet object.
#'
#' The \code{norm_cpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing normalised counts-per-million values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{norm_cpm}{SCESet}(object)
#'
#' \S4method{norm_cpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_cpm
#' @rdname norm_cpm
#' @aliases norm_cpm norm_cpm,SCESet-method norm_cpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#' @aliases norm_cpm norm_cpm,SCESet-method norm_cpm<-,SCESet,matrix-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' norm_cpm(example_sceset)
#'
norm_cpm.SCESet <- function(object) {
    object@assayData$norm_cpm
}

#' @name norm_cpm
#' @rdname norm_cpm
#' @export
#' @aliases norm_cpm,SCESet-method
setMethod("norm_cpm", signature(object = "SCESet"), norm_cpm.SCESet)

#' @name norm_cpm<-
#' @rdname norm_cpm
#' @exportMethod "norm_cpm<-"
#' @aliases norm_cpm<-,SCESet,matrix-method
setReplaceMethod("norm_cpm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "norm_cpm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### fpkm

#' Accessors for the 'fpkm' (fragments per kilobase of exon per million reads mapped) element of an SCESet object.
#'
#' The \code{fpkm} element of the arrayData slot in an SCESet object holds
#' a matrix containing fragments per kilobase of exon per million reads mapped
#' (FPKM) values. It has the same dimensions as the 'exprs' and 'counts'
#' elements, which hold the transformed expression data and count data,
#' respectively.
#'
#' @usage
#' \S4method{fpkm}{SCESet}(object)
#'
#' \S4method{fpkm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name fpkm
#' @rdname fpkm
#' @aliases fpkm fpkm,SCESet-method fpkm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' fpkm(example_sceset)
#'
fpkm.SCESet <- function(object) {
    object@assayData$fpkm
}

#' @name fpkm
#' @rdname fpkm
#' @export
#' @aliases fpkm,SCESet-method
setMethod("fpkm", signature(object = "SCESet"), fpkm.SCESet)

#' @name fpkm<-
#' @rdname fpkm
#' @exportMethod "fpkm<-"
#' @aliases fpkm<-,SCESet,matrix-method
setReplaceMethod("fpkm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "fpkm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### norm_fpkm

#' Accessors for the 'norm_fpkm' (normalised fragments per kilobase of exon per million reads mapped) element of an SCESet object.
#'
#' The \code{norm_fpkm} element of the arrayData slot in an SCESet object holds
#' a matrix containing normalised fragments per kilobase of exon per million
#' reads mapped (FPKM) values. It has the same dimensions as the 'exprs' and
#' 'counts' elements, which hold the transformed expression data and count data,
#' respectively.
#'
#' @usage
#' \S4method{norm_fpkm}{SCESet}(object)
#'
#' \S4method{norm_fpkm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_fpkm
#' @rdname norm_fpkm
#' @aliases norm_fpkm norm_fpkm,SCESet-method norm_fpkm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Davis McCarthy
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' norm_fpkm(example_sceset)
#'
norm_fpkm.SCESet <- function(object) {
    object@assayData$norm_fpkm
}

#' @name norm_fpkm
#' @rdname norm_fpkm
#' @export
#' @aliases norm_fpkm,SCESet-method
setMethod("norm_fpkm", signature(object = "SCESet"), norm_fpkm.SCESet)

#' @name norm_fpkm<-
#' @rdname norm_fpkm
#' @exportMethod "norm_fpkm<-"
#' @aliases norm_fpkm<-,SCESet,matrix-method
setReplaceMethod("norm_fpkm", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     Biobase::assayDataElement(object, "norm_fpkm") <- value
                     validObject(object)
                     object
                 })


################################################################################
### sizeFactors

#' Accessors size factors of an SCESet object.
#'
#' For normalisation, library-specific size factors can be defined. Raw values
#' can be divided by the appropriate size factors to obtain normalised counts,
#' TPM, etc.
#'
#'
#' @usage
#' \S4method{sizeFactors}{SCESet}(object,type)
#'
#' \S4method{sizeFactors}{SCESet,numeric}(object,type)<-value
#' \S4method{sizeFactors}{SCESet,NULL}(object,type)<-value
#'
#' @docType methods
#' @name sizeFactors
#' @rdname sizeFactors
#' @aliases sizeFactors sizeFactors,SCESet-method sizeFactors<-,SCESet,numeric-method sizeFactors<-,SCESet,NULL-method
#'
#' @param object a \code{SCESet} object.
#' @param value a vector of class \code{"numeric"} or \code{NULL}
#' @param type optional character argument providing the type or name of the
#' size factors to be accessed or assigned.
#' @param ... further arguments passed to the function
#'
#' @details The size factors can alternatively be directly accessed from the
#' \code{SCESet} object with \code{object$size_factor_type} (where "type" in the
#' preceding is replaced by the actual type name).
#'
#' @author Davis McCarthy and Aaron Lun
#' @export
#'
#' @importFrom BiocGenerics sizeFactors
#' @importFrom BiocGenerics sizeFactors<-
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' sizeFactors(example_sceset)
#' sizeFactors(example_sceset, NULL) <- 2 ^ rnorm(ncol(example_sceset))
#'
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                                      feature_controls = list(set1 = 1:40))
#' sizeFactors(example_sceset, "set1") <- 2 ^ rnorm(ncol(example_sceset))
#' sizeFactors(example_sceset)
#'
sizeFactors.SCESet <- function(object, type=NULL) {
    ofield <- .construct_sf_field(object, type)
    out <- pData(object)[[ofield]]
    if ( is.null(out) ) {
        wstring <- "'sizeFactors' have not been set"
        if (!is.null(type)) {
            wstring <- paste0(wstring, " for '", type, "'")
        }
        warning(wstring)
        return(NULL)
    }
    names(out) <- colnames(object)
    return(out)
}

.construct_sf_field <- function(object, type) {
    ofield <- "size_factor"
    if (!is.null(type)) {
        fc_available <- .fcontrol_names(object)
        if (length(fc_available) == 0L) {
            stop("no named controls specified in the SCESet object")
        }
        type <- match.arg(type, fc_available)
        ofield <- paste0(ofield, "_", type)
    }
    return(ofield)
}

#' @name sizeFactors
#' @rdname sizeFactors
#' @export
#' @aliases sizeFactors,SCESet-method
setMethod("sizeFactors", signature(object = "SCESet"), sizeFactors.SCESet)

#' @name sizeFactors<-
#' @rdname sizeFactors
#' @exportMethod "sizeFactors<-"
#' @aliases sizeFactors<-,SCESet,numeric-method
setReplaceMethod("sizeFactors", signature(object = "SCESet", value = "numeric"),
                 function(object, type = NULL, ..., value) {
                     ofield <- .construct_sf_field(object, type)
                     pData(object)[[ofield]] <- value
                     validObject(object)
                     object
                 })

#' @name sizeFactors<-
#' @rdname sizeFactors
#' @exportMethod "sizeFactors<-"
#' @aliases sizeFactors<-,SCESet,NULL-method
setReplaceMethod("sizeFactors", signature(object = "SCESet", value = "NULL"),
                 function(object, type = NULL, ..., value) {
                     ofield <- .construct_sf_field(object, type)
                     pData(object)[[ofield]] <- NULL
                     validObject(object)
                     object
                 })


################################################################################
### bootstraps

#' Accessor and replacement for bootstrap results in an SCESet object
#'
#' SCESet objects can contain an of bootstrap expression values (for example, as
#' generated by the kallisto software for quantifying feature abundance). These
#'  functions conveniently access and replace the 'bootstrap' slot with the value
#'  supplied, which must be an matrix of the correct size, namely the same
#'  number of rows and columns as the \code{SCEset} object as a whole.
#'
#' @docType methods
#' @name bootstraps
#' @rdname bootstraps
#' @aliases bootstraps bootstraps,SCESet-method bootstraps<-,SCESet,array-method
#'
#' @param object a \code{SCESet} object.
#' @param value an array of class \code{"numeric"} containing bootstrap
#' expression values
#' @author Davis McCarthy
#'
#' @return If accessing bootstraps slot of an \code{SCESet}, then an array with
#' the bootstrap values, otherwise an \code{SCESet} object containing new
#' bootstrap values.
#'
#' @export
#' @aliases bootstraps bootstraps,SCESet-method bootstraps<-,SCE-Set,array-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' bootstraps(example_sceset)
#'
bootstraps.SCESet <- function(object) {
    object@bootstraps
}

#' @rdname bootstraps
#' @aliases bootstraps
#' @export
setMethod("bootstraps", signature(object = "SCESet"), bootstraps.SCESet)


#' @name bootstraps<-
#' @aliases bootstraps
#' @rdname bootstraps
#' @export "bootstraps<-"
setReplaceMethod("bootstraps", signature(object = "SCESet", value = "array"),
                 function(object, value) {
                     if ( (nrow(value) == nrow(object)) &&
                          (ncol(value) == ncol(object)) ) {
                         object@bootstraps <- value
                         return(object)
                     } else
                         stop("Array supplied is of incorrect size.")
                 } )


################################################################################
### reducedDimension

#' Reduced dimension representation for cells in an SCESet object
#'
#' SCESet objects can contain a matrix of reduced dimension coordinates for
#' cells. These functions conveniently access and replace the reduced dimension
#' coordinates with the value supplied, which must be a matrix of the correct
#' size. The function \code{redDim} is simply shorthand for
#' \code{reducedDimension}.
#'
#' @docType methods
#' @name reducedDimension
#' @rdname reducedDimension
#' @aliases reducedDimension reducedDimension,SCESet-method reducedDimension<-,SCESet,matrix-method redDim,SCESet-method redDim<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing reduced dimension
#' coordinates for cells.
#' @author Davis McCarthy
#'
#' @return If accessing the \code{reducedDimension} slot, then the matrix of
#' reduced dimension coordinates. If replacing the \code{reducedDimension} slot
#' then the new matrix is added to the \code{SCESet} object.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' reducedDimension(example_sceset)
#'
reducedDimension.SCESet <- function(object) {
    object@reducedDimension
}

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
setMethod("reducedDimension", signature(object = "SCESet"),
          reducedDimension.SCESet)

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
redDim.SCESet <- function(object) {
    object@reducedDimension
}

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
setMethod("redDim", signature(object = "SCESet"), redDim.SCESet)

#' @name reducedDimension<-
#' @aliases reducedDimension
#' @rdname reducedDimension
#' @exportMethod "reducedDimension<-"
setReplaceMethod("reducedDimension", signature(object = "SCESet",
                                               value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == ncol(object) ) {
                         if (is.null(rownames(value))) {
                             rownames(value) <- sampleNames(object)
                         } else {
                             if (!identical(rownames(value), sampleNames(object)))
                                 stop("Rownames of reduced dimension must be NULL or equal to sampleNames(SCESet)")
                         }
                         object@reducedDimension <- value
                         return(object)
                     }
                     else
                         stop("Reduced dimension matrix supplied is of incorrect size.
                              Rows of reduced dimension matrix should correspond to cells, i.e. columns of SCESet object.")
                 } )

#' @name redDim<-
#' @aliases reducedDimension
#' @rdname reducedDimension
#' @exportMethod "redDim<-"
setReplaceMethod("redDim", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == ncol(object) ) {
                         if (is.null(rownames(value))) {
                             rownames(value) <- sampleNames(object)
                         } else {
                             if (!identical(rownames(value), sampleNames(object)))
                                 stop("Rownames of reduced dimension must be NULL or equal to sampleNames(SCESet)")
                         }
                         object@reducedDimension <- value
                         return(object)
                     }
                     else
                         stop("Reduced dimension matrix supplied is of incorrect size.
                              Rows of reduced dimension matrix should correspond to cells, i.e. columns of SCESet object.")
                 } )



################################################################################
### cellPairwiseDistances

#' cellPairwiseDistances in an SCESet object
#'
#' SCESet objects can contain a matrix of pairwise distances between cells. These
#'  functions conveniently access and replace the cell pairwise distances with the value
#'  supplied, which must be a matrix of the correct size. The function \code{cellDist}
#'  is simply shorthand for \code{cellPairwiseDistances}.
#'
#' @docType methods
#' @name cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances cellPairwiseDistances,SCESet-method cellPairwiseDistances<-,SCESet,matrix-method cellPairwiseDistances<-,SCESet,dist-method cellDist,SCESet-method cellDist<-,SCESet,matrix-method cellDist<-,SCESet,dist-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing cell pairwise
#' distances
#' @author Davis McCarthy
#'
#' @return An SCESet object containing new cell pairwise distances matrix.
#'
#' @importFrom stats as.dist dist
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' cellPairwiseDistances(example_sceset)
#'
cellPairwiseDistances.SCESet <- function(object) {
    object@cellPairwiseDistances
}

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
setMethod("cellPairwiseDistances", signature(object = "SCESet"),
          cellPairwiseDistances.SCESet)

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
cellDistSCESet <- function(object) {
    object@cellPairwiseDistances
}

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
setMethod("cellDist", signature(object = "SCESet"), cellDistSCESet)

#' @name cellPairwiseDistances<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellPairwiseDistances<-"
setReplaceMethod("cellPairwiseDistances", signature(object = "SCESet",
                                                    value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == ncol(object) ) {
                         object@cellPairwiseDistances <- as.dist(value)
                         return(object)
                     }
                     else
                         stop("Cell pairwise distance matrix supplied is of incorrect size.")
                 } )


#' @name cellPairwiseDistances<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellPairwiseDistances<-"
setReplaceMethod("cellPairwiseDistances", signature(object = "SCESet", value = "dist"),
                 function(object, value) {
                     if ( length(value) == choose(ncol(object), 2) ) {
                         object@cellPairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Cell pairwise dist object supplied is of incorrect size.")
                 } )


#' @name cellDist<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellDist<-"
setReplaceMethod("cellDist", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == ncol(object) ) {
                         object@cellPairwiseDistances <- as.dist(value)
                         return(object)
                     }
                     else
                         stop("Cell pairwise distance matrix supplied is of incorrect size.")
                 } )

#' @name cellDist<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellDist<-"
setReplaceMethod("cellDist", signature(object = "SCESet", value = "dist"),
                 function(object, value) {
                     if ( length(value) == choose(ncol(object), 2) ) {
                         object@cellPairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Cell pairwise dist object supplied is of incorrect size.")
                 } )



################################################################################
### featurePairwiseDistances

#' featurePairwiseDistances in an SCESet object
#'
#' SCESet objects can contain a matrix of pairwise distances between features
#' (e.g. genes, transcripts). These functions conveniently access and replace
#' the gene pairwise distances with the value supplied, which must be a matrix
#' of the correct size. The function \code{featDist} is simply shorthand for
#' \code{featurePairwiseDistances}.
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing feature pairwise
#' distances
#'
#' @docType methods
#' @name featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances featurePairwiseDistances,SCESet-method featurePairwiseDistances<-,SCESet,matrix-method featurePairwiseDistances<-,SCESet,dist-method featDist featDist,SCESet-method featDist<-,SCESet,matrix-method featDist<-,SCESet,dist-method
#'
#' @author Davis McCarthy
#'
#' @return An SCESet object containing new feature pairwise distances matrix.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' featurePairwiseDistances(example_sceset)
#'
featurePairwiseDistancesSCESet <- function(object) {
    object@featurePairwiseDistances
}

#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances
#' @export
setMethod("featurePairwiseDistances", signature(object = "SCESet"),
          featurePairwiseDistancesSCESet)

#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export
featDistSCESet <- function(object) {
    object@featurePairwiseDistances
}

#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export
setMethod("featDist", signature(object = "SCESet"), featDistSCESet)

#' @name featurePairwiseDistances<-
#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export "featurePairwiseDistances<-"
setReplaceMethod("featurePairwiseDistances", signature(object = "SCESet",
                                                       value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == nrow(object) ) {
                         object@featurePairwiseDistances <- as.dist(value)
                         return(object)
                     }
                     else
                         stop("Feature pairwise distance matrix supplied is of incorrect size.")
                 } )

#' @name featurePairwiseDistances<-
#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export "featurePairwiseDistances<-"
setReplaceMethod("featurePairwiseDistances", signature(object = "SCESet",
                                                       value = "dist"),
                 function(object, value) {
                     if ( length(value) == choose(nrow(object)) ) {
                         object@featurePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Feature pairwise dist object supplied is of incorrect size.")
                 } )

#' @name featDist<-
#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances
#' @export "featDist<-"
setReplaceMethod("featDist", signature(object = "SCESet", value = "matrix"),
                 function(object, value) {
                     if ( nrow(value) == nrow(object) ) {
                         object@featurePairwiseDistances <- as.dist(value)
                         return(object)
                     }
                     else
                         stop("Feature pairwise distance matrix supplied is of incorrect size.")
                 } )


#' @name featDist<-
#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances
#' @export "featDist<-"
setReplaceMethod("featDist", signature(object = "SCESet", value = "dist"),
                 function(object, value) {
                     if ( nrow(value) == choose(nrow(object), 2) ) {
                         object@featurePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Feature pairwise distance matrix supplied is of incorrect size.")
                 } )


################################################################################
### featureControlInfo

#' featureControlInfo in an SCESet object
#'
#' Each SCESet object stores optional information about the controls in the
#' \code{featureControlInfo} slot. These functions can be used to access,
#' replace or modify this information.
#'
#' @param object a \code{SCESet} object.
#' @param value an AnnotatedDataFrame object, where each row contains
#' information for a single set of control features.
#'
#' @docType methods
#' @name featureControlInfo
#' @rdname featureControlInfo
#' @aliases featureControlInfo featureControlInfo,SCESet-method featureControlInfo<-,SCESet,AnnotatedDataFrame-method featureControlInfo<-
#'
#' @author Aaron Lun
#'
#' @return An SCESet object containing new feature control information.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                             feature_controls = list(ERCC = 1:40, Mito=41:50))
#' featureControlInfo(example_sceset)
#' featureControlInfo(example_sceset)$IsSpike <- c(TRUE, FALSE)
featureControlInfo.SCESet <- function(object) {
    object@featureControlInfo
}

#' @rdname featureControlInfo
#' @aliases featureControlInfo
#' @export
setMethod("featureControlInfo", signature(object = "SCESet"),
          featureControlInfo.SCESet)

#' @name featureControlInfo<-
#' @aliases featureControlInfo
#' @rdname featureControlInfo
#' @export "featureControlInfo<-"
setReplaceMethod("featureControlInfo", signature(object = "SCESet",
                                                 value = "AnnotatedDataFrame"),
                 function(object, value) {
                     object@featureControlInfo <- value
                     return(object)
                 })

################################################################################
### isSpike

#' Get spike-in features in an SCESet object
#'
#' Get the features in the SCESet object that are spike-in controls, as 
#' specified using \code{\link{setSpike<-}}.
#'
#' @param object a \code{SCESet} object.
#' @param type a character vector specifying the feature control sets to use. All 
#' specified spike-in sets in \code{featureControlInfo(object)} are used by default.
#' @param warning A logical scalar specifying if a warning should be raised 
#' if spike-in controls are unavailable.
#'
#' @docType methods
#' @name isSpike
#' @rdname isSpike
#' @aliases isSpike isSpike,SCESet-method 
#'
#' @author Aaron Lun
#'
#' @return A logical vector specifying if each row is a spike-in feature.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                             feature_controls = list(ERCC = 1:40, Mito=41:50))
#' setSpike(example_sceset) <- "ERCC"
#' summary(isSpike(example_sceset))
setMethod("isSpike", "SCESet", 
          function(object, type=NULL, warning=TRUE) {
              if (is.null(type)) {
                  is.spike <- fData(object)$is_feature_spike
              } else {
                  not.in <- !(type %in% .spike_fcontrol_names(object))
                  if (any(not.in)) { 
                      stop(sprintf("'%s' is not specified as a spike-in control", 
                                   type[not.in][1]))
                  } 
                  
                  # Returning directly if possible.
                  if (length(type)==1L) {
                      is.spike <- fData(object)[[paste0("is_feature_control_", type)]]
                  } else {
                    # Combining the spike-in identities. 
                      is.spike <- logical(nrow(object)) 
                      for (f in type) {
                          is.spike <- is.spike | fData(object)[[paste0("is_feature_control_", f)]]
                      }
                  }
              }

              if (warning && is.null(is.spike)) {
                  if (!is.null(type)) {
                      extra <- sprintf(" for '%s'", type)
                  } else {
                      extra <- ""
                  }
                  warning(sprintf("no spike-ins specified%s, returning NULL", extra)) 
              }
              return(is.spike)
          })

################################################################################
### setSpike

#' Set spike-in features in an SCESet object
#'
#' Specify which feature control sets in the SCESet object are spike-ins, i.e., 
#' RNA of the same type and quantity added to each cell during the scRNA-seq protocol.
#'
#' @param object a \code{SCESet} object.
#' @param value a character vector containing the names of the feature control
#' sets that are spike-ins. If \code{NULL}, all spike-in information is removed. 
#'
#' @details 
#' While it is possible to declare overlapping sets as the spike-in sets with \code{isSpike(x)<-}, this is not advisable.
#' This is because some downstream operations assume that each row belongs to only one set (i.e., one of the spike-in sets, or the set of endogenous genes).
#' For example, normalization will use size factors from only one of the sets, so correspondence to multiple sets will not be honoured.
#' Thus, a warning will be raised if overlapping sets are specified in \code{value}.
#'
#' @docType methods
#' @name setSpike
#' @rdname setSpike
#' @aliases setSpike setSpike,SCESet,NULL-method setSpike,SCESet,character-method
#'
#' @author Aaron Lun
#'
#' @return A \code{SCESet} object containing spike-in information in 
#' \code{featureControlInfo} and an updated \code{is_feature_spike} vector for 
#' extraction with \code{\link{isSpike}}.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                             feature_controls = list(ERCC = 1:40, Mito=41:50))
#' setSpike(example_sceset) <- "ERCC"
#' featureControlInfo(example_sceset)
#' summary(isSpike(example_sceset))
setReplaceMethod("setSpike", signature(object="SCESet", value="NULL"), 
                 function(object, value) {
                     fData(object)$is_feature_spike <- NULL 
                     featureControlInfo(object)$spike <- NULL
                     return(object) 
                 })

setReplaceMethod("setSpike", signature(object="SCESet", value="character"), 
                 function(object, value) {
                     # Recording all those that were listed as spikes.
                     featureControlInfo(object)$spike <- .fcontrol_names(object) %in% value
                     
                     # Setting the default is_feature_spike.
                     fData(object)$is_feature_spike <- isSpike(object, value)
                     
                     # Checking that they don't overlap.
                     if (length(value) > 1L) { 
                         total.hits <- integer(nrow(object))
                         for (v in value) {
                             total.hits <- total.hits + isSpike(object, v)
                         }
                         if (any(total.hits > 1L)) { 
                             warning("overlapping spike-in sets detected")
                         }
                     }
                     
                     return(object) 
                 })

################################################################################
### whichSpike

#' Identify spike-in feature control sets in an SCESet object
#'
#' Get the names of the feature control sets that are spike-ins.
#'
#' @param object a \code{SCESet} object.
#'
#' @docType methods
#' @name whichSpike
#' @rdname whichSpike
#' @aliases whichSpike whichSpike,SCESet-method
#'
#' @author Aaron Lun
#'
#' @return A character vector containing the names of feature control sets
#' that are spike-in sets.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                             feature_controls = list(ERCC = 1:40, Mito=41:50))
#' setSpike(example_sceset) <- "ERCC"
#' whichSpike(example_sceset)
setMethod("whichSpike", signature("SCESet"), 
          function(object) .spike_fcontrol_names(object))

################################################################################
### spikes

#' Extract expression values for spike-in features in an SCESet object
#'
#' Extract a matrix of expression values for features in spike-in control sets.
#'
#' @param object a \code{SCESet} object.
#' @param exprs_values a string specifying the type of expression values to extract.
#' @param type a character vector containing the names of the spike-in control sets to extract. By default,
#' expression values for features in all spike-in control sets are extracted.
#' 
#' @details
#' If \code{exprs_values="exprs"}, users should have run \code{normalize(object)} first,
#' so that spike-in features are normalized with spike-in size factors.
#'
#' @docType methods
#' @name spikes
#' @rdname spikes
#' @aliases spikes spikes,SCESet-method
#'
#' @author Aaron Lun
#'
#' @return A matrix of expression values for features in the specified spike-in control sets.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset,
#'                             feature_controls = list(ERCC = 1:40, Mito=41:50))
#' setSpike(example_sceset) <- "ERCC"
#' head(spikes(example_scesets))
setMethod("spikes", "SCESet", 
          function(object, exprs_values="counts", type=NULL) {
              is.spike <- isSpike(object, type=type)
              get_exprs(object, exprs_values)[is.spike,,drop=FALSE]
          })


################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}
#'
#' @param sce An \code{SCESet} object
#' @param exprs_values What should \code{exprs(cds)} be mapped from in the \code{SCESet}? Should be
#' one of "exprs", "tpm", "fpkm", "counts"
#'
#' @export
#' @rdname toCellDataSet
#' @name toCellDataSet
#' @return An object of class \code{SCESet}
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' if ( requireNamespace("monocle") ) {
#'     toCellDataSet(example_sceset)
#' }
toCellDataSet <- function(sce, exprs_values = "exprs") {
    pkgAvail <- requireNamespace("monocle")
    if (pkgAvail) {
        if (!is(sce,'SCESet')) stop('sce must be of type SCESet')
        exprsData <- NULL
        exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", "counts"))
        exprsData <- switch(exprs_values,
                            exprs = exprs(sce),
                            tpm = tpm(sce),
                            fpkm = fpkm(sce),
                            counts = counts(sce))
        if ( exprs_values == "exprs" && sce@logged )
            exprsData <- 2 ^ exprsData - sce@logExprsOffset
        celldataset <- monocle::newCellDataSet(
            exprsData, phenoData = phenoData(sce),
            featureData = featureData(sce),
            lowerDetectionLimit = sce@lowerDetectionLimit)
        celldataset@reducedDimS <- t(redDim(sce))
        celldataset
    }
    else
        stop("Require package monocle to be installed to use this function.")
}

#' Convert a \code{CellDataSet} to an \code{SCESet}
#'
#' @param cds A \code{CellDataSet} from the \code{monocle} package
#' @param exprs_values What should \code{exprs(cds)} be mapped to in the \code{SCESet}? Should be
#' one of "exprs", "tpm", "fpkm", "counts"
#' @param logged logical, if a value is supplied for the exprsData argument, are
#'  the expression values already on the log2 scale, or not?
#'
#' @export
#' @importFrom Biobase featureData
#' @importFrom Biobase phenoData
#' @rdname fromCellDataSet
#' @name fromCellDataSet
#' @return An object of class \code{SCESet}
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' if ( requireNamespace("monocle") ) {
#'     # cds <- toCellDataSet(example_sceset) # not run
#'     # sceset <- fromCellDataSet(cds) # not run
#' }
fromCellDataSet <- function(cds, exprs_values = "tpm", logged = FALSE) {
    pkgAvail <- requireNamespace("monocle")
    if (pkgAvail) {
        if (!is(cds,'CellDataSet')) stop('cds must be of type CellDataSet from package monocle')
        exprsData <- countData <- NULL
        exprs_values <- match.arg(exprs_values,
                                  c("exprs", "tpm", "fpkm", "counts"))
        exprsData <- countData <- tpmData <- fpkmData <- NULL
        if (exprs_values == "exprs") {
            exprsData <- exprs(cds)
        } else if (exprs_values == "tpm") {
            tpmData <- exprs(cds)
        } else if (exprs_values == "fpkm") {
            fpkmData <- exprs(cds)
        } else {
            countData <- exprs(cds)
        }

        sce <- newSCESet(exprsData = exprsData, tpmData = tpmData,
                         fpkmData = fpkmData, countData = countData,
                         phenoData = phenoData(cds),
                         featureData = featureData(cds),
                         lowerDetectionLimit = cds@lowerDetectionLimit)

        ## now try and preserve a reduced dimension representation
        ## this is really not elegant - KC
        rds <- cds@reducedDimS
        if (length(rds) == 0) {
            rds <- cds@reducedDimA
        }
        if (length(rds) == 2 * ncol(sce)) { # something is there and of the right dimension
            if (nrow(rds) == 2) {
                redDim(sce) <- t(rds)
            } else if (ncol(rds) == 2) {
                redDim(sce) <- rds
            } # else do nothing
        }

        return( sce )
    }
    else
        stop("Require package monocle to be installed to use this function.")
}

################################################################################

#' Retrieve a representation of gene expression
#'
#' Deprecated from scater version 1.3.29.
#'
#' @param object An object of type \code{SCESet}
#' @return A matrix representation of expression corresponding to \code{object@useForExprs}.
#'
getExprs <- function(object) {
    stop("Deprecated from scater 1.3.29")
}


################################################################################

#' Merge SCESet objects
#'
#' Merge two SCESet objects that have the same features but contain different cells/samples.
#'
#' @param x an \code{\link{SCESet}} object
#' @param y an \code{\link{SCESet}} object
#' @param fdata_cols_x a logical or numeric vector indicating which columns of featureData
#' for \code{x} are shared between \code{x} and \code{y} and should feature in the
#' returned merged \code{SCESet}. Default is all columns of \code{fData(x)}.
#' @param fdata_cols_y a logical or numeric vector indicating which columns of featureData
#' for \code{y} are shared between \code{x} and \code{y} and should feature in the
#' returned merged \code{SCESet}. Default is \code{fdata_cols_x}.
#' @param pdata_cols_x a logical or numeric vector indicating which columns of phenoData
#' of \code{x} should be retained.
#' @param pdata_cols_y a logical or numeric vector indicating which columns of phenoData
#' of \code{y} should be retained.
#'
#' @details Existing cell-cell pairwise distances and feature-feature pairwise distances will not be valid for a merged SCESet so these are set to \code{NULL} in the returned object. Similarly \code{experimentData} will need to be added anew to the merged SCESet returned.
#'
#' @return a merged \code{SCESet} object combining data and metadata from \code{x} and \code{y}
#'
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40])
#'
#' ## with specification of columns of fData
#' example_sceset <- calculateQCMetrics(example_sceset)
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40], fdata_cols_x = c(1, 7))
#'
#' ## with specification of columns of pData
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40], pdata_cols_x = 1:6)
#'
#'
mergeSCESet <- function(x, y, fdata_cols_x = 1:ncol(fData(x)), fdata_cols_y = fdata_cols_x,
                        pdata_cols_x = NULL, pdata_cols_y = NULL) {
    if (!is(x,'SCESet')) stop('x must be of type SCESet')
    if (!is(y,'SCESet')) stop('y must be of type SCESet')
    if (!identical(featureNames(x), featureNames(y))) stop("feature names of x and y must be identical")

    if (x@logged != y@logged)
        stop("x and y do not have the same value for the 'logged' slot.")
    if (x@lowerDetectionLimit != y@lowerDetectionLimit)
        stop("x and y do not have the same lowerDetectionLimit.")
    if (x@logExprsOffset != y@logExprsOffset)
        stop("x and y do not have the same logExprsOffset.")

    ## combine fData
    if (ncol(fData(x)) == 0) {
        if (!identical(fData(x), fData(y)))
            stop("featureData do not match for x and y.")
        new_fdata <- as(fData(x), "AnnotatedDataFrame")
    } else {
        fdata1 <- fData(x)[, fdata_cols_x, drop = FALSE]
        fdata2 <- fData(y)[, fdata_cols_y, drop = FALSE]
        if (!identical(fdata1, fdata2))
            stop("featureData columns specified are not identical for x and y.")
        new_fdata <- as(fdata1, "AnnotatedDataFrame")
    }
    ## combine pData
    if (ncol(pData(x)) == 0 || pData(y) == 0)
        stop("phenoData slot is empty for x or y.")
    pdata_x <- pData(x)
    pdata_y <- pData(y)
    if (is.null(pdata_cols_x)) {
        if (is.null(pdata_cols_y)) {
            pdata_cols_x <- which(colnames(pdata_x) %in% colnames(pdata_y))
            pdata_cols_y <- which(colnames(pdata_y) %in% colnames(pdata_x))
        } else
            pdata_cols_x <- which(colnames(pdata_x) %in% colnames(pdata_y)[pdata_cols_y])
    } else {
        if (is.null(pdata_cols_y))
            pdata_cols_y <- which(colnames(pdata_y) %in% colnames(pdata_x)[pdata_cols_x])
    }
    if (length(pdata_cols_x) == 0 | length(pdata_cols_y) == 0)
        stop("no phenoData column names found in common between x and y.")
    ## make sure ordering of columns is correct
    pdata_x <- pdata_x[, pdata_cols_x, drop = FALSE]
    pdata_y <- pdata_y[, pdata_cols_y, drop = FALSE]
    mm <- match(colnames(pdata_x), colnames(pdata_y))
    pdata_y <- pdata_y[, mm]
    if (!identical(colnames(pdata_x), colnames(pdata_y)))
        stop("phenoData columns specified are not identical for x and y.")
    new_pdata <- rbind(pdata_x, pdata_y)
    new_pdata <- as(new_pdata, "AnnotatedDataFrame")
    ## combine exprsData
    new_exprs <- Biobase::combine(exprs(x), exprs(y))
    ## new SCESet
    merged_sceset <- newSCESet(exprsData = new_exprs, featureData = new_fdata,
                               phenoData = new_pdata,
                               logExprsOffset = x@logExprsOffset)
    ## add remaining assayData to merged SCESet
    assay_names <- intersect(names(Biobase::assayData(x)),
                             names(Biobase::assayData(y)))
    for (assaydat in assay_names) {
        new_dat <- Biobase::combine(get_exprs(x, assaydat), get_exprs(y, assaydat))
        set_exprs(merged_sceset, assaydat) <- new_dat
    }
    merged_sceset
}



################################################################################
## writeSCESet

#' Write an SCESet object to an HDF5 file
#'
#' @param object \code{\link{SCESet}} object to be writted to file
#' @param file_path path to written file containing data from SCESet object
#' @param type character string indicating type of output file. Default is "HDF5".
#' @param overwrite_existing logical, if a file of the same name already exists
#' should it be overwritten? Default is \code{FALSE}.
#'
#' @details Currently writing to HDF5 files is supported. The \pkg{\link{rhdf5}}
#' package is used to write data to file and can be used to read data from HDF5
#' files into R. For further details about the HDF5 data format see
#' \url{https://support.hdfgroup.org/HDF5/}.
#'
#' @return Return is \code{NULL}, having written the \code{SCESet} object to file.
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#'
#' \dontrun{
#' writeSCESet(example_sceset, "test.h5")
#' file.remove("test.h5")
#' }
#'
writeSCESet <- function(object, file_path, type = "HDF5", overwrite_existing = FALSE) {
    if (!is(object,'SCESet')) stop('object must be of type SCESet')
    if (file.exists(file_path) && !overwrite_existing)
        stop("To overwrite an existing file use argument overwrite_existing=TRUE")
    if (file.exists(file_path))
        file.remove(file_path)
    if (type == "HDF5") {
        rhdf5::H5close()
        rhdf5::h5createFile(file_path)
        tryCatch({
            rhdf5::h5write(featureNames(object), file = file_path,
                           name = "featureNames")
            rhdf5::h5write(cellNames(object), file = file_path, name = "cellNames")
            rhdf5::h5write(object@logged, file = file_path, name = "logged")
            rhdf5::h5write(object@logExprsOffset, file = file_path,
                           name = "logExprsOffset")
            rhdf5::h5write(object@lowerDetectionLimit, file = file_path,
                           name = "lowerDetectionLimit")
            if (ncol(pData(object)) > 0)
                rhdf5::h5write(pData(object), file = file_path, name = "phenoData",
                               write.attributes = FALSE)
            if (ncol(fData(object)) > 0)
                rhdf5::h5write(fData(object), file = file_path, name = "featureData",
                               write.attributes = FALSE)
            rhdf5::h5createGroup(file_path, "assayData")
            for (assay in names(Biobase::assayData(object))) {
                group_set <- paste0("assayData/", assay)
                rhdf5::h5write(get_exprs(object, assay), file = file_path,
                               name = group_set,
                               write.attributes = FALSE)
            }
            rhdf5::H5close()
        }, finally = rhdf5::H5close())
    } else
        stop("HDF5 is the only format currently supported. See also saveRDS() to write to an object readable with R.")
}
