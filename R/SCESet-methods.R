### Methods for the SCESet class

################################################################################
### constructor function for SCESet class

#' Create a new SCESet object.
#'
#' Scater requires that all data be housed in SCESet objects. SCESet extends
#' Bioconductor's ExpressionSet class, and the same basic interface is
#' supported. newSCESet() expects a matrix of expression values as its first
#' argument, with rows as features (usually genes) and columns as cells.
#' Per-feature and per-cell metadata can be supplied with the featureData and
#' phenoData arguments, respectively. Use of these optional arguments is
#' strongly encouraged. The SCESet also includes a slot 'counts' to store an
#' object containing raw count data.
#'
#' @param exprsData expression data matrix for an experiment
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
#' @param logged logical, if a value is supplied for the exprsData argument, are
#'  the expression values already on the log2 scale, or not?

#' @param useForExprs character string, either 'exprs' (default),'tpm','counts' or
#'    'fpkm' indicating which expression representation both internal methods and
#'    external packages should use when performing analyses.
#' @return a new SCESet object
#'
#' @details
#' SCESet objects store a matrix of expression values. These values are
#' typically transcripts-per-million (tpm), counts-per-million (cpm), fragments
#' per kilobase per million mapped (FPKM) or some other output from a program
#' that calculates expression values from RNA-Seq reads. We recommend that
#' expression values on the log2 scale are used for the 'exprs' slot in the
#' SCESet. For example, you may wish to store raw tpm values in the 'tpm' slot
#' and \code{log2(tpm + 1)} values in the 'exprs' slot. However, expression
#' values could also be values from a single cell qPCR run or some other type of
#'  assay. The newSCESet function can also accept raw count values. In this case
#'  see \code{\link{calculateTPM}} and \code{\link{calculateFPKM}} for computing
#'  TPM and FPKM expression values, respectively, from counts. The function
#'  \code{\link[edgeR]{cpm}} from the package edgeR to can be used to compute
#'  log2(counts-per-million), if desired.
#'
#'  An \code{SCESet} object has to have the \code{'exprs'} slot defined, so if
#'  the \code{exprsData} argument is \code{NULL}, then this function will define
#'  \code{'exprs'} with the following order of precedence: log2(TPM +
#'  logExprsOffset), if \code{tpmData} is defined; log2(FPKM + logExprsOffset)
#'  if \code{fpkmData} is defined; otherwise log2(counts-per-million +
#'  logExprsOffset) are used. The \code{\link[edgeR]{cpm}} function from the
#'  edgeR package is used to compte \code{cpm}. Note that for many analyses
#'  counts-per-million are not recommended, and if possible
#'  transcripts-per-million should be used.
#'
#'  In many downstream functions you will likely find it most convenient if the
#'  \code{'exprs'} values are on the log2-scale, so this is recommended.
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
                      lowerDetectionLimit = 0,
                      logExprsOffset = 1,
                      logged = FALSE,
                      useForExprs = "exprs")
{
    ## Check that we have some expression data
    if ( is.null(exprsData) && is.null(countData) && is.null(tpmData) &
         is.null(fpkmData) && is.null(cpmData))
        stop("Require at least one of exprsData, tpmData, fpkmData or countData arguments.")
    ## Check dimensions of data matrices

    ## Check counts are a matrix; renames is_exprsData if not null
    if ( !is.null(countData) )
        countData <- as.matrix(countData)
    if ( !is.null(is_exprsData) )
        isexprs <- is_exprsData

    ## If no exprsData provided define is_exprs from tpmData, fpkmData or countData
    if ( is.null(exprsData) ) {
        ## Define exprs data if null
        if ( !is.null(tpmData) ) {
            exprsData <- log2(tpmData + logExprsOffset)
            logged <- TRUE
        } else {
            if ( !is.null(fpkmData) ) {
                exprsData <- log2(fpkmData + logExprsOffset)
                logged <- TRUE
            } else {
                if ( !is.null(cpmData) ) {
                    exprsData <- log2(cpmData + logExprsOffset)
                    logged <- TRUE
                }  else {
                    exprsData <- .cpm_default(countData, prior.count = logExprsOffset, log = TRUE)
                    logged <- TRUE
                }
            }
        }
    } else {
        exprsData <- as.matrix(exprsData)
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
            warning("experimentData supplied is not an 'MIAME' object. Thus, experimentData is being set to an empty MIAME object.\n Please supply a valid 'MIAME' class object containing experiment data to experimentData(object).")
        }
    } else {
        expData <- expData_null
    }

    ## Check valid useForExprs
    useForExprs <- match.arg(useForExprs, c("exprs","tpm","counts","fpkm"))

    ## Generate new SCESet object
    if ( !is.null(is_exprsData) )
        assaydata <- assayDataNew("lockedEnvironment", exprs = exprsData,
                                  is_exprs = isexprs)
    else 
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
                   logged = logged,
                   featureControlInfo = AnnotatedDataFrame(),
                   useForExprs = useForExprs)

    ## Add non-null slots to assayData for SCESet object, omitting null slots
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

#' Get cell names from an SCESet object
#'
#' @param object An SCESet object.
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
#' \S4method{get_exprs}{SCESet}(object, exprs_values)
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
get_exprs.SCESet <- function(object, exprs_values = "exprs") {
    exprs_mat <- object@assayData[[exprs_values]]
    if ( is.null(exprs_mat) )
        warning(paste0("The object does not contain ", exprs_values, " expression values. Returning NULL."))
    exprs_mat
}

#' @rdname get_exprs
#' @export
setMethod("get_exprs", signature(object = "SCESet"),
          function(object, exprs_values = "exprs") {
              get_exprs.SCESet(object, exprs_values)
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
#' @aliases set_exprs set_exprs<-,SCESet,ANY,matrix-method
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
#' @name set_exprs
#' @rdname set_exprs
#' @exportMethod "set_exprs<-"
setReplaceMethod("set_exprs", signature(object = "SCESet", value = "matrix"),
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
        fc_available <- .get_feature_control_names(object)
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
                         lowerDetectionLimit = cds@lowerDetectionLimit,
                         logged = logged)

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
#' Gene expression can be summarised in a variety of ways, e.g. as TPM, FPKM or
#' as raw counts. Many internal methods and external packages rely on accessing
#' a generic representation of expression without worrying about the particulars.
#' Scater allows the user to set \code{object@@useForExprs} to the preferred
#' type (either "exprs", "TPM", "fpkm" or "counts") and that particular representation
#' will be returned by calls to \code{getExprs}. Note if such representation is
#' not defined, this method returns \code{NULL}.
#'
#' @param object An object of type \code{SCESet}
#' @return A matrix representation of expression corresponding to \code{object@@useForExprs}.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd, useForExprs = "exprs")
#' all(exprs(example_sceset) == getExprs(example_sceset)) # TRUE
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd, useForExprs = "counts")
#' all(exprs(example_sceset) == getExprs(example_sceset)) # FALSE
#' all(counts(example_sceset) == getExprs(example_sceset)) # TRUE
getExprs <- function(object) {
    if (!is(object,'SCESet')) stop('object must be of type SCESet')

    x <- switch(object@useForExprs,
                exprs = exprs(object),
                tpm = tpm(object),
                fpkm = fpkm(object),
                counts = counts(object))

    if(is.null(x)) warning(paste("Slot for", object@useForExprs, "is empty; returning NULL"))

    return( x )
}
