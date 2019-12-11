#' Developer utilities
#'
#' Various utilities for re-use in packages that happen to depend on \pkg{scater}.
#' These are exported simply to avoid re-writing them in downstream packages, and should not be touched by end-users.
#'
#' @author Aaron Lun
#' @name scater-utils
#' @docType class
#' @aliases .splitRowsByWorkers
#' .splitColsByWorkers
#' .assignIndicesToWorkers 
#' .subset2index
#' .subsetToIndexOrNull
#' .bpNotSharedOrUp
NULL

#' @export
#' @importFrom BiocParallel bpnworkers
.splitRowsByWorkers <- function(x, BPPARAM, subset_row=NULL, subset_col=NULL, assignments=NULL) {
    if (bpnworkers(BPPARAM)==1L) {
        if (!is.null(subset_row)) {
            x <- x[subset_row,,drop=FALSE]
        }
        if (!is.null(subset_col)) {
            x <- x[,subset_col,drop=FALSE]
        }

        list(x)
    } else {
        if (is.null(assignments)) {
            assignments <- .assignIndicesToWorkers(nrow(x), BPPARAM, subset=subset_row)
        }

        for (i in seq_along(assignments)) {
            current <- x[assignments[[i]],,drop=FALSE]
            if (!is.null(subset_col)) {
                current <- current[,subset_col,drop=FALSE]
            }
            assignments[[i]] <- current
        }
    
        assignments
    }
}

#' @export
#' @importFrom BiocParallel bpnworkers
.splitColsByWorkers <- function(x, BPPARAM, subset_row=NULL, subset_col=NULL, assignments=NULL) {
    if (bpnworkers(BPPARAM)==1L) {
        if (!is.null(subset_row)) {
            x <- x[subset_row,,drop=FALSE]
        }
        if (!is.null(subset_col)) {
            x <- x[,subset_col,drop=FALSE]
        }

        list(x)
    } else {
        if (is.null(assignments)) {
            assignments <- .assignIndicesToWorkers(ncol(x), BPPARAM, subset=subset_col)
        }

        for (i in seq_along(assignments)) {
            current <- x[,assignments[[i]],drop=FALSE]
            if (!is.null(subset_row)) {
                current <- current[subset_row,,drop=FALSE]
            }
            assignments[[i]] <- current
        }
    
        assignments
    }
}

#' @export
#' @importFrom BiocParallel bpnworkers
#' @importFrom utils head
.assignIndicesToWorkers <- function(njobs, BPPARAM, subset=NULL) {
    if (!is.null(subset)) {
        subset <- as.vector(subset)
        if (is.logical(subset)) {
            subset <- which(subset)
        }
        njobs <- length(subset)
    }

    n_cores <- bpnworkers(BPPARAM)
    boundaries <- as.integer(seq(from = 0L, to = njobs, length.out = n_cores + 1L))
    per_core <- diff(boundaries)
    work_starts <- head(boundaries, -1L)
    output <- mapply("+", lapply(per_core, seq_len), work_starts, SIMPLIFY=FALSE)

    if (!is.null(subset)) {
        for (i in seq_along(output)) {
            output[[i]] <- subset[output[[i]]]            
        }
    }

    output
}

#' @export
#' @importClassesFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bpisup
.bpNotSharedOrUp  <- function(BPPARAM) !bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")
