#' Parallelization for developers 
#'
#' Parallelization utilities for re-use in packages that happen to depend on \pkg{scater}.
#' These are exported simply to avoid re-writing them in downstream packages, and should not be touched by end-users.
#'
#' @author Aaron Lun
#' @name utils_parallel
#' @docType class
#' @aliases .splitRowsByWorkers
#' .splitColsByWorkers
#' .assignIndicesToWorkers 
#' .subset2index
#' .subsetToIndexOrNull
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

#' @importFrom BiocParallel bpnworkers
#' @importFrom utils head tail
.assign_jobs_to_workers <- function(njobs, BPPARAM) 
# Returns ZERO-INDEXED job start/ends.
{
    n_cores <- bpnworkers(BPPARAM)
    boundaries <- as.integer(seq(from = 0L, to = njobs, length.out = n_cores + 1L)) 
    work_ends <- tail(boundaries, -1L)
    work_starts <- head(boundaries, -1L)
    return(list(start=work_starts, end=work_ends))
}

.split_vector_by_workers <- function(values, BPPARAM) {
    job_indices <- .assign_jobs_to_workers(length(values), BPPARAM)
    out <- vector("list", length(job_indices$start))
    for (i in seq_along(out)) {
        cur_start <- job_indices$start[i]
        cur_end <- job_indices$end[i]
        out[[i]] <- values[cur_start + seq_len(cur_end - cur_start)]
    }
    return(out)
}

.split_subset_by_workers <- function(subset, ..., BPPARAM) {
    if (is.null(subset) && bpnworkers(BPPARAM)==1L) {
        list(NULL)
    } else {
        subset <- .subset2index(subset, ...)
        .split_vector_by_workers(subset, BPPARAM)
    }
}

