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

#' @importFrom BiocParallel bpiterate bpnworkers SerialParam
.iterate_by_chunks <- function(x, FUN, ..., row_sets=list(NULL), col_sets=list(NULL), 
    BPPARAM=SerialParam(), combineROW=c, combineCOL=cbind) 
{
    env <- new.env()
    env$i <- env$j <- 1L

    out <- bpiterate(ITER=function() {
        i <- env$i
        j <- env$j

        if (i==length(row_sets)) {
            env$i <- 1L
            env$j <- env$j + 1L
        } else {
            env$i <- env$i + 1L
        }

        if (j>length(col_sets)) {
            return(NULL)
        } 

        row <- row_sets[[i]]
        col <- col_sets[[j]]
        if (!is.null(row) && !is.null(col)) {
            x[row,col,drop=FALSE]
        } else if (!is.null(row)) {
            x[row,,drop=FALSE]
        } else if (!is.null(col)) {
            x[,col,drop=FALSE]
        } else {
            x
        }
    }, FUN=FUN, ..., BPPARAM=BPPARAM)

    # Combining by row, and then combining by column.
    collected <- list()
    for (i in seq_along(col_sets)) {
        current <- length(row_sets) * (i - 1L) + seq_along(row_sets)
        collected[[i]] <- do.call(combineROW, out[current])
    }

    do.call(combineCOL, collected)
}
