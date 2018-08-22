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
