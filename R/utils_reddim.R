#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
.get_mat_from_sce <- function(x, exprs_values, dimred, n_dimred, assay.type=exprs_values) {
    if (!is.null(dimred)) {
        mat <- reducedDim(x, dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred)==1L) {
                n_dimred <- seq_len(n_dimred)
            }
            mat <- mat[,n_dimred,drop=FALSE]
        }
        mat
    } else {
        assay(x, assay.type)
    }
}

#' @importFrom utils head
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom beachmat realizeFileBackedMatrix
.get_mat_for_reddim <- function(x, subset_row, ntop, scale, get.var=FALSE)
# Picking the 'ntop' most highly variable features or just using a pre-specified set of features.
# Also removing zero-variance columns and scaling the variance of each column.
# Finally, transposing for downstream use (cells are now rows).
{
    use.var <- is.null(subset_row) || scale || get.var
    if (use.var) {
        rv <- rowVars(DelayedArray(x))
    }

    if (is.null(subset_row)) {
        o <- order(rv, decreasing = TRUE)
        subset_row <- head(o, ntop)
    } else if (is.character(subset_row)) {
        subset_row <- .subset2index(subset_row, x, byrow=TRUE)
    }

    x <- x[subset_row,, drop = FALSE]
    if (is.null(rownames(x))) {
        rownames(x) <- subset_row
    }
    if (use.var) {
        rv <- rv[subset_row]
    }

    if (scale) {
        keep <- rv >= 1e-8
        x <- x[keep,,drop=FALSE]/sqrt(rv[keep])
        rv <- rep(1, nrow(x))
    }

    x <- t(x)
    x <- realizeFileBackedMatrix(x)

    if (get.var) {
        list(x=x, v=rv)
    } else {
        x
    }
}

#' @importFrom BiocParallel bpnworkers
#' @importClassesFrom BiocParallel MulticoreParam 
.choose_nthreads <- function(val, BPPARAM) {
    if (is.null(val) && is(BPPARAM, "MulticoreParam")) {
        bpnworkers(BPPARAM)
    } else {
        val
    }
}
