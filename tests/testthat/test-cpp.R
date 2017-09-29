## Testing the C++ code for correctness.
    
library(Matrix)

test_that("margin summaries are correctly computed", {
    set.seed(12938)

    ## For integer base matrices.
    A.int <- matrix(rpois(150, lambda=2), nrow=15)

    # For row sums
    out <- .Call(scater:::cxx_margin_summary, A.int, 0, 1:10 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.int))
    expect_equal(out[[2]], rowSums(A.int > 0))

    out <- .Call(scater:::cxx_margin_summary, A.int, 1, 1:10 - 1L, TRUE) 
    expect_equal(out[[2]], rowSums(A.int > 1))

    out <- .Call(scater:::cxx_margin_summary, A.int, 0, 2:8 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.int[,2:8]))
    expect_equal(out[[2]], rowSums(A.int[,2:8] > 0))

    # For column sums.
    out <- .Call(scater:::cxx_margin_summary, A.int, 0, 1:15 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.int))
    expect_equal(out[[2]], colSums(A.int > 0))

    out <- .Call(scater:::cxx_margin_summary, A.int, 1, 1:15 - 1L, FALSE) 
    expect_equal(out[[2]], colSums(A.int > 1))

    out <- .Call(scater:::cxx_margin_summary, A.int, 0, 2:8 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.int[2:8,]))
    expect_equal(out[[2]], colSums(A.int[2:8,] > 0))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(150, 10, 10), nrow=15)

    # For row sums
    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 1:10 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.dbl))
    expect_equal(out[[2]], rowSums(A.dbl > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 1:10 - 1L, TRUE) 
    expect_equal(out[[2]], rowSums(A.dbl > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 2:8 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.dbl[,2:8]))
    expect_equal(out[[2]], rowSums(A.dbl[,2:8] > 1))

    # For column sums.
    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 1:15 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.dbl))
    expect_equal(out[[2]], colSums(A.dbl > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 1:15 - 1L, FALSE) 
    expect_equal(out[[2]], colSums(A.dbl > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dbl, 1, 2:8 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.dbl[2:8,]))
    expect_equal(out[[2]], colSums(A.dbl[2:8,] > 1))

    ## For sparse matrices
    A.dgc <- rsparsematrix(15, 10, density=0.1)

    # For row sums
    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 1:10 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.dgc))
    expect_equal(out[[2]], rowSums(A.dgc > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 1:10 - 1L, TRUE) 
    expect_equal(out[[2]], rowSums(A.dgc > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 2:8 - 1L, TRUE) 
    expect_equal(out[[1]], rowSums(A.dgc[,2:8]))
    expect_equal(out[[2]], rowSums(A.dgc[,2:8] > 1))

    # For column sums.
    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 1:15 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.dgc))
    expect_equal(out[[2]], colSums(A.dgc > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 1:15 - 1L, FALSE) 
    expect_equal(out[[2]], colSums(A.dgc > 1))

    out <- .Call(scater:::cxx_margin_summary, A.dgc, 1, 2:8 - 1L, FALSE) 
    expect_equal(out[[1]], colSums(A.dgc[2:8,]))
    expect_equal(out[[2]], colSums(A.dgc[2:8,] > 1))
})

test_that("row/column variances are correctly computed", {
    ## For integer base matrices.
    A.int <- matrix(rpois(150, lambda=2), nrow=15)

    # For row variances.
    expect_equal(scater:::.general_rowVars(A.int), apply(A.int, 1, var))
    expect_equal(scater:::.general_rowVars(A.int, cols=2:8), apply(A.int[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.general_colVars(A.int), apply(A.int, 2, var))
    expect_equal(scater:::.general_colVars(A.int, rows=2:8), apply(A.int[2:8,], 2, var))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(150, 10, 10), nrow=15)

    # For row variances.
    expect_equal(scater:::.general_rowVars(A.dbl), apply(A.dbl, 1, var))
    expect_equal(scater:::.general_rowVars(A.dbl, cols=2:8), apply(A.dbl[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.general_colVars(A.dbl), apply(A.dbl, 2, var))
    expect_equal(scater:::.general_colVars(A.dbl, rows=2:8), apply(A.dbl[2:8,], 2, var))

    ## For sparse matrices.
    A.dgc <- rsparsematrix(15, 10, density=0.1)
    A.dgc2 <- as.matrix(A.dgc)

    # For row variances.
    expect_equal(scater:::.general_rowVars(A.dgc), apply(A.dgc2, 1, var))
    expect_equal(scater:::.general_rowVars(A.dgc, cols=2:8), apply(A.dgc2[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.general_colVars(A.dgc), apply(A.dgc2, 2, var))
    expect_equal(scater:::.general_colVars(A.dgc, rows=2:8), apply(A.dgc2[2:8,], 2, var))
})
