## Testing the C++ code for correctness.
library(Matrix)
    
test_that("margin summaries are correctly computed", {
    set.seed(12938)

    ## For integer base matrices.
    A.int <- matrix(rpois(150, lambda=2), nrow=15)

    # For row sums
    expect_equal(scater:::.rowSums(A.int), rowSums(A.int))
    expect_equal(scater:::.rowAbove(A.int), rowSums(A.int > 0))
    expect_equal(scater:::.rowAbove(A.int, value=1), rowSums(A.int > 1))

    expect_equal(scater:::.rowSums(A.int, rows=1:5), rowSums(A.int[1:5,]))
    expect_equal(scater:::.rowAbove(A.int, rows=1:5), rowSums(A.int[1:5,] > 0))
    expect_equal(scater:::.rowAbove(A.int, rows=1:5, value=1), rowSums(A.int[1:5,] > 1))

    expect_equal(scater:::.rowSums(A.int, cols=2:8), rowSums(A.int[,2:8]))
    expect_equal(scater:::.rowAbove(A.int, cols=2:8), rowSums(A.int[,2:8] > 0))
    expect_equal(scater:::.rowAbove(A.int, cols=2:8, value=1), rowSums(A.int[,2:8] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.int), colSums(A.int))
    expect_equal(scater:::.colAbove(A.int), colSums(A.int > 0))
    expect_equal(scater:::.colAbove(A.int, value=1), colSums(A.int > 1))

    expect_equal(scater:::.colSums(A.int, rows=2:8), colSums(A.int[2:8,]))
    expect_equal(scater:::.colAbove(A.int, rows=2:8), colSums(A.int[2:8,] > 0))
    expect_equal(scater:::.colAbove(A.int, rows=2:8, value=1), colSums(A.int[2:8,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=1:5), colSums(A.int[,1:5]))
    expect_equal(scater:::.colAbove(A.int, cols=1:5), colSums(A.int[,1:5] > 0))
    expect_equal(scater:::.colAbove(A.int, cols=1:5, value=1), colSums(A.int[,1:5] > 1))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(150, 10, 10), nrow=15)

    # For row sums
    expect_equal(scater:::.rowSums(A.dbl), rowSums(A.dbl))
    expect_equal(scater:::.rowAbove(A.dbl, value=1), rowSums(A.dbl > 1))

    expect_equal(scater:::.rowSums(A.int, rows=1:5), rowSums(A.int[1:5,]))
    expect_equal(scater:::.rowAbove(A.int, rows=1:5, value=1), rowSums(A.int[1:5,] > 1))

    expect_equal(scater:::.rowSums(A.dbl, cols=2:8), rowSums(A.dbl[,2:8]))
    expect_equal(scater:::.rowAbove(A.dbl, cols=2:8, value=1), rowSums(A.dbl[,2:8] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.dbl), colSums(A.dbl))
    expect_equal(scater:::.colAbove(A.dbl, value=1), colSums(A.dbl > 1))

    expect_equal(scater:::.colSums(A.dbl, rows=2:8), colSums(A.dbl[2:8,]))
    expect_equal(scater:::.colAbove(A.dbl, rows=2:8, value=1), colSums(A.dbl[2:8,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=1:5), colSums(A.int[,1:5]))
    expect_equal(scater:::.colAbove(A.int, cols=1:5, value=1), colSums(A.int[,1:5] > 1))

    ## For sparse matrices
    A.dgc <- rsparsematrix(15, 10, density = 0.1)

    # For row sums
    expect_equal(scater:::.rowSums(A.dgc), Matrix::rowSums(A.dgc))
    expect_equal(scater:::.rowAbove(A.dgc), Matrix::rowSums(A.dgc > 0))
    expect_equal(scater:::.rowAbove(A.dgc, value=1), Matrix::rowSums(A.dgc > 1))

    expect_equal(scater:::.rowSums(A.int, rows=1:5), rowSums(A.int[1:5,]))
    expect_equal(scater:::.rowAbove(A.int, rows=1:5), rowSums(A.int[1:5,] > 0))
    expect_equal(scater:::.rowAbove(A.int, rows=1:5, value=1), rowSums(A.int[1:5,] > 1))

    expect_equal(scater:::.rowSums(A.dgc, cols=2:8), Matrix::rowSums(A.dgc[,2:8]))
    expect_equal(scater:::.rowAbove(A.dgc, cols=2:8), Matrix::rowSums(A.dgc[,2:8] > 0))
    expect_equal(scater:::.rowAbove(A.dgc, cols=2:8, value=1), Matrix::rowSums(A.dgc[,2:8] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.dgc), Matrix::colSums(A.dgc))
    expect_equal(scater:::.colAbove(A.dgc), Matrix::colSums(A.dgc > 0))
    expect_equal(scater:::.colAbove(A.dgc, value=1), Matrix::colSums(A.dgc > 1))

    expect_equal(scater:::.colSums(A.dgc, rows=2:8), Matrix::colSums(A.dgc[2:8,]))
    expect_equal(scater:::.colAbove(A.dgc, rows=2:8), Matrix::colSums(A.dgc[2:8,] > 0))
    expect_equal(scater:::.colAbove(A.dgc, rows=2:8, value=1), Matrix::colSums(A.dgc[2:8,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=1:5), colSums(A.int[,1:5]))
    expect_equal(scater:::.colAbove(A.int, cols=1:5), colSums(A.int[,1:5] > 0))
    expect_equal(scater:::.colAbove(A.int, cols=1:5, value=1), colSums(A.int[,1:5] > 1))

})

test_that("row/column variances are correctly computed", {

    ## For integer base matrices.
    A.int <- matrix(rpois(150, lambda = 2), nrow = 15)

    # For row variances.
    expect_equal(scater:::.rowVars(A.int), apply(A.int, 1, var))
    expect_equal(scater:::.rowVars(A.int, rows = 1:5), apply(A.int[1:5,], 1, var))
    expect_equal(scater:::.rowVars(A.int, cols = 2:8), apply(A.int[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.colVars(A.int), apply(A.int, 2, var))
    expect_equal(scater:::.colVars(A.int, rows = 2:8), apply(A.int[2:8,], 2, var))
    expect_equal(scater:::.colVars(A.int, cols = 1:5), apply(A.int[,1:5], 2, var))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(150, 10, 10), nrow = 15)

    # For row variances.
    expect_equal(scater:::.rowVars(A.dbl), apply(A.dbl, 1, var))
    expect_equal(scater:::.rowVars(A.dbl, rows = 1:5), apply(A.dbl[1:5,], 1, var))
    expect_equal(scater:::.rowVars(A.dbl, cols = 2:8), apply(A.dbl[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.colVars(A.dbl), apply(A.dbl, 2, var))
    expect_equal(scater:::.colVars(A.dbl, rows = 2:8), apply(A.dbl[2:8,], 2, var))
    expect_equal(scater:::.colVars(A.dbl, cols = 1:5), apply(A.dbl[,1:5], 2, var))

    ## For sparse matrices.
    A.dgc <- rsparsematrix(15, 10, density = 0.1)
    A.dgc2 <- as.matrix(A.dgc)

    # For row variances.
    expect_equal(scater:::.rowVars(A.dgc), apply(A.dgc2, 1, var))
    expect_equal(scater:::.rowVars(A.dgc, rows = 1:5), apply(A.dgc2[1:5,], 1, var))
    expect_equal(scater:::.rowVars(A.dgc, cols = 2:8), apply(A.dgc2[,2:8], 1, var))
     
    # For column variances.
    expect_equal(scater:::.colVars(A.dgc), apply(A.dgc2, 2, var))
    expect_equal(scater:::.colVars(A.dgc, rows = 2:8), apply(A.dgc2[2:8,], 2, var))
    expect_equal(scater:::.colVars(A.dgc, cols = 1:5), apply(A.dgc2[,1:5], 2, var))
})
