## Testing the C++ code for correctness.
library(Matrix)

set.seed(12938)
NR <- 15
NC <- 10
rsub1 <- c(2, 12, 15, 3, 9) 
rsub2 <- c(1, 4, 2, 1, 6, 8, 1)
rsub3 <- c(9, 6, 5, 3)
csub1 <- c(2, 3, 7, 3, 9) 
csub2 <- c(4, 6, 4, 10, 5, 10, 7)
csub3 <- c(10, 9, 5, 6)
    
test_that("margin summaries are correctly computed", {

    ## For integer base matrices.
    A.int <- matrix(rpois(NR*NC, lambda=2), nrow=NR, ncol=NC)

    # For row sums
    expect_equal(scater:::.rowSums(A.int), rowSums(A.int))
    expect_equal(scater:::.rowAbove(A.int), rowSums(A.int > 0))
    expect_equal(scater:::.rowAbove(A.int, value=1), rowSums(A.int > 1))

    expect_equal(scater:::.rowSums(A.int, rows=rsub1), rowSums(A.int[rsub1,]))
    expect_equal(scater:::.rowAbove(A.int, rows=rsub2), rowSums(A.int[rsub2,] > 0))
    expect_equal(scater:::.rowAbove(A.int, rows=rsub3, value=1), rowSums(A.int[rsub3,] > 1))

    expect_equal(scater:::.rowSums(A.int, cols=csub1), rowSums(A.int[,csub1]))
    expect_equal(scater:::.rowAbove(A.int, cols=csub2), rowSums(A.int[,csub2] > 0))
    expect_equal(scater:::.rowAbove(A.int, cols=csub3, value=1), rowSums(A.int[,csub3] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.int), colSums(A.int))
    expect_equal(scater:::.colAbove(A.int), colSums(A.int > 0))
    expect_equal(scater:::.colAbove(A.int, value=1), colSums(A.int > 1))

    expect_equal(scater:::.colSums(A.int, rows=rsub1), colSums(A.int[rsub1,]))
    expect_equal(scater:::.colAbove(A.int, rows=rsub2), colSums(A.int[rsub2,] > 0))
    expect_equal(scater:::.colAbove(A.int, rows=rsub3, value=1), colSums(A.int[rsub3,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=csub1), colSums(A.int[,csub1]))
    expect_equal(scater:::.colAbove(A.int, cols=csub2), colSums(A.int[,csub2] > 0))
    expect_equal(scater:::.colAbove(A.int, cols=csub3, value=1), colSums(A.int[,csub3] > 1))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(NR*NC, 10, 10), nrow=NR, ncol=NC)

    # For row sums
    expect_equal(scater:::.rowSums(A.dbl), rowSums(A.dbl))
    expect_equal(scater:::.rowAbove(A.dbl, value=1), rowSums(A.dbl > 1))

    expect_equal(scater:::.rowSums(A.int, rows=rsub1), rowSums(A.int[rsub1,]))
    expect_equal(scater:::.rowAbove(A.int, rows=rsub2, value=1), rowSums(A.int[rsub2,] > 1))

    expect_equal(scater:::.rowSums(A.dbl, cols=csub1), rowSums(A.dbl[,csub1]))
    expect_equal(scater:::.rowAbove(A.dbl, cols=csub2, value=1), rowSums(A.dbl[,csub2] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.dbl), colSums(A.dbl))
    expect_equal(scater:::.colAbove(A.dbl, value=1), colSums(A.dbl > 1))

    expect_equal(scater:::.colSums(A.dbl, rows=rsub1), colSums(A.dbl[rsub1,]))
    expect_equal(scater:::.colAbove(A.dbl, rows=rsub2, value=1), colSums(A.dbl[rsub2,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=csub1), colSums(A.int[,csub1]))
    expect_equal(scater:::.colAbove(A.int, cols=csub2, value=1), colSums(A.int[,csub2] > 1))

    ## For sparse matrices
    A.dgc <- rsparsematrix(NR, NC, density = 0.1)

    # For row sums
    expect_equal(scater:::.rowSums(A.dgc), Matrix::rowSums(A.dgc))
    expect_equal(scater:::.rowAbove(A.dgc), Matrix::rowSums(A.dgc > 0))
    expect_equal(scater:::.rowAbove(A.dgc, value=1), Matrix::rowSums(A.dgc > 1))

    expect_equal(scater:::.rowSums(A.int, rows=rsub1), rowSums(A.int[rsub1,]))
    expect_equal(scater:::.rowAbove(A.int, rows=rsub2), rowSums(A.int[rsub2,] > 0))
    expect_equal(scater:::.rowAbove(A.int, rows=rsub3, value=1), rowSums(A.int[rsub3,] > 1))

    expect_equal(scater:::.rowSums(A.dgc, cols=csub1), Matrix::rowSums(A.dgc[,csub1]))
    expect_equal(scater:::.rowAbove(A.dgc, cols=csub2), Matrix::rowSums(A.dgc[,csub2] > 0))
    expect_equal(scater:::.rowAbove(A.dgc, cols=csub3, value=1), Matrix::rowSums(A.dgc[,csub3] > 1))

    # For column sums.
    expect_equal(scater:::.colSums(A.dgc), Matrix::colSums(A.dgc))
    expect_equal(scater:::.colAbove(A.dgc), Matrix::colSums(A.dgc > 0))
    expect_equal(scater:::.colAbove(A.dgc, value=1), Matrix::colSums(A.dgc > 1))

    expect_equal(scater:::.colSums(A.dgc, rows=rsub1), Matrix::colSums(A.dgc[rsub1,]))
    expect_equal(scater:::.colAbove(A.dgc, rows=rsub2), Matrix::colSums(A.dgc[rsub2,] > 0))
    expect_equal(scater:::.colAbove(A.dgc, rows=rsub3, value=1), Matrix::colSums(A.dgc[rsub3,] > 1))

    expect_equal(scater:::.colSums(A.int, cols=csub1), colSums(A.int[,csub1]))
    expect_equal(scater:::.colAbove(A.int, cols=csub2), colSums(A.int[,csub2] > 0))
    expect_equal(scater:::.colAbove(A.int, cols=csub3, value=1), colSums(A.int[,csub3] > 1))
})

test_that("row/column variances are correctly computed", {

    ## For integer base matrices.
    A.int <- matrix(rpois(NR*NC, lambda = 2), nrow = NR, ncol=NC)

    # For row variances.
    expect_equal(scater:::.rowVars(A.int), apply(A.int, 1, var))
    expect_equal(scater:::.rowVars(A.int, rows = rsub1), apply(A.int[rsub1,], 1, var))
    expect_equal(scater:::.rowVars(A.int, rows = rsub2), apply(A.int[rsub2,], 1, var))
    expect_equal(scater:::.rowVars(A.int, cols = csub1), apply(A.int[,csub1], 1, var))
    expect_equal(scater:::.rowVars(A.int, cols = csub2), apply(A.int[,csub2], 1, var))
     
    # For column variances.
    expect_equal(scater:::.colVars(A.int), apply(A.int, 2, var))
    expect_equal(scater:::.colVars(A.int, rows = rsub1), apply(A.int[rsub1,], 2, var))
    expect_equal(scater:::.colVars(A.int, rows = rsub2), apply(A.int[rsub2,], 2, var))
    expect_equal(scater:::.colVars(A.int, cols = csub1), apply(A.int[,csub1], 2, var))
    expect_equal(scater:::.colVars(A.int, cols = csub2), apply(A.int[,csub2], 2, var))

    ## For double-precision base matrices.
    A.dbl <- matrix(rgamma(NR*NC, 10, 10), nrow = NR, ncol=NC)

    # For row variances.
    expect_equal(scater:::.rowVars(A.dbl), apply(A.dbl, 1, var))
    expect_equal(scater:::.rowVars(A.dbl, rows = rsub1), apply(A.dbl[rsub1,], 1, var))
    expect_equal(scater:::.rowVars(A.dbl, rows = rsub2), apply(A.dbl[rsub2,], 1, var))
    expect_equal(scater:::.rowVars(A.dbl, cols = csub1), apply(A.dbl[,csub1], 1, var))
    expect_equal(scater:::.rowVars(A.dbl, cols = csub2), apply(A.dbl[,csub2], 1, var))
     
    # For column variances.
    expect_equal(scater:::.colVars(A.dbl), apply(A.dbl, 2, var))
    expect_equal(scater:::.colVars(A.dbl, rows = rsub1), apply(A.dbl[rsub1,], 2, var))
    expect_equal(scater:::.colVars(A.dbl, rows = rsub2), apply(A.dbl[rsub2,], 2, var))
    expect_equal(scater:::.colVars(A.dbl, cols = csub1), apply(A.dbl[,csub1], 2, var))
    expect_equal(scater:::.colVars(A.dbl, cols = csub2), apply(A.dbl[,csub2], 2, var))

    ## For sparse matrices.
    A.dgc <- rsparsematrix(NR, NC, density = 0.1)
    A.dgc2 <- as.matrix(A.dgc)

    # For row variances.
    expect_equal(scater:::.rowVars(A.dgc), apply(A.dgc2, 1, var))
    expect_equal(scater:::.rowVars(A.dgc, rows = rsub1), apply(A.dgc[rsub1,], 1, var))
    expect_equal(scater:::.rowVars(A.dgc, rows = rsub2), apply(A.dgc[rsub2,], 1, var))
    expect_equal(scater:::.rowVars(A.dgc, cols = csub1), apply(A.dgc[,csub1], 1, var))
    expect_equal(scater:::.rowVars(A.dgc, cols = csub2), apply(A.dgc[,csub2], 1, var))
   
    # For column variances.
    expect_equal(scater:::.colVars(A.dgc), apply(A.dgc2, 2, var))
    expect_equal(scater:::.colVars(A.dgc, rows = rsub1), apply(A.dgc[rsub1,], 2, var))
    expect_equal(scater:::.colVars(A.dgc, rows = rsub2), apply(A.dgc[rsub2,], 2, var))
    expect_equal(scater:::.colVars(A.dgc, cols = csub1), apply(A.dgc[,csub1], 2, var))
    expect_equal(scater:::.colVars(A.dgc, cols = csub2), apply(A.dgc[,csub2], 2, var))
})
