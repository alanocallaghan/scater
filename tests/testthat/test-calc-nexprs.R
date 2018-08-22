## Test the nexprs() function. 
## library(scater); library(testthat); source("setup-sce.R"); source("test-calc-nexprs.R")

original <- sce

test_that("nexprs works as expected", {
    expect_equal(nexprs(original), colSums(counts(original) > 0))
    expect_equal(nexprs(original), nexprs(counts(original)))

    expect_equal(nexprs(original, byrow = TRUE), rowSums(counts(original) > 0))
    expect_equal(nexprs(original, byrow = TRUE), nexprs(counts(original), byrow = TRUE))
})

test_that("nexprs responds to subsetting", {
    expect_equal(nexprs(original, subset_row = 20:40), colSums(counts(original)[20:40,] > 0))
    expect_equal(nexprs(original, byrow = TRUE, subset_col = 20:40), rowSums(counts(original)[,20:40] > 0))

    expect_equal(nexprs(original, subset_row = 20:40, subset_col=1:10), colSums(counts(original)[20:40,1:10] > 0))
    expect_equal(nexprs(original, byrow = TRUE, subset_row=1:10, subset_col = 20:40), rowSums(counts(original)[1:10,20:40] > 0))
})

test_that("nexprs responds to other options", {    
    expect_equal(nexprs(original, detection_limit=5), colSums(counts(original) > 5))
    expect_equal(nexprs(original, byrow = TRUE, detection_limit=5), rowSums(counts(original) > 5))

    # Handles parallelization.
    expect_equal(nexprs(original), nexprs(original, BPPARAM=MulticoreParam(2)))
    expect_equal(nexprs(original), nexprs(original, BPPARAM=SnowParam(3)))
    expect_equal(nexprs(original, byrow=TRUE), nexprs(original, byrow=TRUE, BPPARAM=MulticoreParam(2)))
    expect_equal(nexprs(original, byrow=TRUE), nexprs(original, byrow=TRUE, BPPARAM=SnowParam(3)))
})

test_that("nexprs works on a sparse matrix", {
    sparsified <- original
    counts(sparsified) <- as(counts(original), "dgCMatrix")
    expect_equal(nexprs(sparsified), Matrix::colSums(counts(sparsified) > 0))
    expect_equal(nexprs(sparsified), nexprs(counts(sparsified)))
})

test_that("nexprs handles silly inputs properly", {
    expect_equivalent(nexprs(original, subset_row=integer(0)), integer(ncol(original)))
    expect_equivalent(nexprs(original, subset_col=integer(0)), integer(0))
    expect_equivalent(nexprs(original, subset_row=integer(0), byrow=TRUE), integer(0))
    expect_equivalent(nexprs(original, subset_col=integer(0), byrow=TRUE), integer(nrow(original)))
})

