## test calculation of average expression
## library(scater); library(testthat); source("setup-sce.R"); source("test-calc-ave.R")

original <- sce

set.seed(10000)
test_that("calculateAverage works as expected", {
    ave_counts <- calculateAverage(original)
    lib.sizes <- colSums(counts(original))
    expected_vals <- colMeans(t(counts(original)) / (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    expect_equal(ave_counts, calculateAverage(counts(original)))

    ## Repeating with subsets.
    sub1 <- calculateAverage(counts(original), subset.row=1:10)
    expect_identical(sub1, calculateAverage(counts(original)[1:10,]))

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calculateAverage(counts(original), subset.row=logi)
    expect_identical(sub2, calculateAverage(counts(original)[logi,]))

    chosen <- sample(rownames(original), 20)
    sub3 <- calculateAverage(counts(original), subset.row=chosen)
    expect_identical(sub3, calculateAverage(counts(original)[chosen,]))
})

test_that("calculateAverage responds to size factor options", {
    ave_counts <- calculateAverage(original)

    ## Responsive to size factors.
    sizeFactors(original) <- runif(ncol(original))
    ave_counts <- calculateAverage(original)
    expect_equal(ave_counts, calculateAverage(counts(original), size.factors=sizeFactors(original)))
    
    sf <- sizeFactors(original) 
    sf <- sf/mean(sf)    
    expected_vals <- colMeans(t(counts(original)) / sf)
    expect_equal(ave_counts, expected_vals)

    # Ignores or overrides the size factors if requested.
    new_sf <- runif(ncol(original))
    ave_counts <- calculateAverage(original, size.factors=new_sf)
    expect_equal(calculateAverage(counts(original), size.factors=new_sf), ave_counts)
})

test_that("calculateAverage responds to other choices", {
    ave_counts <- calculateAverage(original)

    ## Responsive to other assay names.
    whee <- original
    assay(whee, "whee") <- counts(original)*2
    whee_counts <- calculateAverage(whee, assay.type="whee")
    expect_identical(whee_counts, ave_counts*2)

    ## Responsive to parallelization.
    expect_equal(ave_counts, calculateAverage(original, BPPARAM=safeBPParam(2)))
    expect_equal(ave_counts, calculateAverage(original, BPPARAM=safeBPParam(3)))

    ## Repeating with a sparse matrix, to check that the specialized code is correct.
    sparsified <- original
    counts(sparsified) <- as(counts(original), "dgCMatrix")
    expect_equal(ave_counts, calculateAverage(sparsified))
    expect_equal(calculateAverage(original, subset.row=30:20),
        calculateAverage(sparsified, subset.row=30:20))

    unknown <- original
    counts(unknown) <- as(counts(original), "dgTMatrix")
    expect_equal(ave_counts, calculateAverage(unknown))
    expect_equal(calculateAverage(unknown, subset.row=25:15),
         calculateAverage(original, subset.row=25:15))
})

test_that("calculateAverage handles silly inputs", {
    expect_equivalent(calculateAverage(original[0,]), numeric(0)) 
    expect_equivalent(calculateAverage(original[,0]), rep(NaN, nrow(original)))
})
