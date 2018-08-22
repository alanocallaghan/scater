## test calculation of average expression
## library(scater); library(testthat); source("setup-sce.R"); source("test-calc-ave.R")

original <- sce

set.seed(10000)
test_that("calcAverage works as expected", {
    ave_counts <- calcAverage(original)
    lib.sizes <- colSums(counts(original))
    expected_vals <- colMeans(t(counts(original)) / (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    expect_equal(ave_counts, calcAverage(counts(original)))

    ## Repeating with subsets.
    sub1 <- calcAverage(counts(original), subset_row=1:10)
    expect_identical(sub1, calcAverage(counts(original)[1:10,]))

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calcAverage(counts(original), subset_row=logi)
    expect_identical(sub2, calcAverage(counts(original)[logi,]))

    chosen <- sample(rownames(original), 20)
    sub3 <- calcAverage(counts(original), subset_row=chosen)
    expect_identical(sub3, calcAverage(counts(original)[chosen,]))
})

test_that("calcAverage responds to size factor options", {
    ave_counts <- calcAverage(original)

    ## Responsive to size factors.
    sizeFactors(original) <- runif(ncol(original))
    ave_counts <- calcAverage(original)
    expect_equal(ave_counts, calcAverage(counts(original), use_size_factors=sizeFactors(original)))
    
    sf <- sizeFactors(original) 
    sf <- sf/mean(sf)    
    expected_vals <- colMeans(t(counts(original)) / sf)
    expect_equal(ave_counts, expected_vals)

    ## Responsive to multiple size factors.
    spiked <- original
    sizeFactors(spiked, "WHEE") <- runif(ncol(original))
    is_spike <- 10:20
    isSpike(spiked, "WHEE") <- is_spike
    ave_counts <- calcAverage(spiked)
    expect_equal(ave_counts[is_spike], calcAverage(counts(original)[is_spike,], use_size_factors=sizeFactors(spiked, "WHEE")))
    expect_equal(ave_counts[-is_spike], calcAverage(counts(original)[-is_spike,], use_size_factors=sizeFactors(spiked)))

    # Ignores or overrides the size factors if requested.
    expect_equal(calcAverage(counts(original)), calcAverage(original, use_size_factors=FALSE))
    expect_equal(calcAverage(counts(spiked)), calcAverage(spiked, use_size_factors=FALSE))

    new_sf <- runif(ncol(spiked))
    ave_counts <- calcAverage(spiked, use_size_factors=new_sf)
    spiked2 <- spiked
    sizeFactors(spiked2) <- new_sf
    expect_equal(calcAverage(spiked2), ave_counts)

    # Warnings are correctly thrown (or not) when no size factors are available for spike-ins.
    sizeFactors(spiked2, "WHEE") <- NULL
    expect_warning(calcAverage(spiked2), "spike-in set")
    expect_warning(calcAverage(spiked2, use_size_factors=new_sf), "spike-in set")
    expect_warning(calcAverage(spiked2, use_size_factors=FALSE), NA)
})

test_that("calcAverage responds to other choices", {
    ave_counts <- calcAverage(original)

    ## Responsive to other assay names.
    whee <- original
    assayNames(whee) <- "whee"
    whee_counts <- calcAverage(whee, exprs_values="whee")
    expect_identical(whee_counts, ave_counts)

    ## Responsive to paralellization.
    expect_identical(ave_counts, calcAverage(original, BPPARAM=MulticoreParam(2)))
    expect_identical(ave_counts, calcAverage(original, BPPARAM=SnowParam(3)))

    ## Repeating with a sparse matrix.    
    sparsified <- original
    counts(sparsified) <- as(counts(original), "dgCMatrix")
    ave_counts <- calcAverage(sparsified)
    expect_equal(ave_counts, calcAverage(original, use_size_factors=FALSE))  
})

test_that("calcAverage handles silly inputs", {
    expect_equivalent(calcAverage(original[0,]), numeric(0)) 
    expect_equivalent(calcAverage(original[,0]), rep(NaN, nrow(original)))
    expect_equivalent(calcAverage(original, use_size_factors=1), rowMeans(counts(original))) # rep()'s 1 to the number of cells.
})
