# This tests the batchCorrectedAverages function.
# library(testthat); library(gp.sa.solo); source("test-average.R")

library(scater)

test_that("batchCorrectedAverages works correctly in raw mode", {
    # No batches.
    y <- matrix(rnorm(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- batchCorrectedAverages(y, group, block)
    ref <- sumCountsAcrossCells(y, group, average=TRUE)
    expect_equal(averaged, assay(ref))
 
    # Perfectly balanced.
    y <- matrix(rnorm(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep(1:10, 10)
    averaged <- batchCorrectedAverages(y, group, block)
 
    ref <- sumCountsAcrossCells(y, group, average=TRUE)
    ref <- assay(ref)
 
    ref2 <- ref - rowMeans(ref)
    averaged2 <- averaged - rowMeans(averaged)
    expect_equal(ref2, averaged2)
    
    # Effectively ignores batches with only one cluster.
    y <- matrix(rnorm(210), ncol=21)
    group <- c(1:10, 1:10, 1)
    block <- rep(1:3, c(10, 10, 1))
    averaged <- batchCorrectedAverages(y, group, block)
    averaged2 <- batchCorrectedAverages(y[,-21], group[-21], block[-21])
    expect_equal(averaged, averaged2)
 
    # Handles batch-specific clusters.
    y <- matrix(rnorm(220), ncol=22)
    group <- c(1:10, 1:10, c(1,11))
    block <- rep(1:3, c(10, 10, 2))
    averaged <- batchCorrectedAverages(y, group, block)
    expect_identical(sum(!colAnyNAs(averaged)), 11L)
})

test_that("batchCorrectedAverages works correctly in log mode", {
    # Actually has an effect.
    y <- matrix(rexp(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged1 <- batchCorrectedAverages(y, group, block)
    averaged2 <- batchCorrectedAverages(y, group, block, transform="log")
    expect_false(identical(averaged1, averaged2))

    # Survives a round trip with correct untransformation. 
    y <- matrix(rexp(1000), ncol=10)
    group <- 1:10
    block <- rep('all', ncol(y))
    averaged <- batchCorrectedAverages(y, group, block, transform="log")
    expect_equivalent(averaged, y)

    # Handles zeroes.
    y <- matrix(0, ncol=100, nrow=10)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- batchCorrectedAverages(y, group, block, transform="log")
    expect_true(all(averaged==0))
})

test_that("batchCorrectedAverages works correctly in logit mode", {
    # Actually has an effect.
    y <- matrix(runif(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged1 <- batchCorrectedAverages(y, group, block)
    averaged2 <- batchCorrectedAverages(y, group, block, transform="log")
    averaged3 <- batchCorrectedAverages(y, group, block, transform="logit")
    expect_false(identical(averaged1, averaged3))
    expect_false(identical(averaged2, averaged3))

    # Survives a round trip with correct untransformation. 
    y <- matrix(runif(1000), ncol=10)
    group <- 1:10
    block <- rep('all', ncol(y))
    averaged <- batchCorrectedAverages(y, group, block, transform="logit")
    expect_equivalent(averaged, y)

    # Handles boundaries.
    y <- matrix(0, ncol=100, nrow=10)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- batchCorrectedAverages(y, group, block, transform="logit")
    expect_true(all(abs(averaged) < 1e-10))

    y <- matrix(1, ncol=100, nrow=10)
    averaged <- batchCorrectedAverages(y, group, block, transform="logit")
    expect_true(all(averaged==1))
})

