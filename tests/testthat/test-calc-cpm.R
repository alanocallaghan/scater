## Test CPM calculations and related functions.
## library(scater); library(testthat); source("setup-sce.R"); source("test-calc-cpm.R")

library(Matrix)
original <- sce
sparsified <- original
counts(sparsified) <- as(counts(original), "dgCMatrix")

test_that("we can calculate CPM from counts", {
    cpm_out <- calculateCPM(original)
    expect_equal(cpm_out, t(t(counts(original))/(colSums(counts(original)/1e6))))
    expect_identical(cpm_out, calculateCPM(counts(original)))

    ## Repeating with subsets.
    sub1 <- calculateCPM(counts(original), subset.row=1:10)
    expect_identical(sub1, calculateCPM(counts(original)[1:10,]))

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calculateCPM(counts(original), subset.row=logi)
    expect_identical(sub2, calculateCPM(counts(original)[logi,]))

    chosen <- sample(rownames(original), 20)
    sub3 <- calculateCPM(counts(original), subset.row=chosen)
    expect_identical(sub3, calculateCPM(counts(original)[chosen,]))
})

test_that("calculateCPM is responsive to size factors", {
    sizeFactors(original) <- runif(ncol(original))
    cpm_out <- calculateCPM(original)
    expect_equal(cpm_out, calculateCPM(counts(original), size.factors=sizeFactors(original)))
   
    FUN <- function(counts, sf, libsize = colSums(counts)) {
        eff_lib <- sf/mean(sf) * mean(libsize)
        t(t(counts) / (eff_lib/1e6))
    } 
    expect_equal(cpm_out, FUN(counts(original), sizeFactors(original)))

    # Ignores or overrides the size factors if requested.
    new_sf <- runif(ncol(original))
    cpm_out <- calculateCPM(original, size.factors=new_sf)
    expect_equal(calculateCPM(counts(original), new_sf), cpm_out)
})

test_that("calculateCPM works with alternative inputs", {
    # Checking that it works on a sparse matrix. 
    cpm_out <- calculateCPM(sparsified)
    expect_equal(as.matrix(cpm_out), calculateCPM(original))

    # Checking that it works when there are no columns or rows.
    cpm_out <- calculateCPM(original[0,])
    expect_identical(dim(cpm_out), c(0L, ncol(original)))
    cpm_out <- calculateCPM(original[,0])
    expect_identical(dim(cpm_out), c(nrow(original), 0L))
})

test_that("we can calculate FPKM from counts", {
    effective_length <- runif(nrow(original), 1000, 2000)
    fpkms <- calculateFPKM(original, effective_length)

    ref <- counts(original)/(effective_length/1e3)
    ref <- t(t(ref)/(colSums(counts(original))/1e6))
    expect_equal(fpkms, ref)

    # Repeating with subsets.
    out <- calculateFPKM(original, effective_length, subset.row=1:10)
    sub <- calculateFPKM(original[1:10,], effective_length[1:10])
    expect_equal(out, sub)

    # Handles sparse matrices.
    tout2 <- calculateFPKM(sparsified, effective_length)
    expect_s4_class(tout2, "dgCMatrix")
    expect_equal(ref, as.matrix(tout2))
})

test_that("we can calculate TPM from counts", {
    effective_length <- runif(nrow(original), 1000, 2000)
    tout <- calculateTPM(original, effective_length)

    ref <- counts(original)/effective_length
    ref <- t(t(ref)/(colSums(ref)/1e6))
    expect_equal(tout, ref)

    # Behaves when length is not supplied.
    expect_equal(calculateTPM(original, NULL), calculateCPM(original))

    # Repeating with subsets.
    out <- calculateTPM(original, effective_length, subset.row=1:10)
    sub <- calculateTPM(original[1:10,], effective_length[1:10])
    expect_equal(out, sub)

    # Handles sparse matrices.
    tout2 <- calculateTPM(sparsified, effective_length)
    expect_s4_class(tout2, "dgCMatrix")
    expect_equal(tout, as.matrix(tout2))
})
