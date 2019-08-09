# Tests for normalisation methods
# library(scater); library(testthat); source("test-norm.R")

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))

X <- SingleCellExperiment(list(counts=dummy))
ref <- runif(ncells)
sizeFactors(X) <- ref

#######################################################

test_that("normalizeCounts works as expected", {
    out <- normalizeCounts(dummy, ref, center_size_factors=FALSE)
    expect_equal(out, log2(t(t(dummy)/ref)+1))

    # With size factor centering.
    sf <- ref/mean(ref)
    out <- normalizeCounts(dummy, ref)
    expect_equal(out, log2(t(t(dummy)/sf)+1))

    # Without log-transformation.
    out <- normalizeCounts(dummy, ref, log=FALSE)
    expect_equal(out, t(t(dummy)/sf))

    # With subsetting.
    out <- normalizeCounts(dummy, ref, subset_row=1:10)
    sub <- normalizeCounts(dummy[1:10,], ref)
    expect_equal(out, sub)

    out <- normalizeCounts(dummy, subset_row=1:10)
    sub <- normalizeCounts(dummy[1:10,])
    expect_equal(out, sub)

    chosen <- sample(rownames(dummy), 10)
    out <- normalizeCounts(dummy, ref, subset_row=chosen)
    sub <- normalizeCounts(dummy[chosen,], ref)
    expect_equal(out, sub)

    # Handles silly inputs correctly.
    out <- normalizeCounts(dummy[0,,drop=FALSE], ref, log=FALSE)
    expect_identical(dim(out), c(0L, as.integer(ncells)))
    out <- normalizeCounts(dummy[,0,drop=FALSE], ref[0], log=FALSE)
    expect_identical(dim(out), c(as.integer(ngenes), as.integer(0L)))

    expect_error(normalizeCounts(dummy, ref[0]), "does not equal")
    expect_error(normalizeCounts(dummy, rep(0, ncol(dummy))), "should be positive")
    expect_error(normalizeCounts(dummy, rep(NA_real_, ncol(dummy))), "should be positive")
})

test_that("normalizeCounts behaves with sparse inputs", {
    zeroed <- dummy
    zeroed[rbinom(length(zeroed), 1, 0.95)==1] <- 0

    library(Matrix)
    sparsed <- as(zeroed, "dgCMatrix")
   
    expect_s4_class(out <- normalizeCounts(sparsed, ref), "dgCMatrix")
    expect_equal(normalizeCounts(zeroed, ref), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(sparsed, ref, log=FALSE), "dgCMatrix")
    expect_equal(normalizeCounts(zeroed, ref, log=FALSE), as.matrix(out))

    expect_true(is.matrix(out <- normalizeCounts(sparsed, ref, pseudo_count=2)))
    expect_equal(normalizeCounts(zeroed, ref, pseudo_count=2), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(sparsed, ref, subset_row=1:10), "dgCMatrix")
    expect_equal(normalizeCounts(zeroed, ref, subset_row=1:10), as.matrix(out))
})

test_that("normalizeCounts behaves with DelayedArray inputs", {
    library(DelayedArray)
    dadum <- DelayedArray(dummy)      
    
    expect_s4_class(out <- normalizeCounts(dadum, ref), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy, ref), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(dadum, ref, log=FALSE), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy, ref, log=FALSE), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(dadum, ref, pseudo_count=2), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy, ref, pseudo_count=2), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(dadum, ref, subset_row=1:10), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy, ref, subset_row=1:10), as.matrix(out))

    # Library sizes are correctly obtained.
    expect_s4_class(out <- normalizeCounts(dadum), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy), as.matrix(out))

    expect_s4_class(out <- normalizeCounts(dadum, subset_row=1:10), "DelayedMatrix")
    expect_equal(normalizeCounts(dummy, subset_row=1:10), as.matrix(out))
})

test_that("normalizeCounts behaves with S(C)E inputs", {
    expect_equivalent(normalizeCounts(counts(X)), 
        normalizeCounts(as(X, "SummarizedExperiment")))

    sf <- runif(ncol(X))
    expect_equivalent(normalizeCounts(counts(X), sf), 
        normalizeCounts(as(X, "SummarizedExperiment"), sf))

    expect_equal(normalizeCounts(counts(X), sizeFactors(X)), normalizeCounts(X))
    expect_equal(normalizeCounts(counts(X), sf), normalizeCounts(X, sf))
})

#######################################################

test_that("logNormCounts works for SE objects", {
    se <- as(X, "SummarizedExperiment")

    cn <- function(se) assay(se, "counts")
    lc <- function(se) assay(se, "logcounts")
    nc <- function(se) assay(se, "normcounts")

    expect_equal(lc(logNormCounts(se)), normalizeCounts(cn(se)))
    expect_equal(nc(logNormCounts(se, log=FALSE)), normalizeCounts(cn(se), log=FALSE))
    expect_equal(lc(logNormCounts(se, pseudo_count=2)), normalizeCounts(cn(se), pseudo_count=2))

    sf <- runif(ncol(X))
    expect_equal(lc(logNormCounts(se, size_factors=sf)),
        normalizeCounts(cn(se), size_factors=sf))
    expect_equal(lc(logNormCounts(se, size_factors=sf, center_size_factors=FALSE)), 
        normalizeCounts(cn(se), size_factors=sf, center_size_factors=FALSE))
 
    ## Doesn't break on silly inputs.
    expect_equal(unname(dim(logNormCounts(X[,0,drop=FALSE]))), c(ngenes, 0L))
    expect_equal(unname(dim(logNormCounts(X[0,,drop=FALSE]))), c(0L, ncells)) 
})

test_that("logNormCounts works for SCE objects (basic)", {
    expect_equal(logcounts(logNormCounts(X)), normalizeCounts(counts(X), sizeFactors(X)))
    expect_equal(normcounts(logNormCounts(X, log=FALSE)), normalizeCounts(counts(X), sizeFactors(X), log=FALSE))
    expect_equal(logcounts(logNormCounts(X, pseudo_count=2)), normalizeCounts(counts(X), sizeFactors(X), pseudo_count=2))

    # Checking that size factors are correctly reported.
    Y <- X
    sizeFactors(Y) <- NULL
    Y <- logNormCounts(Y)
    expect_identical(sizeFactors(Y), librarySizeFactors(X))

    sf <- runif(ncol(X))
    Y <- logNormCounts(X, size_factors=sf)
    expect_identical(sizeFactors(Y), sf/mean(sf))

    Y <- logNormCounts(X, size_factors=sf, center_size_factors=FALSE)
    expect_identical(sizeFactors(Y), sf)

    # Checking that my pseudo-count appears and does not overwrite other scater stuff.
    expect_identical(int_metadata(Y)$scater$pseudo_count, 1)

    Z <- X
    int_metadata(Z)$scater <- list(whee="YAY")
    Z <- logNormCounts(Z)
    expect_identical(int_metadata(Z)$scater$pseudo_count, 1)
    expect_identical(int_metadata(Z)$scater$whee, "YAY")

    # Diverts to other names.
    Y <- logNormCounts(X, name="blah")
    expect_identical(assay(Y, "blah"), logcounts(logNormCounts(X)))
})

test_that("logNormCounts works for SCE objects (altExp)", {
    sce <- X
    altExp(sce, "BLAH") <- X
    sce1 <- logNormCounts(sce)

    # Do a class round-trip to wipe out metadata added to the int_* fields.
    COMPFUN <- function(left, right) {
        left <- as(left, "SummarizedExperiment")
        left <- as(left, "SingleCellExperiment")
        right <- as(right, "SummarizedExperiment")
        right <- as(right, "SingleCellExperiment")
        expect_equal(left, right)
    }

    COMPFUN(altExp(sce1, "BLAH"), logNormCounts(X))

    # Global size factors are respected.
    sf <- runif(ncol(X))
    sce2 <- logNormCounts(sce, size_factors=sf)
    COMPFUN(altExp(sce2), logNormCounts(X, size_factors=sf))

    # Other parameters are respected.
    sce3a <- logNormCounts(sce, pseudo_count=2)
    COMPFUN(altExp(sce3a), logNormCounts(X, pseudo_count=2))

    sce3b <- logNormCounts(sce, log=FALSE)
    COMPFUN(altExp(sce3b), logNormCounts(X, log=FALSE))

    # Internal size factors do not propagate to alternative experiments.
    sce4 <- sce
    sizeFactors(sce4) <- sf
    sce4 <- logNormCounts(sce4)
    COMPFUN(altExp(sce4), logNormCounts(X))

    # Lack of centering is respected in downstream methods.
    sce5 <- logNormCounts(sce, center_size_factors=FALSE)
    COMPFUN(altExp(sce5), logNormCounts(X, center_size_factors=FALSE))

    # Throws errors with zero-valued size factors.
    sce6 <- sce
    sizeFactors(altExp(sce6)) <- 0
    expect_error(logNormCounts(sce6), 'altExp')
})
