# Tests for normalisation methods
# library(scater); library(testthat); source("test-normalisation.R")

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))

X <- SingleCellExperiment(list(counts=dummy))
ref <- colSums(dummy)
sizeFactors(X) <- ref

#######################################################

test_that("normalizeCounts works as expected", {
    out <- normalizeCounts(dummy, ref)
    expect_equivalent(out, log2(t(t(dummy)/ref)+1))

    # With size factor centering.
    sf <- ref/mean(ref)
    out <- normalizeCounts(dummy, ref, centre_size_factors=TRUE)
    expect_equivalent(out, log2(t(t(dummy)/sf)+1))

    # Without log-transformation.
    out <- normalizeCounts(dummy, ref, return_log=FALSE)
    expect_equivalent(out, t(t(dummy)/ref))

    # With subsetting.
    out <- normalizeCounts(dummy, ref, subset_row=1:10)
    sub <- normalizeCounts(dummy[1:10,], ref)
    expect_equivalent(out, sub)

    chosen <- sample(rownames(dummy), 10)
    out <- normalizeCounts(dummy, ref, subset_row=chosen)
    sub <- normalizeCounts(dummy[chosen,], ref)
    expect_equivalent(out, sub)

    # Handles silly inputs correctly.
    out <- normalizeCounts(dummy[0,,drop=FALSE], ref, return_log=FALSE)
    expect_identical(dim(out), c(0L, as.integer(ncells)))
    out <- normalizeCounts(dummy[,0,drop=FALSE], ref[0], return_log=FALSE)
    expect_identical(dim(out), c(as.integer(ngenes), as.integer(0L)))

    expect_error(normalizeCounts(dummy, ref[0], return_log=FALSE), "does not equal")
})

#######################################################

areSizeFactorsCentred <- function(object, centre=1, tol=1e-6) {
    all.sf.sets <- c(list(NULL), as.list(sizeFactorNames(object)))
    for (sfname in all.sf.sets) {
        sf <- sizeFactors(object, type=sfname)
        if (!is.null(sf) && abs(mean(sf) - centre) > tol) {
            return(FALSE)
        }
    }
    return(TRUE)
}

test_that("scater::normalize works on endogenous genes", {
    out <- normalize(X)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1))

    expect_equivalent(sf, sizeFactors(out)) # checking that size factor centering works properly
    expect_false(areSizeFactorsCentred(X))
    expect_true(areSizeFactorsCentred(out))
    
    ## repeating with different set of size factors
    ref <- runif(ncells, 10, 20)
    sizeFactors(X) <- ref
    out <- normalize(X)
    sf <- ref/mean(ref)

    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1)) 
    expect_equivalent(sf, sizeFactors(out)) # again, centred size factors.
    expect_false(areSizeFactorsCentred(X))
    expect_true(areSizeFactorsCentred(out))
 
    ## Doesn't break on silly inputs.
    expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
    expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 
})

test_that("scater::normalize works with library sizes", {
    sizeFactors(X) <- NULL
    expect_warning(outb <- normalize(X), "using library sizes")

    lib.sizes <- colSums(counts(X))
    lib.sf <- librarySizeFactors(X)
    expect_equivalent(lib.sf, lib.sizes/mean(lib.sizes))

    expect_equivalent(logcounts(outb), log2(t(t(dummy)/lib.sf)+1))

    # Subsetting by row works correctly in librarySizeFactors().
    expect_identical(librarySizeFactors(X, subset_row=20:1), 
        librarySizeFactors(counts(X)[20:1,]))
})

test_that("scater::normalize works on spike-in genes", {
    out <- normalize(X)
    chosen <- rbinom(ngenes, 1, 0.7)==0L
    isSpike(X, "whee") <- chosen

    ## warning if we don't get any size factors for the spike-ins
    expect_warning(X3 <- normalize(X), "spike-in set 'whee'")
    expect_equal(exprs(out), exprs(X3))

    ## checking that it correctly uses the spike-in size factors
    sizeFactors(X, type="whee") <- colSums(counts(X)[chosen,])
    expect_warning(X4 <- normalize(X), NA) # i.e., no warning.
    expect_equivalent(exprs(out)[!chosen,], exprs(X4)[!chosen,])

    ref <- sizeFactors(X, type="whee")
    sf <- ref/mean(ref)
    expect_equivalent(exprs(X4)[chosen,], log2(t(t(dummy[chosen,])/sf)+1))

    # Checking that the spike-in size factors are correctly centered.
    expect_equivalent(sizeFactors(X4, type="whee"), sf)
    expect_false(areSizeFactorsCentred(X))
    expect_true(areSizeFactorsCentred(X4)) 

    # Without centering of the size factors.
    X4b <- normalize(X, centre_size_factors=FALSE)
    expect_equivalent(logcounts(X4b)[!chosen,], log2(t(t(counts(X)[!chosen,])/sizeFactors(X))+1))
    expect_equivalent(logcounts(X4b)[chosen,], log2(t(t(counts(X)[chosen,])/sizeFactors(X, "whee"))+1))
    expect_equivalent(sizeFactors(X4b), sizeFactors(X))
    expect_equivalent(sizeFactors(X4b, type="whee"), sizeFactors(X, type="whee"))
    expect_false(areSizeFactorsCentred(X4b))

    # Defaults to library sizes when no size factors are available.
    X5 <- X
    sizeFactors(X5) <- NULL
    sizeFactors(X5, "whee") <- NULL
    expect_warning(outd <- normalize(X5), "using library sizes")
    expect_equal(logcounts(outd), log2(t(t(counts(X5))/librarySizeFactors(X5)+1)))
})

#######################################################

test_that("scater::normalize responds to changes in the prior count", {
    sf <- ref/mean(ref)

    ## Responds to differences in the prior count.
    out <- normalize(X, log_exprs_offset=3)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))

    Y <- X
    metadata(Y)$log.exprs.offset <- 3
    out2 <- normalize(Y)
    expect_equal(exprs(out), exprs(out2))

    # Preserves sparsity if requested.
    out3 <- normalize(X, log_exprs_offset=3, preserve_zeroes=TRUE)
    sf3 <- ref/mean(ref) * 3
    expect_equivalent(exprs(out3), log2(t(t(dummy)/sf3)+1))
    expect_equivalent(exprs(out3), exprs(out2) - log2(3))
    expect_equal(sizeFactors(out3), sf3)
})

test_that("scater:normalize works with alternative size factor settings", {
    # No centering.
    out <- normalize(X, centre_size_factors=FALSE)
    expect_equal(ref, sizeFactors(out))
    expect_equal(logcounts(out), log2(t(t(counts(X))/ref+1)))

    # Manual centering.
    out <- normalize(X)
    Xb <- centreSizeFactors(X)
    expect_equal(sizeFactors(Xb), sizeFactors(out))
    expect_equal(out, normalize(Xb, centre_size_factors=FALSE))

    # No centering _and_ non-unity pseudo count.
    out <- normalize(X, log_exprs_offset=3, centre_size_factors=FALSE)
    expect_equivalent(exprs(out), log2(t(t(dummy)/ref)+3))

    out <- normalize(X, log_exprs_offset=3, centre_size_factors=FALSE, preserve_zeroes=TRUE)
    expect_equivalent(exprs(out), log2(t(t(counts(X))/(ref*3)) + 1))
})

#######################################################

test_that("scater::normalize can return un-logged values", {
    sf <- ref/mean(ref)

    # Checking return_log=FALSE (prior count should turn off automatically).
    out <- normalize(X, return_log=FALSE)
    expect_equivalent(normcounts(out), t(t(dummy)/sf))

    out2 <- normalize(X, return_log=FALSE, log_exprs_offset=3)
    expect_equal(normcounts(out), normcounts(out2))
})

test_that("scater::normalize preserves sparsity", {
    Y <- X
    library(Matrix)
    counts(Y) <- as(counts(X), "dgCMatrix")
    out <- normalize(Y)
    expect_s4_class(logcounts(out), "dgCMatrix")
    expect_equal(as.matrix(logcounts(out)), logcounts(normalize(X)))
    
    out2 <- normalize(Y, return_log=FALSE)
    expect_s4_class(normcounts(out2), "dgCMatrix")
    expect_equal(as.matrix(normcounts(out2)), normcounts(normalize(X, return_log=FALSE)))
})

test_that("scater::normalize works with other exprs_values", {
    ref <- X
    ref <- normalize(ref)

    Y2 <- X
    assayNames(Y2) <- "whee"
    Y2 <- normalize(Y2, exprs_values="whee")
    expect_identical(logcounts(ref), logcounts(Y2))

    # Checking that exprs_values= gets passed to librarySizeFactors. 
    sizeFactors(Y2) <- NULL
    expect_warning(Y2 <- normalize(Y2, exprs_values="whee"), "library sizes")
    expect_identical(sizeFactors(Y2), librarySizeFactors(ref))
})

