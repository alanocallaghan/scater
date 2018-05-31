# Tests for normalisation methods
# library(scater); library(testthat); source("test-normalisation.R")

context("test expected usage")

#test_that("normaliseExprs does not fail on input with zero-variance features", {
#    data("sc_example_counts")
#    data("sc_example_cell_info")
#    example_sce <- SingleCellExperiment(
#        assays = list(counts = sc_example_counts), 
#        colData = sc_example_cell_info)
#    expect_that(normaliseExprs(example_sce, method = "none", 
#                                feature_set = 1:100), 
#                is_a("SingleCellExperiment"))
#})
#
#test_that("we can compute normalised expression values with TMM method", {
#    data("sc_example_counts")
#    data("sc_example_cell_info")
#    example_sce <- SingleCellExperiment(
#        assays = list(counts = sc_example_counts), 
#        colData = sc_example_cell_info)
#    keep_gene <- rowSums(counts(example_sce)) > 0
#    example_sce <- example_sce[keep_gene,]
#    
#    example_sce <- normaliseExprs(example_sce, method = "TMM", 
#                                     feature_set = 1:100)
#    
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#})
#
#test_that("we can compute normalised expression values with RLE method", {
#    data("sc_example_counts")
#    data("sc_example_cell_info")
#    example_sce <- SingleCellExperiment(
#        assays = list(counts = sc_example_counts), 
#        colData = sc_example_cell_info)
#    keep_gene <- rowSums(counts(example_sce)) > 0
#    example_sce <- example_sce[keep_gene,]
#    
#    example_sce <- normaliseExprs(example_sce, method = "RLE", 
#                                     feature_set = 1:100)
#    
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#})

# test_that("we can compute normalised expression values with upperquartile 
#           method", {
#     data("sc_example_counts")
#     data("sc_example_cell_info")
#     pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#     example_sce <- newSCESet(countData = sc_example_counts + 1, 
#                                 phenoData = pd)
#     keep_gene <- rowSums(counts(example_sce)) > 0
#     example_sce <- example_sce[keep_gene,]
#     example_sce <- example_sce[
#         matrixStats::rowVars(counts(example_sce)) > 0,]
#     
#     example_sce <- normaliseExprs(example_sce, method = "upperquartile", 
#                                      feature_set = 1:200)
#     
#     expect_that(example_sce, is_a("SCESet"))
# })

#test_that("we can compute normalised expression values with none method", {
#    data("sc_example_counts")
#    data("sc_example_cell_info")
#    example_sce <- SingleCellExperiment(
#        assays = list(counts = sc_example_counts), 
#        colData = sc_example_cell_info)
#    keep_gene <- rowSums(counts(example_sce)) > 0
#    example_sce <- example_sce[keep_gene,]
#    
#    example_sce <- normaliseExprs(example_sce, method = "none", 
#                                     feature_set = 1:100)
#    
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#})
#
#test_that("we can compute normalised expression values with a design matrix", {
#    data("sc_example_counts")
#    data("sc_example_cell_info")
#    example_sce <- SingleCellExperiment(
#        assays = list(counts = sc_example_counts), 
#        colData = sc_example_cell_info)
#    exprs(example_sce) <- log2(calculateCPM(example_sce, 
#                                            use_size_factors = FALSE) + 1)
#    keep_gene <- rowSums(counts(example_sce)) > 0
#    example_sce <- calculateQCMetrics(example_sce[keep_gene,], 
#                                         feature_controls = list(set1 = 1:40))
#    design <- model.matrix(~example_sce$Cell_Cycle +
#                               example_sce$pct_counts_top_200_features +
#                               example_sce$total_features)
#    example_sce <- normaliseExprs(example_sce, method = "none", design = design)
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#    
#    example_sce <- normaliseExprs(example_sce, method = "TMM", 
#                                     design = design)
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#    
#    example_sce <- normaliseExprs(example_sce, method = "RLE", 
#                                     design = design)
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#    
#    example_sce <- normaliseExprs(example_sce, exprs_values = "exprs",
#                                  design = design)
#    expect_that(example_sce, is_a("SingleCellExperiment"))
#    
#})

####################################################################################################
# Checking out the behaviour of the normalize() function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))

X <- SingleCellExperiment(list(counts=dummy))
ref <- colSums(dummy)
sizeFactors(X) <- ref

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

    # All size factors are ignored when use_size_factors=FALSE.
    X5 <- X
    sizeFactors(X5) <- NULL
    sizeFactors(X5, "whee") <- NULL
    expect_warning(outd <- normalize(X5), "using library sizes")
    expect_equal(logcounts(outd), log2(t(t(counts(X5))/librarySizeFactors(X5)+1)))
})

test_that("scater::normalize works with different settings", {
    ## Responds to differences in the prior count.
    out <- normalize(X, log_exprs_offset=3)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))

    Y <- X
    metadata(Y)$log.exprs.offset <- 3
    out2 <- normalize(Y)
    expect_equal(exprs(out), exprs(out2))

    # Checking return_log=FALSE (prior count should turn off automatically).
    out <- normalize(X, return_log=FALSE)
    expect_equivalent(normcounts(out), t(t(dummy)/sf))

    out2 <- normalize(X, return_log=FALSE, log_exprs_offset=3)
    expect_equal(normcounts(out), normcounts(out2))

    # Checking that we get sparse matrices out.
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

test_that("scater:normalize works with alternative size factor settings", {
    # No centering.
    out <- normalize(X, centre_size_factors=FALSE)
    expect_equal(ref, sizeFactors(out))

    # Manual centering.
    out <- normalize(X)
    Xb <- centreSizeFactors(X)
    expect_equal(sizeFactors(Xb), sizeFactors(out))
    expect_equal(out, normalize(Xb))

    # Size factor grouping.
    grouping <- sample(5, ncol(X), replace=TRUE)
    out6 <- normalize(X, size_factor_grouping = grouping)

    Xb <- X
    sf.mod <- sizeFactors(Xb) 
    for (x in unique(grouping)) {
        current <- sf.mod[x==grouping]
        sf.mod[x==grouping] <- current/mean(current)
    }
    sizeFactors(Xb) <- sf.mod
    expect_equal(out6, normalize(Xb))

    # Manual grouping.
    Xb <- centreSizeFactors(X, grouping = grouping)
    expect_equal(sizeFactors(Xb), sizeFactors(out6))
    expect_equal(normalize(Xb), out6)
})  


