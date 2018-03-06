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
    
    ## repeating with different set of size factors
    ref <- runif(ncells, 10, 20)
    sizeFactors(X) <- ref
    out <- normalize(X)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1)) 

    ## checking that size factor centering works properly
    expect_equivalent(sf, sizeFactors(out))
    Xb <- X
    sizeFactors(Xb) <- ref
    outb <- normalize(Xb, centre_size_factors=FALSE)
    expect_equivalent(ref, sizeFactors(outb))
    expect_equivalent(exprs(out), exprs(outb))
    
    ## warning if no size factors or other silly inputs.
    sizeFactors(Xb) <- NULL
    expect_warning(normalize(Xb), "using library sizes")

    expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
    expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 
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

    expect_equivalent(sizeFactors(X4, type="whee"), sf)
    X4b <- normalize(X, centre_size_factors=FALSE)
    expect_equivalent(sizeFactors(X4b, type="whee"), sizeFactors(X, type="whee"))
    expect_equivalent(exprs(X4), exprs(X4b))
})

test_that("scater::normalize works with different settings", {
    out <- normalize(X, log_exprs_offset=3)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))

    metadata(X)$log.exprs.offset <- 3
    out2 <- normalize(X)
    expect_equal(exprs(out), exprs(out2))

    # Checking return_log=FALSE (prior count should turn off automatically).
    out <- normalize(X, return_log=FALSE)
    expect_equivalent(normcounts(out), t(t(dummy)/sf))

    out2 <- normalize(X, return_log=FALSE, log_exprs_offset=3)
    expect_equal(normcounts(out), normcounts(out2))

    # Check with no centering.
    expect_equal(sf, sizeFactors(out))
    out <- normalize(X, centre_size_factors=FALSE)
    expect_equal(ref, sizeFactors(out))
})
