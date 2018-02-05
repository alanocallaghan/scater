# Tests the options in the various reduced dimension functions
# library(scater); library(testthat); source("test-reddim.R")

data("sc_example_counts")
data("sc_example_cell_info")
sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)
suppressWarnings(sce <- normalize(sce))

#############################################
# Check PCA.

test_that("runPCA works as expected", {
    sce2 <- runPCA(sce, ncomponents=4)
    expect_identical(reducedDimNames(sce2), "PCA")     
    expect_identical(dim(reducedDim(sce2, "PCA")), c(ncol(sce2), 4L))
    expect_equal(sum(attr(reducedDim(sce2), "percentVar")), 1)

    # Testing that various settings give different results.
    sce3 <- runTSNE(sce, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    # Testing that ntop selection works correctly.
    most_var <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(sce)))
    keep <- head(order(most_var, decreasing=TRUE), 100)
    sce3 <- runPCA(sce, ncomponents=4, feature_set=keep)
    sce4 <- runPCA(sce, ncomponents=4, ntop=100)
    expect_equal(reducedDim(sce3), reducedDim(sce4))    

    # Testing out the scaling (ntop=Inf, otherwise it will pick features based on scaled variance).
    sce_alt <- sce
    logcounts(sce_alt) <- t(scale(t(logcounts(sce_alt)), scale = TRUE))
    sce3 <- runPCA(sce_alt, ncomponents=4, scale_features=FALSE, ntop=Inf)
    sce4 <- runPCA(sce, ncomponents=4, scale_features=TRUE, ntop=Inf)
    expect_equal(reducedDim(sce3), reducedDim(sce4))    

    # Testing out alternative assays.
    sce_alt <- sce
    assay(sce_alt, "whee") <- logcounts(sce)
    logcounts(sce_alt) <- NULL
    sce3 <- runPCA(sce_alt, ncomponents=4, exprs_values="whee")
    expect_equal(reducedDim(sce3), reducedDim(sce2))
})

test_that("runPCA works as expected for QC metrics", {
    sceQ <- calculateQCMetrics(sce, feature_controls=list(ERCC=1:50))
    expect_warning(sceQ <- runPCA(sceQ, use_coldata=TRUE), NA)
    expect_identical(reducedDimNames(sceQ), "PCA_coldata")

    sceQ2 <- calculateQCMetrics(sce, feature_controls=list(ERCC=1:50), compact=TRUE)
    expect_warning(sceQ2 <- runPCA(sceQ2, use_coldata=TRUE), NA)
    expect_identical(reducedDim(sceQ, "PCA_coldata"), reducedDim(sceQ2, "PCA_coldata"))

    # Checking outlier detection works correctly.
    expect_identical(sceQ2$outlier, NULL)
    expect_warning(sceQ2 <- runPCA(sceQ2, use_coldata=TRUE, detect_outliers=TRUE), NA)
    expect_type(sceQ2$outlier, "logical")
    expect_identical(length(sceQ2$outlier), ncol(sceQ2))

    # Check manual specification of variables works properly.
    expect_warning(sceQ3 <- runPCA(sceQ, use_coldata=TRUE, 
        selected_variables = c("total_counts", 
                               "total_features_by_counts",
                               "total_counts_endogenous",
                               "total_counts_feature_control")), NA)
})

#############################################
# Check t-SNE.

test_that("runTSNE works as expected", {
    sce2 <- runTSNE(sce, ncomponents = 3, rand_seed = 100)
    expect_identical(reducedDimNames(sce2), "TSNE")     
    expect_identical(dim(reducedDim(sce2, "TSNE")), c(ncol(sce2), 3L))

    # Testing that rand_seed actually works.
    sce3 <- runTSNE(sce, rand_seed = 100)
    expect_equal(reducedDim(sce2), reducedDim(sce3))

    # Testing that various settings have some effect. 
    sce3 <- runTSNE(sce, scale_features = FALSE, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, ntop = 100, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, exprs_values = "counts", rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, feature_set = 1:100, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))
    
    sce3 <- runTSNE(sce, perplexity = 5, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, pca = FALSE, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runTSNE(sce, initial_dims = 10, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings).
    sceP <- runPCA(sce, ncomponents = 4)
    sce2 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10)

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, ntop = 20)
    expect_identical(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, scale_features = FALSE)
    expect_identical(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, feature_set = 1:20)
    expect_identical(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, pca = FALSE)
    expect_identical(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, initial_dims = 10)
    expect_identical(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))

    sce3 <- runTSNE(sceP, use_dimred = "PCA", rand_seed = 10, n_dimred=3)
    expect_false(isTRUE(all.equal(reducedDim(sce2, "TSNE"), reducedDim(sce3, "TSNE"))))
})

#############################################
# Check MDS.

test_that("runMDS works as expected", {
    sce2 <- runMDS(sce, ncomponents = 3)
    expect_identical(reducedDimNames(sce2), "MDS")     
    expect_identical(dim(reducedDim(sce2, "MDS")), c(ncol(sce2), 3L))

    # Testing that various settings work.
    sce3 <- runMDS(sce, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runMDS(sce, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runMDS(sce, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runMDS(sce, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runMDS(sce, method = "manhattan")
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings).
    sceP <- runPCA(sce, ncomponents = 4)
    sce2 <- runMDS(sceP, use_dimred = "PCA")

    sce3 <- runMDS(sceP, use_dimred = "PCA", ntop = 20)
    expect_identical(reducedDim(sce2, "MDS"), reducedDim(sce3, "MDS"))

    sce3 <- runMDS(sceP, use_dimred = "PCA", scale_features = FALSE)
    expect_identical(reducedDim(sce2, "MDS"), reducedDim(sce3, "MDS"))

    sce3 <- runMDS(sceP, use_dimred = "PCA", feature_set = 1:20)
    expect_identical(reducedDim(sce2, "MDS"), reducedDim(sce3, "MDS"))

    # This does, in fact, happen to be equal, due to the relationship between MDS and PCA.
    sce3 <- runMDS(sceP, use_dimred = "PCA", n_dimred=3)
    expect_equal(reducedDim(sce2, "MDS"), reducedDim(sce3, "MDS"))
})

#############################################
# Check DiffusionMaps, which seems to oscillate the sign of particular components.

SIGNAGNOSTIC <- function(x, y, same = TRUE) {
    ratios <- x/y
    ratios <- t(t(ratios)/colSums(ratios))
    if (same) {
        expect_true(sd(ratios) < 1e-6)
    } else {
        expect_false(sd(ratios) < 1e-6)
    }
}

test_that("runDiffusionMap works as expected", {
    sce2 <- runDiffusionMap(sce, ncomponents = 3, rand_seed = 100)
    expect_identical(reducedDimNames(sce2), "DiffusionMap")     
    expect_identical(dim(reducedDim(sce2, "DiffusionMap")), c(ncol(sce2), 3L))

    # Testing that various settings work.
    sce3 <- runDiffusionMap(sce, scale_features = FALSE, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runDiffusionMap(sce, ntop = 100, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runDiffusionMap(sce, exprs_values = "counts", rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runDiffusionMap(sce, feature_set = 1:100, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))
    
    sce3 <- runDiffusionMap(sce, sigma = 10, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    sce3 <- runDiffusionMap(sce, k = 13, rand_seed = 100)
    expect_false(isTRUE(all.equal(reducedDim(sce2), reducedDim(sce3))))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings)
    sceP <- runPCA(sce, ncomponents = 4)
    sce2 <- runDiffusionMap(sceP, use_dimred = "PCA", rand_seed = 10)

    sce3 <- runDiffusionMap(sceP, use_dimred = "PCA", rand_seed = 10, ntop = 20)
    SIGNAGNOSTIC(reducedDim(sce2, "DiffusionMap"), reducedDim(sce3, "DiffusionMap"))

    sce3 <- runDiffusionMap(sceP, use_dimred = "PCA", rand_seed = 10, scale_features = FALSE)
    SIGNAGNOSTIC(reducedDim(sce2, "DiffusionMap"), reducedDim(sce3, "DiffusionMap"))

    sce3 <- runDiffusionMap(sceP, use_dimred = "PCA", rand_seed = 10, feature_set = 1:20)
    SIGNAGNOSTIC(reducedDim(sce2, "DiffusionMap"), reducedDim(sce3, "DiffusionMap"))

    sce3 <- runDiffusionMap(sceP, use_dimred = "PCA", n_dimred=3)
    SIGNAGNOSTIC(reducedDim(sce2, "DiffusionMap"), reducedDim(sce3, "DiffusionMap"), same = FALSE)

})


