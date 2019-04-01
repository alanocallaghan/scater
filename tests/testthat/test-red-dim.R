# Tests the options in the various reduced dimension functions
# library(scater); library(testthat); source("setup-sce.R"); source("test-red-dim.R")

#############################################
# Check the feature selection and scaling work.

test_that("feature selection is operational", {
    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", feature_set = seq_len(nrow(normed)))
    expect_equal(out, t(logcounts(normed, withDimnames=FALSE)))

    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "counts", feature_set = seq_len(nrow(normed)))
    expect_equal(out, t(counts(normed, withDimnames=FALSE)))

    # Ntop selection works.
    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", ntop = 10, feature_set = NULL)
    rv <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    keep <- head(order(rv, decreasing=TRUE), 10)
    expect_equal(out, t(logcounts(normed, withDimnames=FALSE)[keep,]))

    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", ntop = Inf, feature_set = NULL)
    o <- order(rv, decreasing=TRUE)
    expect_equal(out, t(logcounts(normed, withDimnames=FALSE)[o,]))
    
    # Feature selection works.
    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", feature_set = 10:1)
    expect_equal(out, t(logcounts(normed, withDimnames=FALSE)[10:1,]))

    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", feature_set = rownames(normed)[10:1])
    expect_equal(out, t(logcounts(normed, withDimnames=FALSE)[10:1,]))
})

test_that("scaling by feature variances work correctly", {
    MAT <- t(logcounts(normed, withDimnames=FALSE))
    cv <- DelayedMatrixStats::colVars(DelayedArray(MAT))
    novar <- cv < 1e-8
    expect_true(any(novar))

    NOMAT <- MAT[,!novar]
    XX <- scater:::.scale_columns(MAT)
    scaled <- t(t(NOMAT)/sqrt(cv[!novar]))
    expect_equal(XX, scaled)

    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", feature_set = seq_len(nrow(normed)), scale = TRUE)
    expect_equal(out, scaled)
    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", feature_set = 10:1, scale = TRUE)
    expect_equal(out, scaled[,10:1])

    out <- scater:::.get_mat_for_reddim(normed, exprs_values = "logcounts", ntop = 10, scale = TRUE) # In combination with non-trivial selection.
    rv <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    expect_equal(out, scaled[,head(order(rv[!novar], decreasing=TRUE), 10)])
})

#############################################
# Check PCA.

test_that("runPCA works as expected", {
    normedX <- runPCA(normed, ncomponents=4)
    expect_identical(reducedDimNames(normedX), "PCA")     
    expect_identical(dim(reducedDim(normedX, "PCA")), c(ncol(normedX), 4L))
    expect_true(sum(attr(reducedDim(normedX), "percentVar")) < 1)

    # Testing that various settings give different results.
    normed2 <- runPCA(normed)
    normed3 <- runPCA(normed, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
    expect_true(sum(attr(reducedDim(normed2), "percentVar")) < 1)
    expect_true(sum(attr(reducedDim(normed3), "percentVar")) < 1)

    normed3 <- runPCA(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runPCA(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runPCA(normed, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    # Testing that ntop selection works correctly.
    most_var <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    keep <- head(order(most_var, decreasing=TRUE), 100)
    normed3 <- runPCA(normed, ncomponents=4, feature_set=keep)
    normed4 <- runPCA(normed, ncomponents=4, ntop=100)
    expect_equal(reducedDim(normed3), reducedDim(normed4))    

    # Testing out the scaling (ntop=Inf, otherwise it will pick features based on scaled variance).
    normed_alt <- normed
    rescaled <- t(scale(t(logcounts(normed_alt)), scale = TRUE))
    rescaled[is.na(rescaled)] <- 0
    logcounts(normed_alt) <- rescaled

    normed3 <- runPCA(normed_alt, ncomponents=4, scale_features=FALSE, ntop=Inf)
    normed4 <- runPCA(normed, ncomponents=4, scale_features=TRUE, ntop=Inf)
    expect_equal(reducedDim(normed3), reducedDim(normed4))    

    # Testing out alternative assays.
    normed_alt <- normed
    assay(normed_alt, "whee") <- logcounts(normed)
    logcounts(normed_alt) <- NULL
    normed3 <- runPCA(normed_alt, ncomponents=4, exprs_values="whee")
    expect_equal(reducedDim(normed3), reducedDim(normedX))
})

test_that("runPCA works as expected for QC metrics", {
    normedQ <- calculateQCMetrics(normed, feature_controls=list(ERCC=1:50))
    expect_warning(normedQ <- runPCA(normedQ, use_coldata=TRUE), NA)
    expect_identical(reducedDimNames(normedQ), "PCA_coldata")

    normedQ2 <- calculateQCMetrics(normed, feature_controls=list(ERCC=1:50), compact=TRUE)
    expect_warning(normedQ2 <- runPCA(normedQ2, use_coldata=TRUE), NA)
    expect_identical(reducedDim(normedQ, "PCA_coldata"), reducedDim(normedQ2, "PCA_coldata"))

    # Checking outlier detection works correctly.
    expect_identical(normedQ2$outlier, NULL)
    expect_warning(normedQ2 <- runPCA(normedQ2, use_coldata=TRUE, detect_outliers=TRUE), NA)
    expect_type(normedQ2$outlier, "logical")
    expect_identical(length(normedQ2$outlier), ncol(normedQ2))

    # Check manual specification of variables works properly.
    expect_warning(normedQ3 <- runPCA(normedQ, use_coldata=TRUE, 
        selected_variables = c("total_counts", 
                               "total_features_by_counts",
                               "total_counts_endogenous",
                               "total_counts_feature_control")), NA)
})

test_that("runPCA works with irlba code", {
    set.seed(10)
    normedX <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    expect_identical(reducedDimNames(normedX), "PCA")     
    expect_identical(dim(reducedDim(normedX, "PCA")), c(ncol(normedX), 4L))
    expect_identical(length(attr(reducedDim(normedX), "percentVar")), 4L)
    expect_true(sum(attr(reducedDim(normedX), "percentVar")) < 1)

    # Checking that seed setting works.
    set.seed(100)
    normedX2 <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    set.seed(100)
    normedX3 <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    expect_false(isTRUE(all.equal(reducedDim(normedX), reducedDim(normedX2))))
    expect_equal(reducedDim(normedX2), reducedDim(normedX3))
})

#############################################
# Check t-SNE.

test_that("runTSNE works as expected", {
    set.seed(100)
    normedX <- runTSNE(normed, ncomponents = 3)
    expect_identical(reducedDimNames(normedX), "TSNE")     
    expect_identical(dim(reducedDim(normedX, "TSNE")), c(ncol(normedX), 3L))

    # Testing that setting the seed actually works.
    set.seed(100)
    normed2 <- runTSNE(normed)
    set.seed(100)
    normed3 <- runTSNE(normed)
    expect_equal(reducedDim(normed2), reducedDim(normed3))

    # Testing that various settings have some effect. 
    set.seed(100)
    normed3 <- runTSNE(normed, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
    
    set.seed(100)
    normed3 <- runTSNE(normed, perplexity = 5)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, pca = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, initial_dims = 10)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, normalize=FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, theta=0.1)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
})

test_that("runTSNE on existing reduced dimension results works as expected", {
    # Function should not respond to any feature settings.
    set.seed(10)
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runTSNE(normedP, use_dimred = "PCA")

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", scale_features = FALSE)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", feature_set = 1:20)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", pca = FALSE)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", initial_dims = 10)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, use_dimred = "PCA", n_dimred=3)
    expect_false(isTRUE(all.equal(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))))
})

test_that("runTSNE works with externally computed nearest neighbor results", {
    skip_on_os("windows") # https://github.com/jkrijthe/Rtsne/commit/f3f42504eeac627e4d886b1489ee289f8f9d082b#comments

    normedP <- runPCA(normed, ncomponents = 20)

    # Need to set the random seed to avoid different RNG states after the NN search. 
    set.seed(20) 
    init <- matrix(rnorm(ncol(normedP)*2), ncol=2)

    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init, perplexity=8.6)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, perplexity=8.6, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init, theta=0.1)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, theta=0.1, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init, normalize=FALSE)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, normalize=FALSE, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    # Works with alternative neighbor searching options.
    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, external_neighbors=TRUE, BNPARAM=BiocNeighbors::VptreeParam())
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, use_dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, use_dimred="PCA", Y_init=init, external_neighbors=TRUE, BPPARAM=MulticoreParam(2))
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))
})

#############################################
# Check UMAP.

test_that("runUMAP works as expected", {
    set.seed(100)
    normedX <- runUMAP(normed, ncomponents = 3)
    expect_identical(reducedDimNames(normedX), "UMAP")     
    expect_identical(dim(reducedDim(normedX, "UMAP")), c(ncol(normedX), 3L))

    # Testing that setting the seed actually works.
    set.seed(100)
    normed2 <- runUMAP(normed)
    set.seed(100)
    normed3 <- runUMAP(normed)
    expect_equal(reducedDim(normed2), reducedDim(normed3))

    # Testing that various settings have some effect. 
    set.seed(100)
    normed3 <- runUMAP(normed, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
})

test_that("runUMAP on existing reduced dimension results works as expected", {
    # Function should not respond to any feature settings.
    set.seed(10)
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runUMAP(normedP, use_dimred = "PCA")

    set.seed(10)
    normed3 <- runUMAP(normedP, use_dimred = "PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, use_dimred = "PCA", scale_features = FALSE)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, use_dimred = "PCA", feature_set = 1:20)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, use_dimred = "PCA", n_dimred=3)
    expect_false(isTRUE(all.equal(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))))
})

test_that("runUMAP works with externally computed nearest neighbor results", {
    normedP <- runPCA(normed, ncomponents = 20)

    # Need to cajolethe random seed to avoid different RNG states after the NN search. 
    seedSet <- function(...) invisible(BiocNeighbors::buildIndex(reducedDim(normedP, "PCA"), ...))

    set.seed(20) 
    seedSet()
    ref <- runUMAP(normedP, use_dimred="PCA")
    set.seed(20) 
    alt <- runUMAP(normedP, use_dimred="PCA", external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    set.seed(21) 
    seedSet()
    ref <- runUMAP(normedP, use_dimred="PCA", n_neighbors=10)
    set.seed(21) 
    alt <- runUMAP(normedP, use_dimred="PCA", n_neighbors=10, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    set.seed(22) 
    seedSet()
    ref <- runUMAP(normedP, use_dimred="PCA", bandwidth=1.5)
    set.seed(22) 
    alt <- runUMAP(normedP, use_dimred="PCA", bandwidth=1.5, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    # Works with alternative neighbor searching options.
    set.seed(23) 
    seedSet(BNPARAM=BiocNeighbors::VptreeParam())
    ref <- runUMAP(normedP, use_dimred="PCA")
    set.seed(23) 
    alt <- runUMAP(normedP, use_dimred="PCA", external_neighbors=TRUE, BNPARAM=BiocNeighbors::VptreeParam())
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    set.seed(24) 
    seedSet()
    invisible(MulticoreParam(2)) # more seed-related cajoling!
    ref <- runUMAP(normedP, use_dimred="PCA")
    set.seed(24) 
    alt <- runUMAP(normedP, use_dimred="PCA", external_neighbors=TRUE, BPPARAM=MulticoreParam(2))
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))
})

#############################################
# Check MDS.

test_that("runMDS works as expected", {
    normedX <- runMDS(normed, ncomponents = 3)
    expect_identical(reducedDimNames(normedX), "MDS")     
    expect_identical(dim(reducedDim(normedX, "MDS")), c(ncol(normedX), 3L))

    # Testing that various settings work.
    normed2 <- runMDS(normed)
    normed3 <- runMDS(normed, scale_features = FALSE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, feature_set = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, method = "manhattan")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings).
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runMDS(normedP, use_dimred = "PCA")

    normed3 <- runMDS(normedP, use_dimred = "PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    normed3 <- runMDS(normedP, use_dimred = "PCA", scale_features = FALSE)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    normed3 <- runMDS(normedP, use_dimred = "PCA", feature_set = 1:20)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    # This does, in fact, happen to be equal, due to the relationship between MDS and PCA.
    # This is not identifiable by the sign, hence the finagling.
    normed3 <- runMDS(normedP, use_dimred = "PCA", n_dimred=3)
    fold <- reducedDim(normed2, "MDS")/reducedDim(normed3, "MDS")
    expect_equal(abs(colSums(fold)), rep(nrow(fold), ncol(fold)))
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
    normedX <- runDiffusionMap(normed, ncomponents = 3)
    expect_identical(reducedDimNames(normedX), "DiffusionMap")     
    expect_identical(dim(reducedDim(normedX, "DiffusionMap")), c(ncol(normedX), 3L))

    # Testing that various settings work.
    normed2 <- runDiffusionMap(normed)
    normed3 <- runDiffusionMap(normed, scale_features = FALSE)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, ntop = 100)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, exprs_values = "counts")
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, feature_set = 1:100)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))
    
    normed3 <- runDiffusionMap(normed, sigma = 10)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, k = 13)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings)
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runDiffusionMap(normedP, use_dimred = "PCA")

    normed3 <- runDiffusionMap(normedP, use_dimred = "PCA", ntop = 20)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, use_dimred = "PCA", scale_features = FALSE)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, use_dimred = "PCA", feature_set = 1:20)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, use_dimred = "PCA", n_dimred=3)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"), same = FALSE)

})

#############################################
# Check defences against sparse matrices.

test_that("run* functions work with sparse matrices", {
    library(Matrix)          
    counts(normed) <- as(counts(normed), "dgCMatrix")
    logcounts(normed) <- as(logcounts(normed), "dgCMatrix")

    expect_error(runPCA(normed), NA)
    expect_error(runPCA(normed, BSPARAM=BiocSingular::IrlbaParam()), NA)
    expect_error(runTSNE(normed), NA)
    expect_error(runUMAP(normed), NA)
    expect_error(runDiffusionMap(normed), NA)
    expect_error(runMDS(normed), NA)
})

