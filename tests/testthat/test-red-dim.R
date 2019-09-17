# Tests the options in the various reduced dimension functions
# library(scater); library(testthat); source("setup-sce.R"); source("test-red-dim.R")

#############################################
# Check the feature selection and scaling work.

test_that("feature selection is operational", {
    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row = seq_len(nrow(normed)), scale=FALSE)
    expect_equal(out, t(logcounts(normed)))

    # Ntop selection works.
    out <- scater:::.get_mat_for_reddim(logcounts(normed), ntop = 10, subset_row = NULL, scale=FALSE)
    rv <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    keep <- head(order(rv, decreasing=TRUE), 10)
    expect_equal(out, t(logcounts(normed)[keep,]))

    out <- scater:::.get_mat_for_reddim(logcounts(normed), ntop = Inf, subset_row = NULL, scale=FALSE)
    o <- order(rv, decreasing=TRUE)
    expect_equal(out, t(logcounts(normed)[o,]))

    # Feature selection works.
    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row = 10:1, scale=FALSE)
    expect_equal(out, t(logcounts(normed)[10:1,]))

    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row = rownames(normed)[10:1], scale=FALSE)
    expect_equal(out, t(logcounts(normed)[10:1,]))
})

test_that("scaling by feature variances work correctly", {
    logcounts(normed)[100,] <- 0
    MAT <- t(logcounts(normed))
    cv <- DelayedMatrixStats::colVars(DelayedArray(MAT))
    novar <- cv < 1e-8
    expect_true(any(novar))

    NOMAT <- MAT[,!novar]
    scaled <- t(t(NOMAT)/sqrt(cv[!novar]))

    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row = seq_len(nrow(normed)), scale = TRUE)
    expect_equal(out, scaled)
    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row = 10:1, scale = TRUE)
    expect_equal(out, scaled[,10:1])

    out <- scater:::.get_mat_for_reddim(logcounts(normed), subset_row=NULL, ntop = 10, scale = TRUE) # In combination with non-trivial selection.
    rv <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    expect_equal(out, scaled[,head(order(rv[!novar], decreasing=TRUE), 10)])
})

#############################################
# Check PCA.

test_that("runPCA works as expected", {
    # Using ExactParam() simply to avoid issues due to randomness.
    normedX <- runPCA(normed) 
    expect_identical(reducedDimNames(normedX), "PCA")     
    fullN <- min(dim(normed), 50L)
    expect_identical(dim(reducedDim(normedX, "PCA")), c(ncol(normed), fullN))

    normedX <- runPCA(normed, ncol(normed)) 
    expect_equal(sum(attr(reducedDim(normedX), "percentVar")), 100)

    # Checking that it works with a restricted number of components.
    normedX <- runPCA(normed, ncomponents=4)
    expect_identical(reducedDimNames(normedX), "PCA")     
    expect_identical(dim(reducedDim(normedX, "PCA")), c(ncol(normedX), 4L))
    expect_true(sum(attr(reducedDim(normedX), "percentVar")) < 100)
})

test_that("runPCA responds to changes to various settings", {
    # Testing that various settings give different results.
    normed2 <- runPCA(normed)
    normed3 <- runPCA(normed, scale = TRUE)

    fullN <- min(dim(normed), 50L)
    expect_identical(ncol(reducedDim(normed2)), fullN)
    expect_identical(ncol(reducedDim(normed3)), fullN)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runPCA(normed, ntop = 100)
    expect_identical(ncol(reducedDim(normed3)), fullN)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runPCA(normed, exprs_values = "counts")
    expect_identical(ncol(reducedDim(normed3)), fullN)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runPCA(normed, subset_row = 1:100)
    expect_identical(ncol(reducedDim(normed3)), fullN)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
})

test_that("runPCA handles ntop selection", {
    most_var <- DelayedMatrixStats::rowVars(DelayedArray(logcounts(normed)))
    keep <- head(order(most_var, decreasing=TRUE), 100)
    normed3 <- runPCA(normed, ncomponents=4, subset_row=keep)
    normed4 <- runPCA(normed, ncomponents=4, ntop=100)
    expect_equal(reducedDim(normed3), reducedDim(normed4))    
})

test_that("runPCA handles scaling", {
    # Setting ntop=Inf, otherwise it will pick 'normed_alt' features based on scaled variance.
    normed_alt <- normed
    rescaled <- t(scale(t(logcounts(normed_alt)), scale = TRUE))
    rescaled[is.na(rescaled)] <- 0
    logcounts(normed_alt) <- rescaled

    normed3 <- runPCA(normed_alt, ncomponents=4, scale=FALSE, ntop=Inf)
    normed4 <- runPCA(normed, ncomponents=4, scale=TRUE, ntop=Inf)
    expect_equal(reducedDim(normed3), reducedDim(normed4))

    # Checking that percentVar is computed correctly.
    normed3 <- runPCA(normed, scale=TRUE, ncol(normed)) 
    expect_equal(sum(attr(reducedDim(normed3), "percentVar")), 100)
})

test_that("runPCA behaves with alternative assays", {
    normed_alt <- normed
    assay(normed_alt, "whee") <- logcounts(normed)
    logcounts(normed_alt) <- NULL

    normed3 <- runPCA(normed_alt, ncomponents=4, exprs_values="whee")
    normed4 <- runPCA(normed, ncomponents=4)
    expect_identical(reducedDim(normed3), reducedDim(normed4))
})

test_that("runColDataPCA works as expected for QC metrics", {
    normed$field1 <- runif(ncol(normed))
    normed$field2 <- runif(ncol(normed))
    normed$field3 <- runif(ncol(normed))
    normed$field4 <- runif(ncol(normed))

    vars <- c("field1", "field2", "field3", "field4")
    normed <- runColDataPCA(normed, variables=vars)
    out <- reducedDim(normed, "PCA_coldata")
    ref <- calculatePCA(t(as.matrix(colData(normed)[,vars])), ncomponents=2, scale=TRUE)
    expect_true(all(abs(abs(colMeans(out/ref)) - 1) < 1e-6))
    
    # Checking outlier detection works correctly.
    expect_identical(normed$outlier, NULL)
    expect_warning(normed <- runColDataPCA(normed, variables=vars, outliers=TRUE), NA)
    expect_type(normed$outlier, "logical")
    expect_identical(length(normed$outlier), ncol(normed))
})

test_that("runPCA works with irlba code", {
    set.seed(10)
    normedX <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    expect_identical(reducedDimNames(normedX), "PCA")     
    expect_identical(dim(reducedDim(normedX, "PCA")), c(ncol(normedX), 4L))
    expect_identical(length(attr(reducedDim(normedX), "percentVar")), 4L)
    expect_true(sum(attr(reducedDim(normedX), "percentVar")) < 100)

    # Checking that seed setting works.
    set.seed(100)
    normedX2 <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    set.seed(100)
    normedX3 <- runPCA(normed, ncomponents=4, BSPARAM=BiocSingular::IrlbaParam())
    expect_false(isTRUE(all.equal(reducedDim(normedX), reducedDim(normedX2))))
    expect_identical(reducedDim(normedX2), reducedDim(normedX3))
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
    normed3 <- runTSNE(normed, scale=TRUE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runTSNE(normed, subset_row = 1:100)
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
    normed2 <- runTSNE(normedP, dimred="PCA")

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", scale=TRUE)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", subset_row = 1:20)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", pca = FALSE)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", initial_dims = 10)
    expect_identical(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))

    set.seed(10)
    normed3 <- runTSNE(normedP, dimred="PCA", n_dimred=3)
    expect_false(isTRUE(all.equal(reducedDim(normed2, "TSNE"), reducedDim(normed3, "TSNE"))))
})

test_that("runTSNE works with externally computed nearest neighbor results", {
    skip_on_os("windows") # https://github.com/jkrijthe/Rtsne/commit/f3f42504eeac627e4d886b1489ee289f8f9d082b#comments

    normedP <- runPCA(normed, ncomponents = 20)

    # Need to set the random seed to avoid different RNG states after the NN search. 
    set.seed(20) 
    init <- matrix(rnorm(ncol(normedP)*2), ncol=2)

    ref <- runTSNE(normedP, dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, dimred="PCA", Y_init=init, perplexity=8.6)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, perplexity=8.6, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, dimred="PCA", Y_init=init, theta=0.1)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, theta=0.1, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, dimred="PCA", Y_init=init, normalize=FALSE)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, normalize=FALSE, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    # Works with alternative neighbor searching options.
    ref <- runTSNE(normedP, dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, external_neighbors=TRUE, BNPARAM=BiocNeighbors::VptreeParam())
    expect_identical(reducedDim(ref, "TSNE"), reducedDim(alt, "TSNE"))

    ref <- runTSNE(normedP, dimred="PCA", Y_init=init)
    alt <- runTSNE(normedP, dimred="PCA", Y_init=init, external_neighbors=TRUE, BPPARAM=safeBPParam(2))
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
    normed3 <- runUMAP(normed, scale=TRUE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    set.seed(100)
    normed3 <- runUMAP(normed, subset_row = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))
})

test_that("runUMAP on existing reduced dimension results works as expected", {
    # Function should not respond to any feature settings.
    set.seed(10)
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runUMAP(normedP, dimred="PCA")

    set.seed(10)
    normed3 <- runUMAP(normedP, dimred="PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, dimred="PCA", scale=TRUE)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, dimred="PCA", subset_row = 1:20)
    expect_identical(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))

    set.seed(10)
    normed3 <- runUMAP(normedP, dimred="PCA", n_dimred=3)
    expect_false(isTRUE(all.equal(reducedDim(normed2, "UMAP"), reducedDim(normed3, "UMAP"))))
})

test_that("runUMAP works with externally computed nearest neighbor results", {
    skip_on_os("windows") # Use with VP-tree gives different results from internal NN search on Win32. Why? Who knows. 

    normedP <- runPCA(normed, ncomponents = 20)

    # Need to cajole the random seed to avoid different RNG states after the NN search. 
    seedSet <- function(...) invisible(BiocNeighbors::buildIndex(reducedDim(normedP, "PCA"), ...))

    set.seed(20) 
    seedSet()
    ref <- runUMAP(normedP, dimred="PCA")
    set.seed(20) 
    alt <- runUMAP(normedP, dimred="PCA", external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    set.seed(21) 
    seedSet()
    ref <- runUMAP(normedP, dimred="PCA", n_neighbors=10)
    set.seed(21) 
    alt <- runUMAP(normedP, dimred="PCA", n_neighbors=10, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    set.seed(22) 
    seedSet()
    ref <- runUMAP(normedP, dimred="PCA", bandwidth=1.5)
    set.seed(22) 
    alt <- runUMAP(normedP, dimred="PCA", bandwidth=1.5, external_neighbors=TRUE)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    # Works with alternative neighbor searching options.
    set.seed(23) 
    seedSet(BNPARAM=BiocNeighbors::VptreeParam())
    ref <- runUMAP(normedP, dimred="PCA")
    set.seed(23) 
    alt <- runUMAP(normedP, dimred="PCA", external_neighbors=TRUE, BNPARAM=BiocNeighbors::VptreeParam())
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))

    # Works with parallelization and more seed-related cajoling.
    BPPARAM <- safeBPParam(2)
    BiocParallel::bpstart(BPPARAM)

    set.seed(24) 
    seedSet()
    ref <- runUMAP(normedP, dimred="PCA")
    set.seed(24) 
    alt <- runUMAP(normedP, dimred="PCA", external_neighbors=TRUE, BPPARAM=BPPARAM)
    expect_identical(reducedDim(ref, "UMAP"), reducedDim(alt, "UMAP"))
    
    BiocParallel::bpstop(BPPARAM)
})

#############################################
# Check MDS.

test_that("runMDS works as expected", {
    normedX <- runMDS(normed, ncomponents = 3)
    expect_identical(reducedDimNames(normedX), "MDS")     
    expect_identical(dim(reducedDim(normedX, "MDS")), c(ncol(normedX), 3L))

    # Testing that various settings work.
    normed2 <- runMDS(normed)
    normed3 <- runMDS(normed, scale=TRUE)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, ntop = 100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, exprs_values = "counts")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, subset_row = 1:100)
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    normed3 <- runMDS(normed, method = "manhattan")
    expect_false(isTRUE(all.equal(reducedDim(normed2), reducedDim(normed3))))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings).
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runMDS(normedP, dimred="PCA")

    normed3 <- runMDS(normedP, dimred="PCA", ntop = 20)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    normed3 <- runMDS(normedP, dimred="PCA", scale=TRUE)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    normed3 <- runMDS(normedP, dimred="PCA", subset_row = 1:20)
    expect_identical(reducedDim(normed2, "MDS"), reducedDim(normed3, "MDS"))

    # This does, in fact, happen to be equal, due to the relationship between MDS and PCA.
    # This is not identifiable by the sign, hence the finagling.
    normed3 <- runMDS(normedP, dimred="PCA", n_dimred=3)
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
    normed3 <- runDiffusionMap(normed, scale=TRUE)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, ntop = 100)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, exprs_values = "counts")
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    normed3 <- runDiffusionMap(normed, subset_row = 1:100)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))
    
    normed3 <- runDiffusionMap(normed, k = 13)
    SIGNAGNOSTIC(same=FALSE, reducedDim(normed2), reducedDim(normed3))

    # Testing out the use of existing reduced dimensions (this should not respond to any feature settings)
    normedP <- runPCA(normed, ncomponents = 4)
    normed2 <- runDiffusionMap(normedP, dimred="PCA")

    normed3 <- runDiffusionMap(normedP, dimred="PCA", ntop = 20)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, dimred="PCA", scale=TRUE)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, dimred="PCA", subset_row = 1:20)
    SIGNAGNOSTIC(reducedDim(normed2, "DiffusionMap"), reducedDim(normed3, "DiffusionMap"))

    normed3 <- runDiffusionMap(normedP, dimred="PCA", n_dimred=3)
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
