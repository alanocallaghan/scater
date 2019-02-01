## Tests for dimensionality reduction plotting functions
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-dimred.R")

example_sce <- normed 

test_that("we can produce PCA scatterplots", {
    example_sce <- runPCA(example_sce)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_s3_class(P <- plotPCA(example_sce), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotPCA(example_sce, size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", add_legend = FALSE), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Treatment", by_show_single = TRUE), "ggplot")
    expect_s3_class(plotPCA(example_sce, percentVar = c(19, 5)), "ggplot")
    expect_s3_class(plotPCA(example_sce, text_by="Cell_Cycle"), "ggplot")
    
    # Checking that re-running works, responsive to feature scaling.
    example_sceX <- runPCA(example_sce, scale_features=FALSE)
    expect_s3_class(Px <- plotPCA(example_sceX), "ggplot")
    expect_s3_class(Px2 <- plotPCA(example_sce, rerun=TRUE, run_args=list(scale_features=FALSE)), "ggplot")
    expect_equal(Px, Px2)
    expect_false(isTRUE(all.equal(P, Px2)))

    reducedDim(example_sceX, "PCA") <- NULL
    expect_s3_class(Px3 <- plotPCA(example_sceX), "ggplot")
    expect_equal(P, Px3)
    expect_s3_class(Px4 <- plotPCA(example_sceX, run_args=list(scale_features=FALSE)), "ggplot")
    expect_equal(Px, Px4)
    
    # Checking that specification of multiple ncomponents works.
    expect_s3_class(Pv <- plotPCA(example_sce, ncomponents=1:2), "ggplot")
    expect_equal(P, Pv)
    expect_s3_class(Pv2 <- plotPCA(example_sce, ncomponents=2:1), "ggplot")
    expect_false(isTRUE(all.equal(P, Pv2)))
    expect_error(plotPCA(example_sce, ncomponents=c(3,1)), "larger than")

    example_sceY <- example_sce
    reducedDim(example_sceY, "PCA") <- NULL
    expect_s3_class(Py <- plotPCA(example_sceY, ncomponents=1:2), "ggplot")
    expect_equal(Pv, Py)
    expect_s3_class(Py2 <- plotPCA(example_sceY, ncomponents=2:1), "ggplot")
    expect_equal(Pv2, Py2)
    expect_s3_class(Py3 <- plotPCA(example_sceY, ncomponents=c(3, 1)), "ggplot")
    expect_false(isTRUE(all.equal(Py2, Py3)))
})

test_that("we can produce PCA pairplots", {
    example_sce <- runPCA(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_s3_class(P <- plotPCA(example_sce, ncomponents=4), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", add_legend = FALSE), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Treatment", by_show_single = TRUE), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, percentVar=c(19, 5, 3, 2)), "ggplot")
    
    # Checking that re-running works, responsive to feature scaling.
    expect_s3_class(P2 <- plotPCA(example_sce, ncomponents=4, rerun=TRUE, run_args=list(scale_features=FALSE)), "ggplot")
    expect_false(isTRUE(all.equal(P, P2)))

    example_sceX <- example_sce
    reducedDim(example_sceX, "PCA") <- NULL
    expect_s3_class(Px <- plotPCA(example_sceX, ncomponents=4), "ggplot")
    expect_equal(P, Px)

    expect_s3_class(Px2 <- plotPCA(example_sceX, ncomponents=4, run_args=list(scale_features=FALSE)), "ggplot")
    expect_equal(P2, Px2)
    expect_false(isTRUE(all.equal(Px, Px2)))

    # Checking that specification of multiple ncomponents works.
    expect_s3_class(Pv <- plotPCA(example_sce, ncomponents=1:4), "ggplot")
    expect_equal(P, Pv)
    expect_s3_class(Pv2 <- plotPCA(example_sce, ncomponents=4:1), "ggplot")
    expect_false(isTRUE(all.equal(P, Pv2)))
    expect_error(plotPCA(example_sce, ncomponents=5:1), "larger than")

    example_sceY <- example_sce
    reducedDim(example_sceY, "PCA") <- NULL
    expect_s3_class(Py <- plotPCA(example_sceY, ncomponents=1:4), "ggplot")
    expect_equal(Pv, Py)
    expect_s3_class(Py2 <- plotPCA(example_sceY, ncomponents=4:1), "ggplot")
    expect_equal(Pv2, Py2)
    expect_s3_class(Py3 <- plotPCA(example_sceY, ncomponents=5:1), "ggplot")
    expect_false(isTRUE(all.equal(Py2, Py3)))
})

test_that("we can produce TSNE plots", {
    set.seed(100)
    example_sce <- runTSNE(example_sce)
    expect_identical(reducedDimNames(example_sce), "TSNE")
    expect_s3_class(P <- plotTSNE(example_sce), "ggplot")

    # plotTSNE re-runs it correctly.
    reducedDim(example_sce, "TSNE") <- NULL

    set.seed(100)
    P2 <- plotTSNE(example_sce)
    expect_s3_class(P2, "ggplot")
    expect_equal(P, P2)

    # Responsive to changes in parameters.
    set.seed(100)
    P3 <- plotTSNE(example_sce, run_args=list(perplexity=10))
    expect_s3_class(P3, "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))
    
    # Handles multiple components properly.
    set.seed(20)
    P4 <- plotTSNE(example_sce, ncomponents=3)
    expect_s3_class(P4, "ggplot")

    set.seed(20)
    example_sce <- runTSNE(example_sce, ncomponents=3)
    expect_equal(plotTSNE(example_sce, ncomponents=3), P4)
})

test_that("we can produce UMAP plots", {
    set.seed(100)
    example_sce <- runUMAP(example_sce)
    expect_identical(reducedDimNames(example_sce), "UMAP")
    expect_s3_class(P <- plotUMAP(example_sce), "ggplot")

    # plotUMAP re-runs it correctly.
    reducedDim(example_sce, "UMAP") <- NULL

    set.seed(100)
    P2 <- plotUMAP(example_sce)
    expect_s3_class(P2, "ggplot")
    expect_equal(P, P2)

    # Responsive to changes in parameters.
    set.seed(100)
    P3 <- plotUMAP(example_sce, run_args=list(n_neighbors=10))
    expect_s3_class(P3, "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))
    
    # Handles multiple components properly.
    set.seed(20)
    P4 <- plotUMAP(example_sce, ncomponents=4)
    expect_s3_class(P4, "ggplot")

    set.seed(20)
    example_sce <- runUMAP(example_sce, ncomponents=4)
    expect_equal(plotUMAP(example_sce, ncomponents=4), P4)
})

test_that("we can produce diffusion maps", {
    set.seed(100)        
    example_sce <- runDiffusionMap(example_sce)
    expect_identical(reducedDimNames(example_sce), "DiffusionMap")
    expect_s3_class(P <- plotDiffusionMap(example_sce), "ggplot")

    # plotDiffusionMap re-runs it correctly.
    reducedDim(example_sce, "DiffusionMap") <- NULL

    set.seed(100)
    P2 <- plotDiffusionMap(example_sce)
    expect_s3_class(P2, "ggplot")
#    expect_equal(P, P2) # it seems as if destiny::DiffusionMap does not respond to the seed!

    # Responsive to changes in parameters.
    set.seed(100)
    P3 <- plotDiffusionMap(example_sce, run_args=list(k=13))
    expect_s3_class(P3, "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))

    # Handles multiple components properly.
    set.seed(20)        
    expect_s3_class(P4 <- plotDiffusionMap(example_sce, ncomponents=4), "ggplot")
#    set.seed(20)
#    example_sce <- runDiffusionMap(example_sce, ncomponents=4)
#    expect_equal(plotDiffusionMap(example_sce, ncomponents=4), P4)
})

test_that("we can produce MDS plots", {
    example_sce <- runMDS(example_sce)
    expect_identical(reducedDimNames(example_sce), "MDS")
    expect_s3_class(P <- plotMDS(example_sce), "ggplot")

    # plotMDS re-runs it correctly.
    reducedDim(example_sce, "MDS") <- NULL
    expect_s3_class(P2 <- plotMDS(example_sce), "ggplot")
    expect_equal(P, P2)

    # Responsive to changes in parameters.
    expect_s3_class(P3 <- plotMDS(example_sce, run_args=list(method="manhattan")), "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))

    # Handles multiple components properly.
    expect_s3_class(P4 <- plotMDS(example_sce, ncomponents=4), "ggplot")
    example_sce <- runMDS(example_sce, ncomponents=4)
    expect_equal(plotMDS(example_sce, ncomponents=4), P4)
})
