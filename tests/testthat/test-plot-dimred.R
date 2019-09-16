## Tests for dimensionality reduction plotting functions
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-dimred.R")

example_sce <- normed 

test_that("we can produce PCA scatterplots", {
    example_sce <- runPCA(example_sce, BSPARAM=BiocSingular::ExactParam())
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
    
    # Checking that specification of multiple ncomponents works.
    expect_s3_class(Pv <- plotPCA(example_sce, ncomponents=1:2), "ggplot")
    expect_equal(P, Pv)
    expect_s3_class(Pv2 <- plotPCA(example_sce, ncomponents=2:1), "ggplot")
    expect_false(isTRUE(all.equal(P, Pv2)))
    expect_error(plotPCA(example_sce, ncomponents=c(51,1)), "larger than")
})

test_that("we can produce PCA pairplots", {
    example_sce <- runPCA(example_sce, ncomponents=4, BSPARAM=BiocSingular::ExactParam())
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
    
    # Checking that specification of multiple ncomponents works.
    expect_s3_class(Pv <- plotPCA(example_sce, ncomponents=1:4), "ggplot")
    expect_equal(P, Pv)
    expect_s3_class(Pv2 <- plotPCA(example_sce, ncomponents=4:1), "ggplot")
    expect_false(isTRUE(all.equal(P, Pv2)))
    expect_error(plotPCA(example_sce, ncomponents=5:1), "larger than")
})

test_that("we can produce TSNE plots", {
    set.seed(100)
    example_sce <- runTSNE(example_sce)
    expect_identical(reducedDimNames(example_sce), "TSNE")
    expect_s3_class(P <- plotTSNE(example_sce), "ggplot")

    set.seed(20)
    example_sce <- runTSNE(example_sce, ncomponents=3)
    expect_s3_class(plotTSNE(example_sce, ncomponents=3), "ggplot")
})

test_that("we can produce UMAP plots", {
    set.seed(100)
    example_sce <- runUMAP(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "UMAP")
    expect_s3_class(P <- plotUMAP(example_sce), "ggplot")

    # Handles multiple components properly.
    set.seed(20)
    P4 <- plotUMAP(example_sce, ncomponents=4)
    expect_s3_class(P4, "ggplot")
})

test_that("we can produce diffusion maps", {
    set.seed(100)        
    example_sce <- runDiffusionMap(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "DiffusionMap")
    expect_s3_class(P <- plotDiffusionMap(example_sce), "ggplot")

    # Handles multiple components properly.
    set.seed(20)        
    expect_s3_class(P4 <- plotDiffusionMap(example_sce, ncomponents=4), "ggplot")
})

test_that("we can produce MDS plots", {
    example_sce <- runMDS(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "MDS")
    expect_s3_class(P <- plotMDS(example_sce), "ggplot")

    # Handles multiple components properly.
    expect_s3_class(P4 <- plotMDS(example_sce, ncomponents=4), "ggplot")
})
