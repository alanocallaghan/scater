## Tests for dimensionality reduction plotting functions
## library(scater); library(testthat); source("setup.R"); source("test-plot-dimred.R")

example_sce <- normed
rowData(example_sce)$ENS <- gsub("Gene", "ENS", rownames(example_sce))

test_that("we can produce PCA scatterplots", {
    example_sce <- runPCA(example_sce)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_ggplot(P <- plotPCA(example_sce))

    expect_ggplot(plotPCA(example_sce, colour_by = "Cell_Cycle"))
    expect_ggplot(plotPCA(example_sce, size_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, size_by = "Gene_0001", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, colour_by = "ENS_0001", swap_rownames = "ENS"))

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_ggplot(plotPCA(example_sce, colour_by = "Cell_Cycle", add_legend = FALSE))
    expect_ggplot(plotPCA(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"))
    expect_ggplot(plotPCA(example_sce, percentVar = c(19, 5)))
    expect_ggplot(plotPCA(example_sce, text_by="Cell_Cycle"))

    # Checking that specification of multiple ncomponents works.
    expect_ggplot(Pv <- plotPCA(example_sce, ncomponents=1:2))
    expect_equal(P$data, Pv$data)
    expect_ggplot(Pv2 <- plotPCA(example_sce, ncomponents=2:1))
    expect_false(isTRUE(all.equal(P$data, Pv2$data)))
    expect_error(plotPCA(example_sce, ncomponents=c(51,1)), "larger than")

    # Check that dataframes etc are allowed
    reducedDim(example_sce, "PCA") <- DataFrame(reducedDim(example_sce, "PCA"))
    expect_error(plotReducedDim(example_sce, "PCA"), NA)

    # Use scattermore
    expect_ggplot(plotPCA(example_sce, scattermore = TRUE, point_size = 2))
    expect_ggplot(plotPCA(example_sce, scattermore = TRUE, point_size = 2,
                          colour_by = "Cell_Cycle"))

    # Binning
    expect_ggplot(plotPCA(example_sce, bins = 10))
    expect_ggplot(plotPCA(example_sce, bins = 10, hex = TRUE))
    expect_ggplot(plotPCA(example_sce, bins = 10, colour_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, bins = 10, colour_by = "Gene_0001",
                          hex = TRUE))
})

test_that("we can produce PCA pairplots", {
    example_sce <- runPCA(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_ggplot(P <- plotPCA(example_sce, ncomponents = 4))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Cell_Cycle"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, size_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Cell_Cycle", size_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Cell_Cycle", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, size_by = "sizeFactor", shape_by = "Treatment"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "ENS_0001", swap_rownames = "ENS"))


    # Checking other arguments are passed successfully to plotReducedDim.
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Cell_Cycle", add_legend = FALSE))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, colour_by = "Gene_0001", by_exprs_values = "counts"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, percentVar = c(19, 5, 3, 2)))

    # Checking that specification of multiple ncomponents works.
    expect_ggplot(Pv <- plotPCA(example_sce, ncomponents = 1:4))
    expect_equal(P$data, Pv$data)
    expect_ggplot(Pv2 <- plotPCA(example_sce, ncomponents = 4:1))
    expect_false(isTRUE(all.equal(P$data, Pv2$data)))
    expect_error(plotPCA(example_sce, ncomponents = 5:1), "larger than")

    # Use scattermore
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, scattermore = TRUE, point_size = 3))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, scattermore = TRUE, point_size = 3, colour_by = "Gene_0001"))

    # Binning
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, bins = 10))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, bins = 10, hex = TRUE))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, bins = 10,
                          colour_by = "Gene_0001"))
    expect_ggplot(plotPCA(example_sce, ncomponents = 4, bins = 10,
                          colour_by = "Gene_0001", hex = TRUE))
})

test_that("we can produce TSNE plots", {
    set.seed(100)
    example_sce <- runTSNE(example_sce)
    expect_identical(reducedDimNames(example_sce), "TSNE")
    expect_ggplot(P <- plotTSNE(example_sce))

    set.seed(20)
    example_sce <- runTSNE(example_sce, ncomponents=3)
    expect_ggplot(plotTSNE(example_sce, ncomponents=3))
})

test_that("we can produce UMAP plots", {
    set.seed(100)
    example_sce <- runUMAP(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "UMAP")
    expect_ggplot(P <- plotUMAP(example_sce))

    # Handles multiple components properly.
    set.seed(20)
    P4 <- plotUMAP(example_sce, ncomponents=4)
    expect_ggplot(P4)
})

test_that("we can produce NMF plots", {
    set.seed(100)
    example_sce <- runNMF(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "NMF")
    expect_ggplot(P <- plotNMF(example_sce))

    # Handles multiple components properly.
    set.seed(20)
    P4 <- plotNMF(example_sce, ncomponents=4)
    expect_ggplot(P4)
})

test_that("we can produce MDS plots", {
    example_sce <- runMDS(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "MDS")
    expect_ggplot(P <- plotMDS(example_sce))

    # Handles multiple components properly.
    expect_ggplot(P4 <- plotMDS(example_sce, ncomponents=4))
})


test_that("order by works", {
    set.seed(42)
    example_sce <- mockSCE()
    example_sce <- logNormCounts(example_sce)
    example_sce <- runPCA(example_sce)
    p <- plotReducedDim(example_sce, "PCA", order_by = "Gene_0001")
    expect_equal(order(p$data$order_by), seq_len(nrow(p$data)))
    p <- plotReducedDim(example_sce, "PCA", order_by = "Mutation_Status")
    expect_equal(order(p$data$order_by), seq_len(nrow(p$data)))
})
