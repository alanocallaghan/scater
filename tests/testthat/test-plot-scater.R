## This stress-tests the plotScater-related functions. 
## library(scater); library(testthat); source("setup.R"); source("test-plot-scater.R")

example_sce <- sce

test_that("plotScater works as expected", {
    expect_ggplot(plotScater(example_sce))
    expect_ggplot(plotScater(example_sce, colour_by = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, block1 = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, block2 = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, block1 = "Treatment", block2 = "Cell_Cycle"))

    # Different types of colouring are possible
    expect_ggplot(plotScater(example_sce, colour_by = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, colour_by = "Gene_0001"))

    expect_ggplot(plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, block1 = "Mutation_Status", colour_by = "Gene_0001"))

    expect_ggplot(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Cell_Cycle"))
    expect_ggplot(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Gene_0001"))
    
    expect_ggplot(plotScater(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"))

    # Responds to different type of expression values.
    cpm(example_sce) <- calculateCPM(example_sce)
    expect_ggplot(plotScater(example_sce, exprs_values="cpm"))
    expect_error(plotScater(example_sce, exprs_values="tpm"), "not in names")
})

test_that("plotScater's underlying C++ code works as expected", {
    REFFUN <- function(x, top) {
        prop <- cumsum(sort(x, decreasing=TRUE))/sum(x)
        prop[pmin(top, length(x))] * 100
    }

    out <- scater:::top_cumprop(assay(example_sce), 1:50)
    ref <- apply(assay(example_sce), 2, REFFUN, top=1:50)
    expect_equivalent(out, t(ref))

    out <- scater:::top_cumprop(assay(example_sce), 1:20*5)
    ref <- apply(assay(example_sce), 2, REFFUN, top=1:20*5)
    expect_equivalent(out, t(ref))

    # Handles sparse matrices.
    library(Matrix)
    spmat <- as(assay(example_sce), "dgCMatrix")
    out <- scater:::top_cumprop(spmat, 1:100)
    ref <- apply(spmat, 2, REFFUN, top=1:100)
    expect_equivalent(out, t(ref))

    # Behaves with silly inputs.
    out <- scater:::top_cumprop(assay(example_sce), integer(0))
    expect_identical(dim(out), c(ncol(example_sce), 0L))
    expect_identical(scater:::top_cumprop(assay(example_sce), 5:1), 
        scater:::top_cumprop(assay(example_sce), 1:5))
})
