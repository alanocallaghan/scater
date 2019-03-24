## This stress-tests the plotScater-related functions. 
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-scater.R")

example_sce <- sce

test_that("plotScater works as expected", {
    expect_s3_class(plotScater(example_sce), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block2 = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Treatment", block2 = "Cell_Cycle"), "ggplot")

    # Different types of colouring are possible
    expect_s3_class(plotScater(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Gene_0001"), "ggplot")

    expect_s3_class(plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Mutation_Status", colour_by = "Gene_0001"), "ggplot")

    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Gene_0001"), "ggplot")
    
    expect_s3_class(plotScater(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Treatment", by_show_single = TRUE), "ggplot")

    # Responds to different type of expression values.
    cpm(example_sce) <- calculateCPM(example_sce)
    expect_s3_class(plotScater(example_sce, exprs_values="cpm"), "ggplot")
    expect_error(plotScater(example_sce, exprs_values="tpm"), "not in names")
})

test_that("plotScater's underlying C++ code works as expected", {
    REFFUN <- function(x, top) {
        prop <- cumsum(sort(x, decreasing=TRUE))/sum(x)
        prop[pmin(top, length(x))]
    }

    out <- .Call(scater:::cxx_top_cumprop, assay(example_sce), 1:50)
    ref <- apply(assay(example_sce), 2, REFFUN, top=1:50)
    colnames(out) <- colnames(ref)
    expect_equal(out, ref)

    out <- .Call(scater:::cxx_top_cumprop, assay(example_sce), 1:20*5)
    ref <- apply(assay(example_sce), 2, REFFUN, top=1:20*5)
    colnames(out) <- colnames(ref)
    expect_equal(out, ref)

    # Handles sparse matrices.
    library(Matrix)
    spmat <- as(assay(example_sce), "dgCMatrix")
    out <- .Call(scater:::cxx_top_cumprop, spmat, 1:100)
    ref <- apply(spmat, 2, REFFUN, top=1:100)
    colnames(out) <- colnames(ref)
    expect_equal(out, ref)

    # Behaves with silly inputs.
    out <- .Call(scater:::cxx_top_cumprop, assay(example_sce), integer(0))
    expect_identical(dim(out), c(0L, ncol(example_sce)))
    expect_error(.Call(scater:::cxx_top_cumprop, assay(example_sce), 5:1), "sorted")
})
