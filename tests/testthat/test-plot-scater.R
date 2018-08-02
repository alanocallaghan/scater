## Tests for plotting functions
## This stress-tests the plotting functions for all the different parameter settings. 
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
