## Testing functions for classes used ##

context("tests on inputs")

test_that("example datasets work", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)

    expect_that(example_sce, is_a("SingleCellExperiment"))
})


# test_that("we can update a SingleCellExperiment object", {
#     data("sc_example_counts")
#     data("sc_example_cell_info")
#     pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#     example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#     example_sceset
#     
#     expect_that(updateSCESet(example_sceset), is_a("SCESet"))
# })


context("test manipulations")

test_that("we can subset the example SingleCellExperiment", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(calculateCPM(example_sce) + 1)
    expect_warning(plotPCA(example_sce, return_SCE = TRUE),
                   "deprecated")
    example_sce <- runPCA(example_sce)
    ex_subset <- example_sce[1:200, sample(1:40, 25)]

    expect_equal(as.integer(nrow(ex_subset)), 200L)
    expect_equal(as.integer(ncol(ex_subset)), 25L)
    expect_equal(nrow(reducedDim(ex_subset)), 25L)
})




