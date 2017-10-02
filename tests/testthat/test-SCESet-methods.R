## Testing SCESet methods

test_that("newSCESet fails as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    expect_warning((example_sceset <- 
                       newSCESet(countData = sc_example_counts, phenoData = pd)),
                   "deprecated")
    expect_that(example_sceset, is_a("SingleCellExperiment"))
    
    sceset <- readRDS(system.file("extdata", "test_sceset.rds", package = "scater"))
    expect_that(toSingleCellExperiment(sceset), is_a("SingleCellExperiment"))
    expect_that(updateSCESet(sceset), is_a("SingleCellExperiment"))
    
})
