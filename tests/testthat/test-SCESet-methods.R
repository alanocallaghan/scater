## Testing SCESet methods

test_that("newSCESet fails as expected", {
    sceset <- readRDS(system.file("extdata", "test_sceset.rds", package = "scater"))
    expect_that(toSingleCellExperiment(sceset), is_a("SingleCellExperiment"))
    expect_that(updateSCESet(sceset), is_a("SingleCellExperiment"))
    
})
