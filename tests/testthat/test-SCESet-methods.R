## Testing SCESet methods

library(Biobase)
test_that("newSCESet fails as expected", {
    sceset <- readRDS(system.file("extdata", "test_sceset.rds", package = "scater"))
    featureData(sceset)$strand <- NULL # clashes with protected field in a RSE object.
    expect_that(toSingleCellExperiment(sceset), is_a("SingleCellExperiment"))
    expect_that(updateSCESet(sceset), is_a("SingleCellExperiment"))
    
})
