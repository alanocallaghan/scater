## Testing SCESet methods

test_that("newSCESet fails as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    expect_warning((example_sceset <- 
                       newSCESet(countData = sc_example_counts, phenoData = pd)),
                   "deprecated")
    expect_that(example_sceset, is_a("SingleCellExperiment"))

})
