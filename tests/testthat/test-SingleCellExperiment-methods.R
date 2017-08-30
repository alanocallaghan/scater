## Testing SingleCellExperiment methods

test_that("subsetting SingleCellExperiment works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    
    expect_that(example_sce[1:500,], is_a("SingleCellExperiment"))
    expect_that(example_sce[, 10:35], is_a("SingleCellExperiment"))
    expect_that(example_sce[500:1000, 7:27], is_a("SingleCellExperiment"))
})
