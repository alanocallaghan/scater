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

test_that("accessor functions for SingleCellExperiment work as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    assay(example_sce, "exprs") <- log2(calculateCPM(
        example_sce, use.size.factors = FALSE) + 1)
    expect_that(counts(example_sce), is_a("matrix"))
    expect_that(exprs(example_sce), is_null())
    expect_error(cpm(example_sce), "'i' not in names")
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    expect_that(exprs(example_sce), is_a("matrix"))  
    norm_exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    expect_that(norm_exprs(example_sce), is_a("matrix"))  
    stand_exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    expect_that(stand_exprs(example_sce), is_a("matrix"))  
    fpkm(example_sce) <- log2(calculateFPKM(
        example_sce, effective_length = 1000, use.size.factors = FALSE) + 1)
    expect_that(fpkm(example_sce), is_a("matrix"))  
    
})
