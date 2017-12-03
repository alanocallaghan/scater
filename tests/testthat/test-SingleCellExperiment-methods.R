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
        example_sce, use_size_factors = FALSE) + 1)
    expect_that(counts(example_sce), is_a("matrix"))
    expect_that(exprs(example_sce), is_null())
    expect_error(cpm(example_sce), "'cpm' not in names")
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use_size_factors = FALSE) + 1)
    expect_that(exprs(example_sce), is_a("matrix"))  
    norm_exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use_size_factors = FALSE) + 1)
    expect_that(norm_exprs(example_sce), is_a("matrix"))  
    stand_exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use_size_factors = FALSE) + 1)
    expect_that(stand_exprs(example_sce), is_a("matrix"))  
    fpkm(example_sce) <- log2(calculateFPKM(
        example_sce, effective_length = 1000, use_size_factors = FALSE) + 1)
    expect_that(fpkm(example_sce), is_a("matrix"))  
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    expect_that(exprs(sce10x), is_null())
    expect_that(counts(sce10x), is_a("dgCMatrix"))
    exprs(sce10x) <- log2(calculateCPM(
        sce10x, use_size_factors = FALSE) + 1)
    expect_that(exprs(sce10x), is_a("dgeMatrix"))
    expect_error(cpm(sce10x), "'cpm' not in names")
    norm_exprs(sce10x) <- log2(calculateCPM(sce10x, 
                                                 use_size_factors = FALSE) + 1)
    expect_that(norm_exprs(sce10x), is_a("dgeMatrix"))  
    stand_exprs(sce10x) <- log2(calculateCPM(sce10x, 
                                                  use_size_factors = FALSE) + 1)
    expect_that(stand_exprs(sce10x), is_a("dgeMatrix"))  
    fpkm(sce10x) <- log2(calculateFPKM(
        sce10x, effective_length = 1000, use_size_factors = FALSE) + 1)
    expect_that(fpkm(sce10x), is_a("dgeMatrix"))  
})
