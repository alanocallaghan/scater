## Testing SingleCellExperiment methods

data("sc_example_counts")
data("sc_example_cell_info")
original <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)


test_that("subsetting SingleCellExperiment works as expected", {
    expect_that(original[1:500,], is_a("SingleCellExperiment"))
    expect_that(original[, 10:35], is_a("SingleCellExperiment"))
    expect_that(original[500:1000, 7:27], is_a("SingleCellExperiment"))
})

test_that("accessor functions for SingleCellExperiment work as expected", {
    example_sce <- original 
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

    # Same again for sparse matrices.    
    library(Matrix)
    sparsified <- original
    counts(sparsified) <- as(counts(original), "dgCMatrix")
    
    expect_that(exprs(sparsified), is_null())
    expect_that(counts(sparsified), is_a("dgCMatrix"))
    exprs(sparsified) <- log2(calculateCPM(
        sparsified, use_size_factors = FALSE) + 1)
    expect_that(exprs(sparsified), is_a("dgeMatrix"))
    expect_error(cpm(sparsified), "'cpm' not in names")
    norm_exprs(sparsified) <- log2(calculateCPM(sparsified, 
                                                 use_size_factors = FALSE) + 1)
    expect_that(norm_exprs(sparsified), is_a("dgeMatrix"))  
    stand_exprs(sparsified) <- log2(calculateCPM(sparsified, 
                                                  use_size_factors = FALSE) + 1)
    expect_that(stand_exprs(sparsified), is_a("dgeMatrix"))  
    fpkm(sparsified) <- log2(calculateFPKM(
        sparsified, effective_length = 1000, use_size_factors = FALSE) + 1)
    expect_that(fpkm(sparsified), is_a("dgeMatrix"))  
})
