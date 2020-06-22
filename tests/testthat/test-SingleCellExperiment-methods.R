## Testing SingleCellExperiment methods

original <- sce

test_that("accessor functions for SingleCellExperiment work as expected", {
    example_sce <- original 
    assay(example_sce, "exprs") <- log2(calculateCPM(example_sce) + 1)

    expect_that(counts(example_sce), is_a("matrix"))
    expect_null(exprs(example_sce))
    expect_error(cpm(example_sce), "'cpm' not in names")

    exprs(example_sce) <- log2(calculateCPM(example_sce) + 1)
    expect_that(exprs(example_sce), is_a("matrix"))  

    norm_exprs(example_sce) <- log2(calculateCPM(example_sce) + 1)
    expect_that(norm_exprs(example_sce), is_a("matrix"))  

    stand_exprs(example_sce) <- log2(calculateCPM(example_sce) + 1)
    expect_that(stand_exprs(example_sce), is_a("matrix"))  

    fpkm(example_sce) <- log2(calculateFPKM(example_sce, lengths = 1000) + 1)
    expect_that(fpkm(example_sce), is_a("matrix"))  

    # Same again for sparse matrices.    
    library(Matrix)
    sparsified <- original
    counts(sparsified) <- as(counts(original), "dgCMatrix")
    
    expect_null(exprs(sparsified))
    expect_that(counts(sparsified), is_a("dgCMatrix"))
    exprs(sparsified) <- log2(calculateCPM(sparsified) + 1)
    expect_that(exprs(sparsified), is_a("dgeMatrix"))

    expect_error(cpm(sparsified), "'cpm' not in names")
    norm_exprs(sparsified) <- log2(calculateCPM(sparsified) + 1)
    expect_that(norm_exprs(sparsified), is_a("dgeMatrix"))  

    stand_exprs(sparsified) <- log2(calculateCPM(sparsified) + 1)
    expect_that(stand_exprs(sparsified), is_a("dgeMatrix"))  

    fpkm(sparsified) <- log2(calculateFPKM(sparsified, lengths= 1000) + 1)
    expect_that(fpkm(sparsified), is_a("dgeMatrix"))  
})
