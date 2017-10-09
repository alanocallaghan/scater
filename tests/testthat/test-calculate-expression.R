## test calculate expression

context("test calculation of TPM and FPKM")

test_that("we can calculate TPM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    effective_length <- rep(1000, 2000)
    tpm(example_sce) <- calculateTPM(example_sce, effective_length, 
                                        calc_from = "counts")
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(example_sce)), is_more_than(0))
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    tpm(sce10x) <- calculateTPM(sce10x, effective_length, 
                                     calc_from = "counts")
    expect_that(sce10x, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(sce10x)), is_more_than(0))
    
})


test_that("we can calculate CPM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    cpm(example_sce) <- calculateCPM(example_sce, use.size.factors = FALSE)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_that(sum(cpm(example_sce)), is_more_than(0))
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    cpm(sce10x) <- calculateCPM(sce10x, use.size.factors = FALSE)
    expect_that(sce10x, is_a("SingleCellExperiment"))
    expect_that(sum(cpm(sce10x)), is_more_than(0))
    
})


test_that("we can calculate FPKM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    effective_length <- rep(1000, 2000)
    fpkm(example_sce) <- calculateFPKM(
        example_sce, effective_length, use.size.factors = FALSE)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_that(sum(fpkm(example_sce)), is_more_than(0))
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    fpkm(sce10x) <- calculateFPKM(sce10x, effective_length, 
                                  use.size.factors = FALSE)
    expect_that(sce10x, is_a("SingleCellExperiment"))
    expect_that(sum(fpkm(sce10x)), is_more_than(0))
    
})


test_that("we can calculate TPM from FPKM", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    effective_length <- rep(1000, 2000)
    fpkm(example_sce) <- calculateFPKM(example_sce, effective_length,
                                       use.size.factors = FALSE)
    tpm(example_sce) <- calculateTPM(example_sce, effective_length, 
                                        calc_from = "fpkm")
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(example_sce)), is_more_than(0))
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    fpkm(sce10x) <- calculateFPKM(sce10x, effective_length,
                                       use.size.factors = FALSE)
    tpm(sce10x) <- calculateTPM(sce10x, effective_length, 
                                     calc_from = "fpkm")
    expect_that(sce10x, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(sce10x)), is_more_than(0))
    
})


test_that("nexprs works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    
    ## Testing nexprs on the counts themselves.
    expect_equal(nexprs(example_sce), unname(colSums(counts(example_sce) > 0)))
    expect_equal(nexprs(example_sce, byrow = TRUE), 
                 unname(rowSums(counts(example_sce) > 0)))
    expect_equal(nexprs(example_sce, subset_row = 20:40), 
                 unname(colSums(counts(example_sce)[20:40,] > 0)))
    expect_equal(nexprs(example_sce, byrow = TRUE, subset_col = 20:40), 
                 unname(rowSums(counts(example_sce)[,20:40] > 0)))

    ## Checking what happens when 'is_exprs' is available.
    assay(example_sce, "is_exprs") <- calcIsExprs(example_sce, 
                                                  lowerDetectionLimit = 5)
    expect_identical(assay(example_sce, "is_exprs"), counts(example_sce) > 5)
    expect_equal(nexprs(example_sce, lowerDetectionLimit = 5), 
                 unname(colSums(counts(example_sce) > 5)))
    expect_equal(nexprs(example_sce, lowerDetectionLimit = 5, byrow = TRUE), 
                 unname(rowSums(counts(example_sce) > 5)))
    expect_equal(nexprs(example_sce, lowerDetectionLimit = 5, subset_row = 20:40), 
                 unname(colSums(counts(example_sce)[20:40,] > 5)))
    expect_equal(nexprs(example_sce, lowerDetectionLimit = 5, byrow = TRUE, 
                        subset_col = 20:40), 
                 unname(rowSums(counts(example_sce)[,20:40] > 5)))
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    ## Testing nexprs on the counts themselves.
    expect_equal(nexprs(sce10x), unname(Matrix::colSums(counts(sce10x) > 0)))
    expect_equal(nexprs(sce10x, byrow = TRUE), 
                 unname(Matrix::rowSums(counts(sce10x) > 0)))
    expect_equal(nexprs(sce10x, subset_row = 20:40), 
                 unname(Matrix::colSums(counts(sce10x)[20:40,] > 0)))
    expect_equal(nexprs(sce10x, byrow = TRUE, subset_col = 20:40), 
                 unname(Matrix::rowSums(counts(sce10x)[,20:40] > 0)))
    
    ## Checking what happens when 'is_exprs' is available.
    assay(sce10x, "is_exprs") <- calcIsExprs(sce10x, 
                                                  lowerDetectionLimit = 5)
    expect_identical(assay(sce10x, "is_exprs"), counts(sce10x) > 5)
    expect_equal(nexprs(sce10x, lowerDetectionLimit = 5), 
                 unname(Matrix::colSums(counts(sce10x) > 5)))
    expect_equal(nexprs(sce10x, lowerDetectionLimit = 5, byrow = TRUE), 
                 unname(Matrix::rowSums(counts(sce10x) > 5)))
    expect_equal(nexprs(sce10x, lowerDetectionLimit = 5, subset_row = 20:40), 
                 unname(Matrix::colSums(counts(sce10x)[20:40,] > 5)))
    expect_equal(nexprs(sce10x, lowerDetectionLimit = 5, byrow = TRUE, 
                        subset_col = 20:40), 
                 unname(Matrix::rowSums(counts(sce10x)[,20:40] > 5)))
})

test_that("calcAverage works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    
    ## calculate average counts
    ave_counts <- calcAverage(example_sce)
    lib.sizes <- colSums(counts(example_sce))
    expected_vals <- colMeans(t(counts(example_sce)) / 
                                  (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    
    sce10x <- read10xResults(system.file("extdata", package = "scater"))
    ## calculate average counts
    ave_counts <- calcAverage(sce10x)
    lib.sizes <- Matrix::colSums(counts(sce10x))
    expected_vals <- Matrix::colMeans(Matrix::t(counts(sce10x)) / 
                                  (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    
})

