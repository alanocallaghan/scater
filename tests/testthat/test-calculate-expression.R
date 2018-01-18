## test calculate expression

context("test calculation of TPM and FPKM")

data("sc_example_counts")
data("sc_example_cell_info")
original <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)

library(Matrix)
sparsified <- original
counts(sparsified) <- as(counts(original), "dgCMatrix")

test_that("we can calculate TPM from counts", {
    effective_length <- rep(1000, 2000)
    tpm(original) <- calculateTPM(original, effective_length, 
                                        calc_from = "counts")
    
    expect_that(original, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(original)), is_more_than(0))
    
    tpm(sparsified) <- calculateTPM(sparsified, effective_length, 
                                     calc_from = "counts")
    expect_that(sparsified, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(sparsified)), is_more_than(0))
    
})


test_that("we can calculate CPM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    cpm(original) <- calculateCPM(original, use_size_factors = FALSE)
    
    expect_that(original, is_a("SingleCellExperiment"))
    expect_that(sum(cpm(original)), is_more_than(0))
    
    sparsified <- read10xResults(system.file("extdata", package = "scater"))
    cpm(sparsified) <- calculateCPM(sparsified, use_size_factors = FALSE)
    expect_that(sparsified, is_a("SingleCellExperiment"))
    expect_that(sum(cpm(sparsified)), is_more_than(0))
    
})


test_that("we can calculate FPKM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    effective_length <- rep(1000, 2000)
    fpkm(original) <- calculateFPKM(
        original, effective_length, use_size_factors = FALSE)
    
    expect_that(original, is_a("SingleCellExperiment"))
    expect_that(sum(fpkm(original)), is_more_than(0))
    
    sparsified <- read10xResults(system.file("extdata", package = "scater"))
    fpkm(sparsified) <- calculateFPKM(sparsified, effective_length, 
                                  use_size_factors = FALSE)
    expect_that(sparsified, is_a("SingleCellExperiment"))
    expect_that(sum(fpkm(sparsified)), is_more_than(0))
    
})


test_that("we can calculate TPM from FPKM", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    effective_length <- rep(1000, 2000)
    fpkm(original) <- calculateFPKM(original, effective_length,
                                       use_size_factors = FALSE)
    tpm(original) <- calculateTPM(original, effective_length, 
                                        calc_from = "fpkm")
    expect_that(original, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(original)), is_more_than(0))
    
    sparsified <- read10xResults(system.file("extdata", package = "scater"))
    fpkm(sparsified) <- calculateFPKM(sparsified, effective_length,
                                       use_size_factors = FALSE)
    tpm(sparsified) <- calculateTPM(sparsified, effective_length, 
                                     calc_from = "fpkm")
    expect_that(sparsified, is_a("SingleCellExperiment"))
    expect_that(sum(tpm(sparsified)), is_more_than(0))
    
})


test_that("nexprs works as expected", {
    ## Testing nexprs on the counts themselves.
    expect_equal(nexprs(original), unname(colSums(counts(original) > 0)))
    expect_equal(nexprs(original), nexprs(counts(original)))

    expect_equal(nexprs(original, byrow = TRUE), 
                 unname(rowSums(counts(original) > 0)))
    expect_equal(nexprs(original, subset_row = 20:40), 
                 unname(colSums(counts(original)[20:40,] > 0)))
    expect_equal(nexprs(original, byrow = TRUE, subset_col = 20:40), 
                 unname(rowSums(counts(original)[,20:40] > 0)))
    expect_equal(nexprs(original, byrow = TRUE, detection_limit=5),
                 unname(rowSums(counts(original) > 5)))
    
    ## Testing nexprs on a sparse matrix.
    expect_equal(nexprs(sparsified), unname(Matrix::colSums(counts(sparsified) > 0)))
    expect_equal(nexprs(sparsified), nexprs(counts(sparsified)))
})


test_that("calcAverage works as expected", {
    ## Calculate average counts
    ave_counts <- calcAverage(original)
    lib.sizes <- colSums(counts(original))
    expected_vals <- colMeans(t(counts(original)) / 
                                  (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    expect_equal(ave_counts, calcAverage(counts(original)))

    ## Responsive to size factors.
    sizeFactors(original) <- runif(ncol(original))
    ave_counts <- calcAverage(original)
    expect_equal(ave_counts, calcAverage(counts(original),
                                         size_factors=sizeFactors(original)))
    
    sf <- sizeFactors(original) 
    sf <- sf/mean(sf)    
    expected_vals <- colMeans(t(counts(original)) / sf)
    expect_equal(ave_counts, expected_vals)

    ## Repeating with a sparse matrix.    
    ave_counts <- calcAverage(sparsified)
    lib.sizes <- Matrix::colSums(counts(sparsified))
    expected_vals <- Matrix::colMeans(Matrix::t(counts(sparsified)) / 
                                  (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)  
})

