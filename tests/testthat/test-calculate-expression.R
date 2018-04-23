## test calculate expression
## library(scater); library(testthat); source("test-calculate-expression.R")

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
    cpm_out <- calculateCPM(original)
    expect_equal(cpm_out, t(t(counts(original))/(colSums(counts(original)/1e6))))
    
    ## Responsive to size factors.
    sizeFactors(original) <- runif(ncol(original))
    cpm_out <- calculateCPM(original)
    expect_equal(cpm_out, calculateCPM(counts(original), use_size_factors=sizeFactors(original)))
   
    FUN <- function(counts, sf, libsize = colSums(counts)) {
        eff_lib <- sf/mean(sf) * mean(libsize)
        t(t(counts) / (eff_lib/1e6))
    } 
    expect_equal(cpm_out, FUN(counts(original), sizeFactors(original)))

    ## Responsive to multiple size factors.
    spiked <- original
    sizeFactors(spiked, "WHEE") <- runif(ncol(original))
    is_spike <- 10:20
    isSpike(spiked, "WHEE") <- is_spike
    cpm_out <- calculateCPM(spiked)
    expect_equal(cpm_out[is_spike,], FUN(counts(original)[is_spike,], sizeFactors(spiked, "WHEE"), colSums(counts(original))))
    expect_equal(cpm_out[-is_spike,], FUN(counts(original)[-is_spike,], sizeFactors(spiked), colSums(counts(original))))

    # Ignores or overrides the size factors if requested.
    expect_equal(calculateCPM(counts(original)),
                 calculateCPM(original, use_size_factors=FALSE))
    expect_equal(calculateCPM(counts(spiked)),
                 calculateCPM(spiked, use_size_factors=FALSE))

    new_sf <- runif(ncol(spiked))
    cpm_out <- calculateCPM(spiked, use_size_factors=new_sf)
    spiked2 <- spiked
    sizeFactors(spiked2) <- new_sf
    expect_equal(calculateCPM(spiked2), cpm_out)

    # Checking that it works on a sparse matrix. 
    cpm_out <- calculateCPM(sparsified)
    expect_equal(as.matrix(cpm_out), calculateCPM(original, use_size_factors=FALSE))
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

set.seed(10000)
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
    expect_equal(ave_counts, calcAverage(counts(original), use_size_factors=sizeFactors(original)))
    
    sf <- sizeFactors(original) 
    sf <- sf/mean(sf)    
    expected_vals <- colMeans(t(counts(original)) / sf)
    expect_equal(ave_counts, expected_vals)

    ## Responsive to multiple size factors.
    spiked <- original
    sizeFactors(spiked, "WHEE") <- runif(ncol(original))
    is_spike <- 10:20
    isSpike(spiked, "WHEE") <- is_spike
    ave_counts <- calcAverage(spiked)
    expect_equal(ave_counts[is_spike], calcAverage(counts(original)[is_spike,], use_size_factors=sizeFactors(spiked, "WHEE")))
    expect_equal(ave_counts[-is_spike], calcAverage(counts(original)[-is_spike,], use_size_factors=sizeFactors(spiked)))

    # Ignores or overrides the size factors if requested.
    expect_equal(calcAverage(counts(original)),
                 calcAverage(original, use_size_factors=FALSE))
    expect_equal(calcAverage(counts(spiked)),
                 calcAverage(spiked, use_size_factors=FALSE))

    new_sf <- runif(ncol(spiked))
    ave_counts <- calcAverage(spiked, use_size_factors=new_sf)
    spiked2 <- spiked
    sizeFactors(spiked2) <- new_sf
    expect_equal(calcAverage(spiked2), ave_counts)

    ## Repeating with a sparse matrix.    
    ave_counts <- calcAverage(sparsified)
    expect_equal(ave_counts, calcAverage(original, use_size_factors=FALSE))  

    ## Repeating with subsets.
    ref <- calcAverage(counts(original))
    sub1 <- calcAverage(counts(original), subset_row=1:10)
    expect_identical(ref[1:10], sub1)

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calcAverage(counts(original), subset_row=logi)
    expect_identical(ref[logi], sub2)

    chosen <- sample(rownames(original), 20)
    sub3 <- calcAverage(counts(original), subset_row=chosen)
    expect_identical(ref[chosen], sub3)
})

