## test calculate expression
## library(scater); library(testthat); source("setup-sce.R"); source("test-calculate-expression.R")

library(Matrix)
original <- sce
sparsified <- original
counts(sparsified) <- as(counts(original), "dgCMatrix")

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

    ## Repeating with subsets.
    sub1 <- calculateCPM(counts(original), subset_row=1:10)
    expect_identical(sub1, calculateCPM(counts(original)[1:10,]))

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calculateCPM(counts(original), subset_row=logi)
    expect_identical(sub2, calculateCPM(counts(original)[logi,]))

    chosen <- sample(rownames(original), 20)
    sub3 <- calculateCPM(counts(original), subset_row=chosen)
    expect_identical(sub3, calculateCPM(counts(original)[chosen,]))
})

test_that("we can calculate FPKM from counts", {
    effective_length <- runif(nrow(original), 1000, 2000)
    fpkms <- calculateFPKM(original, effective_length)

    ref <- counts(original)/(effective_length/1e3)
    ref <- t(t(ref)/(colSums(counts(original))/1e6))
    expect_equal(fpkms, ref)

    # Repeating with subsets.
    out <- calculateFPKM(original, effective_length, subset_row=1:10)
    sub <- calculateFPKM(original[1:10,], effective_length[1:10])
    expect_equal(out, sub)

    # Handles sparse matrices.
    tout2 <- calculateFPKM(sparsified, effective_length)
    expect_s4_class(tout2, "dgCMatrix")
    expect_equal(ref, as.matrix(tout2))
})

test_that("we can calculate TPM from counts", {
    effective_length <- runif(nrow(original), 1000, 2000)
    tout <- calculateTPM(original, effective_length)

    ref <- counts(original)/effective_length
    ref <- t(t(ref)/(colSums(ref)/1e6))
    expect_equal(tout, ref)

    # Behaves when length is not supplied.
    expect_equal(calculateTPM(original, NULL), calculateCPM(original))

    # Repeating with subsets.
    out <- calculateTPM(original, effective_length, subset_row=1:10)
    sub <- calculateTPM(original[1:10,], effective_length[1:10])
    expect_equal(out, sub)

    # Handles sparse matrices.
    tout2 <- calculateTPM(sparsified, effective_length)
    expect_s4_class(tout2, "dgCMatrix")
    expect_equal(tout, as.matrix(tout2))
})

test_that("nexprs works as expected", {
    ## Testing nexprs on the counts themselves.
    expect_equal(nexprs(original), unname(colSums(counts(original) > 0)))
    expect_equal(nexprs(original), nexprs(counts(original)))

    expect_equal(nexprs(original, byrow = TRUE), unname(rowSums(counts(original) > 0)))
    expect_equal(nexprs(original, byrow = TRUE), nexprs(counts(original), byrow = TRUE))

    expect_equal(nexprs(original, subset_row = 20:40), 
                 unname(colSums(counts(original)[20:40,] > 0)))
    expect_equal(nexprs(original, detection_limit=5),
                 unname(colSums(counts(original) > 5)))

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
    expected_vals <- colMeans(t(counts(original)) / (lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
    expect_equal(ave_counts, calcAverage(counts(original)))

    ## Responsive to other assay names.
    whee <- original
    assayNames(whee) <- "whee"
    whee_counts <- calcAverage(whee, exprs_values="whee")
    expect_identical(whee_counts, ave_counts)

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
    expect_equal(calcAverage(counts(original)), calcAverage(original, use_size_factors=FALSE))
    expect_equal(calcAverage(counts(spiked)), calcAverage(spiked, use_size_factors=FALSE))

    new_sf <- runif(ncol(spiked))
    ave_counts <- calcAverage(spiked, use_size_factors=new_sf)
    spiked2 <- spiked
    sizeFactors(spiked2) <- new_sf
    expect_equal(calcAverage(spiked2), ave_counts)

    # Warnings are correctly thrown (or not) when no size factors are available for spike-ins.
    sizeFactors(spiked2, "WHEE") <- NULL
    expect_warning(calcAverage(spiked2), "spike-in set")
    expect_warning(calcAverage(spiked2, use_size_factors=new_sf), "spike-in set")
    expect_warning(calcAverage(spiked2, use_size_factors=FALSE), NA)

    ## Repeating with a sparse matrix.    
    ave_counts <- calcAverage(sparsified)
    expect_equal(ave_counts, calcAverage(original, use_size_factors=FALSE))  

    ## Repeating with subsets.
    sub1 <- calcAverage(counts(original), subset_row=1:10)
    expect_identical(sub1, calcAverage(counts(original)[1:10,]))

    logi <- rbinom(nrow(original), 1, 0.5)==1
    sub2 <- calcAverage(counts(original), subset_row=logi)
    expect_identical(sub2, calcAverage(counts(original)[logi,]))

    chosen <- sample(rownames(original), 20)
    sub3 <- calcAverage(counts(original), subset_row=chosen)
    expect_identical(sub3, calcAverage(counts(original)[chosen,]))
})

