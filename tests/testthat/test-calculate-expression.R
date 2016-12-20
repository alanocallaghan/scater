## test calculate expression

context("test calculation of TPM and FPKM")

test_that("we can calculate TPM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    effective_length <- rep(1000, 2000)
    tpm(example_sceset) <- calculateTPM(example_sceset, effective_length, 
                                        calc_from = "counts")
    
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(sum(tpm(example_sceset)),is_more_than(0))
})


test_that("we can calculate FPKM from counts", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    effective_length <- rep(1000, 2000)
    fpkm(example_sceset) <- calculateFPKM(example_sceset, effective_length)
    
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(sum(fpkm(example_sceset)),is_more_than(0))
})


test_that("we can calculate TPM from FPKM", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    effective_length <- rep(1000, 2000)
    fpkm(example_sceset) <- calculateFPKM(example_sceset, effective_length)
    tpm(example_sceset) <- calculateTPM(example_sceset, effective_length, 
                                        calc_from = "fpkm")
    
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(sum(tpm(example_sceset)),is_more_than(0))
})


test_that("nexprs works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## Testing nexprs on the counts themselves.
    expect_equal(nexprs(example_sceset), colSums(counts(example_sceset) > 0))
    expect_equal(nexprs(example_sceset, byrow=TRUE), rowSums(counts(example_sceset) > 0))
    expect_equal(nexprs(example_sceset, subset_row=20:40), colSums(counts(example_sceset)[20:40,] > 0))
    expect_equal(nexprs(example_sceset, byrow=TRUE, subset_col=20:40), rowSums(counts(example_sceset)[,20:40] > 0))

    ## Checking what happens when 'is_exprs' is available.
    is_exprs(example_sceset) <- calcIsExprs(example_sceset, lowerDetectionLimit=5)
    expect_identical(is_exprs(example_sceset), counts(example_sceset) > 5)
    expect_equal(nexprs(example_sceset), colSums(counts(example_sceset) > 5))
    expect_equal(nexprs(example_sceset, byrow=TRUE), rowSums(counts(example_sceset) > 5))
    expect_equal(nexprs(example_sceset, subset_row=20:40), colSums(counts(example_sceset)[20:40,] > 5))
    expect_equal(nexprs(example_sceset, byrow=TRUE, subset_col=20:40), rowSums(counts(example_sceset)[,20:40] > 5))
})

test_that("calcAverage works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## calculate average counts
    ave_counts <- calcAverage(example_sceset)
    lib.sizes <- colSums(counts(example_sceset))
    expected_vals <- colMeans(t(counts(example_sceset))/(lib.sizes/mean(lib.sizes)))
    expect_equal(ave_counts, expected_vals)
})

