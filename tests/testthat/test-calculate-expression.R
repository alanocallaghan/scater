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


test_that("calcAverage works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## calculate average counts
    ave_counts <- calcAverage(example_sceset)
    expected_vals <- c(305.551749, 325.719897, 183.090462, 162.143201, 1.231123)
    expect_true(all(abs(ave_counts[1:5] - expected_vals) < 1e-06))
})

