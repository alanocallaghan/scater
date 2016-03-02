## Testing functions for classes used ##

context("tests on inputs")

test_that("tests for presence of cellData or counts", {
    expect_that(newSCESet(),
                throws_error("Require at least one of exprsData, tpmData, fpkmData or countData argument."))
})

test_that("example datasets work", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
    example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
    example_sceset

    expect_that(example_sceset, is_a("SCESet"))
})


context("test manipulations")

test_that("we can subset the example SCESet", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
    example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
    example_sceset <- plotPCA(example_sceset, return_SCESet=TRUE)
    ex_subset <- example_sceset[1:200, sample(1:40, 25)]

    expect_equal(as.integer(nrow(ex_subset)), 200L)
    expect_equal(as.integer(ncol(ex_subset)), 25L)
    expect_equal(nrow(redDim(ex_subset)), 25L)
})




