# Tests for tximport wrapper

context("test input usage")

test_that("the system works", {
    
    expect_that(read10xResults(system.file("extdata", package="scater")),
                is_a("SingleCellExperiment"))
})
