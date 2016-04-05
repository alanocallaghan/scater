## Test functions for QC

context("test controls functionality")

test_that("we can compute standard QC metrics", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute standard QC metrics with feature controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:20)
    
    expect_that(example_sceset, is_a("SCESet"))
})


test_that("we can compute standard QC metrics with multiple sets of feature and 
          cell controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(
        example_sceset, feature_controls = list(controls1 = 1:20, 
                                                controls2 = 500:1000),
        cell_controls = list(set_1 = 1:5, set_2 = 31:40))
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("failure is as expected for misspecified arg to plotExplanatoryVariables()", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    expect_error(plotExplanatoryVariables(example_sceset, "expl"))
})


test_that("failure is as expected for input with zero-variance features", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    err_string <- "Some features have zero variance"
    expect_error(plotExplanatoryVariables(example_sceset, "density"), err_string)
})





