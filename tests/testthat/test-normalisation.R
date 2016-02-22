# Tests for normalisation methods

context("test expected usage")

test_that("we can compute normalised expression values with TMM method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    example_sceset <- normaliseExprs(example_sceset, method = "TMM", 
                                     feature_set = 1:100)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with RLE method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    example_sceset <- normaliseExprs(example_sceset, method = "RLE", 
                                     feature_set = 1:100)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with upperquartile 
          method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    example_sceset <- normaliseExprs(example_sceset, method = "upperquartile", 
                                     feature_set = 1:200)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with none method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    example_sceset <- normaliseExprs(example_sceset, method = "none", 
                                     feature_set = 1:100)
    
    expect_that(example_sceset, is_a("SCESet"))
})


