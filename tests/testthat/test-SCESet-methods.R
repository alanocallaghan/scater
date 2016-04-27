## Testing SCESet methods

test_that("sizeFactors() works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## expect NULL is returnd if size factors haven't been defined
    expect_null(sizeFactors(example_sceset))
    
    sizeFactors(example_sceset) <- rep(1, ncol(example_sceset))
    
    ## expect output still an SCESet
    expect_that(example_sceset, is_a("SCESet"))
    
    ## size factors should return a vector
    expect_that(sizeFactors(example_sceset), is_a("numeric"))
    
    ## reset size factors with NULL
    sizeFactors(example_sceset) <- NULL
    expect_null(sizeFactors(example_sceset))
})
