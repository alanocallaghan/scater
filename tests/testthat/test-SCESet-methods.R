## Testing SCESet methods

test_that("sizeFactors() works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## expect NULL is returnd if size factors haven't been defined
    expect_null(sizeFactors(example_sceset))
    
    sizeFactors(example_sceset, NULL) <- rep(1, ncol(example_sceset))
    
    ## expect output still an SCESet
    expect_that(example_sceset, is_a("SCESet"))
    
    ## size factors should return a vector
    expect_that(sizeFactors(example_sceset), is_a("numeric"))
    
    ## check with a named set of control features
    example_sceset <- calculateQCMetrics(example_sceset, 
                                         feature_controls = list(set1 = 1:40))
    sizeFactors(example_sceset, "set1") <- 2 ^ rnorm(ncol(example_sceset))
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(sizeFactors(example_sceset), is_a("numeric"))
    
    ## reset size factors with NULL
    sizeFactors(example_sceset) <- NULL
    expect_null(sizeFactors(example_sceset))
    
})

test_that("set_exprs and get_exprs work as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## expect get_exprs to return a matrix of exprs values
    expect_that(exprs(example_sceset), is_a("matrix"))
    expect_that(get_exprs(example_sceset, "exprs"), is_a("matrix"))
        
    ## expect NULL is returned if an assayData element hasn't been defined
    expect_warning(tmp <- get_exprs(example_sceset, "norm_exprs"), 
                   "The object does not contain")
    expect_null(tmp)
    
    ## expect we still get an SCESet after using set_exprs
    set_exprs(example_sceset, "new_exprs") <- exprs(example_sceset) / 2
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(example_sceset@assayData[["new_exprs"]], is_a("matrix"))
    expect_that(get_exprs(example_sceset, "new_exprs"), is_a("matrix"))
    
})

