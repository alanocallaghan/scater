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
    exprs(example_sceset)[1:5,] <- 0
    err_string <- "Some features have zero variance"
    expect_error(plotExplanatoryVariables(example_sceset, "density"), err_string)
})


test_that("plotHighestExprs works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:500)
    expect_that(
        plotHighestExprs(example_sceset, col_by_variable = "Mutation_Status"), 
        is_a("ggplot"))
})

test_that("plotExplanatoryVariables works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:500)
    drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
    example_sceset <- example_sceset[!drop_genes, ]
    example_sceset <- calculateQCMetrics(example_sceset)
    vars <- names(pData(example_sceset))[c(2:3, 5:14)]
    expect_that(
        plotExplanatoryVariables(example_sceset, variables = vars), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sceset, variables = vars[1]), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sceset, variables = vars, 
                                 method = "pairs"), 
        is_a("ggplot"))
    err_string <- "Only one variable"
    expect_error(plotExplanatoryVariables(example_sceset, variables = vars[1], 
                                          method = "pairs"), err_string)
})


test_that("plotExprsFreqVsMean works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
    rownames(pd) <- pd$Cell
    ex_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
    ex_sceset <- calculateQCMetrics(ex_sceset)
    expect_that(plotExprsFreqVsMean(ex_sceset), is_a("ggplot"))
    
    ex_sceset <- calculateQCMetrics(
        ex_sceset, feature_controls = list(controls1 = 1:20,
                                           controls2 = 500:1000),
        cell_controls = list(set_1 = 1:5,
                             set_2 = 31:40))
    expect_that(plotExprsFreqVsMean(ex_sceset), is_a("ggplot"))
})

