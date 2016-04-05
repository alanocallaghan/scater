# Tests for normalisation methods

context("test expected usage")

test_that("failure is as expected for input with zero-variance features", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    err_string <- "Some features have zero variance"
    expect_error(normaliseExprs(example_sceset, method = "none", 
                                feature_set = 1:100), err_string)
})

test_that("we can compute normalised expression values with TMM method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    keep_gene <- rowSums(counts(example_sceset)) > 0
    example_sceset <- example_sceset[keep_gene,]
    
    example_sceset <- normaliseExprs(example_sceset, method = "TMM", 
                                     feature_set = 1:100)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with RLE method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    keep_gene <- rowSums(counts(example_sceset)) > 0
    example_sceset <- example_sceset[keep_gene,]
    
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
    keep_gene <- rowSums(counts(example_sceset)) > 0
    example_sceset <- example_sceset[keep_gene,]
    
    example_sceset <- normaliseExprs(example_sceset, method = "upperquartile", 
                                     feature_set = 1:200)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with none method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    keep_gene <- rowSums(counts(example_sceset)) > 0
    example_sceset <- example_sceset[keep_gene,]
    
    example_sceset <- normaliseExprs(example_sceset, method = "none", 
                                     feature_set = 1:100)
    
    expect_that(example_sceset, is_a("SCESet"))
})

test_that("we can compute normalised expression values with a design matrix", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    keep_gene <- rowSums(counts(example_sceset)) > 0
    example_sceset <- calculateQCMetrics(example_sceset[keep_gene,], 
                                         feature_controls = 1:40)
    design <- model.matrix(~example_sceset$Cell_Cycle +
                               example_sceset$pct_exprs_top_200_features +
                               example_sceset$pct_dropout +
                               example_sceset$total_features)
    example_sceset <- normaliseExprs(example_sceset, method = "none", 
                                     design = design)
    expect_that(example_sceset, is_a("SCESet"))
    
    example_sceset <- normaliseExprs(example_sceset, method = "TMM", 
                                     design = design)
    expect_that(example_sceset, is_a("SCESet"))
    
    example_sceset <- normaliseExprs(example_sceset, method = "RLE", 
                                     design = design)
    expect_that(example_sceset, is_a("SCESet"))
})

