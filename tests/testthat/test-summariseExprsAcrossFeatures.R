## Tests for summariseExprsAcrossFeatures()

context("test expected usage")

test_that("we can compute standard QC metrics", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    fd <- new("AnnotatedDataFrame", data = 
    data.frame(gene_id = featureNames(example_sceset), 
    feature_id = paste("feature", rep(1:500, each = 4, sep = "_"))))
    fData(example_sceset) <- fd
    example_sceset_summarised <- 
        summariseExprsAcrossFeatures(example_sceset, exprs_values = "counts")
    
    expect_that(example_sceset_summarised, is_a("SCESet"))
    expect_equivalent(nrow(example_sceset_summarised), 500)
})


