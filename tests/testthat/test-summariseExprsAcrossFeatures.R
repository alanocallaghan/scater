## Tests for summariseExprsAcrossFeatures()

context("test expected usage")

test_that("we can compute standard QC metrics", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    fd <- data.frame(gene_id = rownames(example_sce), 
                     feature_id = paste(
                         "feature", rep(1:500, each = 4, sep = "_")))
    rowData(example_sce) <- fd
    example_sceset_summarised <- 
        summariseExprsAcrossFeatures(example_sce, exprs_values = "counts")
    
    expect_that(example_sceset_summarised, is_a("SingleCellExperiment"))
    expect_equivalent(nrow(example_sceset_summarised), 500)
})


