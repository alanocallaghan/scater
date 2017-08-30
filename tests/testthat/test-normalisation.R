# Tests for normalisation methods

context("test expected usage")

test_that("normaliseExprs does not fail on input with zero-variance features", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    expect_that(normaliseExprs(example_sce, method = "none", 
                                feature_set = 1:100), 
                is_a("SingleCellExperiment"))
})

test_that("we can compute normalised expression values with TMM method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    keep_gene <- rowSums(counts(example_sce)) > 0
    example_sce <- example_sce[keep_gene,]
    
    example_sce <- normaliseExprs(example_sce, method = "TMM", 
                                     feature_set = 1:100)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
})

test_that("we can compute normalised expression values with RLE method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    keep_gene <- rowSums(counts(example_sce)) > 0
    example_sce <- example_sce[keep_gene,]
    
    example_sce <- normaliseExprs(example_sce, method = "RLE", 
                                     feature_set = 1:100)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
})

# test_that("we can compute normalised expression values with upperquartile 
#           method", {
#     data("sc_example_counts")
#     data("sc_example_cell_info")
#     pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#     example_sce <- newSCESet(countData = sc_example_counts + 1, 
#                                 phenoData = pd)
#     keep_gene <- rowSums(counts(example_sce)) > 0
#     example_sce <- example_sce[keep_gene,]
#     example_sce <- example_sce[
#         matrixStats::rowVars(counts(example_sce)) > 0,]
#     
#     example_sce <- normaliseExprs(example_sce, method = "upperquartile", 
#                                      feature_set = 1:200)
#     
#     expect_that(example_sce, is_a("SCESet"))
# })

test_that("we can compute normalised expression values with none method", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    keep_gene <- rowSums(counts(example_sce)) > 0
    example_sce <- example_sce[keep_gene,]
    
    example_sce <- normaliseExprs(example_sce, method = "none", 
                                     feature_set = 1:100)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
})

test_that("we can compute normalised expression values with a design matrix", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    keep_gene <- rowSums(counts(example_sce)) > 0
    example_sce <- calculateQCMetrics(example_sce[keep_gene,], 
                                         feature_controls = list(set1 = 1:40))
    design <- model.matrix(~example_sce$Cell_Cycle +
                               example_sce$pct_counts_top_200_features +
                               example_sce$total_features)
    example_sce <- normaliseExprs(example_sce, method = "none", design = design)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
    example_sce <- normaliseExprs(example_sce, method = "TMM", 
                                     design = design)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
    example_sce <- normaliseExprs(example_sce, method = "RLE", 
                                     design = design)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
    example_sce <- normaliseExprs(example_sce, exprs_values = "exprs",
                                  design = design)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
})

test_that("we can compute normalise the object", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    keep_gene <- rowSums(counts(example_sce)) > 0
    example_sce <- example_sce[keep_gene,]
    
    example_sce <- normaliseExprs(example_sce, method = "none", 
                                     feature_set = 1:100)
    ## normalize
    example_sce <- normalize(example_sce)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
    ## normalise
    example_sce <- normalise(example_sce)    
    expect_that(example_sce, is_a("SingleCellExperiment"))
    
    ## check error if no size factors
    sizeFactors(example_sce) <- NULL

    expect_warning(normalize(example_sce), "using library sizes")
    
})



