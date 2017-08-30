## Test functions for QC

context("test controls functionality")

test_that("we can compute standard QC metrics", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))
})


test_that("we can compute standard QC metrics with feature controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    expect_error(
        example_sce <- calculateQCMetrics(example_sce, feature_controls = 1:20),
        "feature_controls should be named")
    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:20))
    expect_that(example_sce, is_a("SingleCellExperiment"))
    example_sce <- calculateQCMetrics(
        example_sce, 
        feature_controls = list(set1 = c(rep(TRUE, 20), 
                                         rep(FALSE, nrow(example_sce) - 20))))
    expect_that(example_sce, is_a("SingleCellExperiment"))
})


test_that("we can compute standard QC metrics with multiple sets of feature and 
          cell controls", {
              data("sc_example_counts")
              data("sc_example_cell_info")
              example_sce <- SingleCellExperiment(
                  assays = list(counts = sc_example_counts), 
                  colData = sc_example_cell_info)
              example_sce <- calculateQCMetrics(
                  example_sce, feature_controls = list(controls1 = 1:20, 
                                                          controls2 = 500:1000),
                  cell_controls = list(set_1 = 1:5, set_2 = 31:40))
              
              expect_that(example_sce, is_a("SingleCellExperiment"))
          })


test_that("computing standard QC metrics with FPKM data fails as expected", {
    gene_df <- data.frame(Gene = rownames(sc_example_counts))
    rownames(gene_df) <- gene_df$Gene
    example_sce <- SingleCellExperiment(
        assays = list(fpkm = sc_example_counts), 
        colData = sc_example_cell_info, rowData = gene_df)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_error(example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(set1 = 1:20)),
        "not in names")
})

test_that("computing standard QC metrics with TPM data fails as expected", {
    gene_df <- data.frame(Gene = rownames(sc_example_counts))
    rownames(gene_df) <- gene_df$Gene
    example_sce <- SingleCellExperiment(
        assays = list(tpm = sc_example_counts), 
        colData = sc_example_cell_info, rowData = gene_df)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_error(example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(set1 = 1:20)),
        "not in names")
})

test_that("failure is as expected for misspecified arg to plotExplanatoryVariables()", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    expect_error(plotExplanatoryVariables(example_sce, "expl"))
})


test_that("failure is as expected for input with zero-variance features", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    exprs(example_sce)[1:5,] <- 0
    expect_that(
        plotExplanatoryVariables(example_sce, "density"), is_a("ggplot"))
})


test_that("plotHighestExprs works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:500))
    expect_that(
        plotHighestExprs(example_sce, col_by_variable = "Mutation_Status"), 
        is_a("ggplot"))
})


test_that("plotExplanatoryVariables works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:500))
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    drop_genes <- apply(exprs(example_sce), 1, function(x) { var(x) == 0 })
    example_sce <- example_sce[!drop_genes, ]
    example_sce <- calculateQCMetrics(example_sce)
    vars <- colnames(colData(example_sce))[c(2:3, 5:14)]
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars[1]), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars, 
                                 method = "pairs"), 
        is_a("ggplot"))
    err_string <- "Only one variable"
    expect_error(plotExplanatoryVariables(example_sce, variables = vars[1], 
                                          method = "pairs"), err_string)
})


test_that("plotExprsFreqVsMean works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce)
    expect_that(plotExprsFreqVsMean(example_sce), is_a("ggplot"))
    
    example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(controls1 = 1:20,
                                           controls2 = 500:1000),
        cell_controls = list(set_1 = 1:5,
                             set_2 = 31:40))
    expect_that(plotExprsFreqVsMean(example_sce), is_a("ggplot"))
})



test_that("plotRLE works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)

    p <- plotRLE(example_sce, list(exprs = "logcounts", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "minimal")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, 
                 list(exprs = "exprs", counts = counts(example_sce)), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    expect_error(plotRLE(example_sce, 
                         list("exprs", counts = counts(example_sce)), 
                         c(TRUE, FALSE)), 
                 regexp = "exprs_mats must be a named list")
    
    expect_error(plotRLE(example_sce, 
                         list(exprs = "exprs", counts = counts(example_sce)[, 1:30]), 
                         c(TRUE, FALSE)), 
                 regexp = "Number of cells")
    
    expect_error(plotRLE(example_sce, 
                         list(exprs = "exprs"), style = "blah", 
                         c(TRUE, FALSE)), 
                 regexp = "should be one of")
    
})

