## tests for plotting functions

context("test plotPCA, plotTSNE and plotDiffusionMap")

test_that("we can produce PCA plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotPCA(example_sceset), is_a("ggplot"))
})

test_that("we can produce t-SNE plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotTSNE(example_sceset), is_a("ggplot"))
})

test_that("we can produce Diffusion Map plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
              example_sceset <- newSCESet(countData = sc_example_counts,
                                          phenoData = pd)
              example_sceset <- calculateQCMetrics(example_sceset)

              expect_that(plotDiffusionMap(example_sceset), is_a("ggplot"))
          })

test_that("we can produce MDS plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
              example_sceset <- newSCESet(countData = sc_example_counts,
                                          phenoData = pd)
              example_sceset <- calculateQCMetrics(example_sceset)
              cellDist(example_sceset) <- as.matrix(
                  dist(t(exprs(example_sceset))))

              expect_that(plotMDS(example_sceset), is_a("ggplot"))
              expect_that(plotMDS(example_sceset,
                                          colour_by = "Cell_Cycle",
                                          shape_by = "Treatment",
                                          size_by = "Mutation_Status"),
                            is_a("ggplot"))
          })

context("test plotExpression")

test_that("we can produce expression plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotExpression(example_sceset, 1:4, "Cell_Cycle"), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4, "Gene_0004"), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4, scales = "fixed"), 
                is_a("ggplot"))
})

