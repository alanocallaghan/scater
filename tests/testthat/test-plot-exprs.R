## Tests for plotExpression().
## library(scater); library(testthat); source("setup.R"); source("test-plot-exprs.R")

example_sce <- normed
rowData(example_sce)$ENS <- gsub("Gene", "ENS", rownames(example_sce))

test_that("plotExpression works for various 'x'", {
    for (gene_set in list("Gene_0001", rownames(example_sce)[1:5])) { # different numbers of genes, types of specification.
        for (x in list(NULL, "Cell_Cycle", "Gene_0100")) { # nothing, categorical, or continuous.
            expect_ggplot(plotExpression(example_sce, gene_set, x = x))
            expect_ggplot(plotExpression(example_sce, gene_set, x = x, colour_by = "Cell_Cycle"))
            expect_ggplot(plotExpression(example_sce, gene_set, x = x, size_by = "Gene_0001"))
            expect_ggplot(plotExpression(example_sce, gene_set, x = x, shape_by = "Treatment"))
        }
    }
})

test_that("plotExpression works for various aesthetics", {
    gene_set <- rownames(example_sce)[1:20]
    expect_ggplot(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", size_by = "Gene_0001"))
    expect_ggplot(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", shape_by = "Treatment"))
    expect_ggplot(plotExpression(example_sce, gene_set, size_by = "Gene_0001", shape_by = "Treatment"))
    expect_ggplot(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"))

    expect_ggplot(plotExpression(example_sce, gene_set, size_by = "Gene_0001", shape_by = "Treatment", by_exprs_values = "counts"))

    expect_ggplot(plotExpression(example_sce, rowData(example_sce)[1:10, "ENS"], colour_by = "ENS_0001", swap_rownames="ENS"))

    # Testing options when dealing with many genes and no 'x' specified.
    expect_ggplot(plotExpression(example_sce, gene_set, one_facet=FALSE))
    expect_ggplot(plotExpression(example_sce, gene_set, feature_colours=FALSE))
})

test_that("plotExpression works with semi-analysis options", {
    gene_set <- sample(rownames(example_sce)[5:15])
    expect_ggplot(plotExpression(example_sce, gene_set, show_violin=FALSE))
    expect_ggplot(plotExpression(example_sce, gene_set, show_median=TRUE))
    expect_ggplot(plotExpression(example_sce, gene_set, jitter_type="jitter"))

    expect_ggplot(plotExpression(example_sce, gene_set, x="Gene_0001", show_smooth=TRUE))
    expect_ggplot(plotExpression(example_sce, gene_set, x="Gene_0001", show_smooth=TRUE, show_se=FALSE))
})

test_that("plotExpression works for different exprs_values", {
    gene_set <- tail(rownames(example_sce))
    expect_ggplot(plotExpression(example_sce, gene_set, x = "Mutation_Status", exprs_values = "counts"))
    expect_ggplot(plotExpression(example_sce, gene_set, x = "Mutation_Status", exprs_values = "counts", log2_values = TRUE))
    expect_error(plotExpression(example_sce, rownames(example_sce)[1:6], exprs_values = "silly"), "not in names")

    # And on sparse matrices.
    sparsified <- example_sce
    logcounts(sparsified) <- as(logcounts(sparsified), "dgCMatrix")
    sparse <- plotExpression(sparsified, "Gene_0001")
    ref <- plotExpression(example_sce, "Gene_0001")
    expect_equal(sparse$data, ref$data)
})

test_that("plotExpression works for other fields", {
    gg <- plotExpression(example_sce, head(rownames(example_sce)),
        other_fields=c("Cell_Cycle", "Mutation_Status"))

    expect_true("Cell_Cycle" %in% colnames(gg$data))
    expect_true("Mutation_Status" %in% colnames(gg$data))

    expect_identical(gg$data$Cell_Cycle, rep(example_sce$Cell_Cycle, 6)) # default 'n' in head().

    # This should throw a warning.
    example_sce$colour_by <- example_sce$Cell_Cycle
    expect_warning(gg <- plotExpression(example_sce, head(rownames(example_sce)),
        colour_by="Cell_Cycle", other_fields=c("colour_by")), "duplicated")
})

test_that("plotExpression with scattermore", {
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100", scattermore = TRUE, point_size = 2))
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100", scattermore = TRUE, point_size = 2, color_by = "Cell_Cycle"))
})

test_that("plotExpression with binning", {
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100", bins = 10))
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100", bins = 10, hex = TRUE))
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100",
                                 bins = 10, colour_by = "Gene_0002", summary_fun = "mean"))
    expect_ggplot(plotExpression(example_sce, "Gene_0001", x = "Gene_0100",
                                 bins = 10, colour_by = "Gene_0002", hex = TRUE,
                                 summary_fun = "mean"))
})
