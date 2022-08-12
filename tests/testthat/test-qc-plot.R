## Test functions for QC plotting.
## library(scater); library(testthat); source("setup-sce.R"); source("test-qc-plot.R")

wt_qc <- sce
colData(wt_qc) <- cbind(colData(wt_qc), perCellQCMetrics(wt_qc))

library(Matrix)
sparsified <- wt_qc
counts(sparsified) <- as(as(as(counts(sparsified), "dMatrix"), "generalMatrix"), "CsparseMatrix")

#######################################################################

test_that("plotHighestExprs works on vanilla cases", {
    expect_s3_class(plotHighestExprs(wt_qc), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, as_percentage = FALSE), "ggplot")
})

test_that("plotHighestExprs' aesthetics choices work", {
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "sum"), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "Mutation_Status"), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = NULL), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
})

test_that("plotHighestExprs works with different feature selections", {
    expect_s3_class(plotHighestExprs(wt_qc, n=Inf), "ggplot")

    # Responds to other sources for row names.
    rowData(wt_qc)$Whee <- paste("Feature", seq_len(nrow(wt_qc)))
    expect_s3_class(plotHighestExprs(wt_qc, feature_names_to_plot = "Whee"), "ggplot")

    dummy <- wt_qc
    rownames(dummy) <- NULL
    expect_s3_class(plotHighestExprs(dummy), "ggplot")
    
    # Discarding and subset exclusion yield the same results.
    discard <- rbinom(nrow(wt_qc), 1, 0.5)==1
    p <- plotHighestExprs(wt_qc, drop_features=discard, as_percentage=FALSE)
    ref <- plotHighestExprs(wt_qc[!discard,], as_percentage=FALSE)
    expect_equal(p$data, ref$data)

    # Replacement rownames respect original ordering when discarding.
    discard <- rbinom(nrow(dummy), 1, 0.5)==1
    p <- plotHighestExprs(dummy, drop_features=discard, as_percentage=FALSE)
    dummy2 <- dummy
    rownames(dummy2) <- sprintf("Feature %i", seq_len(nrow(dummy2)))
    ref <- plotHighestExprs(dummy2[!discard,], as_percentage=FALSE)
    expect_equal(p$data, ref$data)
})

test_that("plotHighestExprs works on alternative exprs", {
    alt_sce <- wt_qc 
    assayNames(alt_sce) <- "whee"
    expect_s3_class(plotHighestExprs(alt_sce, exprs_values="whee"), "ggplot")

    # Works for sparse matrices.
    expect_s3_class(plotHighestExprs(sparsified), "ggplot")
})
