## Test functions for QC plotting.
## library(scater); library(testthat); source("setup-sce.R"); source("test-qc-plot.R")

wt_qc <- sce
colData(wt_qc) <- cbind(colData(wt_qc), perCellQCMetrics(wt_qc))

library(Matrix)
sparsified <- wt_qc
counts(sparsified) <- as(counts(sparsified), "dgCMatrix")

#######################################################################

test_that("plotHighestExprs works on vanilla cases", {
    expect_ggplot(plotHighestExprs(wt_qc))
    expect_ggplot(plotHighestExprs(wt_qc))
    expect_ggplot(plotHighestExprs(wt_qc, as_percentage = FALSE))
})

test_that("plotHighestExprs' aesthetics choices work", {
    expect_ggplot(plotHighestExprs(wt_qc, colour_cells_by = "sum"))
    expect_ggplot(plotHighestExprs(wt_qc, colour_cells_by = "Mutation_Status"))
    expect_ggplot(plotHighestExprs(wt_qc, colour_cells_by = NULL))
    expect_ggplot(plotHighestExprs(wt_qc, colour_cells_by = "Gene_0001", by_exprs_values = "counts"))
})

test_that("plotHighestExprs works with different feature selections", {
    expect_ggplot(plotHighestExprs(wt_qc, n=Inf))

    # Responds to other sources for row names.
    rowData(wt_qc)$Whee <- paste("Feature", seq_len(nrow(wt_qc)))
    expect_ggplot(plotHighestExprs(wt_qc, feature_names_to_plot = "Whee"))

    dummy <- wt_qc
    rownames(dummy) <- NULL
    expect_ggplot(plotHighestExprs(dummy))
    
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
    expect_ggplot(plotHighestExprs(alt_sce, exprs_values="whee"))

    # Works for sparse matrices.
    expect_ggplot(plotHighestExprs(sparsified))
})
