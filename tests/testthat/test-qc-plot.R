## Test functions for QC plotting.
## library(scater); library(testthat); source("setup-sce.R"); source("test-qc-plot.R")

wo_qc <- sce

wt_qc <- calculateQCMetrics(wo_qc, 
    feature_controls = list(set1 = 1:500),
    cell_controls = list(whee = 1:10))

wt_qc_compact <- calculateQCMetrics(wo_qc, 
    feature_controls = list(set1 = 1:500),
    cell_controls = list(whee = 1:10),
    compact=TRUE)

library(Matrix)
sparsified <- wt_qc
counts(sparsified) <- as(counts(sparsified), "dgCMatrix")

#######################################################################

test_that("the QC hunter works as expected on columns", {
    for (x in c("total_counts", "total_features_by_counts",
                "total_counts_endogenous", "total_counts_feature_control",
                "total_counts_set1", "pct_counts_set1")) {

        expect_error(scater:::.qc_hunter(wo_qc, x, mode = "column"), "failed")
        expect_warning(out <- scater:::.qc_hunter(wo_qc, x, mode = "column", error = FALSE), "failed")
        expect_identical(out, NULL)
        expect_identical(scater:::.qc_hunter(wt_qc, x, mode = "column"), x)
    }

    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts", mode = "column"), 
                     c("scater_qc", "all", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_features_by_counts", mode = "column"), 
                     c("scater_qc", "all", "total_features_by_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_endogenous", mode = "column"), 
                     c("scater_qc", "endogenous", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_feature_control", mode = "column"), 
                     c("scater_qc", "feature_control", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_set1", mode = "column"), 
                     c("scater_qc", "feature_control_set1", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "pct_counts_set1", mode = "column"), 
                     c("scater_qc", "feature_control_set1", "pct_counts"))
})

test_that("the QC hunter works as expected on rows", {
    for (x in c("total_counts", "n_cells_by_counts",
                "total_counts_non_control", "total_counts_cell_control",
                "total_counts_whee", "pct_counts_whee")) {

        expect_error(scater:::.qc_hunter(wo_qc, x, mode = "row"), "failed")
        expect_warning(out <- scater:::.qc_hunter(wo_qc, x, mode = "row", error = FALSE), "failed")
        expect_identical(out, NULL)
        expect_identical(scater:::.qc_hunter(wt_qc, x, mode = "row"), x)
    }

    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts", mode = "row"), 
                     c("scater_qc", "all", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "n_cells_by_counts", mode = "row"), 
                     c("scater_qc", "all", "n_cells_by_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_non_control", mode = "row"), 
                     c("scater_qc", "non_control", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_cell_control", mode = "row"), 
                     c("scater_qc", "cell_control", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "total_counts_whee", mode = "row"), 
                     c("scater_qc", "cell_control_whee", "total_counts"))
    expect_identical(scater:::.qc_hunter(wt_qc_compact, "pct_counts_whee", mode = "row"), 
                     c("scater_qc", "cell_control_whee", "pct_counts"))
})

#######################################################################

test_that("plotHighestExprs works as expected", {
    expect_s3_class(plotHighestExprs(wt_qc), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc_compact), "ggplot")

    # Checking out the error messages.
    expect_error(plotHighestExprs(wo_qc, colour_cells_by = NULL), "failed to find")
    expect_error(plotHighestExprs(wo_qc, controls = NULL), "failed to find")
    expect_s3_class(plotHighestExprs(wo_qc, controls = NULL, colour_cells_by = NULL), "ggplot")

    # Checking out the options.
    expect_s3_class(plotHighestExprs(wt_qc, n=Inf), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, drop_features=1:20), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, as_percentage = FALSE), "ggplot")
    
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "Mutation_Status"), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = NULL), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, controls = "is_feature_control_set1"), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, controls = NULL), "ggplot")

    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "Mutation_Status", by_show_single = FALSE), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc, colour_cells_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")

    rowData(wt_qc)$Whee <- paste("Feature", seq_len(nrow(wt_qc)))
    expect_s3_class(plotHighestExprs(wt_qc, feature_names_to_plot = "Whee"), "ggplot")

    # Checking that the variable pickers work.
    expect_s3_class(plotHighestExprs(wt_qc_compact, controls = c("scater_qc", "is_feature_control_set1")), "ggplot")
    expect_s3_class(plotHighestExprs(wt_qc_compact, colour_cells_by = c("scater_qc", "all", "total_counts")), "ggplot")

    # Recognizes alternative exprs_values.
    alt_sce <- wt_qc 
    assayNames(alt_sce) <- "whee"
    expect_error(plotHighestExprs(alt_sce, exprs_values="whee"), "failed to find")
    alt_sce$total_features_by_whee <- runif(ncol(alt_sce))
    expect_s3_class(plotHighestExprs(alt_sce, exprs_values="whee"), "ggplot")

    # Works for sparse matrices.
    expect_s3_class(plotHighestExprs(sparsified), "ggplot")
})

#######################################################################

test_that("plotExprsFreqVsMean works as expected", {
    expect_s3_class(plotExprsFreqVsMean(wt_qc), "ggplot")

    # Checking arguments are passed to plotRowData.
    expect_s3_class(plotExprsFreqVsMean(wt_qc, legend = FALSE), "ggplot")
    expect_s3_class(plotExprsFreqVsMean(wt_qc, size_by = 'n_cells_by_counts'), "ggplot")

    # Checking we can turn off control colouring and other settings.
    expect_s3_class(plotExprsFreqVsMean(wt_qc, controls=NULL), "ggplot")
    expect_s3_class(plotExprsFreqVsMean(wt_qc, show_smooth=FALSE), "ggplot")
    expect_s3_class(plotExprsFreqVsMean(wt_qc, show_se=FALSE), "ggplot")

    # Avoid computing metrics if there are no controls.
    rowData(wt_qc)$is_feature_control <- FALSE
    expect_s3_class(plotExprsFreqVsMean(wt_qc), "ggplot")

    # Recognizes alternative exprs_values.
    alt_sce <- wt_qc 
    assayNames(alt_sce) <- "whee"
    expect_error(plotExprsFreqVsMean(alt_sce, exprs_values="whee"), "failed to find")
    rowData(alt_sce)$n_cells_by_whee <- runif(ncol(alt_sce))
    rowData(alt_sce)$mean_whee <- runif(ncol(alt_sce))
    expect_s3_class(plotExprsFreqVsMean(alt_sce, exprs_values="whee"), "ggplot")

    # Checking errors.
    rowData(wo_qc)$whee <- runif(nrow(wo_qc))
    rowData(wo_qc)$stuff <- runif(nrow(wo_qc))
    expect_error(plotExprsFreqVsMean(wo_qc), "failed to find")
    expect_error(plotExprsFreqVsMean(wo_qc, freq_exprs="whee"), "failed to find")
    expect_s3_class(plotExprsFreqVsMean(wo_qc, freq_exprs="whee", mean_exprs="stuff", controls=NULL), "ggplot")
})
