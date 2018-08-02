## Tests for plotting functions
## This stress-tests the plotting utilities.
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-utils.R")

example_sce <- calculateQCMetrics(normed, exprs_values="counts")

test_that("visualization variable picker works properly: columns", {
    out <- scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column")
    expect_identical(out$val, example_sce$Mutation_Status)
    expect_identical(out$name, "Mutation_Status")

    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    # Responsive to exprs_values.
    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    # Responsive to search mode.
    example_sce$Gene_0002 <- seq_len(ncol(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "metadata")
    expect_identical(out_m$val, example_sce$Gene_0002)
    expect_identical(out_m$name, "Gene_0002")

    out_f <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out_f$name, "Gene_0002")

    out_a <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "any")
    expect_identical(out_a, out_m)

    # Handles named strings properly.
    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "metadata"), mode="column")
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column", search="metadata") # overriden by 'search'.
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    # Handles nested strings properly.
    example_sce$blah <- DataFrame(X=runif(ncol(example_sce)), Y=sample(LETTERS, ncol(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="column", search="any")
    expect_identical(out$val, example_sce$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="metadata")
    expect_identical(out$val, example_sce$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="exprs"), "character vector")

    # Works for input data.frames.
    thing <- data.frame(B=runif(ncol(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "column")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "column"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "column"), "number of rows")

    # Handles errors properly.
    expect_error(scater:::.choose_vis_values(example_sce, "whee", mode="column", search = "any"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", search = "metadata"), "cannot find .* in metadata")
})

test_that("visualization variable picker works properly: rows", {
    out <- scater:::.choose_vis_values(example_sce, "is_feature_control", mode="row")
    expect_identical(out$val, rowData(example_sce)$is_feature_control)
    expect_identical(out$name, "is_feature_control")

    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    # Responsive to exprs_values.
    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    # Responsive to search mode.
    rowData(example_sce)$Cell_002 <- seq_len(nrow(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "metadata")
    expect_identical(out_m$val, rowData(example_sce)$Cell_002)
    expect_identical(out_m$name, "Cell_002")

    out_f <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out_f$name, "Cell_002")

    out_a <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "any")
    expect_identical(out_a, out_m)

    # Handles named strings properly.
    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "metadata"), mode="row")
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row", search="metadata") # overriden by 'search'.
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    # Handles nested strings properly.
    rowData(example_sce)$blah <- DataFrame(X=runif(nrow(example_sce)), Y=sample(LETTERS, nrow(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="row", search="any")
    expect_identical(out$val, rowData(example_sce)$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="metadata")
    expect_identical(out$val, rowData(example_sce)$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="exprs"), "character vector")

    # Works for input data.frames.
    thing <- data.frame(B=runif(nrow(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "row")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "row"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "row"), "number of rows")

    # Handles errors properly.
    expect_error(scater:::.choose_vis_values(example_sce, "whee", mode = "row", search = "any"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode = "row", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode = "row", search = "metadata"), "cannot find .* in metadata")
})

test_that("visualization variable picker works properly: misc", {
    example_sce$BLAH <- 0
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column")
    expect_identical(out$val, example_sce$BLAH)
    expect_identical(out$name, "BLAH")
                 
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", discard_solo = TRUE)
    expect_identical(out$val, NULL)
    expect_identical(out$name, NULL) 

    # Checking what happens with coerced factors.
    example_sce$BLAH <- rep(LETTERS, length.out=ncol(example_sce))

    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column")
    expect_identical(out$val, example_sce$BLAH)
    expect_identical(out$name, "BLAH")

    out2 <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE)
    expect_identical(out2$val, as.factor(out$val))                 
    expect_identical(out2$name, "BLAH")

    expect_error(scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE, level_limit = 10), 
                 "number of unique levels")

    out3 <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE, level_limit = 100)
    expect_identical(out2, out3)
})
