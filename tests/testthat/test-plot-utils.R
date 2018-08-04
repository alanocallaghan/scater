## Tests for the visualization variable picker.
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-utils.R")

example_sce <- normed

test_that("visualization variable picker works for columns with unnamed strings", {
    # Known metadata.
    out <- scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column")
    expect_identical(out$val, example_sce$Mutation_Status)
    expect_identical(out$name, "Mutation_Status")

    out <- scater:::.choose_vis_values(example_sce, "Cell_Cycle", mode="column")
    expect_identical(out$val, example_sce$Cell_Cycle)
    expect_identical(out$name, "Cell_Cycle")

    # Known gene exprs. 
    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    out <- scater:::.choose_vis_values(example_sce, "Gene_0100", mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0100",])
    expect_identical(out$name, "Gene_0100")

    # Known not to be either.
    expect_error(scater:::.choose_vis_values(example_sce, "WHEE", mode="column"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", search = "metadata"), "cannot find .* in metadata")

    # Responsive to search mode.
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "metadata"), "cannot find")

    example_sce$Gene_0002 <- seq_len(ncol(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "metadata")
    expect_identical(out_m$val, example_sce$Gene_0002)
    expect_identical(out_m$name, "Gene_0002")

    out_f <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out_f$name, "Gene_0002")

    out_a <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "any")
    expect_identical(out_a, out_m)
})

test_that("visualization variable picker works for columns with named strings", {
    example_sce$Gene_0002 <- seq_len(ncol(example_sce))

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "metadata"), mode="column")
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out$name, "Gene_0002")

    # The 'search' specifier will override any naming.
    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column", search="metadata") 
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "metadata"), mode="column", search="exprs") 
    expect_identical(out$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out$name, "Gene_0002")
})

test_that("visualization variable picker works fo columns with nested strings", {
    example_sce$blah <- DataFrame(X=runif(ncol(example_sce)), Y=sample(LETTERS, ncol(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="column", search="any")
    expect_identical(out$val, example_sce$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="metadata")
    expect_identical(out$val, example_sce$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    # Check error condition with search='exprs'.
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="exprs"), "character vector")

    # Works for internal fields.
    sizeFactors(example_sce) <- runif(ncol(example_sce))
    out <- scater:::.choose_vis_values(example_sce, c(NA, "size_factor"), mode="column", search="any")
    expect_identical(out$val, sizeFactors(example_sce))
    expect_identical(out$name, "size_factor")

    out <- scater:::.choose_vis_values(example_sce, c(NA, "size_factor"), mode="column", search="metadata")
    expect_identical(out$val, sizeFactors(example_sce))
    expect_identical(out$name, "size_factor")

    expect_error(scater:::.choose_vis_values(example_sce, "size_factor", mode="column", search="any"), "cannot find .* in any fields")
})

test_that("visualization variable picker works for columns with data.frames", {
    thing <- data.frame(B=runif(ncol(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "column")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    # Check errors.
    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "column"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "column"), "number of rows")
})

###################################################

set.seed(1312313)
rowData(example_sce) <- DataFrame(HAPPY=runif(nrow(example_sce)), SAD=rbinom(nrow(example_sce), 1, 0.5)==1)

test_that("visualization variable picker works for rows with strings", {
    # Known metadata.
    out <- scater:::.choose_vis_values(example_sce, "HAPPY", mode="row")
    expect_identical(out$val, rowData(example_sce)$HAPPY)
    expect_identical(out$name, "HAPPY")

    out <- scater:::.choose_vis_values(example_sce, "SAD", mode="row")
    expect_identical(out$val, rowData(example_sce)$SAD)
    expect_identical(out$name, "SAD")

    # Known exprs.
    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    out <- scater:::.choose_vis_values(example_sce, "Cell_010", mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_010"])
    expect_identical(out$name, "Cell_010")

    # Handles errors properly.
    expect_error(scater:::.choose_vis_values(example_sce, "whee", mode = "row", search = "any"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode = "row", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode = "row", search = "metadata"), "cannot find .* in metadata")

    # Responsive to search mode.
    expect_error(scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "metadata"), "cannot find")

    rowData(example_sce)$Cell_002 <- seq_len(nrow(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "metadata")
    expect_identical(out_m$val, rowData(example_sce)$Cell_002)
    expect_identical(out_m$name, "Cell_002")

    out_f <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out_f$name, "Cell_002")

    out_a <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "any")
    expect_identical(out_a, out_m)
})

test_that("visualization variable picker works for rows with named strings", {
    rowData(example_sce)$Cell_002 <- seq_len(nrow(example_sce))

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "metadata"), mode="row")
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out$name, "Cell_002")

    # Overriden by 'search'.
    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row", search="metadata") 
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "metadata"), mode="row", search="exprs") 
    expect_identical(out$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out$name, "Cell_002")
})

test_that("visualization variable picker works for rows with nested strings", {
    rowData(example_sce)$blah <- DataFrame(X=runif(nrow(example_sce)), Y=sample(LETTERS, nrow(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="row", search="any")
    expect_identical(out$val, rowData(example_sce)$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="metadata")
    expect_identical(out$val, rowData(example_sce)$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    # Throws errors if search='exprs'.
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="exprs"), "character vector")

    # Works for internal fields.
    isSpike(example_sce, "ERCC") <- 1:10
    out <- scater:::.choose_vis_values(example_sce, c(NA, "is_spike_ERCC"), mode="row", search="any")
    expect_identical(out$val, isSpike(example_sce, "ERCC"))
    expect_identical(out$name, "is_spike_ERCC")

    out <- scater:::.choose_vis_values(example_sce, c(NA, "is_spike"), mode="row", search="metadata")
    expect_identical(out$val, isSpike(example_sce))
    expect_identical(out$name, "is_spike")

    expect_error(scater:::.choose_vis_values(example_sce, "is_spike_ERCC", mode="row", search="any"), "cannot find .* in any fields")
})
     
test_that("visualization variable picker works for rows with data.frames", {
    thing <- data.frame(B=runif(nrow(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "row")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "row"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "row"), "number of rows")
})

###################################################

test_that("visualization variable picker discards one-level factors if requested", {
    example_sce$BLAH <- 0
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column")
    expect_identical(out$val, example_sce$BLAH)
    expect_identical(out$name, "BLAH")
                 
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", discard_solo = TRUE)
    expect_identical(out$val, NULL)
    expect_identical(out$name, NULL) 
})

test_that("visualization variable picker is responsive to exprs_values", {
    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")
})

test_that("visualization variable picker correctly coerces into factors", {
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
