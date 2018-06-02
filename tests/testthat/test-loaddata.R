# Tests for various data input methods
# library(scater); library(testthat); source("test-loaddata.R")

context("test expected usage")

library(Matrix)
test_that("readSparseCounts works as expected", {
    a <- matrix(rpois(10000, lambda=1), ncol=50)
    rownames(a) <- paste0("Gene", seq_len(nrow(a)))
    colnames(a) <- paste0("Cell", seq_len(ncol(a)))

    ofile <- tempfile()
    write.table(a, file=ofile, sep="\t", quote=FALSE, col.names=NA) 

    ref <- as(a, "dgCMatrix")
    out <- readSparseCounts(ofile)
    expect_identical(ref, out)

    # Unaffected by chunk size.
    out2 <- readSparseCounts(ofile, chunk=23L)
    expect_identical(ref, out2)
    out2 <- readSparseCounts(ofile, chunk=51L)
    expect_identical(ref, out2)

    # Avoids row names if requested.
    out <- readSparseCounts(ofile, row.names=FALSE, ignore.col=1L)
    ref2 <- ref
    rownames(ref2) <- NULL
    expect_identical(ref2, out)

    expect_error(readSparseCounts(ofile, row.names=FALSE), "expected")

    # Avoids column names if requested.
    out <- readSparseCounts(ofile, col.names=FALSE, ignore.row=1L)
    ref2 <- ref
    colnames(ref2) <- NULL
    expect_identical(ref2, out)

    expect_error(readSparseCounts(ofile, col.names=FALSE), "invalid")

    # Skipping works correctly.
    out <- readSparseCounts(ofile, skip.row=10L)
    expect_identical(ref[11:nrow(ref),], out)

    out <- readSparseCounts(ofile, skip.col=10L)
    expect_identical(ref[,11:ncol(ref)], out)
})
