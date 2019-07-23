# Tests for various data input methods
# library(scater); library(testthat); source("test-load-data.R")

context("test expected usage")

a <- matrix(rpois(10000, lambda=1), ncol=50)
rownames(a) <- paste0("Gene", seq_len(nrow(a)))
colnames(a) <- paste0("Cell", seq_len(ncol(a)))

ofile <- tempfile()
write.table(a, file=ofile, sep="\t", quote=FALSE, col.names=NA) 
ref <- as(a, "dgCMatrix")

library(Matrix)
test_that("readSparseCounts works as expected in basic cases", {
    out <- readSparseCounts(ofile)
    expect_identical(ref, out)

    # Unaffected by chunk size.
    out2 <- readSparseCounts(ofile, chunk=23L)
    expect_identical(ref, out2)
    out2 <- readSparseCounts(ofile, chunk=51L)
    expect_identical(ref, out2)
})

test_that("readSparseCounts avoids row/column names if requested", {
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
})

test_that("Skipping works correctly", {
    out <- readSparseCounts(ofile, skip.row=10L)
    expect_identical(ref[11:nrow(ref),], out)

    out <- readSparseCounts(ofile, skip.col=10L)
    expect_identical(ref[,11:ncol(ref)], out)
})

test_that("Behaves correctly with wonky quotes and comments in the row/column names", {
    a2 <- a
    rownames(a2) <- paste0('"', rownames(a), "#")
    colnames(a2) <- paste0('#"', colnames(a), "'")
    ref2 <- as(a2, "dgCMatrix")

    ofile2 <- tempfile()
    write.table(a2, file=ofile2, sep="\t", quote=FALSE, col.names=NA) 
    out2 <- readSparseCounts(ofile2)
    expect_identical(ref2, out2)
})

test_that("Behaves properly with file handle input", {
    ofile3 <- tempfile(fileext=".gz")
    XHANDLE <- gzfile(ofile3, open='wb')
    write.table(a, file=XHANDLE, sep="\t", quote=FALSE, col.names=NA) 
    close(XHANDLE)

    fhandle <- file(ofile3)
    expect_error(readSparseCounts(fhandle), "read mode")
    close(fhandle)
    expect_error(readSparseCounts(DataFrame(ofile3)), "connection")
    
    outX <- readSparseCounts(ofile3)
    fhandle <- file(ofile3, open='rt')
    outY <- readSparseCounts(fhandle)
    close(fhandle)
    expect_identical(outX, ref)
    expect_identical(outY, ref)
})
