# Tests for tximport wrapper

context("test input usage")

test_that("the system works", {
    
    expect_that(read10xResults(system.file("extdata", package="scater")),
                is_a("SingleCellExperiment"))
})

test_that("downsampling is correct", {
    set.seed(0)
    ncells <- 100
    u <- matrix(rpois(20000, 5), ncol=ncells)

    set.seed(100)
    out <- scater:::downsampleCounts(u, 0.1)
    set.seed(100)
    ref <- rbinom(length(u), u, 0.1)
    dim(ref) <- dim(out)
    expect_identical(out, ref)

    # Checking double-precision values.
    v <- u
    storage.mode(v) <- "double"

    set.seed(200)
    out <- scater:::downsampleCounts(v, 0.1)
    set.seed(200)
    ref <- as.numeric(rbinom(length(u), u, 0.1))
    dim(ref) <- dim(out)
    expect_equal(out, ref)

    # Checking vectors of proportions.
    props <- runif(ncells)

    set.seed(300)
    out <- scater:::downsampleCounts(u, props)
    set.seed(300)
    ref <- rbinom(length(u), u, matrix(props, nrow(u), ncol(u), byrow=TRUE))
    dim(ref) <- dim(out)
    expect_identical(out, ref)

    # Checking sparse matrix inputs.
    w <- as(v, "dgCMatrix")
    
    set.seed(300)
    out <- scater:::downsampleCounts(w, 0.2)
    set.seed(300)
    ref <- rbinom(length(u), u, 0.2)
    dim(ref) <- dim(out)

    coerced <- as.matrix(out)
    dimnames(coerced) <- NULL    
    expect_equal(coerced, ref)
})

