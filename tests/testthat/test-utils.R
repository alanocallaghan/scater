# This tests the various internal utilities.
# library(testthat); library(scater); source("setup.R"); source("test-utils.R")

test_that("converting subset vectors to indices is correct", {
    M <- matrix(rnorm(1000), ncol=20)
    colnames(M) <- head(LETTERS, ncol(M))
    rownames(M) <- sprintf("GENE_%s", seq_len(nrow(M)))

    # For matrix rows:
    expect_identical(.subset2index(NULL, M, byrow=TRUE), seq_len(nrow(M)))
    expect_identical(.subset2index(1:5, M, byrow=TRUE), 1:5)
    expect_identical(.subset2index(rownames(M)[10:5], M, byrow=TRUE), 10:5)
    expect_identical(.subset2index(seq_len(nrow(M)) %in% 20:50, M, byrow=TRUE), 20:50)
    expect_error(.subset2index(1000:2000, M, byrow=TRUE), "invalid")

    # For matrix columns:
    expect_identical(.subset2index(NULL, M, byrow=FALSE), seq_len(ncol(M)))
    expect_identical(.subset2index(1:5, M, byrow=FALSE), 1:5)
    expect_identical(.subset2index(colnames(M)[10:5], M, byrow=FALSE), 10:5)
    expect_identical(.subset2index(seq_len(ncol(M)) %in% 10:15, M, byrow=FALSE), 10:15)
    expect_error(.subset2index(1000:2000, M, byrow=FALSE), "invalid")

    # For vectors:
    V <- M[,1]
    expect_identical(.subset2index(NULL, V, byrow=NA), seq_along(V))
    expect_identical(.subset2index(1:5, V, byrow=NA), 1:5)
    expect_identical(.subset2index(names(V)[10:5], V, byrow=NA), 10:5)
    expect_identical(.subset2index(seq_along(V) %in% 10:15, V, byrow=NA), 10:15)
    expect_error(.subset2index(1000:2000, V, byrow=NA), "invalid")

    # Edge cases.
    expect_identical(.subset2index(numeric(0), M, byrow=FALSE), integer(0))
    expect_identical(.subset2index(character(0), M, byrow=FALSE), integer(0))
    expect_identical(.subset2index(logical(0), M, byrow=FALSE), integer(0))
})

test_that("job assignment to workers is correct", {
    out <- .assignIndicesToWorkers(100, safeBPParam(3))
    expect_identical(length(out), 3L)
    expect_true(all(lengths(out) >= floor(100/3)))
    expect_identical(unlist(out), seq_len(100))

    out <- .assignIndicesToWorkers(100, safeBPParam(6))
    expect_identical(length(out), 6L)
    expect_true(all(lengths(out) >= floor(100/6)))
    expect_identical(unlist(out), seq_len(100))

    out <- .assignIndicesToWorkers(100, safeBPParam(11))
    expect_identical(length(out), 11L)
    expect_true(all(lengths(out) >= floor(100/11)))
    expect_identical(unlist(out), seq_len(100))

    # Works with subsetting.
    chosen <- sample(100, 50)
    out <- .assignIndicesToWorkers(NULL, safeBPParam(11), subset=chosen)
    expect_identical(length(out), 11L)
    expect_true(all(lengths(out) >= floor(length(chosen)/11)))
    expect_identical(unlist(out), chosen)

    chosen <- sample(LETTERS)
    out <- .assignIndicesToWorkers(NULL, safeBPParam(5), subset=chosen)
    expect_identical(length(out), 5L)
    expect_true(all(lengths(out) >= floor(length(chosen)/5)))
    expect_identical(unlist(out), chosen)

    chosen <- rbinom(25, 1, 0.5)==1
    out <- .assignIndicesToWorkers(NULL, safeBPParam(3), subset=chosen)
    expect_identical(length(out), 3L)
    expect_true(all(lengths(out) >= floor(sum(chosen)/3)))
    expect_identical(unlist(out), which(chosen))
})


test_that("splitting a vector to workers is correct", {
    X <- runif(99)
    out <- .splitVectorByWorkers(X, safeBPParam(7))
    expect_identical(unlist(out), X)
    expect_identical(length(out), 7L)

    out <- .splitVectorByWorkers(X, safeBPParam(1))
    expect_identical(out[[1]], X)

    # Behaves with subsetting of all flavors.
    i <- 1:50
    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(7), subset=i),
        .splitVectorByWorkers(X[i], safeBPParam(7))
    )

    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(1), subset=i),
        .splitVectorByWorkers(X[i], safeBPParam(1))
    )

    j <- rbinom(99, 1, 0.2)==1
    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(1), subset=j),
        .splitVectorByWorkers(X[j], safeBPParam(1))
    )
})

test_that("splitting a matrix by row is correct", {
    M <- matrix(rnorm(1000), ncol=10)

    out <- .splitRowsByWorkers(M, safeBPParam(7))
    expect_identical(do.call(rbind, out), M)
    expect_identical(length(out), 7L)

    out <- .splitRowsByWorkers(M, safeBPParam(1))
    expect_identical(out[[1]], M)
    expect_identical(length(out), 1L)

    # Behaves with subsetting by row:
    i <- 1:50
    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(7), subset_row=i),
        .splitRowsByWorkers(M[i,], safeBPParam(7))
    )

    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(1), subset_row=i),
        .splitRowsByWorkers(M[i,], safeBPParam(1))
    )

    # Behaves with column subsetting:
    i <- 5:1
    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(7), subset_col=i),
        .splitRowsByWorkers(M[,i], safeBPParam(7))
    )

    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(1), subset_col=i),
        .splitRowsByWorkers(M[,i], safeBPParam(1))
    )
})

test_that("splitting a matrix by column is correct", {
    M <- matrix(rnorm(1000), nrow=20)

    out <- .splitColsByWorkers(M, safeBPParam(8))
    expect_identical(do.call(cbind, out), M)
    expect_identical(length(out), 8L)

    out <- .splitColsByWorkers(M, safeBPParam(1))
    expect_identical(out[[1]], M)
    expect_identical(length(out), 1L)

    # Behaves with subsetting by row:
    i <- 1:10
    expect_identical(
        .splitColsByWorkers(M, safeBPParam(4), subset_row=i),
        .splitColsByWorkers(M[i,], safeBPParam(4))
    )

    expect_identical(
        .splitColsByWorkers(M, safeBPParam(1), subset_row=i),
        .splitColsByWorkers(M[i,], safeBPParam(1))
    )

    # Behaves with column subsetting:
    i <- sample(ncol(M), 20)
    expect_identical(
        .splitColsByWorkers(M, safeBPParam(4), subset_col=i),
        .splitColsByWorkers(M[,i], safeBPParam(4))
    )

    expect_identical(
        .splitColsByWorkers(M, safeBPParam(1), subset_col=i),
        .splitColsByWorkers(M[,i], safeBPParam(1))
    )
})
