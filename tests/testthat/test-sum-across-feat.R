# tests for feature pre-processing functions.
# library(scater); library(testthat); source("setup-sce.R"); source("test-sum-across-feat.R")

library(Matrix)
library(DelayedArray)

##########################################################

set.seed(10001)
test_that("we can summarise counts at feature set level", {
    ids <- sample(nrow(sce)/2, nrow(sce), replace=TRUE)
    out <- sumCountsAcrossFeatures(sce, ids)
    expect_identical(out, rowsum(counts(sce), ids))
    expect_identical(rownames(out), as.character(sort(unique(ids))))

    out2 <- sumCountsAcrossFeatures(counts(sce), ids)
    expect_identical(out, out2)

    # Handles averaging correctly.
    out2 <- sumCountsAcrossFeatures(sce, ids, average=TRUE)
    expect_identical(out2, rowsum(counts(sce), ids)/as.integer(table(ids)))

    # exprs_values= works correctly.
    alt <- sce
    assayNames(alt) <- "whee"
    out2 <- sumCountsAcrossFeatures(alt, ids, exprs_values="whee")
    expect_identical(out, out2)

    # Respects levels properly.
    fids <- factor(ids, levels=rev(sort(unique(ids))))
    fout <- sumCountsAcrossFeatures(sce, fids)
    expect_identical(out, fout[nrow(fout):1,])
})

set.seed(10001)
test_that("count summarization at feature set level respects NAs", {
    ids2 <- sample(LETTERS, nrow(sce), replace=TRUE)
    out2 <- sumCountsAcrossFeatures(sce, ids2)

    ids3 <- ids2
    ids3[ids3=="A"] <- NA
    out3 <- sumCountsAcrossFeatures(sce, ids3)
    expect_identical(out2[setdiff(rownames(out2), "A"),], out3)

    ids4 <- ids2
    ids4[1:10] <- NA
    out4a <- sumCountsAcrossFeatures(sce, ids4)
    out4b <- sumCountsAcrossFeatures(sce[-(1:10),], ids4[-(1:10)])
    expect_identical(out4a, out4b)
})

set.seed(10002)
test_that("by-feature count summarization behaves with lists", {
    idl <- list(10:1, sample(nrow(sce), 100), nrow(sce) - 1:10)
    outl <- sumCountsAcrossFeatures(sce, idl)

    manual <- list()
    for (i in seq_along(idl)) {
        manual[[i]] <- colSums(counts(sce)[idl[[i]],])
    }
    expect_identical(outl, do.call(rbind, manual))

    expect_identical(outl, sumCountsAcrossFeatures(sce, lapply(idl, function(i) rownames(sce)[i])))
    
    expect_identical(outl/lengths(idl), sumCountsAcrossFeatures(sce, idl, average=TRUE))
})

set.seed(100021)
test_that("by-feature count summarization responds to subsetting", {
    ids <- sample(LETTERS, nrow(sce), replace=TRUE)

    keep <- rbinom(nrow(sce), 1, 0.5)>0
    ref <- sumCountsAcrossFeatures(sce[keep,], ids[keep])
    out <- sumCountsAcrossFeatures(sce, ids, subset_row=keep)
    expect_identical(out, ref)

    keep2 <- rbinom(ncol(sce), 1, 0.5)>0
    ref <- sumCountsAcrossFeatures(sce[,keep2], ids)
    out <- sumCountsAcrossFeatures(sce, ids, subset_col=keep2)
    expect_identical(out, ref)
})

##########################################################

set.seed(10003)
test_that("by-feature count summarization behaves with different classes", {
    ids <- sample(nrow(sce)/2, nrow(sce), replace=TRUE)
    ref <- sumCountsAcrossFeatures(sce, ids)

    # Handles sparse matrices properly.
    library(Matrix)
    sparsified <- sce
    counts(sparsified) <- as(counts(sparsified), "dgCMatrix")
    spack <- sumCountsAcrossFeatures(sparsified, ids)
    expect_equal(ref, as.matrix(spack))

    unknown <- sce
    counts(unknown) <- as(counts(unknown), "dgTMatrix")
    spack <- sumCountsAcrossFeatures(unknown, ids)
    expect_equivalent(ref, as.matrix(spack))

    # Handles DelayedArrays properly.
    delayed <- sce
    counts(delayed) <- DelayedArray(counts(delayed))
    dack <- sumCountsAcrossFeatures(delayed, ids)
    expect_equivalent(ref, as.matrix(dack))
})

set.seed(100031)
test_that("by-feature count summarization parallelizes properly", {
    ids <- sample(nrow(sce)/2, nrow(sce), replace=TRUE)
    ref <- sumCountsAcrossFeatures(sce, ids)

    # Handles parallelization properly.
    alt <- sumCountsAcrossFeatures(sce, ids, BPPARAM=safeBPParam(2))
    expect_identical(alt, ref)

    alt <- sumCountsAcrossFeatures(sce, ids, BPPARAM=safeBPParam(3))
    expect_identical(alt, ref)
})

##########################################################

set.seed(10004)
test_that("Aggregation across features works correctly", {
    ids <- paste0("GENE_", sample(nrow(sce)/2, nrow(sce), replace=TRUE))
    alt <- aggregateAcrossFeatures(sce, ids)

    expect_identical(rownames(alt), sort(unique(ids)))
    expect_identical(counts(alt), sumCountsAcrossFeatures(counts(sce), ids))

    # Behaves in the presence of multiple assays.
    normcounts(sce) <- normalizeCounts(sce, log=FALSE)
    alt2 <- aggregateAcrossFeatures(sce, ids)
    expect_identical(alt, alt2)

    alt3 <- aggregateAcrossFeatures(sce, ids, use_exprs_values=c("counts", "normcounts"))
    expect_identical(counts(alt), counts(alt3))
    expect_identical(normcounts(alt3), sumCountsAcrossFeatures(sce, ids, exprs_values="normcounts"))
})

############################################

test_that("numDetectedAcrossFeatures works as expected", {
    ids <- sample(LETTERS[1:5], nrow(sce), replace=TRUE)

    expect_equal(numDetectedAcrossFeatures(counts(sce), ids),
        rowsum((counts(sce) > 0)+0, ids)) 
    expect_identical(numDetectedAcrossFeatures(counts(sce), ids, average=TRUE),
        rowsum((counts(sce) > 0)+0, ids)/as.integer(table(ids)))

    # Checking that it works direclty with SCEs.
    expect_equal(numDetectedAcrossFeatures(counts(sce), ids),
        numDetectedAcrossFeatures(sce, ids))
    expect_equal(numDetectedAcrossFeatures(counts(sce), ids, average=TRUE),
        numDetectedAcrossFeatures(sce, ids, average=TRUE))

    # Checking that subsetting works.
    expect_identical(numDetectedAcrossFeatures(counts(sce), ids, subset_col=10:1),
        numDetectedAcrossFeatures(counts(sce), ids)[,10:1])

    expect_identical(numDetectedAcrossFeatures(counts(sce), ids, subset_row=2:15),
        numDetectedAcrossFeatures(counts(sce)[2:15,], ids[2:15]))

    ids[c(1,3,5,6)] <- NA
    expect_identical(numDetectedAcrossFeatures(counts(sce), ids),
        numDetectedAcrossFeatures(counts(sce)[!is.na(ids),], ids[!is.na(ids)]))

    # Comparing to sumCountsAcrossFeatures.
    expect_equal(numDetectedAcrossFeatures(counts(sce), ids),
        sumCountsAcrossFeatures((counts(sce) > 0)+0, ids))
    expect_equal(numDetectedAcrossFeatures(counts(sce), ids, average=TRUE),
        sumCountsAcrossFeatures((counts(sce) > 0)+0, ids, average=TRUE))
})

test_that("numDetectedAcrossFeatures handles other matrix classes", {
    thing <- matrix(rpois(2000, lambda=0.5), ncol=100, nrow=20)
    ids <- sample(LETTERS[1:6], nrow(thing), replace=TRUE)

    ref <- numDetectedAcrossFeatures(thing, ids)
    expect_equal(rowSums(ref), rowSums(thing > 0))

    sparse <- as(thing, 'dgCMatrix')
    expect_equal(numDetectedAcrossFeatures(sparse, ids), ref)

    delayed <- DelayedArray(thing)
    expect_equal(numDetectedAcrossFeatures(delayed, ids), ref)
})
