# tests for feature pre-processing functions.
# library(scater); library(testthat); source("setup-sce.R"); source("test-feat-proc.R")

context("test feature pre-processing functions")

##########################################################

set.seed(10001)
test_that("we can summarise counts at feature set level", {
    ids <- sample(nrow(sce)/2, nrow(sce), replace=TRUE)
    out <- sumCountsAcrossFeatures(sce, ids)
    expect_identical(out, rowsum(counts(sce), ids))
    expect_identical(rownames(out), as.character(sort(unique(ids))))

    out2 <- sumCountsAcrossFeatures(counts(sce), ids)
    expect_identical(out, out2)

    # assay.type= works correctly.
    alt <- sce
    assayNames(alt) <- "whee"
    out2 <- sumCountsAcrossFeatures(alt, ids, assay.type="whee")
    expect_identical(out, out2)

    # Respects levels properly.
    fids <- factor(ids, levels=rev(sort(unique(ids))))
    fout <- sumCountsAcrossFeatures(sce, fids)
    expect_identical(out, fout[nrow(fout):1,])

    # Handles NA's correctly.
    ids2 <- sample(LETTERS, nrow(sce), replace=TRUE)
    out2 <- sumCountsAcrossFeatures(sce, ids2)

    ids3 <- ids2
    ids3[ids3=="A"] <- NA
    out3 <- sumCountsAcrossFeatures(sce, ids3)
    expect_identical(out2[setdiff(rownames(out2), "A"),], out3)
})

set.seed(10002)
test_that("by-feature count summarization behaves with odd inputs", {
    ids <- sample(nrow(sce)/2, nrow(sce), replace=TRUE)
    ref <- sumCountsAcrossFeatures(sce, ids)

    # Handles sparse matrices properly.
    library(Matrix)
    sparsified <- sce
    counts(sparsified) <- as(counts(sparsified), "dgCMatrix")
    spack <- sumCountsAcrossFeatures(sparsified, ids)
    expect_s4_class(spack, "dgCMatrix")
    expect_equal(ref, as.matrix(spack))

    unknown <- sce
    counts(unknown) <- as(counts(unknown), "dgTMatrix")
    spack <- sumCountsAcrossFeatures(unknown, ids)
    expect_true(is.matrix(spack))
    expect_equivalent(ref, as.matrix(spack))

    # Handles parallelization properly.
    alt <- sumCountsAcrossFeatures(sce, ids, BPPARAM=safeBPParam(2))
    expect_identical(alt, ref)

    alt <- sumCountsAcrossFeatures(sce, ids, BPPARAM=safeBPParam(3))
    expect_identical(alt, ref)
})

set.seed(100021)
test_that("Aggregation across features works correctly", {
    ids <- paste0("GENE_", sample(nrow(sce)/2, nrow(sce), replace=TRUE))
    alt <- aggregateAcrossFeatures(sce, ids)

    expect_identical(rownames(alt), sort(unique(ids)))
    expect_identical(counts(alt), sumCountsAcrossFeatures(counts(sce), ids))

    # Behaves in the presence of multiple assays.
    normcounts(sce) <- normalizeCounts(sce, log=FALSE)
    alt2 <- aggregateAcrossFeatures(sce, ids)
    expect_identical(alt, alt2)

    alt3 <- aggregateAcrossFeatures(sce, ids, use.assay.types=c("counts", "normcounts"))
    expect_identical(counts(alt), counts(alt3))
    expect_identical(normcounts(alt3), sumCountsAcrossFeatures(sce, ids, assay.type="normcounts"))
})

##########################################################

set.seed(10003)
test_that("we can summarise counts at cell cluster level", {
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    out <- sumCountsAcrossCells(sce, ids)
    expect_identical(out, colsum(counts(sce), ids))

    out2 <- sumCountsAcrossCells(counts(sce), ids)
    expect_identical(out, out2)

    # assay.type= works correctly.
    alt <- sce
    assayNames(alt) <- "whee"
    out2 <- sumCountsAcrossCells(alt, ids, assay.type="whee")
    expect_identical(out, out2)

    # Respects levels properly.
    fids <- factor(ids, levels=rev(sort(unique(ids))))
    fout <- sumCountsAcrossCells(sce, fids)
    expect_identical(out, fout[,ncol(fout):1])

    # Handles NA's correctly.
    ids2 <- sample(LETTERS, ncol(sce), replace=TRUE)
    out2 <- sumCountsAcrossCells(sce, ids2)

    ids3 <- ids2
    ids3[ids3=="A"] <- NA
    out3 <- sumCountsAcrossCells(sce, ids3)

    expect_identical(out2[,setdiff(colnames(out2), "A")], out3)
})

set.seed(10004)
test_that("by-cell count summarization behaves with odd inputs", {
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    ref <- sumCountsAcrossCells(sce, ids)

    # Handles sparse matrices properly.
    library(Matrix)
    sparsified <- sce
    counts(sparsified) <- as(counts(sparsified), "dgCMatrix")
    spack <- sumCountsAcrossCells(sparsified, ids)
    expect_s4_class(spack, "dgCMatrix")
    expect_equal(ref, as.matrix(spack))

    unknown <- sce
    counts(unknown) <- as(counts(unknown), "dgTMatrix")
    spack <- sumCountsAcrossCells(unknown, ids)
    expect_true(is.matrix(spack))
    expect_equivalent(ref, as.matrix(spack))

    # Handles parallelization properly.
    alt <- sumCountsAcrossCells(sce, ids, BPPARAM=safeBPParam(2))
    expect_identical(alt, ref)

    alt <- sumCountsAcrossCells(sce, ids, BPPARAM=safeBPParam(3))
    expect_identical(alt, ref)
})

set.seed(100041)
test_that("Aggregation across cells works correctly for SCEs", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    alt <- aggregateAcrossCells(sce, ids)

    expect_identical(colnames(alt), sort(unique(ids)))
    expect_identical(counts(alt), sumCountsAcrossCells(counts(sce), ids))

    # Behaves in the presence of multiple assays.
    normcounts(sce) <- normalizeCounts(sce, log=FALSE)
    alt2 <- aggregateAcrossCells(sce, ids)
    expect_identical(alt, alt2)

    alt3 <- aggregateAcrossCells(sce, ids, use.assay.types=c("counts", "normcounts"))
    expect_identical(counts(alt), counts(alt3))
    expect_identical(normcounts(alt3), sumCountsAcrossCells(sce, ids, assay.type="normcounts"))

    # Behaves for alternative experiments.
    copy <- sce
    altExp(copy, 1) <- sce
    copy <- aggregateAcrossCells(copy, ids)
    expect_identical(altExp(copy, 1), alt)
})

set.seed(100042)
test_that("Aggregation across cells works correctly for SEs", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    alt <- aggregateAcrossCells(sce, ids)
    expect_identical(colnames(alt), sort(unique(ids)))
    expect_identical(counts(alt), sumCountsAcrossCells(counts(sce), ids))
})

##########################################################

test_that("we can uniquify the feature names", {
    all.genes <- sample(c(LETTERS, LETTERS[1:5], NA, NA))
    all.ids <- paste0("GENE", seq_along(all.genes))
    out <- uniquifyFeatureNames(all.ids, all.genes)

    lost <- is.na(all.genes)
    expect_identical(out[lost], all.ids[lost])
    dup <- all.genes %in% all.genes[duplicated(all.genes)]
    expect_identical(out[!dup & !lost], all.genes[!dup & !lost])
    expect_identical(out[dup & !lost], paste0(all.genes, "_", all.ids)[dup & !lost])
})

