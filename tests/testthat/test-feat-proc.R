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

    alt3 <- aggregateAcrossFeatures(sce, ids, use_exprs_values=c("counts", "normcounts"))
    expect_identical(counts(alt), counts(alt3))
    expect_identical(normcounts(alt3), sumCountsAcrossFeatures(sce, ids, exprs_values="normcounts"))
})

##########################################################

test_that("we can uniquify the feature names", {
    all.genes <- sample(c(LETTERS, LETTERS[1:5], NA, NA))
    all.ids <- paste0("GENE", seq_along(all.genes))
    out <- uniquifyFeatureNames(all.ids, all.genes)
    out.factor.names <- uniquifyFeatureNames(all.ids, factor(all.genes))
    out.factor.id <- uniquifyFeatureNames(factor(all.ids), all.genes)
    
    lost <- is.na(all.genes)
    expect_identical(out[lost], all.ids[lost])
    dup <- all.genes %in% all.genes[duplicated(all.genes)]
    expect_identical(out[!dup & !lost], all.genes[!dup & !lost])
    expect_identical(out[dup & !lost], paste0(all.genes, "_", all.ids)[dup & !lost])
    
    expect_identical(out, out.factor.names)
    expect_identical(out, out.factor.id)
})

