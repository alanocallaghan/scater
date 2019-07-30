## Test functions for QC calculation.
## library(scater); library(testthat); source("setup-sce.R"); source("test-qc-calc.R")

original <- sce

test_that("we can compute standard per-cell QC metrics", {
    df <- perCellQCMetrics(original, flatten=FALSE)
    expect_identical(rownames(df), colnames(original))

    # Testing total metrics for cells.
    expect_equal(df$sum, unname(colSums(counts(original))))
    expect_equal(df$detected, unname(colSums(counts(original) > 0)))

    # Testing percentage metrics for cells.
    for (i in seq_len(ncol(original))) { 
        cur_counts <- counts(original)[,i]
        o <- order(cur_counts, decreasing=TRUE)
        lib_size <- sum(cur_counts)

        for (x in c(50, 100, 200, 500)) { 
            chosen <- o[seq_len(x)]
            expect_equivalent(df$percent_top[i,as.character(x)], sum(cur_counts[chosen])/lib_size * 100) 
        }
    }
})

test_that("we can compute standard QC metrics with subsets", {
    ref <- perCellQCMetrics(original, flatten=FALSE)
    df <- perCellQCMetrics(original, subsets = list(set1 = 1:20), flatten=FALSE)
    expect_identical(df[,1:3], ref[,1:3])

    expect_equivalent(df$subsets$set1$sum, colSums(counts(original[1:20,])))
    expect_equivalent(df$subsets$set1$detected, colSums(counts(original[1:20,])> 0))
    expect_equivalent(df$subsets$set1$percent, df$subsets$set1$sum/df$sum * 100)

    # Testing behaviour with multiple feature controls.
    multi_controls <- list(controls1 = 1:20, controls2 = rownames(original)[500:1000])
    df2 <- perCellQCMetrics(original, subsets = multi_controls, flatten=FALSE)

    expect_equivalent(df$subsets$set1, df2$subsets$controls1)
    expect_equivalent(df2$subsets$controls2$sum, colSums(counts(original[500:1000,])))
    expect_equivalent(df2$subsets$controls2$detected, colSums(counts(original[500:1000,])> 0))
    expect_equivalent(df2$subsets$controls2$percent, df2$subsets$controls2$sum/df2$sum * 100)

    # Flattening works as expected.
    flat <- perCellQCMetrics(original, subsets = list(set1=1:20))
    expect_identical(flat$subsets_set1_sum, df$subsets$set1$sum)
    expect_identical(flat$percent_top_50, df$percent_top[,"50"])
})

test_that("perCellQCMetrics works with alternative experiments", {
    sce <- original
    altExp(sce, "alpha") <- original[1:10,]
    altExp(sce, "bravo") <- original[10:20,]

    ref <- perCellQCMetrics(original, flatten=FALSE)
    df <- perCellQCMetrics(sce, flatten=FALSE)
    expect_identical(df[,1:3], ref[,1:3])

    for (x in altExpNames(sce)) {
        current <- perCellQCMetrics(altExp(sce, x), flatten=FALSE)
        expect_identical(df$altexps[[x]]$sum, current$sum)
        expect_identical(df$altexps[[x]]$detected, current$detected)
        expect_equal(df$altexps[[x]]$percent, current$sum/df$total*100)
    }
    
    expect_identical(df$total, df$sum + df$altexps$alpha$sum + df$altexps$bravo$sum)

    # Flattening works as expected.
    flat <- perCellQCMetrics(sce)
    expect_identical(flat$altexps_alpha_sum, df$altexps$alpha$sum)
    expect_identical(flat$altexps_bravo_detected, df$altexps$bravo$detected)
})

test_that("perCellQCMetrics handles silly inputs", {
    expect_error(perCellQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    # Doesn't choke with no entries.
    thing <- perCellQCMetrics(original[0,], flatten=FALSE)
    expect_identical(rownames(thing), colnames(original))
    expect_true(all(thing$sum==0L))

    thing2 <- perCellQCMetrics(original[,0], flatten=FALSE)
    expect_identical(nrow(thing2), 0L)
    expect_identical(colnames(thing), colnames(thing2))

    # Percentage holds at the limit.
    df <- perCellQCMetrics(original[1:10,], flatten=FALSE)
    expect_true(all(df$percent_top==100))

    df <- perCellQCMetrics(original, percent_top=integer(0), flatten=FALSE)
    expect_identical(ncol(df$percent_top), 0L)

    # Responds to alternative inputs.
    blah <- sce
    assayNames(blah) <- "whee"
    expect_error(perCellQCMetrics(blah, flatten=FALSE), "counts")
    expect_error(perCellQCMetrics(blah, exprs_values="whee", flatten=FALSE), NA)
})

#######################################################################
# Works for per-feature metrics.

test_that("perFeatureQCMetrics works correctly", {
    out <- perFeatureQCMetrics(original, flatten=FALSE)
    expect_equal(out$mean, unname(rowMeans(counts(original))))
    expect_equal(out$detected, unname(rowMeans(counts(original) > 0))*100)
})

test_that("we can compute standard QC metrics with cell controls", {
    expect_error(perFeatureQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    df <- perFeatureQCMetrics(original, subsets = list(set1 = 1:20), flatten=FALSE)
    sub_counts <- counts(original)[,1:20]

    expect_equal(df$subsets$set1$mean, unname(rowMeans(sub_counts)))
    expect_equal(df$subsets$set1$detected, unname(rowMeans(sub_counts > 0) * 100))
    expect_equal(df$subsets$set1$ratio, df$subsets$set1$mean/df$mean)

    # Testing behaviour with multiple cell controls.
    multi_controls <- list(controls1 = 1:5, controls2 = 10:20)
    df2 <- perFeatureQCMetrics(original, subsets = multi_controls, flatten=FALSE)

    expect_equivalent(df2$subsets$controls2$mean, rowMeans(counts(original[,10:20])))
    expect_equivalent(df2$subsets$controls2$detected, rowMeans(counts(original[,10:20])> 0)*100)
    expect_equivalent(df2$subsets$controls2$ratio, df2$subsets$controls2$mean/df2$mean)

    # Flattening works as expected.
    flat <- perFeatureQCMetrics(original, subsets = list(set1 = 1:20))
    expect_identical(flat$subsets_set1_mean, df$subsets$set1$mean)
    expect_identical(flat$subsets_set1_ratio, df$subsets$set1$ratio)
})

test_that("perFeatureQCmetrics handles silly inputs", {
    expect_error(perFeatureQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    # Doesn't choke with no entries.
    thing <- perFeatureQCMetrics(original[,0], flatten=FALSE)
    expect_identical(rownames(thing), rownames(original))
    expect_true(all(thing$sum==0L))

    thing2 <- perFeatureQCMetrics(original[0,], flatten=FALSE)
    expect_identical(nrow(thing2), 0L)
    expect_identical(colnames(thing), colnames(thing2))

    # Responds to alternative inputs.
    blah <- sce
    assayNames(blah) <- "whee"
    expect_error(perFeatureQCMetrics(blah, flatten=FALSE), "counts")
    expect_error(perFeatureQCMetrics(blah, exprs_values="whee", flatten=FALSE), NA)
})

#######################################################################
# Responds to special settings: 

test_that("we can compute standard QC metrics on sparse counts matrix", {
    alt <- original

    library(Matrix)
    counts(alt) <- as(counts(alt), "dgCMatrix")

    expect_equal(perCellQCMetrics(alt, flatten=FALSE), perCellQCMetrics(original, flatten=FALSE))
    expect_equal(perFeatureQCMetrics(alt, flatten=FALSE), perFeatureQCMetrics(original, flatten=FALSE))

    expect_equal(perCellQCMetrics(alt, subset=list(set=1:10), flatten=FALSE), 
        perCellQCMetrics(original, subset=list(set=1:10), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(alt, subset=list(set=1:10), flatten=FALSE), 
        perFeatureQCMetrics(original, subset=list(set=1:10), flatten=FALSE))
})

test_that("we can compute standard QC metrics across multiple cores", {
    expect_equal(perCellQCMetrics(original, flatten=FALSE), 
        perCellQCMetrics(original, BPPARAM=safeBPParam(3), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(original, flatten=FALSE), 
        perFeatureQCMetrics(original, BPPARAM=safeBPParam(3), flatten=FALSE))

    expect_equal(perCellQCMetrics(original, subset=list(set=1:10), flatten=FALSE), 
        perCellQCMetrics(original, subset=list(set=1:10), BPPARAM=safeBPParam(3), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(original, subset=list(set=1:10), flatten=FALSE), 
        perFeatureQCMetrics(original, subset=list(set=1:10), BPPARAM=safeBPParam(3), flatten=FALSE))
})
