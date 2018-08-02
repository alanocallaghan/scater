## Test functions related to calculation of the variance explained.
## library(scater); library(testthat); source("setup-sce.R"); source("test-expvar.R")

#############################################################
# getVarianceExplained() tests:

varexp <- getVarianceExplained(normed) 
  
test_that("getVarianceExplained matches with a reference function", {
    expect_identical(colnames(varexp), colnames(colData(normed)))
    for (v in colnames(varexp)) {
        X <- colData(normed)[[v]]
        Y <- logcounts(normed, withDimnames=FALSE)

        # Using lm() per gene, which is a bit slow but ensures we get a reference R2.
        output <- numeric(nrow(normed))
        for (g in seq_len(nrow(normed))) { 
            y <- Y[g,]
            fit <- lm(y ~ X)
            output[g] <- suppressWarnings(summary(fit)$r.squared)
        }

        expect_equal(output, unname(varexp[,v]))
    }
})

test_that("getVarianceExplained responds to the options", {
    # Responds to differences in the assay name.
    blah <- normed
    assayNames(blah) <- c("yay", "whee")
    expect_error(getVarianceExplained(blah), "logcounts")
    expect_identical(varexp, getVarianceExplained(blah, exprs_values="whee"))
    
    # Responds to choice of variable.
    expect_identical(varexp[,1,drop=FALSE], getVarianceExplained(normed, variables=colnames(varexp)[1]))
    expect_identical(varexp[,c(3,2),drop=FALSE], getVarianceExplained(normed, variables=colnames(varexp)[c(3,2)]))
    expect_identical(varexp, getVarianceExplained(normed, variables=colnames(varexp)))

    # Unaffected by chunk size.
    expect_identical(varexp, getVarianceExplained(normed, chunk=10))
    expect_identical(varexp, getVarianceExplained(normed, chunk=1))
    expect_identical(varexp, getVarianceExplained(normed, chunk=100000))
})

test_that("getVarianceExplained handles sparse inputs", {
    varexp_sparse <- getVarianceExplained(normed_sparse)
    expect_equal(varexp, varexp_sparse)
})

test_that("getVarianceExplained handles silly inputs correctly", {
    # Misspecified variables.
    expect_error(getVarianceExplained(normed, variables="yay"), "invalid names")

    # Empty inputs.
    out <- getVarianceExplained(normed[0,])
    expect_identical(colnames(out), colnames(varexp))
    expect_identical(nrow(out), 0L)

    expect_warning(out <- getVarianceExplained(normed[,0]), "2 unique levels")
    expect_identical(dimnames(out), dimnames(varexp))
    expect_true(all(is.na(out)))
    
    # Inputs with only one unique level.
    normed$whee <- 1
    expect_warning(out2 <- getVarianceExplained(normed, variables="whee"), "2 unique levels")
    expect_identical(rownames(out2), rownames(varexp))
    expect_identical(colnames(out2), "whee")
    expect_true(all(is.na(out)))
})

#############################################################
# plotExplanatoryVariables() tests:

test_that("plotExplanatoryVariables works as expected", {
    out <- plotExplanatoryVariables(normed)
    ref <- plotExplanatoryVariables(varexp)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    medians <- apply(varexp, 2, median, na.rm=TRUE)

    # Responds to choice of number of variables
    out <- plotExplanatoryVariables(normed, nvars_to_plot=2)
    ref <- plotExplanatoryVariables(varexp[,order(medians, decreasing=TRUE)[1:2]])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, nvars_to_plot=Inf)
    ref <- plotExplanatoryVariables(varexp)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, nvars_to_plot=0)
    ref <- plotExplanatoryVariables(varexp[,0,drop=FALSE])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    # Responds to choice of minimum marginal R2.
    out <- plotExplanatoryVariables(normed, min_marginal_r2=0.5)
    ref <- plotExplanatoryVariables(varexp[,medians >= 0.5,drop=FALSE])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, min_marginal_r2=0.05)
    ref <- plotExplanatoryVariables(varexp[,medians >= 0.05,drop=FALSE])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    # Handles silly inputs.
    expect_s3_class(plotExplanatoryVariables(varexp[0,,drop=FALSE]), "ggplot")
    expect_s3_class(plotExplanatoryVariables(varexp[,0,drop=FALSE]), "ggplot")
})

#############################################################
# getExplanatoryPCs() tests:

exppcs <- getExplanatoryPCs(normed, ncomponents=10)

test_that("getExplanatoryPCs matches with a reference function", {
    expect_identical(colnames(exppcs), colnames(colData(normed)))
    normed <- runPCA(normed, ncomponents=nrow(exppcs))

    for (v in colnames(varexp)) {
        X <- colData(normed)[[v]]
        Y <- reducedDim(normed, "PCA")

        # Using lm() per gene, which is a bit slow but ensures we get a reference R2.
        output <- numeric(ncol(Y))
        for (g in seq_along(output)) {
            y <- Y[,g]
            fit <- lm(y ~ X)
            output[g] <- suppressWarnings(summary(fit)$r.squared)
        }

        expect_equal(output, unname(exppcs[,v]))
    }
})

test_that("getVarianceExplained responds to the options", {
    # Responds to differences in the reduced dimension slot.
    blah <- normed
    normed <- runPCA(normed, ncomponents=10)
    reducedDim(blah, "WHEE") <- reducedDim(normed, "PCA")
    expect_identical(getExplanatoryPCs(normed), getExplanatoryPCs(blah, use_dimred="WHEE"))

    reducedDim(blah, "WHEE") <- reducedDim(normed, "PCA")[,1:2]
    expect_identical(getExplanatoryPCs(normed)[1:2,], getExplanatoryPCs(blah, use_dimred="WHEE"))
    reducedDims(normed) <- List()
    expect_identical(getExplanatoryPCs(normed, ncomponents=2), getExplanatoryPCs(blah, use_dimred="WHEE"))

    # Responds to choice of variable.
    expect_identical(exppcs[1:2,1,drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)[1]))
    expect_identical(exppcs[1:2,c(3,2),drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)[c(3,2)]))
    expect_identical(exppcs[1:2,,drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)))

    # Unaffected by chunk size.
    expect_identical(exppcs, getExplanatoryPCs(normed, ncomponents=nrow(exppcs), chunk=10))
    expect_identical(exppcs, getExplanatoryPCs(normed, ncomponents=nrow(exppcs), chunk=1))
    expect_identical(exppcs, getExplanatoryPCs(normed, ncomponents=nrow(exppcs), chunk=100000))
})

#############################################################
# plotExplanatoryPCs() tests:

test_that("plotExplanatoryVariables works as expected", {
    out <- plotExplanatoryPCs(normed, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    # Handles situations where many more PCs are requested.
    out <- plotExplanatoryPCs(normed, npcs=5)
    allpcs <- runPCA(normed, ncomponents=5)
    ref <- plotExplanatoryPCs(allpcs)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryPCs(normed, npcs=Inf)
    allpcs <- runPCA(normed, ncomponents=Inf)
    ref <- plotExplanatoryPCs(allpcs)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    # Responds to choice of number of variables
    maxes <- apply(exppcs, 2, max, na.rm=TRUE)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=2, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs[,order(maxes, decreasing=TRUE)[1:2]])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=Inf, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs)
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=0, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs[,0,drop=FALSE])
    expect_s3_class(out, "ggplot")
    expect_identical(out$data, ref$data)

    # Handles silly inputs.
    expect_s3_class(suppressWarnings(plotExplanatoryPCs(exppcs[0,,drop=FALSE])), "ggplot")
    expect_s3_class(plotExplanatoryPCs(exppcs[,0,drop=FALSE]), "ggplot")
})
