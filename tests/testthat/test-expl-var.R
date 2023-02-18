## Test functions related to calculation of the variance explained.
## library(scater); library(testthat); source("setup.R"); source("test-expl-var.R")

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

        expect_equal(output * 100, unname(varexp[,v]))
    }
})

test_that("getVarianceExplained responds to the options", {
    # Responds to differences in the assay name.
    blah <- normed
    assayNames(blah) <- c("yay", "whee")
    expect_error(getVarianceExplained(blah), "logcounts")
    expect_identical(varexp, getVarianceExplained(blah, assay.type="whee"))
    
    # Responds to choice of variable.
    expect_identical(varexp[,1,drop=FALSE], getVarianceExplained(normed, variables=colnames(varexp)[1]))
    expect_identical(varexp, getVarianceExplained(normed, variables=colnames(varexp)))
    expect_identical(varexp, getVarianceExplained(logcounts(normed), variables=colData(normed)))

    # Responds to subsetting.
    expect_identical(varexp[1:10,], getVarianceExplained(normed, subset_row=1:10))
    expect_identical(varexp[2:20,], getVarianceExplained(normed, subset_row=rownames(normed)[2:20]))
})

test_that("getVarianceExplained handles sparse inputs", {
    normed_sparse <- normed
    library(Matrix)
    counts(normed_sparse) <- as(counts(normed), "dgCMatrix")
    logcounts(normed_sparse) <- as(logcounts(normed), "dgCMatrix")

    varexp_sparse <- getVarianceExplained(normed_sparse)
    expect_equal(varexp, varexp_sparse)
})

test_that("getVarianceExplained handles NA values in the metadata", {
    whee <- runif(ncol(normed))
    whee[sample(ncol(normed), ncol(normed)/3)] <- NA
    normed$whee <- whee

    out <- getVarianceExplained(normed, variables="whee")
    ref <- getVarianceExplained(normed[,!is.na(whee)], variables="whee")
    expect_identical(out, ref)
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
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    medians <- apply(varexp, 2, median, na.rm=TRUE)

    # Responds to choice of number of variables
    out <- plotExplanatoryVariables(normed, nvars_to_plot=2)
    ref <- plotExplanatoryVariables(varexp[,order(medians, decreasing=TRUE)[1:2]])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, nvars_to_plot=Inf)
    ref <- plotExplanatoryVariables(varexp)
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, nvars_to_plot=0)
    ref <- plotExplanatoryVariables(varexp[,0,drop=FALSE])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    # Responds to choice of minimum marginal R2.
    out <- plotExplanatoryVariables(normed, min_marginal_r2=0.5)
    ref <- plotExplanatoryVariables(varexp[,medians >= 0.5,drop=FALSE])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryVariables(normed, min_marginal_r2=0.05)
    ref <- plotExplanatoryVariables(varexp[,medians >= 0.05,drop=FALSE])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    # Handles silly inputs.
    expect_ggplot(plotExplanatoryVariables(varexp[0,,drop=FALSE]))
    expect_ggplot(plotExplanatoryVariables(varexp[,0,drop=FALSE]))
})

#############################################################
# getExplanatoryPCs() tests:

normed <- runPCA(normed, ncomponents=20)
exppcs <- getExplanatoryPCs(normed, n_dimred=10)

test_that("getExplanatoryPCs matches with a reference function", {
    expect_identical(nrow(exppcs), 10L)
    expect_identical(colnames(exppcs), colnames(colData(normed)))

    for (v in colnames(varexp)) {
        X <- colData(normed)[[v]]
        Y <- reducedDim(normed, "PCA")[,seq_len(nrow(exppcs))]

        # Using lm() per PC, which is a bit slow but ensures we get a reference R2.
        output <- numeric(ncol(Y))
        for (g in seq_along(output)) {
            y <- Y[,g]
            fit <- lm(y ~ X)
            output[g] <- suppressWarnings(summary(fit)$r.squared)
        }

        expect_equal(output * 100, unname(exppcs[,v]))
    }
})

test_that("getExplanatoryPCs responds to PC-specific options", {
    # Responds to differences in the reduced dimension slot.
    blah <- normed
    normed2 <- runPCA(normed, ncomponents=10)
    reducedDim(blah, "WHEE") <- reducedDim(normed2, "PCA")
    expect_identical(res <- getExplanatoryPCs(normed2), getExplanatoryPCs(blah, dimred="WHEE"))
    expect_identical(nrow(res), 10L)

    reducedDim(blah, "WHEE") <- reducedDim(normed2, "PCA")[,1:2]
    expect_identical(getExplanatoryPCs(normed2)[1:2,], getExplanatoryPCs(blah, dimred="WHEE"))
    
    # Correctly truncates existing PCs.
    expect_identical(res <- getExplanatoryPCs(normed2, n_dimred=2), getExplanatoryPCs(blah, dimred="WHEE"))
    expect_identical(nrow(res), 2L)

    expect_identical(getExplanatoryPCs(normed2, n_dimred=Inf), getExplanatoryPCs(normed)) # ignores Inf, as it's maxed out.
})

test_that("getExplanatoryPCs responds to getVarianceExplained options", {
    # Responds to choice of variable.
    expect_identical(exppcs[,1,drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)[1]))
    expect_identical(exppcs[,c(3,2),drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)[c(3,2)]))
    expect_identical(exppcs[,,drop=FALSE], getExplanatoryPCs(normed, variables=colnames(varexp)))
})

#############################################################
# plotExplanatoryPCs() tests:

test_that("plotExplanatoryPCs works with PC choice options", {
    out <- plotExplanatoryPCs(normed, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs)
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    # Handles situations where different numbers of PCs are requested.
    out <- plotExplanatoryPCs(normed, npcs=5)
    allpcs <- runPCA(normed, ncomponents=5)
    ref <- plotExplanatoryPCs(allpcs)
    expect_ggplot(out)
    expect_identical(out$data, ref$data)
})

test_that("plotExplanatoryPCs responds to choice of number of variables", {
    maxes <- apply(exppcs, 2, max, na.rm=TRUE)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=2, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs[,order(maxes, decreasing=TRUE)[1:2]])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=Inf, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs)
    expect_ggplot(out)
    expect_identical(out$data, ref$data)

    out <- plotExplanatoryPCs(normed, nvars_to_plot=0, npcs=nrow(exppcs))
    ref <- plotExplanatoryPCs(exppcs[,0,drop=FALSE])
    expect_ggplot(out)
    expect_identical(out$data, ref$data)
})

test_that("plotExplanatoryPCs handles silly inputs.", {
    expect_ggplot(suppressWarnings(plotExplanatoryPCs(exppcs[0,,drop=FALSE])))
    expect_ggplot(plotExplanatoryPCs(exppcs[,0,drop=FALSE]))
})
