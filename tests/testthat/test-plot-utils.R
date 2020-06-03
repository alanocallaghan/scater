## Tests for the visualization variable picker.
## library(scater); library(testthat); source("setup.R"); source("test-plot-utils.R")

example_sce <- normed

altExp(example_sce, "thing") <- normed[1:10,]
rownames(altExp(example_sce)) <- paste0(rownames(altExp(example_sce)), ".0")

altExp(example_sce, "other") <- normed[1:5,]
rownames(altExp(example_sce, 2)) <- paste0(rownames(altExp(example_sce, 2)), "-R")

test_that("retrieveCellInfo works in the basic case", {
    out <- retrieveCellInfo(example_sce, "Mutation_Status")
    expect_identical(out$val, example_sce$Mutation_Status)
    expect_identical(out$name, "Mutation_Status")

    out <- retrieveCellInfo(example_sce, "Cell_Cycle")
    expect_identical(out$val, example_sce$Cell_Cycle)
    expect_identical(out$name, "Cell_Cycle")

    # Known gene exprs. 
    out <- retrieveCellInfo(example_sce, "Gene_0001")
    expect_identical(out$val, logcounts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    out <- retrieveCellInfo(example_sce, "Gene_0100")
    expect_identical(out$val, logcounts(example_sce)["Gene_0100",])
    expect_identical(out$name, "Gene_0100")

    # Alternative experiments.
    out <- retrieveCellInfo(example_sce, "Gene_0005.0")
    expect_identical(out$val, logcounts(altExp(example_sce))["Gene_0005.0",])
    expect_identical(out$name, "Gene_0005.0")

    out <- retrieveCellInfo(example_sce, "Gene_0005-R")
    expect_identical(out$val, logcounts(altExp(example_sce, 2))["Gene_0005-R",])
    expect_identical(out$name, "Gene_0005-R")

    # Known not to be either.
    expect_error(retrieveCellInfo(example_sce, "WHEE"), "cannot find")
    expect_error(retrieveCellInfo(example_sce, "Mutation_Status", search = "assays"), "cannot find")
    expect_error(retrieveCellInfo(example_sce, "Gene_0002", search = "colData"), "cannot find")

    expect_identical(retrieveCellInfo(example_sce, NULL), list(name=NULL, value=NULL))
})

test_that("retrieveCellInfo handles clashes correctly", {
    example_sce$Gene_0002 <- seq_len(ncol(example_sce))
    out_m <- retrieveCellInfo(example_sce, "Gene_0002", search = "colData")
    expect_identical(out_m$val, example_sce$Gene_0002)
    expect_identical(out_m$name, "Gene_0002")

    out_f <- retrieveCellInfo(example_sce, "Gene_0002", search = "assays")
    expect_identical(out_f$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out_f$name, "Gene_0002")

    # Respects ordering of inputs.
    expect_identical(out_f, retrieveCellInfo(example_sce, "Gene_0002", search = c("assays", "colData")))
    expect_identical(out_m, retrieveCellInfo(example_sce, "Gene_0002", search = c("colData", "assays")))
})

test_that("retrieveCellInfo handles wrapped elements and columns with data.frames", {
    thing <- data.frame(B=runif(ncol(example_sce)))
    out <- retrieveCellInfo(example_sce, thing)
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    thing <- DataFrame(B=runif(ncol(example_sce)))
    out <- retrieveCellInfo(example_sce, thing)
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    thing <- I(runif(ncol(example_sce)))
    out <- retrieveCellInfo(example_sce, thing)
    expect_identical(out$val, as.numeric(thing))
    expect_identical(out$name, "")

    # Check errors.
    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(retrieveCellInfo(example_sce, rething), "one column")
    expect_error(retrieveCellInfo(example_sce, rething[1:10,1,drop=FALSE]), "number of rows")
})

###################################################

set.seed(1312313)
rowData(example_sce) <- DataFrame(HAPPY=runif(nrow(example_sce)), SAD=rbinom(nrow(example_sce), 1, 0.5)==1)

test_that("retrieveFeatureInfo works for rows with strings", {
    # Known metadata.
    out <- retrieveFeatureInfo(example_sce, "HAPPY")
    expect_identical(out$val, rowData(example_sce)$HAPPY)
    expect_identical(out$name, "HAPPY")

    out <- retrieveFeatureInfo(example_sce, "SAD")
    expect_identical(out$val, rowData(example_sce)$SAD)
    expect_identical(out$name, "SAD")

    # Known exprs.
    out <- retrieveFeatureInfo(example_sce, "Cell_001")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    out <- retrieveFeatureInfo(example_sce, "Cell_010")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_010"])
    expect_identical(out$name, "Cell_010")

    # Handles errors properly.
    expect_error(retrieveFeatureInfo(example_sce, "whee"), "cannot find")
    expect_error(retrieveFeatureInfo(example_sce, "HAPPY", search = "assays"), "cannot find")
    expect_error(retrieveFeatureInfo(example_sce, "Cell_010", search = "rowData"), "cannot find")

    expect_identical(retrieveFeatureInfo(example_sce, NULL), list(name=NULL, value=NULL))
})

test_that("retrieveFeatureInfo is responsive to search mode", {
    rowData(example_sce)$Cell_002 <- seq_len(nrow(example_sce))
    out_m <- retrieveFeatureInfo(example_sce, "Cell_002", search = "rowData")
    expect_identical(out_m$val, rowData(example_sce)$Cell_002)
    expect_identical(out_m$name, "Cell_002")

    out_f <- retrieveFeatureInfo(example_sce, "Cell_002", , search = "assays")
    expect_identical(out_f$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out_f$name, "Cell_002")

    expect_identical(out_m, retrieveFeatureInfo(example_sce, "Cell_002", search=c("rowData", "assays")))
    expect_identical(out_f, retrieveFeatureInfo(example_sce, "Cell_002", search=c("assays", "rowData")))
})
     
test_that("retrieveFeatureInfo works for rows with data.frames", {
    thing <- data.frame(B=runif(nrow(example_sce)))
    out <- retrieveFeatureInfo(example_sce, thing, )
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    thing <- DataFrame(B=runif(nrow(example_sce)))
    out <- retrieveFeatureInfo(example_sce, thing)
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    thing <- I(runif(nrow(example_sce)))
    out <- retrieveFeatureInfo(example_sce, thing)
    expect_identical(out$val, as.numeric(thing))
    expect_identical(out$name, "")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(retrieveFeatureInfo(example_sce, rething, ), "one column")
    expect_error(retrieveFeatureInfo(example_sce, rething[1:10,1,drop=FALSE], ), "number of rows")
})

test_that("gg functions work as expected", {
    gg <- ggcells(example_sce, mapping=aes(x=Gene_0001, y=Gene_0002))
    expect_s3_class(gg, "ggplot")
    
    gg <- ggfeatures(example_sce, mapping=aes(x=Cell_001, y=Cell_002))
    expect_s3_class(gg, "ggplot")
})

