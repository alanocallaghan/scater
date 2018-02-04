## Tests for plotting functions
## This stress-tests the plotting functions for all the different parameter settings. 
## library(scater); library(testthat); source("test-plotting.R")

data("sc_example_counts")
data("sc_example_cell_info")
example_sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)
example_sce <- calculateQCMetrics(example_sce)
example_sce <- normalize(example_sce)

#################################################
# Testing the baseline visualization picker.

test_that("visualization variable picker works properly: columns", {
    out <- scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column")
    expect_identical(out$val, example_sce$Mutation_Status)
    expect_identical(out$name, "Mutation_Status")

    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    # Responsive to exprs_values.
    out <- scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)["Gene_0001",])
    expect_identical(out$name, "Gene_0001")

    # Responsive to search mode.
    example_sce$Gene_0002 <- seq_len(ncol(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "metadata")
    expect_identical(out_m$val, example_sce$Gene_0002)
    expect_identical(out_m$name, "Gene_0002")

    out_f <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out_f$name, "Gene_0002")

    out_a <- scater:::.choose_vis_values(example_sce, "Gene_0002", mode="column", search = "any")
    expect_identical(out_a, out_m)

    # Handles named strings properly.
    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "metadata"), mode="column")
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column")
    expect_identical(out$val, logcounts(example_sce)["Gene_0002",])
    expect_identical(out$name, "Gene_0002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Gene_0002", "exprs"), mode="column", search="metadata") # overriden by 'search'.
    expect_identical(out$val, example_sce$Gene_0002)
    expect_identical(out$name, "Gene_0002")

    # Handles nested strings properly.
    example_sce$blah <- DataFrame(X=runif(ncol(example_sce)), Y=sample(LETTERS, ncol(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="column", search="any")
    expect_identical(out$val, example_sce$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="metadata")
    expect_identical(out$val, example_sce$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="column", search="exprs"), "character vector")

    # Works for input data.frames.
    thing <- data.frame(B=runif(ncol(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "column")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "column"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "column"), "number of rows")

    # Handles errors properly.
    expect_error(scater:::.choose_vis_values(example_sce, "whee", mode="column", search = "any"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode="column", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode="column", search = "metadata"), "cannot find .* in metadata")
})

test_that("visualization variable picker works properly: rows", {
    out <- scater:::.choose_vis_values(example_sce, "is_feature_control", mode="row")
    expect_identical(out$val, rowData(example_sce)$is_feature_control)
    expect_identical(out$name, "is_feature_control")

    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    # Responsive to exprs_values.
    out <- scater:::.choose_vis_values(example_sce, "Cell_001", mode="row", exprs_values="counts")
    expect_identical(out$val, counts(example_sce)[,"Cell_001"])
    expect_identical(out$name, "Cell_001")

    # Responsive to search mode.
    rowData(example_sce)$Cell_002 <- seq_len(nrow(example_sce))
    out_m <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "metadata")
    expect_identical(out_m$val, rowData(example_sce)$Cell_002)
    expect_identical(out_m$name, "Cell_002")

    out_f <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "exprs")
    expect_identical(out_f$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out_f$name, "Cell_002")

    out_a <- scater:::.choose_vis_values(example_sce, "Cell_002", mode="row", search = "any")
    expect_identical(out_a, out_m)

    # Handles named strings properly.
    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "metadata"), mode="row")
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row")
    expect_identical(out$val, logcounts(example_sce)[,"Cell_002"])
    expect_identical(out$name, "Cell_002")

    out <- scater:::.choose_vis_values(example_sce, setNames("Cell_002", "exprs"), mode="row", search="metadata") # overriden by 'search'.
    expect_identical(out$val, rowData(example_sce)$Cell_002)
    expect_identical(out$name, "Cell_002")

    # Handles nested strings properly.
    rowData(example_sce)$blah <- DataFrame(X=runif(nrow(example_sce)), Y=sample(LETTERS, nrow(example_sce), replace=TRUE))
    out <- scater:::.choose_vis_values(example_sce, c("blah", "X"), mode="row", search="any")
    expect_identical(out$val, rowData(example_sce)$blah$X)
    expect_identical(out$name, "blah:X")

    out <- scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="metadata")
    expect_identical(out$val, rowData(example_sce)$blah$Y)
    expect_identical(out$name, "blah:Y")
    
    expect_error(scater:::.choose_vis_values(example_sce, c("blah", "Y"), mode="row", search="exprs"), "character vector")

    # Works for input data.frames.
    thing <- data.frame(B=runif(nrow(example_sce)))
    out <- scater:::.choose_vis_values(example_sce, thing, mode = "row")
    expect_identical(out$val, thing$B)
    expect_identical(out$name, "B")

    rething <- data.frame(B=runif(ncol(example_sce)), C=2)
    expect_error(scater:::.choose_vis_values(example_sce, rething, mode = "row"), "one column")
    expect_error(scater:::.choose_vis_values(example_sce, thing[1:10,,drop=FALSE], mode = "row"), "number of rows")

    # Handles errors properly.
    expect_error(scater:::.choose_vis_values(example_sce, "whee", mode = "row", search = "any"), "cannot find .* in any fields")
    expect_error(scater:::.choose_vis_values(example_sce, "Mutation_Status", mode = "row", search = "exprs"), "cannot find .* in exprs")
    expect_error(scater:::.choose_vis_values(example_sce, "Gene_0001", mode = "row", search = "metadata"), "cannot find .* in metadata")
})

test_that("visualization variable picker works properly: misc", {
    example_sce$BLAH <- 0
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column")
    expect_identical(out$val, example_sce$BLAH)
    expect_identical(out$name, "BLAH")
                 
    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", discard_solo = TRUE)
    expect_identical(out$val, NULL)
    expect_identical(out$name, NULL) 

    # Checking what happens with coerced factors.
    example_sce$BLAH <- rep(LETTERS, length.out=ncol(example_sce))

    out <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column")
    expect_identical(out$val, example_sce$BLAH)
    expect_identical(out$name, "BLAH")

    out2 <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE)
    expect_identical(out2$val, as.factor(out$val))                 
    expect_identical(out2$name, "BLAH")

    expect_error(scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE, level_limit = 10), 
                 "number of unique levels")

    out3 <- scater:::.choose_vis_values(example_sce, "BLAH", mode="column", coerce_factor = TRUE, level_limit = 100)
    expect_identical(out2, out3)
})

#################################################
# Checking the specific plotting functions.

test_that("plotScater works as expected", {
    expect_s3_class(plotScater(example_sce), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block2 = "Cell_Cycle"), "ggplot")

    expect_s3_class(plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Mutation_Status", block2 = "Cell_Cycle", colour_by = "Gene_0001"), "ggplot")

    cpm(example_sce) <- calculateCPM(example_sce)
    expect_s3_class(plotScater(example_sce, exprs_values="cpm"), "ggplot")
    expect_error(plotScater(example_sce, exprs_values="tpm"), "not in names")
})

#################################################
# Checking the reduced dimension wrappers.

test_that("we can produce PCA scatterplots", {
    example_sce <- runPCA(example_sce)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_s3_class(P <- plotPCA(example_sce), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotPCA(example_sce, size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", legend="none"), "ggplot")
    expect_s3_class(plotPCA(example_sce, exprs_values="counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, percentVar=c(19, 5)), "ggplot")
    
    # Checking that re-running works, responsive to feature scaling.
    expect_s3_class(P2 <- plotPCA(example_sce, rerun=TRUE, run_args=list(scale_features=FALSE)), "ggplot")
    expect_false(isTRUE(all.equal(P, P2)))

    reducedDim(example_sce, "PCA") <- NULL
    expect_s3_class(P3 <- plotPCA(example_sce), "ggplot")
    expect_equal(P, P3)

    expect_s3_class(P4 <- plotPCA(example_sce, run_args=list(scale_features=FALSE)), "ggplot")
    expect_equal(P2, P4)
    expect_false(isTRUE(all.equal(P, P4)))
})

test_that("we can produce PCA pairplots", {
    example_sce <- runPCA(example_sce, ncomponents=4)
    expect_identical(reducedDimNames(example_sce), "PCA")

    # Checking that visual parameters work.
    expect_s3_class(P <- plotPCA(example_sce, ncomponents=4), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", legend="none"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, exprs_values="counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, percentVar=c(19, 5, 3, 2)), "ggplot")
    
    # Checking that re-running works, responsive to feature scaling.
    expect_s3_class(P2 <- plotPCA(example_sce, ncomponents=4, rerun=TRUE, run_args=list(scale_features=FALSE)), "ggplot")
    expect_false(isTRUE(all.equal(P, P2)))

    reducedDim(example_sce, "PCA") <- NULL
    expect_s3_class(P3 <- plotPCA(example_sce, ncomponents=4), "ggplot")
    expect_equal(P, P3)

    expect_s3_class(P4 <- plotPCA(example_sce, ncomponents=4, run_args=list(scale_features=FALSE)), "ggplot")
    expect_equal(P2, P4)
    expect_false(isTRUE(all.equal(P, P4)))
})

test_that("we can produce TSNE plots", {
    example_sce <- runTSNE(example_sce, rand_seed=100)
    expect_identical(reducedDimNames(example_sce), "TSNE")
    expect_s3_class(P <- plotTSNE(example_sce), "ggplot")

    # plotTSNE re-runs it correctly.
    reducedDim(example_sce, "TSNE") <- NULL
    expect_s3_class(P2 <- plotTSNE(example_sce, run_args=list(rand_seed=100)), "ggplot")
    expect_equal(P, P2)

    # Responsive to changes in parameters.
    expect_s3_class(P3 <- plotTSNE(example_sce, run_args=list(perplexity=10)), "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))

    # Handles multiple components properly.
    expect_s3_class(P4 <- plotTSNE(example_sce, ncomponents=4, run_args=list(rand_seed=20)), "ggplot")
    example_sce <- runTSNE(example_sce, ncomponents=4, rand_seed=20)
    expect_equal(plotTSNE(example_sce, ncomponents=4), P4)
})

test_that("we can produce diffusion maps", {
    example_sce <- runDiffusionMap(example_sce, rand_seed=100)
    expect_identical(reducedDimNames(example_sce), "DiffusionMap")
    expect_s3_class(P <- plotDiffusionMap(example_sce), "ggplot")

    # plotDiffusionMap re-runs it correctly.
    reducedDim(example_sce, "DiffusionMap") <- NULL
    expect_s3_class(P2 <- plotDiffusionMap(example_sce, run_args=list(rand_seed=100)), "ggplot")
#    expect_equal(P, P2) # it seems as if destiny::DiffusionMap does not respond to the seed!

    # Responsive to changes in parameters.
    expect_s3_class(P3 <- plotDiffusionMap(example_sce, run_args=list(k=13)), "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))

    # Handles multiple components properly.
    expect_s3_class(P4 <- plotDiffusionMap(example_sce, ncomponents=4, run_args=list(rand_seed=20)), "ggplot")
#    example_sce <- runDiffusionMap(example_sce, ncomponents=4, rand_seed=20)
#    expect_equal(plotDiffusionMap(example_sce, ncomponents=4), P4)
})

test_that("we can produce MDS plots", {
    example_sce <- runMDS(example_sce)
    expect_identical(reducedDimNames(example_sce), "MDS")
    expect_s3_class(P <- plotMDS(example_sce), "ggplot")

    # plotMDS re-runs it correctly.
    reducedDim(example_sce, "MDS") <- NULL
    expect_s3_class(P2 <- plotMDS(example_sce), "ggplot")
    expect_equal(P, P2)

    # Responsive to changes in parameters.
    expect_s3_class(P3 <- plotMDS(example_sce, run_args=list(method="manhattan")), "ggplot")
    expect_false(isTRUE(all.equal(P, P3)))

    # Handles multiple components properly.
    expect_s3_class(P4 <- plotMDS(example_sce, ncomponents=4), "ggplot")
    example_sce <- runMDS(example_sce, ncomponents=4)
    expect_equal(plotMDS(example_sce, ncomponents=4), P4)
})

#################################################
# Testing plotExpression

test_that("we can produce expression plots with different expression values", {
    # Testing various 'by x' scenarios.
    for (gene_set in list("Gene_0001", rownames(example_sce)[1:5], 10, 6:20)) { # different numbers of genes, types of specification.
        for (x in list(NULL, "Cell_Cycle", "Gene_0100")) { # nothing, categorical, or continuous.
            expect_s3_class(plotExpression(example_sce, gene_set, x = x), "ggplot")
            expect_s3_class(plotExpression(example_sce, gene_set, x = x, colour_by = "Cell_Cycle"), "ggplot")
            expect_s3_class(plotExpression(example_sce, gene_set, x = x, size_by = "Gene_0001"), "ggplot")
            expect_s3_class(plotExpression(example_sce, gene_set, x = x, shape_by = "Treatment"), "ggplot")
        }
    }

    # Testing different visualization schemes.
    gene_set <- rownames(example_sce)[1:20]
    expect_s3_class(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Testing options when dealing with many genes and no 'x' specified.
    expect_s3_class(plotExpression(example_sce, gene_set, one_facet=FALSE), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, feature_colours=FALSE), "ggplot")

    # Fiddling with all the semi-analysis options.
    expect_s3_class(plotExpression(example_sce, gene_set, show_violin=FALSE), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, show_median=TRUE), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, jitter="jitter"), "ggplot")

    expect_s3_class(plotExpression(example_sce, gene_set, x="Gene_0001", show_smooth=TRUE), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, x="Gene_0001", show_smooth=TRUE, show_se=FALSE), "ggplot")

    # Checking for behaviour with different values.
    expect_s3_class(plotExpression(example_sce, gene_set, x = "Mutation_Status", exprs_values = "counts"),  "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, x = "Mutation_Status", exprs_values = "counts", log2_values = TRUE),  "ggplot")
    expect_error(plotExpression(example_sce, rownames(example_sce)[1:6], exprs_values = "silly"), "not in names")
})

#################################################
# Testing plotPlatePosition

test_that("we can produce plots showing cells in plate position", {
    ## Define plate positions
    row <- rep(LETTERS[1:5], each = 8)
    col <- rep(1:8, 5)
    plate_position <- paste0(row, col)

    # Different types of inputs are accepted.    
    expect_s3_class(plotPlatePosition(example_sce, list(row=row, column=col)), "ggplot")
    expect_s3_class(P <- plotPlatePosition(example_sce, plate_position), "ggplot")

    alt <- example_sce
    alt$plate_position <- plate_position
    expect_s3_class(plotPlatePosition(alt, colour_by="Cell_Cycle"), "ggplot")

    # Different types of colouring, shaping and sizing are possible.
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking that other arguments are passed through.
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment", legend="none"), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", exprs_values = "counts"), "ggplot")

    # Checking that an error is thrown,
    expect_error(plotPlatePosition(example_sce, paste0(col, row)), "invalid format") 
})

test_that("we can produce plots for metadata", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))

    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))

    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))
    
    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))

    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))

    expect_that(plotColData(example_sce, x = "total_counts",
        y = "total_features_by_counts", colour = "Mutation_Status"), is_a("ggplot"))

    expect_that(plotRowData(example_sce, x = "n_cells_by_counts", 
        y = "log10_total_counts"), is_a("ggplot"))

    expect_that(plotRowData(example_sce, x = "n_cells_by_counts", 
        y = "log10_total_counts"), is_a("ggplot"))
})

test_that("plotExprsVsTxLength works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    fd <- data.frame(
        gene_id = rownames(sc_example_counts), 
        feature_id = paste("feature", rep(1:500, each = 4), sep = "_"),
        median_tx_length = rnorm(2000, mean = 5000, sd = 500))
    rownames(fd) <- rownames(sc_example_counts)
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info, rowData = fd)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)
    rowData(example_sce)$group <- rep(1:4, each = 500)
    
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length")
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length", 
                              show_smooth = TRUE)
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length", 
                              show_smooth = TRUE, show_exprs_sd = TRUE)
    expect_that(p1, is_a("ggplot"))
    
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length",
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, colour_by = "group")
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length", 
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, size_by = "group")
    expect_that(p1, is_a("ggplot"))
    
    rowData(example_sce)$group <- rep(letters[1:4], each = 500)
    p1 <- plotExprsVsTxLength(example_sce, "median_tx_length", 
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, shape_by = "group")
    expect_that(p1, is_a("ggplot"))
    
    
    ## using matrix of tx length values in assayData(object)
    mat <- matrix(rnorm(ncol(example_sce) * nrow(example_sce), mean = 5000,
                        sd = 500), nrow = nrow(example_sce))
  
    dimnames(mat) <- dimnames(example_sce)
    assay(example_sce, "tx_len") <- mat
    p1 <-  plotExprsVsTxLength(example_sce, "tx_len", length_is_assay = TRUE,
                               show_smooth = TRUE, show_exprs_sd = TRUE)
    expect_that(p1, is_a("ggplot"))
    
    ## using a vector of tx length values
    p1 <- plotExprsVsTxLength(example_sce, data.frame(Length=rnorm(2000, mean = 5000, sd = 500)))
    expect_that(p1, is_a("ggplot"))
})

