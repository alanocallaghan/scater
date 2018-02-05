## Tests for plotting functions
## This stress-tests the plotting functions for all the different parameter settings. 
## library(scater); library(testthat); source("test-plotting.R")

data("sc_example_counts")
data("sc_example_cell_info")
example_sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)
example_sce <- calculateQCMetrics(example_sce)
suppressWarnings(example_sce <- normalize(example_sce))

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
# Checking plotScater.

test_that("plotScater works as expected", {
    expect_s3_class(plotScater(example_sce), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block2 = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Treatment", block2 = "Cell_Cycle"), "ggplot")

    # Different types of colouring are possible
    expect_s3_class(plotScater(example_sce, colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Gene_0001"), "ggplot")

    expect_s3_class(plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Mutation_Status", colour_by = "Gene_0001"), "ggplot")

    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Cell_Cycle"), "ggplot")
    expect_s3_class(plotScater(example_sce, block1 = "Cell_Cycle", block2 = "Treatment", colour_by = "Gene_0001"), "ggplot")
    
    expect_s3_class(plotScater(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotScater(example_sce, colour_by = "Treatment", by_show_single = TRUE), "ggplot")

    # Responds to different type of expression values.
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
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")

    # Checking other arguments are passed successfully to plotReducedDim.
    expect_s3_class(plotPCA(example_sce, colour_by = "Cell_Cycle", legend = FALSE), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, colour_by = "Treatment", by_show_single = TRUE), "ggplot")
    expect_s3_class(plotPCA(example_sce, percentVar = c(19, 5)), "ggplot")
    
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
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Cell_Cycle", legend = FALSE), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotPCA(example_sce, ncomponents=4, colour_by = "Treatment", by_show_single = TRUE), "ggplot")
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
    
    # Throws a warning if we try to specify this locally.
    expect_warning(plotTSNE(example_sce, perplexity=10), "non-plotting arguments")

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

    # Throws a warning if we try to specify this locally.
    expect_warning(plotDiffusionMap(example_sce, k=10), "non-plotting arguments")

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

    # Throws a warning if we try to specify this locally.
    expect_warning(plotMDS(example_sce, method="manhattan"), "non-plotting arguments")

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
    
    expect_s3_class(plotExpression(example_sce, gene_set, size_by = "Gene_0001", shape_by = "Treatment", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotExpression(example_sce, gene_set, size_by = "Gene_0001", shape_by = "Treatment", by_show_single = TRUE), "ggplot")

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
# Testing plotHeatmap

test_that("we can produce heatmaps", {
    # Testing out the options.          
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10])
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], columns = 1:20)
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], exprs_values='counts')

    # Colour parameters for the expression values.
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], zlim=c(0, 2))
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], color=viridis::viridis(20))
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], center=TRUE, symmetric=TRUE)
         
    # Testing out the column colouring. 
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                colour_columns_by=c("Mutation_Status", "Cell_Cycle"))
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                colour_columns_by=c("Mutation_Status", "Gene_0001"), 
                by_exprs_values = "logcounts", by_show_single = TRUE)

    # Testing out passing arguments to pheatmap.
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], fontsize = 20, legend = FALSE)
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
    
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", shape_by = "Treatment", by_show_single = TRUE), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", shape_by = "Treatment", by_exprs_values = "counts"), "ggplot")

    # Checking that other arguments are passed through.
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment", legend = FALSE), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")

    # Checking that an error is thrown,
    expect_error(plotPlatePosition(example_sce, paste0(col, row)), "invalid format") 
})

#################################################
# Testing plotColData and plotRowData

test_that("we can produce plots for column metadata", {
    for (y in c("total_features_by_counts", "Mutation_Status")) { # discrete or continuous.
        for (x in list(NULL, "total_counts", "Cell_Cycle")) { # nothing, discrete or continuous. 
              expect_s3_class(plotColData(example_sce, x = x, y = y, colour_by = "Treatment"), "ggplot")
              expect_s3_class(plotColData(example_sce, x = x, y = y, size_by = "Gene_0001"), "ggplot")
              expect_s3_class(plotColData(example_sce, x = x, y = y, shape_by = "Treatment"), "ggplot")
        }
    }

    # Testing more visualization schemes.
    expect_s3_class(plotColData(example_sce, "total_counts", colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    
    # Testing that other arguments are passed through.
    expect_s3_class(plotColData(example_sce, "total_counts", colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment", legend = FALSE), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", size_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", colour_by = "Treatment", by_show_single = TRUE), "ggplot")

    # Fiddling with all the semi-analysis options.
    expect_s3_class(plotColData(example_sce, "total_counts", show_violin=FALSE), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", show_median=TRUE), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", jitter="jitter"), "ggplot")

    expect_s3_class(plotColData(example_sce, "total_counts", x="total_features_by_counts", show_smooth=TRUE), "ggplot")
    expect_s3_class(plotColData(example_sce, "total_counts", x="total_features_by_counts", show_smooth=TRUE, show_se=FALSE), "ggplot")

    # Checking that it doesn't try to retrieve expression data.
    expect_error(plotColData(example_sce, "Gene_0001", exprs_values = "counts"), "cannot find .* in metadata fields")
    expect_error(plotColData(example_sce, "total_counts", x="Gene_0001", exprs_values = "counts"), "cannot find .* in metadata fields")
})

test_that("we can produce plots for row metadata", {
    rowData(example_sce)$WHEE <- rep(LETTERS[1:10], length.out=nrow(example_sce))

    for (y in c("mean_counts", "is_feature_control")) { # discrete or continuous.
        for (x in list(NULL, "n_cells_by_counts", "WHEE")) { # nothing, discrete or continuous. 
              expect_s3_class(plotRowData(example_sce, x = x, y = y, colour_by = "WHEE"), "ggplot")
              expect_s3_class(plotRowData(example_sce, x = x, y = y, size_by = "Cell_001"), "ggplot")
              expect_s3_class(plotRowData(example_sce, x = x, y = y, shape_by = "WHEE"), "ggplot")
        }
    }

    # Testing more visualization schemes.
    expect_s3_class(plotRowData(example_sce, "total_counts", colour_by = "is_feature_control", size_by = "Cell_002"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", colour_by = "is_feature_control", shape_by = "WHEE"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", size_by = "Cell_002", shape_by = "WHEE"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", colour_by = "is_feature_control", size_by = "Cell_002", shape_by = "WHEE"), "ggplot")
    
    # Testing that other arguments are passed through.
    expect_s3_class(plotRowData(example_sce, "total_counts", colour_by = "is_feature_control", size_by = "Cell_002", shape_by = "WHEE", legend = FALSE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", size_by = "Cell_002", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", colour_by = "is_feature_control", by_show_single = TRUE), "ggplot")

    # Fiddling with all the semi-analysis options.
    expect_s3_class(plotRowData(example_sce, "total_counts", show_violin=FALSE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", show_median=TRUE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", jitter="jitter"), "ggplot")

    expect_s3_class(plotRowData(example_sce, "total_counts", x="n_cells_by_counts", show_smooth=TRUE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "total_counts", x="n_cells_by_counts", show_smooth=TRUE, show_se=FALSE), "ggplot")

    # Checking that it doesn't try to retrieve expression data.
    expect_error(plotRowData(example_sce, "Cell_002", exprs_values = "counts"), "cannot find .* in metadata fields")
    expect_error(plotRowData(example_sce, "total_counts", x="Cell_002", exprs_values = "counts"), "cannot find .* in metadata fields")
})

#################################################
# Testing plotExprsVsTxLength

test_that("plotExprsVsTxLength works as expected", {
    rowData(example_sce)$median_tx_length <- rnorm(2000, mean = 5000, sd = 500)
    rowData(example_sce)$group <- rep(1:4, each = 500)
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", exprs_values="counts"), "ggplot")
   
    # Testing more visualization schemes.
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "is_feature_control", size_by = "Cell_002"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "is_feature_control", shape_by = "group"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", size_by = "Cell_002", shape_by = "group"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "is_feature_control", size_by = "Cell_002", shape_by = "group"), "ggplot")

    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "Cell_002", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "is_feature_control", by_show_single = TRUE), "ggplot")
 
    # Testing various semi-analysis options.
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_exprs_sd = TRUE), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE, show_exprs_sd = TRUE), "ggplot")
    
    # Testing visualization options.
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", colour_by = "group"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_exprs_sd = TRUE, colour_by = "group"), "ggplot") # checking proper interaction with geom_pointrange.
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_exprs_sd = TRUE, shape_by= "group"), "ggplot")
    expect_s3_class(plotExprsVsTxLength(example_sce, "median_tx_length", show_exprs_sd = TRUE, size_by= "n_cells_by_counts"), "ggplot")
    
    ## using matrix of tx length values in assayData(object)
    mat <- matrix(rnorm(ncol(example_sce) * nrow(example_sce), mean = 5000, sd = 500), nrow = nrow(example_sce))
    dimnames(mat) <- dimnames(example_sce)
    assay(example_sce, "tx_len") <- mat

    expect_s3_class(plotExprsVsTxLength(example_sce, "tx_len", length_is_assay = TRUE, show_smooth = TRUE, show_exprs_sd = TRUE), "ggplot")
    
    ## using a vector of tx length values
    expect_s3_class(plotExprsVsTxLength(example_sce, data.frame(Length=rnorm(2000, mean = 5000, sd = 500))), "ggplot")
})

