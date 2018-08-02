## Tests for plotting functions
## This stress-tests the plotting functions for all the different parameter settings. 
## library(scater); library(testthat); source("setup-sce.R"); source("test-plotting.R")

example_sce <- calculateQCMetrics(normed, exprs_values="counts")

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
    # Testing out the options (need an expect_* clause to avoid skipping).
    expect_error(plotHeatmap(example_sce, features=rownames(example_sce)[1:10]), NA)
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

#################################################
# Testing plotRLE

test_that("plotRLE works as expected", {
    cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)

    # With minimal or full.
    for (style in c("minimal", "full")) {
        p <- plotRLE(example_sce, colour_by = "Mutation_Status", style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, legend=FALSE, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, exprs_values="cpm", exprs_logged=FALSE, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, colour_by = "Gene_0004", style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, colour_by = "Gene_0004", by_exprs_values = "cpm", style=style)
        expect_s3_class(p, "ggplot")
    }
})

