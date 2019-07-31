## Tests for other plotting functions
## library(scater); library(testthat); source("setup-sce.R"); source("test-plot-other.R")

example_sce <- normed 
colData(example_sce) <- cbind(colData(example_sce), perCellQCMetrics(example_sce))
rowData(example_sce) <- cbind(rowData(example_sce), perFeatureQCMetrics(example_sce))

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

    # Testing that column colouring still works when columns have no names.
    unnamed <- example_sce
    colnames(unnamed) <- NULL
    plotHeatmap(unnamed, features=rownames(unnamed)[1:10],
                colour_columns_by=c("Mutation_Status", "Gene_0001")) 

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
    expect_s3_class(plotPlatePosition(alt, colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment", add_legend = FALSE), "ggplot")
    expect_s3_class(plotPlatePosition(alt, size_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")

    # Checking that an error is thrown,
    expect_error(plotPlatePosition(example_sce, paste0(col, row)), "invalid format") 
})

#################################################
# Testing plotColData and plotRowData

test_that("we can produce plots for column metadata", {
    example_sce <- addQCPerCell(example_sce)

    for (y in c("detected", "Mutation_Status")) { # discrete or continuous.
        for (x in list(NULL, "sum", "Cell_Cycle")) { # nothing, discrete or continuous. 
              expect_s3_class(plotColData(example_sce, x = x, y = y, colour_by = "Treatment"), "ggplot")
              expect_s3_class(plotColData(example_sce, x = x, y = y, size_by = "Gene_0001"), "ggplot")
              expect_s3_class(plotColData(example_sce, x = x, y = y, shape_by = "Treatment"), "ggplot")
        }
    }

    # Testing more visualization schemes.
    expect_s3_class(plotColData(example_sce, "sum", colour_by = "Cell_Cycle", size_by = "Gene_0001"), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", colour_by = "Cell_Cycle", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment"), "ggplot")
    
    # Testing that other arguments are passed through.
    expect_s3_class(plotColData(example_sce, "sum", colour_by = "Cell_Cycle", size_by = "Gene_0001", shape_by = "Treatment", add_legend = FALSE), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", size_by = "Gene_0001", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", colour_by = "Treatment", by_show_single = TRUE), "ggplot")

    # Fiddling with all the semi-analysis options.
    expect_s3_class(plotColData(example_sce, "sum", show_violin=FALSE), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", show_median=TRUE), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", jitter_type="jitter"), "ggplot")

    expect_s3_class(plotColData(example_sce, "sum", x="detected", show_smooth=TRUE), "ggplot")
    expect_s3_class(plotColData(example_sce, "sum", x="detected", show_smooth=TRUE, show_se=FALSE), "ggplot")

    # Checking that it doesn't try to retrieve expression data.
    expect_error(plotColData(example_sce, "Gene_0001", exprs_values = "counts"), "cannot find")
    expect_error(plotColData(example_sce, "sum", x="Gene_0001", exprs_values = "counts"), "cannot find")
})

test_that("we can produce plots for row metadata", {
    rowData(example_sce)$WHEE <- rep(LETTERS[1:10], length.out=nrow(example_sce))
    rowData(example_sce)$is_feature_control <- rbinom(nrow(example_sce), 1, 0.5)
    example_sce <- addQCPerFeature(example_sce)

    for (y in c("mean", "is_feature_control")) { # discrete or continuous.
        for (x in list(NULL, "detected", "WHEE")) { # nothing, discrete or continuous. 
              expect_s3_class(plotRowData(example_sce, x = x, y = y, colour_by = "WHEE"), "ggplot")
              expect_s3_class(plotRowData(example_sce, x = x, y = y, size_by = "Cell_001"), "ggplot")
              expect_s3_class(plotRowData(example_sce, x = x, y = y, shape_by = "WHEE"), "ggplot")
        }
    }

    # Testing more visualization schemes.
    expect_s3_class(plotRowData(example_sce, "mean", colour_by = "is_feature_control", size_by = "Cell_002"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", colour_by = "is_feature_control", shape_by = "WHEE"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", size_by = "Cell_002", shape_by = "WHEE"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", colour_by = "is_feature_control", size_by = "Cell_002", shape_by = "WHEE"), "ggplot")
    
    # Testing that other arguments are passed through.
    expect_s3_class(plotRowData(example_sce, "mean", colour_by = "is_feature_control", size_by = "Cell_002", shape_by = "WHEE", add_legend = FALSE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", size_by = "Cell_002", by_exprs_values = "counts"), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", colour_by = "is_feature_control", by_show_single = TRUE), "ggplot")

    # Fiddling with all the semi-analysis options.
    expect_s3_class(plotRowData(example_sce, "mean", show_violin=FALSE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", show_median=TRUE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", jitter_type="jitter"), "ggplot")

    expect_s3_class(plotRowData(example_sce, "mean", x="detected", show_smooth=TRUE), "ggplot")
    expect_s3_class(plotRowData(example_sce, "mean", x="detected", show_smooth=TRUE, show_se=FALSE), "ggplot")

    # Checking that it doesn't try to retrieve expression data.
    expect_error(plotRowData(example_sce, "Cell_002", exprs_values = "counts"), "cannot find")
    expect_error(plotRowData(example_sce, "mean", x="Cell_002", exprs_values = "counts"), "cannot find")
})

test_that("plotRowData works for other fields", {
    rowData(example_sce)$WHEE <- rep(LETTERS[1:10], length.out=nrow(example_sce))
    rowData(example_sce)$STUFF <- rep(letters[1:10], length.out=nrow(example_sce))

    gg <- plotRowData(example_sce, "mean",
        other_fields=c("WHEE", "STUFF"))
    expect_true("WHEE" %in% colnames(gg$data))
    expect_true("STUFF" %in% colnames(gg$data))

    # This should throw a warning.
    rowData(example_sce)$colour_by <- rowData(example_sce)$STUFF
    expect_warning(gg <- plotRowData(example_sce, "mean", 
        colour_by="STUFF", other_fields=c("colour_by")), "duplicated")
})

#################################################
# Testing plotRLE

test_that("plotRLE works as expected", {
    cpm(example_sce) <- calculateCPM(example_sce)

    # With minimal or full.
    for (style in c("minimal", "full")) {
        p <- plotRLE(example_sce, colour_by = "Mutation_Status", style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, add_legend=FALSE, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, exprs_values="cpm", exprs_logged=FALSE, style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, colour_by = "Gene_0004", style=style)
        expect_s3_class(p, "ggplot")
    
        p <- plotRLE(example_sce, colour_by = "Gene_0004", by_exprs_values = "cpm", style=style)
        expect_s3_class(p, "ggplot")
    }
})

