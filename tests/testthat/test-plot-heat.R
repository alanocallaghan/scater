## Tests for the heatmap-related plotting functions
## library(scater); library(testthat); source("setup.R"); source("test-plot-heat.R")

example_sce <- normed 
colData(example_sce) <- cbind(colData(example_sce), perCellQCMetrics(example_sce))
rowData(example_sce) <- cbind(rowData(example_sce), perFeatureQCMetrics(example_sce))
rowData(example_sce)$ENS <- gsub("Gene", "ENS", rownames(example_sce))
rowData(example_sce)$ENS_e1 <- rowData(example_sce)$ENS
rowData(example_sce)$ENS_e1[1:10] <- NA
rowData(example_sce)$ENS_e2 <- "constant"

#################################################
# Testing plotHeatmap

test_that("we can produce heatmaps", {
    # Testing out the options (need an expect_* clause to avoid skipping).
    expect_error(plotHeatmap(example_sce, features=rownames(example_sce)[1:10]), NA)
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], columns = 1:20)
    plotHeatmap(example_sce, features=rowData(example_sce)[1:10, "ENS"], 
        swap_rownames = "ENS", columns = 1:20)
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], exprs_values='counts')

    # Colour parameters for the expression values.
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], zlim=c(0, 2))
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], color=viridis::viridis(20))
    expect_warning(
        plotHeatmap(example_sce, features=rownames(example_sce)[1:10], center=TRUE, symmetric=TRUE)
    )
         
    # Testing out the column colouring. 
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                colour_columns_by=c("Mutation_Status", "Cell_Cycle"))

    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                colour_columns_by=c("Mutation_Status", "Gene_0001"), 
                by_exprs_values = "logcounts")

    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                colour_columns_by=list(I(example_sce$Mutation_Status), "Gene_0001"), 
                by_exprs_values = "logcounts")

    # Testing out the column ordering + colouring. 
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10],
                order_columns_by=c("Mutation_Status", "Gene_0001"), 
                by_exprs_values = "logcounts")

    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], columns=1:10,
                order_columns_by=c("Mutation_Status", "Gene_0001"), 
                by_exprs_values = "logcounts")

    # Testing that column colouring still works when columns have no names.
    unnamed <- example_sce
    colnames(unnamed) <- NULL
    plotHeatmap(unnamed, features=rownames(unnamed)[1:10],
                colour_columns_by=c("Mutation_Status", "Gene_0001")) 

    # Testing out passing arguments to pheatmap.
    plotHeatmap(example_sce, features=rownames(example_sce)[1:10], fontsize = 20, legend = FALSE)

    plotHeatmap(example_sce, features = rowData(example_sce)[1:10, "ENS"], 
        swap_rownames = "ENS", columns = 1:20)
    expect_error(
        plotHeatmap(example_sce, features = NA, swap_rownames = "ENS_e1",
            columns = 1:20),
        "must have n >= 2 objects to cluster"
    )

    expect_error(
        plotHeatmap(example_sce, features = "constant", swap_rownames = "ENS_e2",
            columns = 1:20),
        "must have n >= 2 objects to cluster"
    )

    expect_error(
        plotHeatmap(example_sce, features = NA, swap_rownames = "sdfsf",
            columns = 1:20),
        "Cannot find column sdfsf in rowData"
    )
})

#################################################
# Testing plotGroupedHeatmap

test_that("we can produce grouped heatmaps", {
    example_sce$Group <- paste0(example_sce$Treatment, "+", example_sce$Mutation_Status)

    # Testing out the options (need an expect_* clause to avoid skipping)
    expect_error(plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10], group="Group"), NA)
    plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10], columns=1:20, group="Group")
    plotGroupedHeatmap(example_sce, features=1:10, group="Group")

    # Works with blocking.
    plotGroupedHeatmap(example_sce, features=1:10, group="Group", block="Cell_Cycle")
    plotGroupedHeatmap(example_sce, features=1:10, group="Group", block="Cell_Cycle", columns=20:30)

    # Works with the various color options.
    expect_warning(
        plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
            group="Group", center=TRUE, symmetric=TRUE)
    )

    expect_warning(
        plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
            group="Group", center=TRUE, symmetric=TRUE)
    )

    plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10], group="Group", color=viridis::viridis(20))
    plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10], group="Group", zlim=c(0, 2))

    # Works with rownames swapping.
    plotGroupedHeatmap(example_sce, features = rowData(example_sce)[1:10, "ENS"], 
        group="Group", swap_rownames = "ENS", columns = 1:20)
})
