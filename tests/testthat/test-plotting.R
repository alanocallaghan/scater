## tests for plotting functions

context("test plotScater, plotPCA, plotTSNE, plotDiffusionMap, plotMDS, plotReducedDim")

test_that("we can produce default plots for SingleCellExperiment objects", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    
    p <- plotScater(example_sce)
    print(p)
    expect_that(p, is_a("ggplot"))
    expect_that(
        plotScater(example_sce, exprs_values = "counts", colour_by = "Cell_Cycle"),
        is_a("ggplot"))
    expect_that(
        plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle"),
        is_a("ggplot"))

    cpm(example_sce) <- calculateCPM(example_sce, use.size.factors = FALSE)
    expect_that(
        plotScater(example_sce, exprs_values = "cpm", block1 = "Treatment",
             block2 = "Mutation_Status", colour_by = "Cell_Cycle"),
        is_a("ggplot"))
    # What happens if chosen expression values are not available?
    expect_error(
        plotScater(example_sce, exprs_values = "tpm", block1 = "Treatment",
             colour_by = "Cell_Cycle"),
        "not in names")
})

test_that("we can produce PCA plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(plotPCA(example_sce), is_a("ggplot"))
})

test_that("we can produce t-SNE plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(calculateCPM(example_sce, 
                                            use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(plotTSNE(example_sce), is_a("ggplot"))
})

test_that("we can produce Diffusion Map plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              example_sce <- SingleCellExperiment(
                  assays = list(counts = sc_example_counts), 
                  colData = sc_example_cell_info)
              exprs(example_sce) <- log2(
                  calculateCPM(example_sce, use.size.factors = FALSE) + 1)
              example_sce <- calculateQCMetrics(example_sce)
              
              #expect_that(plotDiffusionMap(example_sce), is_a("ggplot"))
          })

test_that("we can produce MDS plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              example_sce <- SingleCellExperiment(
                  assays = list(counts = sc_example_counts), 
                  colData = sc_example_cell_info)
              exprs(example_sce) <- log2(
                  calculateCPM(example_sce, use.size.factors = FALSE) + 1)
              example_sce <- calculateQCMetrics(example_sce)

              expect_that(plotMDS(example_sce), is_a("ggplot"))
              expect_that(
                  plotMDS(example_sce, colour_by = "Cell_Cycle",
                          shape_by = "Treatment", size_by = "Mutation_Status"),
                  is_a("ggplot"))
          })


test_that("plotReducedDim works as expexted",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              example_sce <- SingleCellExperiment(
                  assays = list(counts = sc_example_counts), 
                  colData = sc_example_cell_info)
              exprs(example_sce) <- log2(calculateCPM(
                  example_sce, use.size.factors = FALSE) + 1)
              drop_genes <- apply(exprs(example_sce), 1, 
                                  function(x) {var(x) == 0})
              example_sce <- example_sce[!drop_genes, ]
              
              reducedDim(example_sce, "PCA") <- 
                  prcomp(t(exprs(example_sce)), scale. = TRUE)$x
              expect_that(plotReducedDim(example_sce, "PCA"), is_a("ggplot"))
              expect_that(plotReducedDim(
                  example_sce, "PCA", colour_by = "Cell_Cycle"), is_a("ggplot"))
              expect_that(plotReducedDim(
                  example_sce, "PCA", colour_by = "Cell_Cycle", 
                  shape_by = "Treatment"), is_a("ggplot"))
              expect_that(plotReducedDim(
                  example_sce, "PCA", colour_by = "Cell_Cycle", 
                  size_by = "Treatment"), is_a("ggplot"))
              expect_that(plotReducedDim(example_sce, "PCA", ncomponents = 5), 
                          is_a("ggplot"))
              expect_that(plotReducedDim(
                  example_sce, "PCA", ncomponents = 5, colour_by = "Cell_Cycle",
                  shape_by = "Treatment"), is_a("ggplot"))
              expect_that(plotReducedDim(
                  example_sce, "PCA", colour_by = "Gene_0001"), is_a("ggplot"))
              
          })

context("test plotExpression")

test_that("we can produce expression plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)

    expect_that(plotExpression(example_sce, "Gene_0001", x = "Mutation_Status"),
                is_a("ggplot"))
    expect_that(plotExpression(example_sce, c("Gene_0001", "Gene_0004"), 
                               x = "Mutation_Status"), is_a("ggplot"))
    expect_that(plotExpression(example_sce, "Gene_0001", x = "Gene_0002"),
                is_a("ggplot"))
    expect_that(plotExpression(example_sce, c("Gene_0001", "Gene_0004"), 
                               x = "Gene_0002"), is_a("ggplot"))
    expect_that(plotExpression(example_sce, 1:4, "Cell_Cycle"), is_a("ggplot"))
    expect_that(plotExpression(example_sce, 1:4, "Gene_0004"), is_a("ggplot"))
    expect_that(plotExpression(example_sce, 1:4), is_a("ggplot"))
    expect_that(plotExpression(example_sce, 1:4, scales = "fixed"), 
                is_a("ggplot"))
})

test_that("we can plot expression for named genes", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    geneset <- rownames(example_sce)[1:6]
    expect_that(plotExpression(example_sce, geneset, "Cell_Cycle"),
                is_a("ggplot"))
    expect_that(plotExpression(example_sce, geneset, "Gene_0004"), 
                is_a("ggplot"))
})

test_that("plotting expression for an object with non-NULL is_exprs() ", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    geneset <- rownames(example_sce)[1:6]
    assay(example_sce, "is_exprs") <- counts(example_sce) > 0.5
    expect_that(plotExpression(example_sce, geneset, "Cell_Cycle"),
                is_a("ggplot"))
    expect_that(plotExpression(example_sce, geneset, "Gene_0004"), 
                is_a("ggplot"))
})

test_that("we can produce plots showing cells in plate position", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    ## define plate positions
    example_sce$plate_position <- paste0(
        rep(LETTERS[1:5], each = 8), rep(formatC(1:8, width = 2, flag = "0"), 5))
    
    expect_that(plotPlatePosition(example_sce, colour_by = "Cell_Cycle"), 
                is_a("ggplot"))
    expect_that(plotPlatePosition(example_sce, colour_by = "Gene_0004"), 
                is_a("ggplot"))
    ppos <- example_sce$plate_position
    expect_that(plotPlatePosition(example_sce, plate_position = ppos, 
                                  colour_by = "Gene_0004"), is_a("ggplot"))
    example_sce$plate_position <- NULL
    expect_that(plotPlatePosition(example_sce, 
                                  y_position = rep(1:5, each = 8), 
                                  x_position = rep(1:8, 5), 
                                  colour_by = "Gene_0004"), is_a("ggplot"))
    
})

test_that("we can produce plots for metadata", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(plotPhenoData(example_sce, 
                              aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))

    expect_that(plotColData(example_sce, aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))

    expect_that(plotCellData(example_sce, aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))
    
    expect_that(plotPhenoData(example_sce, aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))

    expect_that(plotColData(example_sce, aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))

    expect_that(plotCellData(example_sce, aesth = aes_string(x = "log10(total_counts)",
    y = "total_features", colour = "Mutation_Status")), is_a("ggplot"))

    expect_that(plotFeatureData(
        example_sce, aesth = aes(x = n_cells_counts, y = log10_total_counts)),
        is_a("ggplot"))

    expect_that(plotRowData(
        example_sce, aesth = aes(x = n_cells_counts, y = log10_total_counts)),
        is_a("ggplot"))
    
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
        calculateCPM(example_sce, use.size.factors = FALSE) + 1)
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
    p1 <-  plotExprsVsTxLength(example_sce, "tx_len", show_smooth = TRUE,
                               show_exprs_sd = TRUE)
    expect_that(p1, is_a("ggplot"))
    
    ## using a vector of tx length values
    p1 <- plotExprsVsTxLength(example_sce, rnorm(2000, mean = 5000, sd = 500))
    expect_that(p1, is_a("ggplot"))
    
    ## test errors
    expect_error(plotExprsVsTxLength(example_sce, "foot"), 
                 "should specify a column")
    expect_error(plotExprsVsTxLength(example_sce, 1:10), 
                 "must have length equal")
})

