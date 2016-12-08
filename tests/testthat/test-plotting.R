## tests for plotting functions

context("test plotPCA, plotTSNE and plotDiffusionMap")

test_that("we can produce PCA plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotPCA(example_sceset), is_a("ggplot"))
})

test_that("we can produce t-SNE plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotTSNE(example_sceset), is_a("ggplot"))
})

test_that("we can produce Diffusion Map plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
              example_sceset <- newSCESet(countData = sc_example_counts,
                                          phenoData = pd)
              example_sceset <- calculateQCMetrics(example_sceset)

              expect_that(plotDiffusionMap(example_sceset), is_a("ggplot"))
          })

test_that("we can produce MDS plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
              example_sceset <- newSCESet(countData = sc_example_counts,
                                          phenoData = pd)
              example_sceset <- calculateQCMetrics(example_sceset)
              cellDist(example_sceset) <- as.matrix(
                  dist(t(exprs(example_sceset))))

              expect_that(plotMDS(example_sceset), is_a("ggplot"))
              expect_that(plotMDS(example_sceset,
                                          colour_by = "Cell_Cycle",
                                          shape_by = "Treatment",
                                          size_by = "Mutation_Status"),
                            is_a("ggplot"))
          })

context("test plotExpression")

test_that("we can produce expression plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)

    expect_that(plotExpression(example_sceset, 1:4, "Cell_Cycle"), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4, "Gene_0004"), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4), is_a("ggplot"))
    expect_that(plotExpression(example_sceset, 1:4, scales = "fixed"), 
                is_a("ggplot"))
})

test_that("we can produce plots showing cells in plate position", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    ## define plate positions
    example_sceset$plate_position <- paste0(
        rep(LETTERS[1:5], each = 8), rep(formatC(1:8, width = 2, flag = "0"), 5))
    
    expect_that(plotPlatePosition(example_sceset, colour_by = "Cell_Cycle"), 
                is_a("ggplot"))
    expect_that(plotPlatePosition(example_sceset, colour_by = "Gene_0004"), 
                is_a("ggplot"))
    ppos <- example_sceset$plate_position
    expect_that(plotPlatePosition(example_sceset, plate_position = ppos, 
                                  colour_by = "Gene_0004"), is_a("ggplot"))
    example_sceset$plate_position <- NULL
    expect_that(plotPlatePosition(example_sceset, 
                                  y_position = rep(1:5, each = 8), 
                                  x_position = rep(1:8, 5), 
                                  colour_by = "Gene_0004"), is_a("ggplot"))
    

})

test_that("plotExprsVsTxLength works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    fd <- new("AnnotatedDataFrame", data = 
                  data.frame(gene_id = rownames(sc_example_counts), 
                             feature_id = paste("feature", rep(1:500, each = 4), 
                                                sep = "_"),
                             median_tx_length = rnorm(2000, mean = 5000, sd = 500)))
    rownames(fd) <- rownames(sc_example_counts)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd,
                                featureData = fd)
    fData(example_sceset)$group <- rep(1:4, each = 500)
    
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length")
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length", 
                              show_smooth = TRUE)
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length", 
                              show_smooth = TRUE, show_exprs_sd = TRUE)
    expect_that(p1, is_a("ggplot"))
    
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length",
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, colour_by = "group")
    expect_that(p1, is_a("ggplot"))
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length", 
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, size_by = "group")
    expect_that(p1, is_a("ggplot"))
    
    fData(example_sceset)$group <- rep(letters[1:4], each = 500)
    p1 <- plotExprsVsTxLength(example_sceset, "median_tx_length", 
                              show_smooth = TRUE,
                              show_exprs_sd = FALSE, shape_by = "group")
    expect_that(p1, is_a("ggplot"))
    
    
    ## using matrix of tx length values in assayData(object)
    mat <- matrix(rnorm(ncol(example_sceset) * nrow(example_sceset), mean = 5000,
                        sd = 500), nrow = nrow(example_sceset))
    expect_error(set_exprs(example_sceset, "tx_len") <- mat, 
                 "dimnames for new expression matrix")
    
    dimnames(mat) <- dimnames(example_sceset)
    set_exprs(example_sceset, "tx_len") <- mat
    p1 <-  plotExprsVsTxLength(example_sceset, "tx_len", show_smooth = TRUE,
                               show_exprs_sd = TRUE)
    expect_that(p1, is_a("ggplot"))
    
    ## using a vector of tx length values
    p1 <- plotExprsVsTxLength(example_sceset, rnorm(2000, mean = 5000, sd = 500))
    expect_that(p1, is_a("ggplot"))
    
    ## test errors
    expect_error(plotExprsVsTxLength(example_sceset, "foot"), 
                 "should specify a column")
    expect_error(plotExprsVsTxLength(example_sceset, 1:10), 
                 "must have length equal")
})

