## Test functions for QC

context("test controls functionality")

test_that("we can compute standard QC metrics", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce)
    
    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Testing total metrics for cells.
    expect_equal(example_sce$total_counts, unname(colSums(counts(example_sce))))
    expect_equal(example_sce$total_features, unname(colSums(counts(example_sce) > 0)))
    expect_equal(example_sce$log10_total_counts, log10(example_sce$total_counts+1))
    expect_equal(example_sce$log10_total_features, log10(example_sce$total_features+1))

    # Testing percentage metrics for cells.
    for (i in seq_len(ncol(example_sce))) { 
        cur_counts <- counts(example_sce)[,i]
        o <- order(cur_counts, decreasing=TRUE)
        lib_size <- sum(cur_counts)

        for (x in c(50, 100, 200, 500)) { 
            chosen <- o[seq_len(x)]
            expect_equal(example_sce[[paste0("pct_counts_top_", x, "_features")]][i],
                         sum(cur_counts[chosen])/lib_size * 100) 
        }
    }

    # Testing mean metrics for genes.
    expect_equal(rowData(example_sce)$mean_counts, unname(rowMeans(counts(example_sce))))
    expect_equal(rowData(example_sce)$total_counts, unname(rowSums(counts(example_sce))))
    expect_equal(rowData(example_sce)$n_cells_counts, unname(rowSums(counts(example_sce) > 0)))
    expect_equal(rowData(example_sce)$pct_dropout_counts, unname(100 * rowMeans(counts(example_sce)==0)))

    expect_equal(rowData(example_sce)$rank_counts, unname(rank(-rowData(example_sce)$mean_counts)))
    expect_equal(rowData(example_sce)$log10_mean_counts, unname(log10(rowData(example_sce)$mean_counts+1)))
    expect_equal(rowData(example_sce)$log10_total_counts, unname(log10(rowData(example_sce)$total_counts+1)))
})


test_that("we can compute standard QC metrics with feature controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)

    example_sce <- original
    expect_error(
        example_sce <- calculateQCMetrics(example_sce, feature_controls = 1:20),
        "feature_controls should be named")

    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:20))
    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Testing total metrics for cells.
    expect_identical(which(rowData(example_sce)$is_feature_control_set1), 1:20)
    expect_equal(unname(colSums(counts(example_sce))), example_sce$total_counts)

    sub_counts <- counts(example_sce)[rowData(example_sce)$is_feature_control_set1,]
    expect_equal(example_sce$total_counts_set1, unname(colSums(sub_counts)))
    expect_equal(example_sce$total_features_set1, unname(colSums(sub_counts > 0)))
    expect_equal(example_sce$pct_counts_set1, example_sce$total_counts_set1/example_sce$total_counts * 100)

    expect_equal(example_sce$log10_total_counts_set1, log10(example_sce$total_counts_set1+1))
    expect_equal(example_sce$log10_total_features_set1, log10(example_sce$total_features_set1+1))

    # Same output with logical vectors.
    example_sce2 <- calculateQCMetrics(
        example_sce, 
        feature_controls = list(set1 = c(rep(c(TRUE, FALSE),
                                             c(20, nrow(example_sce) - 20)))))
    expect_equal(example_sce2, example_sce)

    # Testing behaviour with multiple feature controls.
    example_sce <- original
    multi_controls <- list(controls1 = 1:20, controls2 = 500:1000)
    example_sce <- calculateQCMetrics(example_sce, feature_controls = multi_controls)
    expect_equal(unname(colSums(counts(example_sce))), example_sce$total_counts)

    all_hits <- c(multi_controls, 
            list(feature_control=unlist(multi_controls),
                 endogenous=setdiff(seq_len(nrow(original)), unlist(multi_controls))))
    expect_identical(which(rowData(example_sce)$is_feature_control), unname(all_hits$feature_control))

    for (set in names(all_hits)) {
        sub_counts <- counts(example_sce)[all_hits[[set]],]

        NAMER <- function(x) { paste0(x, "_", set) }
        expect_equal(example_sce[[NAMER("total_counts")]], unname(colSums(sub_counts)))
        expect_equal(example_sce[[NAMER("total_features")]], unname(colSums(sub_counts > 0)))
        expect_equal(example_sce[[NAMER("pct_counts")]], example_sce[[NAMER("total_counts")]]/example_sce$total_counts * 100)
    
        expect_equal(example_sce[[NAMER("log10_total_counts")]], log10(example_sce[[NAMER("total_counts")]]+1))
        expect_equal(example_sce[[NAMER("log10_total_features")]], log10(example_sce[[NAMER("total_features")]]+1))
    }
})

test_that("we can compute standard QC metrics with cell controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)

    example_sce <- original
    expect_error(
        example_sce <- calculateQCMetrics(example_sce, cell_controls = 1:20),
        "cell_controls should be named")

    example_sce <- calculateQCMetrics(example_sce, 
                                      cell_controls = list(set1 = 1:20))
    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Testing total metrics for cells.
    expect_identical(which(example_sce$is_cell_control_set1), 1:20)

    sub_counts <- counts(example_sce)[,example_sce$is_cell_control_set1]
    expect_equal(rowData(example_sce)$mean_counts_set1, unname(rowMeans(sub_counts)))
    expect_equal(rowData(example_sce)$total_counts_set1, unname(rowSums(sub_counts)))
    expect_equal(rowData(example_sce)$n_cells_counts_set1, unname(rowSums(sub_counts > 0)))
    expect_equal(rowData(example_sce)$pct_dropout_counts_set1, unname(100 * rowMeans(sub_counts==0)))

    expect_equal(rowData(example_sce)$rank_counts_set1, unname(rank(-rowData(example_sce)$mean_counts_set1)))
    expect_equal(rowData(example_sce)$log10_mean_counts_set1, unname(log10(rowData(example_sce)$mean_counts_set1+1)))
    expect_equal(rowData(example_sce)$log10_total_counts_set1, unname(log10(rowData(example_sce)$total_counts_set1+1)))

     #Same output with logical vectors.
    example_sce2 <- calculateQCMetrics(
        example_sce, 
        cell_controls = list(set1 = c(rep(c(TRUE, FALSE),
                                             c(20, nrow(example_sce) - 20)))))
    expect_equal(example_sce2, example_sce)

    # Testing behaviour with multiple feature controls.
    example_sce <- original
    multi_controls <- list(controls1 = 1:5, controls2 = 10:20)
    example_sce <- calculateQCMetrics(example_sce, cell_controls = multi_controls)

    all_hits <- c(multi_controls, 
            list(cell_control=unlist(multi_controls),
                 non_control=setdiff(seq_len(ncol(original)), unlist(multi_controls))))
    expect_identical(which(example_sce$is_cell_control), unname(all_hits$cell_control))

    for (set in names(all_hits)) {
        sub_counts <- counts(example_sce)[,all_hits[[set]]]

        NAMER <- function(x) { paste0(x, "_", set) }
        expect_equal(rowData(example_sce)[[NAMER("mean_counts")]], unname(rowMeans(sub_counts)))
        expect_equal(rowData(example_sce)[[NAMER("total_counts")]], unname(rowSums(sub_counts)))
        expect_equal(rowData(example_sce)[[NAMER("n_cells_counts")]], unname(rowSums(sub_counts > 0)))
        expect_equal(rowData(example_sce)[[NAMER("pct_dropout_counts")]], unname(100 * rowMeans(sub_counts==0)))

        expect_equal(rowData(example_sce)[[NAMER("rank_counts")]], unname(rank(-rowData(example_sce)[[NAMER("mean_counts")]])))
        expect_equal(rowData(example_sce)[[NAMER("log10_mean_counts")]], unname(log10(rowData(example_sce)[[NAMER("mean_counts")]]+1)))
        expect_equal(rowData(example_sce)[[NAMER("log10_total_counts")]], unname(log10(rowData(example_sce)[[NAMER("total_counts")]]+1)))
    }
})

test_that("we can compute standard QC metrics on sparse counts matrix", {
    original <- read10xResults(system.file("extdata", package = "scater"))

    sce10x <- original
    expect_error(
        example_sce <- calculateQCMetrics(sce10x, feature_controls = 1:20),
        "feature_controls should be named")

    sce10x <- calculateQCMetrics(sce10x, feature_controls = list(set1 = 1:20))
    expect_that(sce10x, is_a("SingleCellExperiment"))

    sce10x <- calculateQCMetrics(
        sce10x, 
        feature_controls = list(set1 = c(rep(TRUE, 20), 
                                         rep(FALSE, nrow(sce10x) - 20))))
    expect_that(sce10x, is_a("SingleCellExperiment"))

    # Checking the values against a reference.
    library(Matrix)
    ref <- original
    counts(ref) <- as.matrix(counts(original))
    ref <- calculateQCMetrics(ref, feature_controls=list(set1=1:20))
    expect_equal(rowData(ref), rowData(sce10x))
    expect_equal(colData(ref), colData(sce10x))
})

test_that("we can compute standard QC metrics with multiple sets of feature and cell controls", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)

    example_sce <- calculateQCMetrics(
        original, feature_controls = list(controls1 = 1:20, 
            controls2 = 500:1000),
        cell_controls = list(set1 = 1:5, set2 = 10:20))

    ref1 <- calculateQCMetrics(
        original, feature_controls = list(controls1 = 1:20, 
            controls2 = 500:1000))
    expect_equal(colData(ref1)[,colnames(colData(ref1))!="is_cell_control"],
        colData(example_sce)[,!grepl("is_cell_control", colnames(colData(example_sce)))])

    ref2 <- calculateQCMetrics(
        original, cell_controls = list(set1 = 1:5, set2 = 10:20))
    expect_equal(rowData(ref2)[,colnames(rowData(ref2))!="is_feature_control"],
        rowData(example_sce)[,!grepl("is_feature_control", colnames(rowData(example_sce)))])

    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Checking for correct overwriting of elements.
    reref <- calculateQCMetrics(example_sce, feature_controls = list(controls1 = 20:50),
                                    cell_controls = list(set2 = 5:10))

    expect_identical(which(rowData(reref)$is_feature_control_controls1), 20:50)
    expect_equal(reref$total_counts_controls1, 
                 unname(colSums(counts(reref)[20:50,])))

    expect_identical(which(reref$is_cell_control_set2), 5:10)
    expect_equal(rowData(reref)$total_counts_set2, 
                 unname(rowSums(counts(reref)[,5:10])))
})

test_that("we can compute standard QC metrics with the compactness format", {
    data("sc_example_counts")
    data("sc_example_cell_info")

    original <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)

    ref <- calculateQCMetrics(
        original, feature_controls = list(controls1 = 1:20, 
            controls2 = 500:1000),
        cell_controls = list(set1 = 1:5, set2 = 10:20))
    
    compact <- calculateQCMetrics(
        original, feature_controls = list(controls1 = 1:20, 
            controls2 = 500:1000),
        cell_controls = list(set1 = 1:5, set2 = 10:20), compact=TRUE)

    # Checking the column data.
    expect_equal(ref$total_counts, compact$scater_qc$counts$all$total_counts)  
    expect_equal(ref$total_features, compact$scater_qc$counts$all$total_features)  
    expect_equal(ref$total_counts_controls1, compact$scater_qc$counts$feature_control_controls1$total_counts)  
    expect_equal(ref$total_counts_controls2, compact$scater_qc$counts$feature_control_controls2$total_counts)  

    expect_identical(ref$is_cell_control, compact$scater_qc$is_cell_control)
    expect_identical(ref$is_cell_control_controls1, compact$scater_qc$is_cell_control_controls1)
    expect_identical(ref$is_cell_control_controls2, compact$scater_qc$is_cell_control_controls2)

    # Checking the row data.
    expect_equal(rowData(ref)$mean_counts, rowData(compact)$scater_qc$counts$all$mean_counts) 
    expect_equal(rowData(ref)$mean_counts_set1, rowData(compact)$scater_qc$counts$cell_control_set1$mean_counts)  
    expect_equal(rowData(ref)$mean_counts_set2, rowData(compact)$scater_qc$counts$cell_control_set2$mean_counts) 

    expect_identical(rowData(ref)$is_feature_control, rowData(compact)$scater_qc$is_feature_control)
    expect_identical(rowData(ref)$is_feature_control_set1, rowData(compact)$scater_qc$is_feature_control_set1)
    expect_identical(rowData(ref)$is_feature_control_set2, rowData(compact)$scater_qc$is_feature_control_set2)

    # Checking for correct overwriting of elements.
    recompact <- calculateQCMetrics(compact, feature_controls = list(controls1 = 20:50),
                                    cell_controls = list(set2 = 5:10), compact=TRUE)

    expect_identical(which(rowData(recompact)$scater_qc$is_feature_control_controls1), 20:50)
    expect_equal(recompact$scater_qc$counts$feature_control_controls1$total_counts, 
                 unname(colSums(counts(recompact)[20:50,])))

    expect_identical(which(recompact$scater_qc$is_cell_control_set2), 5:10)
    expect_equal(rowData(recompact)$scater_qc$counts$cell_control_set2$total_counts, 
                 unname(rowSums(counts(recompact)[,5:10])))
})


test_that("computing standard QC metrics with FPKM data fails as expected", {
    gene_df <- data.frame(Gene = rownames(sc_example_counts))
    rownames(gene_df) <- gene_df$Gene
    example_sce <- SingleCellExperiment(
        assays = list(fpkm = sc_example_counts), 
        colData = sc_example_cell_info, rowData = gene_df)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_error(example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(set1 = 1:20)),
        "not in names")
})

test_that("computing standard QC metrics with TPM data fails as expected", {
    gene_df <- data.frame(Gene = rownames(sc_example_counts))
    rownames(gene_df) <- gene_df$Gene
    example_sce <- SingleCellExperiment(
        assays = list(tpm = sc_example_counts), 
        colData = sc_example_cell_info, rowData = gene_df)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_error(example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(set1 = 1:20)),
        "not in names")
})

test_that("failure is as expected for misspecified arg to plotExplanatoryVariables()", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    expect_error(plotExplanatoryVariables(example_sce, "expl"))
})


test_that("failure is as expected for input with zero-variance features", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)
    exprs(example_sce)[1:5,] <- 0
    expect_that(
        plotExplanatoryVariables(example_sce, "density"), is_a("ggplot"))
})


test_that("plotHighestExprs works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)
    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:500))
    expect_that(
        plotHighestExprs(example_sce), 
        is_a("ggplot"))
    expect_that(
        plotHighestExprs(example_sce, col_by_variable = "Mutation_Status"), 
        is_a("ggplot"))
    
    sce.blank <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    expect_warning(
        plotHighestExprs(sce.blank), 
        "not found")
    sce.blank <- calculateQCMetrics(sce.blank)
    expect_that(
        plotHighestExprs(sce.blank), 
        is_a("ggplot"))
})


test_that("plotExplanatoryVariables works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce, 
                                      feature_controls = list(set1 = 1:500))
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)
    drop_genes <- apply(exprs(example_sce), 1, function(x) { var(x) == 0 })
    example_sce <- example_sce[!drop_genes, ]
    example_sce <- calculateQCMetrics(example_sce)
    vars <- colnames(colData(example_sce))[c(2:3, 5:14)]
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars[1]), 
        is_a("ggplot"))
    expect_that(
        plotExplanatoryVariables(example_sce, variables = vars, 
                                 method = "pairs"), 
        is_a("ggplot"))
    err_string <- "Only one variable"
    expect_error(plotExplanatoryVariables(example_sce, variables = vars[1], 
                                          method = "pairs"), err_string)
})


test_that("plotExprsFreqVsMean works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    example_sce <- calculateQCMetrics(example_sce)
    expect_that(plotExprsFreqVsMean(example_sce), is_a("ggplot"))
    
    example_sce <- calculateQCMetrics(
        example_sce, feature_controls = list(controls1 = 1:20,
                                           controls2 = 500:1000),
        cell_controls = list(set_1 = 1:5,
                             set_2 = 31:40))
    expect_that(plotExprsFreqVsMean(example_sce), is_a("ggplot"))
})



test_that("plotRLE works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    example_sce <- SingleCellExperiment(
        assays = list(counts = sc_example_counts), 
        colData = sc_example_cell_info)
    exprs(example_sce) <- log2(
        calculateCPM(example_sce, use_size_factors = FALSE) + 1)

    p <- plotRLE(example_sce, list(exprs = "logcounts", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "minimal")
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Mutation_Status", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, list(exprs = "exprs", counts = "counts"), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    p <- plotRLE(example_sce, 
                 list(exprs = "exprs", counts = counts(example_sce)), 
                 c(TRUE, FALSE), colour_by = "Gene_0004", style = "full",
                 outlier.alpha = 0.1, outlier.shape = NULL, outlier.size = 0)
    expect_that(p, is_a("ggplot"))
    
    expect_error(plotRLE(example_sce, 
                         list("exprs", counts = counts(example_sce)), 
                         c(TRUE, FALSE)), 
                 regexp = "exprs_mats must be a named list")
    
    expect_error(plotRLE(example_sce, 
                         list(exprs = "exprs", counts = counts(example_sce)[, 1:30]), 
                         c(TRUE, FALSE)), 
                 regexp = "Number of cells")
    
    expect_error(plotRLE(example_sce, 
                         list(exprs = "exprs"), style = "blah", 
                         c(TRUE, FALSE)), 
                 regexp = "should be one of")
    
})

