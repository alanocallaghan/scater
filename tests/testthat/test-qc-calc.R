## Test functions for QC calculation.
## library(scater); library(testthat); source("setup-sce.R"); source("test-qc-calc.R")

original <- sce

test_that("we can compute standard QC metrics", {
    example_sce <- calculateQCMetrics(original)
    expect_that(example_sce, is_a("SingleCellExperiment"))
    expect_identical(counts(example_sce), counts(original))
    expect_identical(example_sce$Cell, original$Cell)
    expect_identical(colData(original), colData(example_sce)[,colnames(colData(original))])

    # Testing total metrics for cells.
    expect_equal(example_sce$total_counts, unname(colSums(counts(example_sce))))
    expect_equal(example_sce$total_features_by_counts, unname(colSums(counts(example_sce) > 0)))
    expect_equal(example_sce$log10_total_counts, log10(example_sce$total_counts+1))
    expect_equal(example_sce$log10_total_features_by_counts, log10(example_sce$total_features_by_counts+1))

    # Testing percentage metrics for cells.
    for (i in seq_len(ncol(example_sce))) { 
        cur_counts <- counts(example_sce)[,i]
        o <- order(cur_counts, decreasing=TRUE)
        lib_size <- sum(cur_counts)

        for (x in c(50, 100, 200, 500)) { 
            chosen <- o[seq_len(x)]
            expect_equal(example_sce[[paste0("pct_counts_in_top_", x, "_features")]][i],
                sum(cur_counts[chosen])/lib_size * 100) 
        }
    }

    # Testing mean metrics for genes.
    expect_equal(rowData(example_sce)$mean_counts, unname(rowMeans(counts(example_sce))))
    expect_equal(rowData(example_sce)$total_counts, unname(rowSums(counts(example_sce))))
    expect_equal(rowData(example_sce)$n_cells_by_counts, unname(rowSums(counts(example_sce) > 0)))
    expect_equal(rowData(example_sce)$pct_dropout_by_counts, unname(100 * rowMeans(counts(example_sce)==0)))

    expect_equal(rowData(example_sce)$log10_mean_counts, unname(log10(rowData(example_sce)$mean_counts+1)))
    expect_equal(rowData(example_sce)$log10_total_counts, unname(log10(rowData(example_sce)$total_counts+1)))
})

test_that("we can compute standard QC metrics with feature controls", {
    expect_error(example_sce <- calculateQCMetrics(original, feature_controls = 1:20), "feature_controls should be named")

    example_sce <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20))
    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Testing total metrics for cells.
    expect_identical(which(rowData(example_sce)$is_feature_control_set1), 1:20)
    expect_equal(unname(colSums(counts(example_sce))), example_sce$total_counts)

    sub_counts <- counts(example_sce)[rowData(example_sce)$is_feature_control_set1,]
    expect_equal(example_sce$total_counts_set1, unname(colSums(sub_counts)))
    expect_equal(example_sce$total_features_by_counts_set1, unname(colSums(sub_counts > 0)))
    expect_equal(example_sce$pct_counts_set1, example_sce$total_counts_set1/example_sce$total_counts * 100)

    expect_equal(example_sce$log10_total_counts_set1, log10(example_sce$total_counts_set1+1))
    expect_equal(example_sce$log10_total_features_by_counts_set1, log10(example_sce$total_features_by_counts_set1+1))

    # Same output with logical vectors.
    set1_con <- logical(nrow(original))
    set1_con[1:20] <- TRUE
    example_sce2 <- calculateQCMetrics(original, feature_controls = list(set1 = set1_con))
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
        expect_equal(example_sce[[NAMER("total_features_by_counts")]], unname(colSums(sub_counts > 0)))
        expect_equal(example_sce[[NAMER("pct_counts")]], example_sce[[NAMER("total_counts")]]/example_sce$total_counts * 100)
    
        expect_equal(example_sce[[NAMER("log10_total_counts")]], log10(example_sce[[NAMER("total_counts")]]+1))
        expect_equal(example_sce[[NAMER("log10_total_features_by_counts")]], log10(example_sce[[NAMER("total_features_by_counts")]]+1))
    }
})

test_that("we can compute standard QC metrics with cell controls", {
    expect_error(example_sce <- calculateQCMetrics(original, cell_controls = 1:20), "cell_controls should be named")

    example_sce <- calculateQCMetrics(original, cell_controls = list(set1 = 1:20))
    expect_that(example_sce, is_a("SingleCellExperiment"))

    # Testing total metrics for cells.
    expect_identical(which(example_sce$is_cell_control_set1), 1:20)

    sub_counts <- counts(example_sce)[,example_sce$is_cell_control_set1]
    expect_equal(rowData(example_sce)$mean_counts_set1, unname(rowMeans(sub_counts)))
    expect_equal(rowData(example_sce)$total_counts_set1, unname(rowSums(sub_counts)))
    expect_equal(rowData(example_sce)$n_cells_by_counts_set1, unname(rowSums(sub_counts > 0)))
    expect_equal(rowData(example_sce)$pct_dropout_by_counts_set1, unname(100 * rowMeans(sub_counts==0)))

    expect_equal(rowData(example_sce)$log10_mean_counts_set1, unname(log10(rowData(example_sce)$mean_counts_set1+1)))
    expect_equal(rowData(example_sce)$log10_total_counts_set1, unname(log10(rowData(example_sce)$total_counts_set1+1)))

    # Same output with logical vectors.
    set1_con <- logical(nrow(original))
    set1_con[1:20] <- TRUE
    example_sce2 <- calculateQCMetrics(original, cell_controls = list(set1 = set1_con))
    expect_equal(example_sce2, example_sce)

    # Testing behaviour with multiple cell controls.
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
        expect_equal(rowData(example_sce)[[NAMER("n_cells_by_counts")]], unname(rowSums(sub_counts > 0)))
        expect_equal(rowData(example_sce)[[NAMER("pct_dropout_by_counts")]], unname(100 * rowMeans(sub_counts==0)))

        expect_equal(rowData(example_sce)[[NAMER("log10_mean_counts")]], unname(log10(rowData(example_sce)[[NAMER("mean_counts")]]+1)))
        expect_equal(rowData(example_sce)[[NAMER("log10_total_counts")]], unname(log10(rowData(example_sce)[[NAMER("total_counts")]]+1)))
    }
})

#######################################################################
# Responds to special settings: 

library(Matrix)
test_that("we can compute standard QC metrics on sparse counts matrix", {
    counts(original) <- as(counts(original), "dgCMatrix")

    expect_error(example_sce <- calculateQCMetrics(original, feature_controls = 1:20),
        "feature_controls should be named")

    sparsified <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20))
    expect_that(sparsified, is_a("SingleCellExperiment"))

    # Checking the values against a reference.
    ref <- original 
    counts(ref) <- as.matrix(counts(original))
    ref <- calculateQCMetrics(ref, feature_controls=list(set1=1:20))
    expect_equal(rowData(ref), rowData(sparsified))
    expect_equal(colData(ref), colData(sparsified))
})

test_that("we can compute standard QC metrics across multiple cores", {
    out <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20), cell_controls=list(whee=5:10))
    out2 <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20), cell_controls=list(whee=5:10), BPPARAM=MulticoreParam(2))
    expect_identical(out, out2)
    out3 <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20), cell_controls=list(whee=5:10), BPPARAM=SnowParam(3))
    expect_identical(out, out3)

    out <- calculateQCMetrics(original, cell_controls=list(whee=5:10, yippe=20:30))
    out2 <- calculateQCMetrics(original, cell_controls=list(whee=5:10, yippe=20:30), BPPARAM=MulticoreParam(2))
    expect_identical(out, out2)
    out3 <- calculateQCMetrics(original, cell_controls=list(whee=5:10, yippe=20:30), BPPARAM=SnowParam(3))
    expect_identical(out, out3)

    out <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20, whee=1000:1010))
    out2 <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20, whee=1000:1010), BPPARAM=MulticoreParam(2))
    expect_identical(out, out2)
    out3 <- calculateQCMetrics(original, feature_controls = list(set1 = 1:20, whee=1000:1010), BPPARAM=SnowParam(3))
    expect_identical(out, out3)
})

#######################################################################
# Repeating with multiple feature and cell controls:

multi_qc <- calculateQCMetrics(original, 
    feature_controls = list(controls1 = 1:20, controls2 = 500:1000),
    cell_controls = list(set1 = 1:5, set2 = 10:20))

test_that("we can compute standard QC metrics with multiple sets of feature and cell controls", {
    expect_that(multi_qc, is_a("SingleCellExperiment"))
   
    ref1 <- calculateQCMetrics(original, feature_controls = list(controls1 = 1:20, controls2 = 500:1000))
    expect_equal(colData(ref1)[,colnames(colData(ref1))!="is_cell_control"],
        colData(multi_qc)[,!grepl("is_cell_control", colnames(colData(multi_qc)))])

    ref2 <- calculateQCMetrics(original, cell_controls = list(set1 = 1:5, set2 = 10:20))
    expect_equal(rowData(ref2)[,colnames(rowData(ref2))!="is_feature_control"],
        rowData(multi_qc)[,!grepl("is_feature_control", colnames(rowData(multi_qc)))])

    # Checking for correct overwriting of elements.
    reref <- calculateQCMetrics(multi_qc, feature_controls = list(controls1 = 20:50),
        cell_controls = list(set2 = 5:10))

    expect_identical(which(rowData(reref)$is_feature_control_controls1), 20:50)
    expect_equal(reref$total_counts_controls1, unname(colSums(counts(reref)[20:50,])))

    expect_identical(which(reref$is_cell_control_set2), 5:10)
    expect_equal(rowData(reref)$total_counts_set2, unname(rowSums(counts(reref)[,5:10])))
})

test_that("we can compute standard QC metrics with spike-in sets", {
    alt_object <- original
    isSpike(alt_object, "controls1") <- 1:20
    isSpike(alt_object, "controls2") <- 500:1000
    
    # Should get the same results without specifying the spike-in set.
    another_qc <- calculateQCMetrics(alt_object, cell_controls = list(set1 = 1:5, set2 = 10:20))
    expect_equal(colData(another_qc), colData(multi_qc))
    expect_equal(rowData(another_qc), rowData(multi_qc))

    # Specified feature controls should override internal spike-ins.
    isSpike(alt_object, "controls1") <- 500:600
    expect_warning(another_qc <- calculateQCMetrics(alt_object, 
        feature_controls = list(controls1 = 1:20),
        cell_controls = list(set1 = 1:5, set2 = 10:20)), "spike-in set")
    expect_equal(colData(another_qc), colData(multi_qc))
    expect_equal(rowData(another_qc), rowData(multi_qc))

    expect_warning(another_qc <- calculateQCMetrics(alt_object, 
        feature_controls = list(controls1 = 1:20, controls2 = 500:1000),
        cell_controls = list(set1 = 1:5, set2 = 10:20),
        use_spike=FALSE), NA)
})

test_that("we can compute standard QC metrics with the compact format", {
    compact <- calculateQCMetrics(original, 
        feature_controls = list(controls1 = 1:20, controls2 = 500:1000),
        cell_controls = list(set1 = 1:5, set2 = 10:20), compact=TRUE)

    # Checking the column data.
    expect_equal(multi_qc$total_counts, compact$scater_qc$all$total_counts)  
    expect_equal(multi_qc$total_features_by_counts, compact$scater_qc$all$total_features_by_counts)  
    expect_equal(multi_qc$total_counts_controls1, compact$scater_qc$feature_control_controls1$total_counts)  
    expect_equal(multi_qc$total_counts_controls2, compact$scater_qc$feature_control_controls2$total_counts)  

    expect_identical(multi_qc$is_cell_control, compact$scater_qc$is_cell_control)
    expect_identical(multi_qc$is_cell_control_controls1, compact$scater_qc$is_cell_control_controls1)
    expect_identical(multi_qc$is_cell_control_controls2, compact$scater_qc$is_cell_control_controls2)

    # Checking the row data.
    expect_equal(rowData(multi_qc)$mean_counts, rowData(compact)$scater_qc$all$mean_counts) 
    expect_equal(rowData(multi_qc)$mean_counts_set1, rowData(compact)$scater_qc$cell_control_set1$mean_counts)  
    expect_equal(rowData(multi_qc)$mean_counts_set2, rowData(compact)$scater_qc$cell_control_set2$mean_counts) 

    expect_identical(rowData(multi_qc)$is_feature_control, rowData(compact)$scater_qc$is_feature_control)
    expect_identical(rowData(multi_qc)$is_feature_control_set1, rowData(compact)$scater_qc$is_feature_control_set1)
    expect_identical(rowData(multi_qc)$is_feature_control_set2, rowData(compact)$scater_qc$is_feature_control_set2)

    # Checking for correct overwriting of elements.
    recompact <- calculateQCMetrics(compact, feature_controls = list(controls1 = 20:50),
        cell_controls = list(set2 = 5:10), compact=TRUE)

    expect_identical(which(rowData(recompact)$scater_qc$is_feature_control_controls1), 20:50)
    expect_equal(recompact$scater_qc$feature_control_controls1$total_counts, 
        unname(colSums(counts(recompact)[20:50,])))

    expect_identical(which(recompact$scater_qc$is_cell_control_set2), 5:10)
    expect_equal(rowData(recompact)$scater_qc$cell_control_set2$total_counts, 
        unname(rowSums(counts(recompact)[,5:10])))
})

#######################################################################
# Checking for other options and silly inputs.

test_that("calculateQCMetrics works with silly inputs", {
    # Responds to alternative inputs.
    blah <- sce
    assayNames(blah) <- "whee"
    expect_error(calculateQCMetrics(blah), "counts")
    expect_error(calculateQCMetrics(blah, exprs_values="whee"), NA)

    # Doesn't throw up with empty objects.
    out <- calculateQCMetrics(sce[0,])
    expect_identical(dim(out), c(0L, ncol(out)))
    out <- calculateQCMetrics(sce[,0])
    expect_identical(dim(out), c(nrow(out), 0L))

    # Handles empty feature/cell control sets.
    out <- calculateQCMetrics(sce, feature_controls=list(whee=integer(0)), cell_controls=list(blah=integer(0)))
    expect_true(!any(rowData(out)$is_feature_control))
    expect_true(!any(out$is_cell_control))

    # Handles empty 'top' specification.
    ref <- calculateQCMetrics(sce)
    expect_true(any(grepl("pct_counts_in_top", colnames(colData(ref)))))
    out <- calculateQCMetrics(sce, percent_top=numeric(0))
    expect_false(any(grepl("pct_counts_in_top", colnames(colData(out)))))
})


