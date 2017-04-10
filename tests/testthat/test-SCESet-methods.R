## Testing SCESet methods

test_that("newSCESet works as expected with data type hierarchy", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    ## new SCESet with count data
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with tpm data
    tpm <- calculateCPM(example_sceset, use.size.factors = FALSE)
    expect_warning(example_sceset <- newSCESet(tpmData = tpm, phenoData = pd),
                   "'logExprsOffset' should be set manually")
    expect_that(example_sceset, is_a("SCESet"))
    example_sceset <- newSCESet(tpmData = tpm, phenoData = pd, 
                                logExprsOffset = 1)
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with cpm data
    expect_warning(example_sceset <- newSCESet(cpmData = tpm, phenoData = pd),
                   "'logExprsOffset' should be set manually")
    example_sceset <- newSCESet(cpmData = tpm, phenoData = pd, 
                                logExprsOffset = 1)
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with fpkm data
    expect_warning(example_sceset <- newSCESet(fpkmData = tpm, phenoData = pd),
                   "'logExprsOffset' should be set manually")
    example_sceset <- newSCESet(fpkmData = tpm, phenoData = pd, 
                                logExprsOffset = 1)
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with cpm data and counts
    expect_warning(example_sceset <- 
                       newSCESet(cpmData = tpm, countData = sc_example_counts,
                                phenoData = pd),  "'cpmData' will be ignored")
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with fpkm data and counts
    expect_warning(example_sceset <- 
                       newSCESet(fpkmData = tpm, countData = sc_example_counts,
                                 phenoData = pd),  "'fpkmData' will be ignored")
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with tpm and cpm data
    expect_warning(example_sceset <- 
                       newSCESet(tpmData = tpm, cpmData = tpm,
                                 phenoData = pd),  "'cpmData' will be ignored")
    expect_that(example_sceset, is_a("SCESet"))
    ## new SCESet with tpm and cpm data
    expect_warning(example_sceset <- 
                       newSCESet(tpmData = tpm, fpkmData = tpm,
                                 phenoData = pd),  "'fpkmData' will be ignored")
    expect_that(example_sceset, is_a("SCESet"))

    ## new SCESet with tpm data and counts
    example_sceset <- newSCESet(tpmData = tpm, countData = sc_example_counts,
                                phenoData = pd)
    expect_that(example_sceset, is_a("SCESet"))
    
    
    
})


test_that("subsetting SCESet works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    expect_that(example_sceset[1:500,], is_a("SCESet"))
    expect_that(example_sceset[, 10:35], is_a("SCESet"))
    expect_that(example_sceset[500:1000, 7:27], is_a("SCESet"))
})


test_that("sizeFactors() works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## expect NULL is returnd if size factors haven't been defined
    expect_null(sizeFactors(example_sceset))
    
    sizeFactors(example_sceset, NULL) <- rep(1, ncol(example_sceset))
    
    ## expect output still an SCESet
    expect_that(example_sceset, is_a("SCESet"))
    
    ## size factors should return a vector
    expect_that(sizeFactors(example_sceset), is_a("numeric"))
    
    ## check with a named set of control features
    example_sceset <- calculateQCMetrics(example_sceset, 
                                         feature_controls = list(set1 = 1:40))
    sizeFactors(example_sceset, "set1") <- 2 ^ rnorm(ncol(example_sceset))
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(sizeFactors(example_sceset), is_a("numeric"))
    
    ## reset size factors with NULL
    sizeFactors(example_sceset) <- NULL
    expect_null(sizeFactors(example_sceset))
    
})

test_that("set_exprs and get_exprs work as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    ## expect get_exprs to return a matrix of exprs values
    expect_that(exprs(example_sceset), is_a("matrix"))
    expect_that(get_exprs(example_sceset, "exprs"), is_a("matrix"))
        
    ## expect NULL is returned if an assayData element hasn't been defined
    expect_warning(tmp <- get_exprs(example_sceset, "norm_exprs"), 
                   "'object' does not contain")
    expect_null(tmp)
    
    ## expect we still get an SCESet after using set_exprs
    set_exprs(example_sceset, "new_exprs") <- exprs(example_sceset) / 2
    expect_that(example_sceset, is_a("SCESet"))
    expect_that(example_sceset@assayData[["new_exprs"]], is_a("matrix"))
    expect_that(get_exprs(example_sceset, "new_exprs"), is_a("matrix"))
    set_exprs(example_sceset, "counts") <- NULL
    expect_that(example_sceset, is_a("SCESet"))
    expect_null(counts(example_sceset))
    ## test the check on dimnames
    tmp <- exprs(example_sceset) / 2
    rownames(tmp) <- rev(rownames(tmp))
    expect_error(set_exprs(example_sceset, "new_exprs") <- tmp,
                 "dimnames for new expression matrix")
    colnames(tmp) <- rev(colnames(tmp))
    expect_error(set_exprs(example_sceset, "new_exprs") <- tmp,
                 "dimnames for new expression matrix")
    rownames(tmp) <- NULL
    expect_error(set_exprs(example_sceset, "new_exprs") <- tmp,
                 "dimnames for new expression matrix")
    dimnames(tmp) <- NULL
    expect_error(set_exprs(example_sceset, "new_exprs") <- tmp,
                 "dimnames for new expression matrix")
    
})



test_that("mergeSCESet works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    fd <- new("AnnotatedDataFrame",
              data = data.frame(ID = rownames(sc_example_counts),
                                row.names = rownames(sc_example_counts)))
    example_sceset <- newSCESet(countData = sc_example_counts)
    
    merged_sceset <- mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40])
    expect_that(merged_sceset, is_a("SCESet"))
    expect_identical(pData(merged_sceset), pData(example_sceset))
    
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    merged_sceset <- mergeSCESet(example_sceset[, 1:20],
                                 example_sceset[, 21:40])
    expect_that(merged_sceset, is_a("SCESet"))
    expect_identical(pData(merged_sceset), pData(example_sceset))

    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd,
                                featureData = fd)
    merged_sceset <- mergeSCESet(example_sceset[, 1:20],
                                 example_sceset[, 21:40])
    expect_that(merged_sceset, is_a("SCESet"))

    tmp_sceset <- example_sceset[, 21:40]
    fData(tmp_sceset)$ID <- "nope"
    expect_error(mergeSCESet(example_sceset[, 1:20], tmp_sceset),
                 "not identical")
    
    # checking that selecting 1 pData column works OK
    expect_that(mergeSCESet(example_sceset[, 1:20], example_sceset[, 40], 
                            pdata_cols = 3), is_a("SCESet"))

    # Checking that featureControlInfo is handled properly.
    new_example_sceset <- calculateQCMetrics(example_sceset, 
                                         feature_controls = list(ERCC = 1:40))
    merged_sceset <- mergeSCESet(new_example_sceset[, 1:20],
                                 new_example_sceset[, 21:40])
    expect_identical(featureControlInfo(new_example_sceset), featureControlInfo(merged_sceset))
    expect_warning(merged_sceset <- mergeSCESet(new_example_sceset[, 1:20],
                                                new_example_sceset[, 21:40], fdata_cols=0),
                   "removing undefined feature control set 'ERCC'")
})


test_that("writeSCESet works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    
    expect_error(writeSCESet(example_sceset, "test.h5", type = "nope",
                             overwrite_existing = TRUE), 
                             "HDF5 is the only format")
    
    # writeSCESet(example_sceset, "test.h5", overwrite_existing = TRUE)
    # expect_true(file.exists("test.h5"))
    # file.remove("test.h5")
})


test_that("cellNames works as expected", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    expect_true(is.vector(cellNames(example_sceset)))
    cellNames(example_sceset) <- 1:ncol(example_sceset)
    expect_true(is.vector(cellNames(example_sceset)))
})

