# tests for feature pre-processing functions.

context("test feature pre-processing functions")

test_that("we can summarise expression at feature level", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    effective_length <- rep(c(1000, 2000), times = 1000)
    tpm(example_sceset) <- calculateTPM(example_sceset, effective_length, 
                                        calc_from = "counts")
    
    fd <- new("AnnotatedDataFrame", 
              data = data.frame(gene_id = featureNames(example_sceset),
                                feature_id = paste("feature", 
                                                   rep(1:500, each = 4), sep = "_")))
    rownames(fd) <- featureNames(example_sceset)
    fData(example_sceset) <- fd
    
    ## tpm with scaled tpm counts
    example_sceset_summarised <-
        summariseExprsAcrossFeatures(example_sceset, exprs_values = "tpm")
    expect_that(example_sceset_summarised, is_a("SCESet"))    
    example_sceset_summarised <-
        summariseExprsAcrossFeatures(example_sceset, exprs_values = "tpm",
                                     scaled_tpm_counts = FALSE)
    expect_that(example_sceset_summarised, is_a("SCESet"))    

    ## counts 
    example_sceset_summarised <-
        summariseExprsAcrossFeatures(example_sceset, exprs_values = "counts")
    expect_that(example_sceset_summarised, is_a("SCESet"))
    
    ## exprs
    example_sceset_summarised <-
        summariseExprsAcrossFeatures(example_sceset, exprs_values = "exprs")
    expect_that(example_sceset_summarised, is_a("SCESet"))
    
    ## errors
    example_sceset2 <- newSCESet(tpmData = tpm(example_sceset), 
                                 phenoData = pd, featureData = fd)
    expect_error(summariseExprsAcrossFeatures(example_sceset2, exprs_values = "tpm"),
                 "lib_size argument")
    expect_error(summariseExprsAcrossFeatures(example_sceset2, exprs_values = "tpm",
                                              lib_size = 1:10),
                 "lib_size argument must have length equal")

})
