# context("test magic")
# 
# test_that("we can apply magic to example data", {
#     data("sc_example_counts")
#     data("sc_example_cell_info")
#     pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#     example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#     example_sceset <- example_sceset[rowSums(counts(example_sceset)) > 0.5, ]
#     mgc <- .magic(example_sceset, exprs_values = "counts", power = 6,
#                  logged_data = FALSE, rescale = 0.9, k = 30, n_eigs = 20,
#                  n_local = 10)
#     mgc_norescale <- .magic(example_sceset, power = 6, k = 30, n_eigs = 20,
#                            rescale = 0, n_local = 10)
# 
#     expect_that(mgc, is_a("matrix"))
#     expect_that(mgc_norescale, is_a("matrix"))
#     expect_true(all(!is.na(mgc)))
#     expect_true(all(!is.na(mgc_norescale)))
#     expect_true(all(mgc >= 0))
# })
# 
# 
