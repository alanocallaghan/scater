Sys.setenv("R_TESTS" = "")
library(testthat)
library(scater)

test_check("scater")
