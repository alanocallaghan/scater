# Tests for dplyr-style verbs
# library(scater); library(testthat); source("setup-sce.R"); source("test-verbs.R")

example_sce <- sce
snames <- paste0("sample", 1:ncol(example_sce))
colnames(example_sce) <- snames

# Test mutate -------------------------------------------------------------

test_that("Mutate correctly adds column name to pData", {
 
  example_sce <- mutate(example_sce, New_CC = Cell_Cycle)
  
  expect_that(example_sce, is_a("SingleCellExperiment"))
  expect_true("New_CC" %in% colnames(colData(example_sce)))  
  expect_equal(example_sce$Cell_Cycle, example_sce$New_CC)
  expect_equal(colnames(example_sce), snames)
})


# Test rename ------------------------------------------------------------

test_that("Rename correctly renames columns", {
  
  old_cc <- example_sce$Cell_Cycle
  example_sce <- rename(example_sce, Cell_Phase = Cell_Cycle)
  
  expect_that(example_sce, is_a("SingleCellExperiment"))
  expect_true("Cell_Phase" %in% colnames(colData(example_sce)))
  expect_false("Cell_Cycle" %in% colnames(colData(example_sce)))
  expect_equal(example_sce$Cell_Phase, old_cc)
  expect_equal(colnames(example_sce), snames)
})


# Test filter -------------------------------------------------------------

test_that("Filter correctly chooses cells", {
  example_sce <- filter(example_sce, Cell_Cycle == "G0")
  
  expect_that(example_sce, is_a("SingleCellExperiment"))
  expect_equal(example_sce$Cell_Cycle, rep("G0", ncol(example_sce)))
})


# Test arrange ------------------------------------------------------------

test_that("Arrange correctly orders", {
  treatment <- example_sce$Treatment
  n_each <- table(treatment)
  arranged <- c(rep("treat1", n_each[1]), rep("treat2", n_each[2]))
  
  example_sce <- arrange(example_sce, Treatment)

  expect_that(example_sce, is_a("SingleCellExperiment"))
  expect_equal(example_sce$Treatment, arranged)
})
