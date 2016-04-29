# Tests for dplyr-style verbs


# Test mutate -------------------------------------------------------------

test_that("Mutate correctly adds column name to pData", {
  data("sc_example_counts")
  data("sc_example_cell_info")
  pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
  example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
  example_sceset <- mutate(example_sceset, New_CC = Cell_Cycle)
  
  expect_that(example_sceset, is_a("SCESet"))
  expect_true("New_CC" %in% varLabels(example_sceset))  
  expect_equal(example_sceset$Cell_Cycle, example_sceset$New_CC)
})


# Test rename ------------------------------------------------------------

test_that("Rename correctly renames columns", {
  data("sc_example_counts")
  data("sc_example_cell_info")
  pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
  example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
  
  old_cc <- example_sceset$Cell_Cycle
  example_sceset <- rename(example_sceset, Cell_Phase = Cell_Cycle)
  
  expect_that(example_sceset, is_a("SCESet"))
  expect_true("Cell_Phase" %in% varLabels(example_sceset))
  expect_false("Cell_Cycle" %in% varLabels(example_sceset))
  expect_equal(example_sceset$Cell_Phase, old_cc)
})


# Test filter -------------------------------------------------------------

test_that("Filter correctly chooses cells", {
  data("sc_example_counts")
  data("sc_example_cell_info")
  pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
  example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
  
  example_sceset <- filter(example_sceset, Cell_Cycle == "G0")
  
  expect_that(example_sceset, is_a("SCESet"))
  expect_equal(example_sceset$Cell_Cycle, rep("G0", ncol(example_sceset)))
})


# Test arrange ------------------------------------------------------------

test_that("Arrange correctly orders", {
  data("sc_example_counts")
  data("sc_example_cell_info")
  pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
  example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
  
  treatment <- example_sceset$Treatment
  n_each <- table(treatment)
  arranged <- c(rep("treat1", n_each[1]), rep("treat2", n_each[2]))
  
  example_sceset <- arrange(example_sceset, Treatment)

  expect_that(example_sceset, is_a("SCESet"))
  expect_equal(example_sceset$Treatment, arranged)
})
