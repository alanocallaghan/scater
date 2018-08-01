# Creating a common SCE for further operations.
data("sc_example_counts")
data("sc_example_cell_info")
sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)

# Generating a normalized SCE for specific methods.
suppressWarnings(normed <- normalize(sce))

# Generating sparse SCEs.
library(Matrix)
sce_sparse <- sce
counts(sce_sparse) <- as(counts(sce_sparse), "dgCMatrix")
suppressWarnings(normed_sparse <- normalize(sce_sparse))
