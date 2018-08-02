# Creating a common SCE for further operations.
data("sc_example_counts")
data("sc_example_cell_info")
sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)

# Generating a normalized SCE for specific methods.
suppressWarnings(normed <- normalize(sce))

