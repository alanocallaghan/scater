# Creating a common SCE for further operations.
data("sc_example_counts")
data("sc_example_cell_info")
sce <- SingleCellExperiment(
    assays = list(counts = sc_example_counts), 
    colData = sc_example_cell_info)

# Generating a normalized SCE for specific methods.
normed <- logNormCounts(sce)

# Because SnowParam() is too slow, yet MulticoreParam() fails on Windows.
# See discussion at https://github.com/Bioconductor/BiocParallel/issues/98.
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
