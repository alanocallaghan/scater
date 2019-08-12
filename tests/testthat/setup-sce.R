# Creating a common SCE for further operations.
set.seed(100)
sce <- mockSCE()

# Killing alternative Experiments for simplicity.
altExps(sce) <- NULL

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
