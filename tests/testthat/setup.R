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
        BiocParallel::SnowParam(nworkers)
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}

# Using an exact algorithm to avoid needing to set the seed for reproducibility.
options(BiocSingularParam.default=BiocSingular::ExactParam())

# Adding a test to flush out any uncontrolled parallelization.
library(BiocParallel)
failgen <- setRefClass("FailParam",
    contains="BiocParallelParam",
    fields=list(),
    methods=list())

FAIL <- failgen()
# register(FAIL) # TODO: once DelayedArray's %*% fix gets in.

library(DelayedArray)
setAutoBPPARAM(FAIL)
