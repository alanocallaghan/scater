.Deprecated(msg="'data(sc_example_counts)' is deprecated.\nUse mockSCE() instead.")
sce_example_counts <- (function() {
    bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
    bpath <- BiocFileCache::bfcrpath(bfc, "https://github.com/davismcc/scater/raw/a3d87aa61ebc050374e2ab1ec29a428d39744281/data/sc_example_counts.RData")
    env <- new.env()
    load(bpath, envir=env)
    env$sc_example_counts
})()
