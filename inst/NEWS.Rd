\name{NEWS}
\title{News for Package \pkg{scater}}

\section{Changes in version 1.7.16, Bioconductor 3.7 Release}{
  \itemize{
    \item Refactored calculateQCMetrics() to ignore potential non-linearity, rank genes by highest expression, 
        rename automatically generated union sets, allow for output of a compact format. 
     \item Refactored all plotting functions to allow access to nested fields in the colData() or rowData(), by supplying a character vector.
     \item Refactored plotTSNE(), plotPCA(), etc. to dispatch to the calculation functions (e.g., runTSNE(), runPCA()), with argument checks.
     \item Refactored plotColData() and plotRowData() to use the same argument types as other functions rather than aes= input.
     \item Removed all plotting functions that do not operate on SingleCellExperiment objects.
     \item Deprecated read10xResults(), downsampleCounts() in favour of methods from the DropletUtils package.
     \item Deprecated scater_gui() in favour of methods from the iSEE package.
     \item Deprecated normalizeExprs() as this function made very little sense.
     \item Added plotHeatmap() function, for easy plotting of heatmaps. 
     \item Added librarySizeFactors() function, to compute size factors from library sizes.
     \item Added by_exprs_values= argument to many plotting functions, to distinguish direct plotting of expression values from their use in aesthetics.
     \item Renamed arguments in plotHighestExprs(), plotExprsVsTxLength(), plotExprsFreqVsMean() for greater clarity. 
  }
}

\section{Changes in version 1.5.11, Bioconductor 3.6 Release (September 2017)}{
  \itemize{
    \item Complete refactoring of the package to use the SingleCellExperiment
    class
  }
}


\section{Changes in version 1.3.49, Bioconductor 3.5 Release (April 2017)}{
  \itemize{
    \item plotRLE() function to make relative log expression plots to assess and
    compare normalizations
    \item Refactored newSCESet() with defined hierarchy of data types
    \item read10XResults() to read in results from 10x Chromium CellRanger output
    \item Refined QC metrics
    \item Bug fixes, efficiency improvements and more tests
  }
}


\section{Changes in version 1.1.26, Bioconductor 3.4 Release (October 2016)}{
  \itemize{
    \item Key code ported to C++ for greater computational and memory efficiency
    \item Added support/integration for SC3 package for single-cell clustering
    \item writeSCESet() function added to write SCESets to HDF5 format on disk
    \item mergeSCESet() function to merge two SCESet objects that share features
    \item plotPlatePosition() function to visualise gene expression and cell metadata for cells in their position on a plate
    \item Adding plotExprsVsTxLength() to plot expression against transcript length
    \item Added fitted line and some tweaks to plotExprsFreqVsMean().
    \item Adding support for scaled TPM counts when summarising expression at feature level. 
    \item Adding NULL method for set_exprs(). Adding tests.
    \item Adding import of effective feature lengths with readKallistoResults()
    \item runSalmon() function for running Salmon from within R, following suggestions from Rob Patro.
    \item Added cellNames<- assignment function
    \item Added extra QC metrics
    \item Numerous other bug fixes and minor improvements
  }
}


\section{Changes in version 0.99.3 (2016-02-29)}{
  \itemize{
    \item Package added to Bioconductor
    \item Bioc-submission branch merged with master
  }
}

\section{Changes in version 0.99.2 (2016-02-21)}{
  \itemize{
    \item Package prepared for Bioconductor submission.
  }
}
