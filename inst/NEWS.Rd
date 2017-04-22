\name{NEWS}
\title{News for Package \pkg{scater}}

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
