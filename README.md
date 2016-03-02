# scater: single-cell analysis toolkit for expression with R

[//]: # ([![Linux Build Status](https://travis-ci.org/davismcc/scater.svg?branch=master)](https://travis-ci.org/davismcc/scater))
[![Linux Build Status](https://semaphoreci.com/api/v1/davismcc/scater/branches/master/badge.svg)](https://semaphoreci.com/davismcc/scater)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/davismcc/scater?svg=true)](https://ci.appveyor.com/project/davismcc/scater)
[![Coverage Status](https://img.shields.io/codecov/c/github/davismcc/scater/master.svg)](https://codecov.io/github/davismcc/scater?branch=master)

This package contains useful tools for the analysis of single-cell
gene expression data using the statistical software R. The package places an
emphasis on tools for quality control, visualisation and pre-processing of data
before further downstream analysis.

We hope that `scater` fills a useful niche between raw RNA-sequencing
count or transcripts-per-million data and more focused downstream
modelling tools such as
[monocle](http://www.bioconductor.org/packages/release/bioc/html/monocle.html),
[scLVM](http://github.com/PMBio/scLVM),
[SCDE](http://pklab.med.harvard.edu/scde/index.html),
[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html),
[limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
and so on.

Briefly, `scater` enables the following:

1. Automated computation of QC metrics
1. Transcript quantification from read data with pseudo-alignment
2. Data format standardisation
3. Rich visualisations for exploratory analysis
4. Seamless integration into the Bioconductor universe
5. Simple normalisation methods

See below for information about installation, getting started and highlights of the package.

## Installation
The `scater` package has been accepted into Bioconductor!
Thus, the most reliable way to install the package is to use the usual
Bioconductor method:

```{R}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
```

Currently, only the "devel" (i.e. development) version of `scater` is
available through Bioconductor. This means that you will need to be
using Bioconductor devel and the development version of R (R 3.3) in
order to install `scater` from Bioconductor.

The `scater` package will become available as a "release" version in
the next Bioconductor release in April 2016. At this point the release
version of `scater` will work with the release version of R and
Bioconductor, and development will continue in the devel version of
the package.

Alternatively, `scater` can be installed directly from GitHub as
described below. In this case, package that `scater` uses ("depends
on" in R parlance) will not be automatically installed, so you will
have to install the required packages as shown below.

I recommend using Hadley Wickham's `devtools` package to install
`scater` directly from GitHub. If you don't have `devtools` installed,
then install that from CRAN (as shown below) and then run the call to
install `scater`:

**If you are using the development version of R, 3.3:**

```{r}
install.packages("devtools")
devtools::install_github("davismcc/scater", build_vignettes = TRUE)
```

**If you are using the current release version of R, 3.2.3:**
```{r}
devtools::install_github("davismcc/scater", ref = "release-R-3.2", build_vignettes = TRUE)
```

If you find that the above will not install on Linux systems, please
try with the option `build_vignettes = FALSE`. This is a known issue
that we are working to resolve.

Using the most recent version of R is strongly recommended (R 3.2.3 at the time
of writing). Effort has been made to ensure the package works with R >3.0, but
the package has not been tested with R <3.1.1.

There are several other packages from CRAN and Bioconductor that `scater` uses,
so you will need to have these packages installed as well. The CRAN packages
should install automatically when `scater` is installed, but you will need to
install the Bioconductor packages manually.

Not all of the following are strictly necessary, but they enhance the
functionality of `scater` and are good packages in their own right. The commands
below should help with package installations.

CRAN packages:

```{r}
install.packages(c("data.table", "ggplot2", "knitr", "matrixStats", "MASS",
                "plyr", "reshape2", "rjson", "testthat", "viridis"))
```

Bioconductor packages:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase", "biomaRt", "edgeR", "limma", "rhdf5"))
```

Optional packages that are not strictly required but enhance the functionality of `scater`:

```{r}
install.packages(c("cowplot", "cluster", "mvoutlier", "parallel", "Rtsne"))
biocLite(c("destiny", "monocle"))
```

You might also like to install `dplyr` for convenient data manipulation:

```{r}
install.packages("dplyr")
```


## Getting started

<!---
The best place to start is the [vignette](http://htmlpreview.github.io/?http://github.com/davismcc/scater/blob/master/vignettes/vignette.html).
-->

The best place to start is the vignette. From inside an R session, load `scater`
and then browse the vignettes:

```{r}
library(scater)
browseVignettes("scater")
```

There is a detailed HTML document available that introduces the main features
and functionality of `scater`.

## `scater` workflow

The diagram below provised an overview of the pre-processing and QC workflow possible in `scater`, listing the functions that can be used at various stages.

![Diagram outlining the scater workflow](inst/scater_qc_workflow.png)


## Highlights

The `scater` package allows you to do some neat things relatively quickly. Some highlights are shown below with example code and screenshots.

1. Automated computation of QC metrics
1. Transcript quantification from read data with pseudo-alignment approaches
2. Data format standardisation
3. Rich visualisations for QC and exploratory analysis
4. Seamless integration into the Bioconductor universe
5. Simple normalisation methods

For details of how to use these functions, please consult the **vignette** and **package documentation**.  The plots shown use the example data included with the package (for which there is no interesting structure) and as shown require only one or two lines of code to generate.

### Automatic computation of QC metrics

Use the `calculateQCMetrics` function to compute many metrics useful for gene/transcript-level and cell-level QC. Metrics computed include number of genes expressed per cell, percentage of expression from control genes (e.g. ERCC spike-ins) and many more.

### Transcript quantification with `kallisto`

The `runKallisto` function provides a wrapper to the [`kallisto`](http://pachterlab.github.io/kallisto) software for quantifying transcript abundance from FASTQ files using a pseudo-alignment approach. This new approach is extremely fast. With `readKallisto`, transcript quantities can be read into a data object in `R`.

### Plotting functions

Default `plot` for an SCESet object gives cumulative expression for the
most-expressed features (genes or transcripts)

The `plotTSNE` function produces a t-distributed stochastic neighbour embedding
plot for the cells.

The `plotPCA` function produces a principal components analysis plot for the
cells.

The `plotDiffusionMap` function produces a diffusion map plot for the cells.

The `plotExpression` function plots the expression values for a selection of
features.

The `plotQC` function produces a variety of QC plots useful for diagnostics and
feature and cell filtering. It can be used to plot the most highly-expressed
genes (or features) in the data set or create density plots to assess the
relative importance of explanatory variables, as well as many other
visualisations useful for QC.

The `plotPhenoData` function plots two phenotype metadata variables (such as QC
metrics).

See also `plotFeatureData` to plot feature (gene) metadata variables, including QC metrics.

Plus many, many more possibilities. Please consult the vignette and documentation for details.

## Acknowledgements and disclaimer

The package leans heavily on previously published work and packages, namely
[edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) and
[limma](http://bioconductor.org/packages/release/bioc/html/limma.html). The
`SCESet` class is inspired by the `CellDataSet` class from [monocle](http://www.bioconductor.org/packages/release/bioc/html/monocle.html),
and `SCESet` objects in `scater` can be easily converted to and from `monocle's`
`CellDataSet` objects.


<!---
It also uses and extends code for an approximate rank-product test by [Heskes et al (2014)](http://dx.doi.org/10.1186/s12859-014-0367-1).
-->


The package is currently in an Beta state. The major functionality of the
package is settled, but it is still under development so may change from time
to time. Please do try it and contact me with bug reports, feedback, feature
requests, questions and suggestions to improve the package.

Davis McCarthy, February 2016
