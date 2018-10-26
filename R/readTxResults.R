#' Read transcript quantification data
#' 
#' Create a \linkS4class{SingleCellExperiment} object from pseudo-aligner results via \pkg{tximport}.
#'
#' @param ... Arguments to be passed to \code{\link[tximport]{tximport}}.
#' @param full_length Logical scalar indicating whether the sequencing data is full-length (e.g., Smart-seq2) or end-biased, e.g., UMI-based protocols.
#' 
#' @details 
#' If \code{full_length=TRUE}, counts are computed from the length-scaled TPMs.
#' Otherwise, counts are not computed from the abundances.
#'
#' \code{readKallistoResults} and \code{readSalmonResults} are simply wrappers around \code{readTxResults} with \code{type} pre-specified.
#'
#' This function has now been deprecated in favour of \code{\link[tximeta]{tximeta}}.
#' The latter produces a SummarizedExperiment that is easily coerced into a SingleCellExperiment.
#' 
#' @return 
#' A SingleCellExperiment containing the abundance, count and feature length information from the supplied samples.
#'
#' @author Davis McCarthy and Aaron Lun
#'
#' @export
#' @rdname readTxResults
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @references 
#' Soneson C, Love MI, Robinson MD (2015).
#' Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. 
#' \emph{F1000Res.} 4, 1521.
#' 
#' @examples 
#' library(tximportData)
#' dir <- system.file("extdata", package="tximportData")
#' samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
#' files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
#' names(files) <- paste0("sample",1:6)
#'                         
#' # tx2gene links transcript IDs to gene IDs for summarization
#' tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
#' 
#' sce <- readTxResults(files, type="salmon", tx2gene=tx2gene)
readTxResults <- function(..., full_length=TRUE) {
    .Deprecated(new="tximeta::tximeta")
    tx_out <- tximport::tximport(..., countsFromAbundance=if (full_length) "lengthScaledTPM" else "no")
    SingleCellExperiment(tx_out[c("counts", "abundance", "length")])
}

#' @rdname readTxResults
#' @export
readKallistoResults <- function(..., full_length=TRUE) {
    readTxResults(..., full_length=full_length, type="kallisto")
}

#' @rdname readTxResults
#' @export
readSalmonResults <- function(..., full_length=TRUE) {
    readTxResults(..., full_length=full_length, type="salmon")
}
