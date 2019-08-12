#' Mock up a SingleCellExperiment
#'
#' Mock up a \linkS4class{SingleCellExperiment} containing simulated data,
#' for use in documentation examples.
#'
#' @param ncells Integer scalar, number of cells to simulate.
#' @param ngenes Integer scalar, number of genes to simulate.
#' @param nspikes Integer scalar, number of spike-in transcripts to simulate.
#'
#' @return A SingleCellExperiment object containing a count matrix in the \code{"counts"} assay,
#' a set of simulated \code{\link{colData}} fields,
#' and spike-in data in the \code{"Spikes"} field of \code{\link{altExps}}. 
#'
#' @author Aaron Lun
#' 
#' @details
#' Users should set a seed to obtain reproducible results from this function.
#'
#' @seealso
#' \code{\link{SingleCellExperiment}}, for the constructor.
#'
#' @examples
#' set.seed(1000)
#' sce <- mockSCE()
#' sce
#' 
#' @export
#' @importFrom stats rnbinom runif
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<-
mockSCE <- function(ncells=200, ngenes=2000, nspikes=100) {
    spike.means <- 2^runif(nspikes, 3, 8)
    spike.disp <- 100/spike.means + 0.5
    spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)
    rownames(spike.data) <- sprintf("Spike_%s", formatC(seq_len(nspikes), width=4, flag=0))

    cell.means <- 2^runif(ngenes, 2, 10)
    cell.disp <- 100/cell.means + 0.5
    cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)
    rownames(cell.data) <- sprintf("Gene_%s", formatC(seq_len(ngenes), width=4, flag=0))
    colnames(cell.data) <- sprintf("Cell_%s", formatC(seq_len(ncells), width=3, flag=0))

    sce <- SingleCellExperiment(list(counts=cell.data))
    colData(sce) <- cbind(colData(sce), DataFrame(
        Mutation_Status=sample(c("positive", "negative"), ncells, replace=TRUE),
        Cell_Cycle=sample(c("S", "G0", "G1", "G2M"), ncells, replace=TRUE),
        Treatment=sample(c("treat1", "treat2"), ncells, replace=TRUE)
    ))

    altExp(sce, "Spikes") <- SingleCellExperiment(list(counts=spike.data))
    sce
}
