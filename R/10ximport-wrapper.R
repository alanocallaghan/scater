#' Load in data from 10x experiment
#' 
#' Creates a full or sparse matrix from a sparse data matrix provided by 10X 
#' genomics.
#' 
#' @param data_dir Directory containing the matrix.mtx, genes.tsv, and 
#' barcodes.tsv files provided by 10x. A vector or named vector can be given in order to load 
#' several data directories. If a named vector is given, the cell barcode names 
#' will be prefixed with the name.
#' @param min_total_cell_counts integer(1) threshold such that cells (barcodes)
#' with total counts below the threshold are filtered out
#' @param min_mean_gene_counts numeric(1) threshold such that genes with mean 
#' counts below the threshold are filtered out.
#' @param ... passed arguments
#' 
#' @details This function was developed from the \code{Read10X} function from 
#' the \pkg{Seurat} package.
#' 
#' @return Returns an SingleCellExperiment object with counts data stored as a
#' sparse matrix. Rows are named with the gene name and columns are named with 
#' the cell barcode (if \code{data_dir} contains one element; otherwise the 
#' columns are unnamed to avoid problems with non-unique barcodes).
#' 
#' @importFrom Matrix readMM
#' @import Matrix
#' @rdname read10xResults
#' @aliases read10xResults read10XResults
#' @export
#' @examples 
#' sce10x <- read10xResults(system.file("extdata", package="scater"))
#' 
read10xResults <- function(data_dir, min_total_cell_counts = NULL, 
                           min_mean_gene_counts = NULL) { 
    
    nsets <- length(data_dir)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)
    
    for (i in seq_len(nsets)) { 
        run <- data_dir[i]
        barcode.loc <- file.path(run, "barcodes.tsv")
        gene.loc <- file.path(run, "genes.tsv")
        matrix.loc <- file.path(run, "matrix.mtx")
        
        ## read sparse count matrix and cell barcodes.
        data_mat <- Matrix::readMM(matrix.loc)
        data_mat <- as(data_mat, "dgCMatrix")
        cell.names <- utils::read.table(barcode.loc, header = FALSE, 
                                        colClasses = "character")[[1]]

        ## define filters
        if (!is.null(min_total_cell_counts)) { 
            keep_barcode <- .general_colSums(data_mat) >= min_total_cell_counts
            data_mat <- data_mat[, keep_barcode]
            cell.names <- cell.names[keep_barcode]
        }

        dataset <- i
        if (!is.null(names(data_dir))) {
            dataset <- names(data_dir)[i]
        }
        
        full_data[[i]] <- data_mat
        gene_info_list[[i]] <- utils::read.table(gene.loc, header = FALSE, 
                                                 colClasses = "character")
        cell_info_list[[i]] <- DataFrame(dataset = dataset, 
                                         barcode = cell.names)
    }
    
    # Checking gene uniqueness. 
    if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
        stop("gene information differs between runs")
    }
    gene_info <- gene_info_list[[1]]
    colnames(gene_info) <- c("id", "symbol")
    rownames(gene_info) <- gene_info$id
    
    # Forming the full data matrix.
    full_data <- do.call(cbind, full_data)
    rownames(full_data) <- gene_info$id

    # Applying some filtering if requested.
    if (!is.null(min_mean_gene_counts)) {
        keep_gene <- .general_rowSums(data_mat) >= min_mean_gene_counts
        full_data <- full_data[keep_gene,]
        gene_info <- gene_info[keep_gene,]
    }
    
    # Adding the cell data (only using as colnames if there is only 1 set - guaranteed unique).
    cell_info <- do.call(rbind, cell_info_list)
    if (nsets==1L) {
        colnames(full_data) <- cell_info$barcode
    }
    SingleCellExperiment(list(counts = full_data), rowData = gene_info, 
                         colData = cell_info)
}


#' @rdname read10xResults
#' @export
read10XResults <- function(...) {
    read10xResults(...)
}


#' Downsample a count matrix
#' 
#' Downsample a count matrix to a desired proportion.
#' 
#' @param x matrix of counts
#' @param prop numeric scalar or vector of length \code{ncol(x)} in [0, 1] 
#' indicating the downsampling proportion
#' 
#' @details Given multiple 10X batches of very different sequencing depths, it
#' can be beneficial to downsample the deepest batches to match the coverage of
#' the shallowest batches. This avoids differences in technical noise that can
#' drive clustering by batch. 
#' 
#' Downsampling without replacement is performed on the counts in each cell to
#' generate the output matrix. Each count in the returned matrix is guaranteed
#' to be smaller than the original value in \code{x}. This provides an 
#' alternative to downsampling in the CellRanger \code{aggr} function.
#' 
#' @return an integer matrix of downsampled counts
#' 
#' @export
#' @examples 
#' sce10x <- read10xResults(system.file("extdata", package="scater"))
#' downsampled <- downsampleCounts(counts(sce10x), prop = 0.5)
#' 
downsampleCounts <- function(x, prop) {
    prop <- rep(prop, length.out = ncol(x))
    .Call(cxx_downsample_matrix, x, prop)
}

