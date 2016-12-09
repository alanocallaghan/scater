#' Load in data from 10X experiment
#' 
#' Creates a full or sparse matrix from a sparse data matrix provided by 10X 
#' genomics.
#' 
#' @param data_dir Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load 
#' several data directories. If a named vector is given, the cell barcode names 
#' will be prefixed with the name.
#' @param min_total_cell_counts integer(1) threshold such that cells (barcodes)
#' with total counts below the threshold are filtered out
#' @param min_mean_gene_counts numeric(1) threshold such that genes with mean 
#' counts below the threshold are filtered out.
#' @param expand logical(1), should the sparse count matrix be expanded 
#' into an SCESet object with dense matrices for expression data (default), or
#' should the sparse count matrix be returned?
#' @param logExprsOffset numeric(1) offset value to apply when computing 
#' expression values  as log2(cpm + offset) for the SCESet. Ignored if 
#' \code{expand = FALSE}.
#' 
#' @details This function was developed from the \code{Read10X} function from 
#' the \code{Seurat} package.
#' 
#' @return If \code{expand} is TRUE, returns an SCESet object with counts data 
#' and log2(cpm + offset) as expression data; else returns a sparse matrix with 
#' rows and columns labeled.
#' 
#' @importFrom Matrix readMM
#' @export
#' @examples 
#' \dontrun{
#' sce10x <- read10XResults("path/to/data/directory")
#' count_matrix_10x <- read10XResults("path/to/data/directory", expand = FALSE)
#' }
read10XResults <- function(data_dir = NULL, min_total_cell_counts = 1000L, 
                           min_mean_gene_counts = NULL, expand = TRUE, 
                           logExprsOffset = 1) {
    full_data <- list()
    gene_info_list <- list()
    
    for (i in seq_along(data_dir)) {
        run <- data_dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist")
        }
        
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        
        barcode.loc <- paste(run, "barcodes.tsv", sep = "")
        gene.loc <- paste(run, "genes.tsv", sep = "")
        matrix.loc <- paste(run, "matrix.mtx", sep = "")
        
        if (!file.exists(barcode.loc)) 
            stop("Barcode file missing")
        if (!file.exists(gene.loc))
            stop("Gene name file missing")
        if (!file.exists(matrix.loc))
            stop("Expression matrix file missing")
        
        ## read sparse count matrix
        data_mat <- Matrix::readMM(matrix.loc)
        
        ## define filters
        keep_barcode <- (Matrix::colSums(data_mat) >= min_total_cell_counts)
        data_mat <- data_mat[, keep_barcode]
        
        
        cell.names <- data.table::fread(barcode.loc, header = FALSE)[[1]]
        # cell.names <- readLines(barcode.loc)
        if (all(grepl("\\-1$", cell.names)) == TRUE) {
            cell.names <- gsub("-.*$", "", cell.names)
            cell.names <- cell.names[keep_barcode]
        } else
            cell.names <- sprintf("cell_%09s", seq_len(sum(keep_barcode)))
        
        
        gene_info_list[[i]] <- data.table::fread(gene.loc, header = FALSE)

        if (is.null(names(data_dir))) {
            if (i < 2) {
                colnames(data_mat) <- cell.names
            }
            else {
                colnames(data_mat) <- paste0(i, "_", cell.names, sep = "") 
            }
        } else {
            colnames(data_mat) <- paste0(names(data_dir)[i],"_",cell.names) 
        }
        full_data <- append(full_data, data_mat)
    }
    
    if (length(unique(gene_info_list)) != 1L)
        stop("gene information differs between runs")
    gene_info <- as(gene_info_list[[1]], "AnnotatedDataFrame")
    
    full_data <- do.call(cbind, full_data)
    if (ncol(full_data) < 1L)
        stop("no cells passing min_total_cell_counts filter")
    
    if (!is.null(min_mean_gene_counts)) {
        keep_gene <- (Matrix::rowSums(data_mat) >= min_mean_gene_counts)
        full_data <- full_data[keep_gene,]
        gene_info <- gene_info[keep_gene,]
    }
    colnames(gene_info) <- c("id", "symbol")
    rownames(gene_info) <- gene_info[[1]]
    rownames(full_data) <- gene_info[[1]]
    
    if (expand)
        out <- newSCESet(countData = as.matrix(full_data), 
                         featureData = gene_info, 
                         logExprsOffset = logExprsOffset)
    else 
        out <- full_data
    out
}


