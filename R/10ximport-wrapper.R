#' Load in data from 10X experiment
#' 
#' Creates a full matrix from a sparse data matrix provided by 10X genomics.
#' 
#' This function is identical to the \code{Read10X} from \code{Seurat} package.
#' 
#' @param data_dir Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load 
#' several data directories. If a named vector is given, the cell barcode names 
#' will be prefixed with the name.
#' @return Returns a full matrix with rows and columns labeled 
#' @importFrom Matrix readMM
#' @export
read10XResults <- function(data_dir = NULL){
    full_data <- list()
    for(i in seq_along(data_dir)){
        run <- data_dir[i]
        if (!dir.exists(run)){
            stop("Directory provided does not exist")
        }
        
        if(!grepl("\\/$", run)){
            run <- paste(run, "/", sep = "")
        }
        
        barcode.loc <- paste(run, "barcodes.tsv", sep ="")
        gene.loc <- paste(run, "genes.tsv", sep ="")
        matrix.loc <- paste(run, "matrix.mtx", sep ="")
        
        if (!file.exists(barcode.loc)){
            stop("Barcode file missing")
        }
        if (!file.exists(gene.loc)){
            stop("Gene name file missing")
        }
        if (!file.exists(matrix.loc)){
            stop("Expression matrix file missing")
        }
        
        data <- readMM(matrix.loc)
        cell.names <- readLines(barcode.loc)
        gene.names <- readLines(gene.loc)
        if(all(grepl("\\-1$", cell.names)) == TRUE) {
            cell.names <- as.vector(as.character(sapply(cell.names, extract_field, 1, delim = "-")))
        }
        rownames(data) <- make.unique(as.character(sapply(gene.names, extract_field, 2, delim = "\\t"))) 
        
        if(is.null(names(data_dir))){
            if(i < 2){
                colnames(data) <- cell.names
            }
            else {
                colnames(data) <- paste0(i, "_", cell.names, sep = "") 
            }
        } else {
            colnames(data) <- paste0(names(data_dir)[i],"_",cell.names) 
        }
        full_data <- append(full_data, data)
    }
    full_data <- do.call(cbind, full_data)
    return(as.matrix(full_data))
}
