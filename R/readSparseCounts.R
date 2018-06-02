#' Read sparse count matrix from file
#'
#' Reads a sparse count matrix from file containing a dense tabular format.
#'
#' @param file A string containing a file path to a count table, or a connection object opened in read-only text mode.
#' @param sep A string specifying the delimiter between fields in \code{file}.
#' @param quote A string specifying the quote character, e.g., in column or row names.
#' @param comment.char A string specifying the comment character after which values are ignored.
#' @param row.names A logical scalar specifying whether row names are present. 
#' @param col.names A logical scalar specifying whether column names are present. 
#' @param ignore.row An integer scalar specifying the number of rows to ignore at the start of the file, \emph{before} the column names.
#' @param skip.row An integer scalar specifying the number of rows to ignore at the start of the file, \emph{after} the column names.
#' @param ignore.col An integer scalar specifying the number of columns to ignore at the start of the file, \emph{before} the column names.
#' @param skip.col An integer scalar specifying the number of columns to ignore at the start of the file, \emph{after} the column names.
#' @param chunk A integer scalar indicating the chunk size to use, i.e., number of rows to read at any one time.
#'
#' @details
#' This function provides a convenient method for reading dense arrays from flat files into a sparse matrix in memory.
#' Memory usage can be further improved by setting \code{chunk} to a smaller positive value.
#'
#' The \code{ignore.*} and \code{skip.*} parameters allow irrelevant rows or columns to be skipped.
#' Note that the distinction between the two parameters is only relevant when \code{row.names=FALSE} (for skipping/ignoring columns) or \code{col.names=FALSE} (for rows).
#' 
#' @return
#' A dgCMatrix containing double-precision values (usually counts) for each row (gene) and column (cell).
#' 
#' @author Aaron Lun
#' 
#' @seealso 
#' \code{\link{read.table}},
#' \code{\link{readMM}}
#'
#' @examples
#' outfile <- tempfile()
#' write.table(data.frame(A=1:5, B=0, C=0:4, row.names=letters[1:5]), 
#'     file=outfile, col.names=NA, sep="\t", quote=FALSE)
#'
#' readSparseCounts(outfile)
#'
#' @export
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#' @importFrom utils tail read.table
readSparseCounts <- function(file, sep="\t", quote=NULL, comment.char="", row.names=TRUE, col.names=TRUE, 
    ignore.row=0L, skip.row=0L, ignore.col=0L, skip.col=0L, chunk=1000L)
{
    if (!is(file, "connection")) {
        fhandle <- file(file, open="r")
        on.exit(close(fhandle))
    } else {
        fhandle <- file
    }

    # Scanning through rows.
    if (ignore.row) {
        readLines(fhandle, n=ignore.row)
    }
    if (col.names) {
        cell.names <- read.table(fhandle, sep=sep, quote=quote, comment.char=comment.char, nrows=1L, 
            stringsAsFactors=FALSE, header=FALSE)
        cell.names <- as.character(cell.names)
    } else {
        cell.names <- NULL
    }
    if (skip.row) {
        readLines(fhandle, n=skip.row)
    }

    # Figuring out how to extract the columns.
    first <- read.table(fhandle, sep=sep, quote=quote, comment.char=comment.char, nrows=1L, 
        stringsAsFactors=FALSE, header=FALSE)

    nentries <- ncol(first)
    what <- vector("list", nentries)
    if (row.names) {
        row.name.col <- ignore.col + 1L
        what[[row.name.col]] <- "character"
        ignore.col <- row.name.col
    } 

    skip.col <- skip.col + ignore.col
    ncells <- nentries - skip.col
    cell.cols <- skip.col + seq_len(ncells)
    what[cell.cols] <- 0

    # Processing the first element.
    gene.names <- NULL
    if (row.names) { 
        gene.names <- first[[row.name.col]]
    }
    output <- list(as.matrix(first[,cell.cols]))
    colnames(output[[1]]) <- NULL
    it <- 2L

    # Reading it in, chunk by chunk.
    repeat {
        current <- scan(fhandle, what=what, sep=sep, quote=quote, comment.char=comment.char, nmax=chunk, quiet=TRUE)
        if (row.names) {
            gene.names <- c(gene.names, current[[row.name.col]])
        }
        output[[it]] <- as(do.call(cbind, current[cell.cols]), "dgCMatrix")
        it <- it + 1L
        if (chunk<0 || length(current[[1]]) < chunk) {
            break
        }
    }
    output <- do.call(rbind, output)

    # Adding row and column names, if available. 
    if (row.names) {
        rownames(output) <- gene.names
    }
    if (col.names) {
        colnames(output) <- tail(cell.names, ncells)
    }
    return(output)
}
