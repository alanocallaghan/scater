#' Make feature names unique
#' 
#' Combine a user-interpretable feature name (e.g., gene symbol) with a standard identifier that is guaranteed to be unique and valid (e.g., Ensembl) for use as row names.
#' 
#' @param ID A character vector of unique identifiers.
#' @param names A character vector of feature names.
#'
#' @details This function will attempt to use \code{names} if it is unique.
#' If not, it will append the \code{_ID} to any non-unique value of \code{names}.
#' Missing \code{names} will be replaced entirely by \code{ID}. 
#' 
#' The output is guaranteed to be unique, assuming that \code{ID} is also unique.
#' This can be directly used as the row names of a SingleCellExperiment object.
#'
#' @return A character vector of unique-ified feature names.
#'
#' @author Aaron Lun
#'
#' @export
#' @examples
#' uniquifyFeatureNames(
#'   ID=paste0("ENSG0000000", 1:5),
#'   names=c("A", NA, "B", "C", "A")
#' )
uniquifyFeatureNames <- function(ID, names) {
    if (length(ID)!=length(names)) {
        stop("lengths of 'ID' and 'names' must be equal")
    }
    if (is.factor(names)) {
        names <- as.character(names)
    }
    missing.name <- is.na(names)
    names[missing.name] <- ID[missing.name]
    dup.name <- names %in% names[duplicated(names)]
    names[dup.name] <- paste0(names[dup.name], "_", ID[dup.name])
    return(names)
}

