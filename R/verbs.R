

# Mutate verb -------------------------------------------------------------

#' Add new variables to \code{colData(object)} (deprecated).
#'
#' @param object a \code{SingleCellExperiment} object.
#' @param ... Additional arguments to be passed to \code{dplyr::mutate} to
#' act on \code{colData(object)}.
#'
#' @return An SingleCellExperiment object.
#'
#' @export
#' @rdname mutate
#' @name mutate
#' @aliases mutate mutate,SingleCellExperiment-method
#'
#' @details 
#' Refer to \url{https://github.com/sa-lee/plyexperiment} for replacement functionality.
#'
#' @rdname mutate
#' @export
setMethod("mutate", "SingleCellExperiment", function(object, ...) {
    .Deprecated()
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, ...)
    rownames(pd) <- colnames(object)
    colData(object) <- DataFrame(pd)
    return(object)
})


# Rename verb -------------------------------------------------------------

#' Rename variables of \code{colData(object)} (deprecated).
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param ... Additional arguments to be passed to \code{dplyr::rename} to
#' act on \code{colData(object)}.
#'
#' @return An SingleCellExperiment object.
#'
#' @export
#' @rdname rename
#' @name rename
#' @aliases rename rename,SingleCellExperiment-method
#'
#' @details 
#' Refer to \url{https://github.com/sa-lee/plyexperiment} for replacement functionality.
#' 
setMethod("rename", "SingleCellExperiment", function(object, ...) {
    .Deprecated()
    pd <- as.data.frame(colData(object))
    pd <- dplyr::rename(pd, ...)
    rownames(pd) <- colnames(object)
    colData(object) <- DataFrame(pd)
    return(object)
})

# Filter verb -------------------------------------------------------------

#' Return \code{SingleCellExperiment} with cells matching conditions (deprecated).
#'
#' Subsets the columns (cells) of a \code{SingleCellExperiment} based on
#' matching conditions in the rows of \code{colData(object)}.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param ... Additional arguments to be passed to \code{dplyr::filter} to
#' act on \code{colData(object)}.
#'
#' @return An SingleCellExperiment object.
#'
#' @details 
#' Refer to \url{https://github.com/sa-lee/plyexperiment} for replacement functionality.
#'
#' @export
#' @rdname filter
#' @name filter
#' @aliases filter filter,SingleCellExperiment-method
setMethod("filter", "SingleCellExperiment", function(object, ...) {
    .Deprecated()
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, scater_placeholder_index = seq_len(nrow(pd)))
    pd <- dplyr::filter(pd, ...) # filter not exported by dplyr
    object[, pd$scater_placeholder_index]
})

# Arrange verb ------------------------------------------------------------

#' Arrange columns (cells) of a SingleCellExperiment object (deprecated).
#'
#' The \code{SingleCellExperiment} returned will have cells ordered by the corresponding
#' variable in \code{colData(object)}.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param ... Additional arguments to be passed to \code{dplyr::arrange} to
#' act on \code{colData(object)}.
#'
#' @return An SingleCellExperiment object.
#'
#' @details 
#' Refer to \url{https://github.com/sa-lee/plyexperiment} for replacement functionality.
#'
#' @export
#' @rdname arrange
#' @name arrange
#' @aliases arrange arrange,SingleCellExperiment-method
setMethod("arrange", "SingleCellExperiment", function(object, ...) {
    .Deprecated()
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, scater_placeholder_index = seq_len(nrow(pd)))
    pd <- dplyr::arrange(pd, ...) # arrange not exported by dplyr
    object[, pd$scater_placeholder_index]
})

