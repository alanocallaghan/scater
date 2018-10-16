

# Mutate verb -------------------------------------------------------------

#' Add new variables to \code{colData(object)}.
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
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- mutate(example_sce, is_quiescent = Cell_Cycle == "G0")
#' @rdname mutate
#' @export
setMethod("mutate", "SingleCellExperiment", function(object, ...) {
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, ...)
    rownames(pd) <- colnames(object)
    colData(object) <- DataFrame(pd)
    return(object)
})


# Rename verb -------------------------------------------------------------

#' Rename variables of \code{colData(object)}.
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
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- rename(example_sce, Cell_Phase = Cell_Cycle)
#' 
setMethod("rename", "SingleCellExperiment", function(object, ...) {
    pd <- as.data.frame(colData(object))
    pd <- dplyr::rename(pd, ...)
    rownames(pd) <- colnames(object)
    colData(object) <- DataFrame(pd)
    return(object)
})

# Filter verb -------------------------------------------------------------

#' Return \code{SingleCellExperiment} with cells matching conditions.
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
#' @export
#' @rdname filter
#' @name filter
#' @aliases filter filter,SingleCellExperiment-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce_treat1 <- filter(example_sce, Treatment == "treat1")
#' 
setMethod("filter", "SingleCellExperiment", function(object, ...) {
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, scater_placeholder_index = seq_len(nrow(pd)))
    pd <- dplyr::filter(pd, ...) # filter not exported by dplyr
    object[, pd$scater_placeholder_index]
})

# Arrange verb ------------------------------------------------------------

#' Arrange columns (cells) of a SingleCellExperiment object
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
#' @export
#' @rdname arrange
#' @name arrange
#' @aliases arrange arrange,SingleCellExperiment-method
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), 
#' colData = sc_example_cell_info)
#' example_sce <- arrange(example_sce, Cell_Cycle)
#' 
setMethod("arrange", "SingleCellExperiment", function(object, ...) {
    pd <- as.data.frame(colData(object))
    pd <- dplyr::mutate(pd, scater_placeholder_index = seq_len(nrow(pd)))
    pd <- dplyr::arrange(pd, ...) # arrange not exported by dplyr
    object[, pd$scater_placeholder_index]
})

