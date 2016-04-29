

# Mutate verb -------------------------------------------------------------

#' @rdname mutate
#' @export
setMethod("mutate", signature(object = "SCESet"),
          function(object, ...) {
            mutate.SCESet(object, ...)
          })

#' Add new variables to \code{pData(object)}.
#' 
#' Adds new columns to \code{pData(object)} preserving existing variables.
#' 
#' 
#' @param object A \code{SCESet} object.
#' @param ... Additional arguments to be passed to \code{dplyr::mutate} to
#' act on \code{pData(object)}.
#' 
#' @export
#' @rdname mutate
#' @name mutate
#' @aliases mutate mutate,SCESet-method
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- mutate(example_sceset, is_quiescent = Cell_Cycle == "G0")
mutate.SCESet <- function(object, ...) {
  pData(object) <- dplyr::mutate(pData(object), ...)
  return( object )
}

# Rename verb -------------------------------------------------------------

#' @rdname rename
#' @export
setMethod("rename", signature(object = "SCESet"),
          function(object, ...) {
            rename.SCESet(object, ...)
          })

#' Rename variables of \code{pData(object)}.
#' 
#' @param object A \code{SCESet} object.
#' @param ... Additional arguments to be passed to \code{dplyr::rename} to
#' act on \code{pData(object)}.
#' 
#' @export
#' @rdname rename
#' @name rename
#' @aliases rename rename,SCESet-method
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- rename(example_sceset, Cell_Phase = Cell_Cycle)
rename.SCESet <- function(object, ...) {
  pData(object) <- dplyr::rename(pData(object), ...)
  return( object )
}


# Filter verb -------------------------------------------------------------

#' @rdname filter
#' @export
setMethod("filter", signature(object = "SCESet"),
          function(object, ...) {
            filter.SCESet(object, ...)
          })

#' Return \code{SCESet} with cells matching conditions.
#' 
#' Subsets the columns (cells) of a \code{SCESet} based on 
#' matching conditions in the rows of \code{pData(object)}.
#' 
#' @param object A \code{SCESet} object.
#' @param ... Additional arguments to be passed to \code{dplyr::filter} to
#' act on \code{pData(object)}.
#' 
#' @export
#' @rdname filter
#' @name filter
#' @aliases filter filter,SCESet-method
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset_treat1 <- filter(example_sceset, Treatment == "treat1")
filter.SCESet <- function(object, ...) {
  pd <- pData(object)
  pd <- dplyr::mutate(pd, scater_placeholder_index = 1:nrow(pd))
  pd <- dplyr::filter(pd, ...) # filter not exported by dplyr
  object <- object[, pd$scater_placeholder_index]
  return( object )
}

# Arrange verb ------------------------------------------------------------

#' @rdname arrange
#' @export
setMethod("arrange", signature(object = "SCESet"),
          function(object, ...) {
            arrange.SCESet(object, ...)
          })

#' Arrange rows of \code{pData(object)} by variables.
#' 
#' The \code{SCESet} returned will have cells ordered by the corresponding
#' variable in \code{pData(object)}.
#' 
#' @param object A \code{SCESet} object.
#' @param ... Additional arguments to be passed to \code{dplyr::arrange} to
#' act on \code{pData(object)}.
#' 
#' @export
#' @rdname arrange
#' @name arrange
#' @aliases arrange arrange,SCESet-method
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- arrange(example_sceset, Cell_Cycle)
arrange.SCESet <- function(object, ...) {
  pd <- pData(object)
  pd <- dplyr::mutate(pd, scater_placeholder_index = 1:nrow(pd))
  pd <- dplyr::arrange(pd, ...) # arrange not exported by dplyr
  object <- object[, pd$scater_placeholder_index]
  return( object )
}
