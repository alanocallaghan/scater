#' @name norm_exprs
#' @aliases norm_exprs
#' @export
#' @docType methods
#' @return a matrix of normalised expression data
#' @rdname accessors
setGeneric("norm_exprs", function(object) standardGeneric("norm_exprs"))

#' @name norm_exprs<-
#' @aliases norm_exprs<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("norm_exprs<-", function(object, value) standardGeneric("norm_exprs<-"))

#' @name stand_exprs
#' @aliases stand_exprs
#' @export
#' @docType methods
#' @return a matrix of standardised expressiond data
#' @rdname accessors
setGeneric("stand_exprs", function(object) standardGeneric("stand_exprs"))

#' @name stand_exprs<-
#' @aliases stand_exprs<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("stand_exprs<-", function(object, value) standardGeneric("stand_exprs<-"))

#' @name bootstraps
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps", function(object) standardGeneric("bootstraps"))

#' @name bootstraps<-
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps<-", function(object, value) standardGeneric("bootstraps<-"))

#' @name fpkm
#' @aliases  fpkm
#' @export
#' @docType methods
#' @return a matrix of FPKM values
#' @rdname accessors
setGeneric("fpkm", function(object) standardGeneric("fpkm"))

#' @name fpkm<-
#' @aliases  fpkm<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("fpkm<-", function(object, value) standardGeneric("fpkm<-"))


## dplyr-style verb generics

#' @name mutate
#' @rdname mutate
#' @docType methods
#' @export
setGeneric("mutate", function(object, ...) standardGeneric("mutate"))

#' @name rename
#' @rdname rename
#' @docType methods
#' @export
setGeneric("rename", function(object, ...) standardGeneric("rename"))

#' @name filter
#' @rdname filter
#' @docType methods
#' @export
setGeneric("filter", function(object, ...) standardGeneric("filter"))

#' @name arrange
#' @rdname arrange
#' @docType methods
#' @export
setGeneric("arrange", function(object, ...) standardGeneric("arrange"))

