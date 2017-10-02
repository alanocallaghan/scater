####################################################################
# Internal utilies that are placed here to make them easier to find.

.subset2index <- function(subset, target, byrow=TRUE) {
    ## Converts a subsetting vector into a integer equivalent.
    ## Requires some care to handle logical/character vectors.

    if (is.na(byrow)) {
        dummy <- seq_along(target)
        names(dummy) <- names(target)
    } else if (byrow) {
        dummy <- seq_len(nrow(target))
        names(dummy) <- rownames(target)
    } else {
        dummy <- seq_len(ncol(target))
        names(dummy) <- colnames(target)
    }

    if (!is.null(subset)) {
        subset <- dummy[subset]
        if (any(is.na(subset))) {
            stop("invalid subset indices specified")
        }
    } else {
        subset <- dummy
    }
    return(unname(subset))
}

.get_all_sf_sets <- function(object) {
    fcontrols <- spikeNames(object)
    
    # Storing the default size factors.
    sf.list <- vector("list", length(fcontrols)+1)
    sf.list[[1]] <- sizeFactors(object)
    to.use <- rep(1L, nrow(object))

    # Filling up the controls.
    counter <- 1L
    okay <- character(length(fcontrols))
    for (fc in fcontrols) {
        specific_sf <- sizeFactors(object, type=fc)
        if (!is.null(specific_sf)) {
            okay[counter] <- fc
            counter <- counter+1L
            which.current <- isSpike(object, type=fc)
            to.use[which.current] <- counter # after increment, as 1 is the NULL sizeFactors.
            sf.list[[counter]] <- specific_sf
        }
    }

    # Returning the output.
    return(list(size.factors=sf.list[seq_len(counter)], index=to.use, 
                available=okay[seq_len(counter-1)]))
}

.compute_exprs <- function(exprs_mat, size_factors, sf_to_use=NULL, log = TRUE,
                           sum = FALSE, logExprsOffset = 1,
                           subset_row = NULL) {

    if (!is.list(size_factors)) { 
        size_factors <- list(size_factors)
        sf_to_use <- rep(1L, nrow(exprs_mat))
    }

    ## Mean centers all the size factors.
    for (s in seq_along(size_factors)) {
        sf <- size_factors[[s]]
        sf <- sf/mean(sf)
        size_factors[[s]] <- sf
    }

    ## Specify the rows to be subsetted.
    subset_row <- .subset2index(subset_row, exprs_mat, byrow=TRUE)
    
    ## computes normalized expression values.
    .Call(cxx_calc_exprs, exprs_mat, size_factors, sf_to_use,
          as.double(logExprsOffset), as.logical(log),
          as.logical(sum), subset_row - 1L)
}

.general_rowVars <- function(x, cols=NULL) 
# Computes variance of values from each row of 'x', using only
# columns specified in 'cols' if not NULL.
{
    if (is.null(cols)) { cols <- seq_len(ncol(x)) }
    .Call(cxx_calc_variance, x, cols - 1L, TRUE)
}

.general_colVars <- function(x, rows=NULL) 
# Computes variance of values from each column of 'x', using only
# rows specified in 'rows' if not NULL.
{
    if (is.null(rows)) { rows <- seq_len(nrow(x)) }
    .Call(cxx_calc_variance, x, rows - 1L, FALSE)
}

.general_colSums <- function(mat) {
    margin.stats <- .Call(cxx_margin_summary, mat, 0, 
                          seq_len(nrow(mat)) - 1L, FALSE)
    margin.stats[[1]]
}

.general_rowSums <- function(mat) {
    margin.stats <- .Call(cxx_margin_summary, mat, 0, 
                          seq_len(ncol(mat)) - 1L, TRUE)
    margin.stats[[1]]
}
