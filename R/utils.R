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
    
    subset <- dummy[subset]
    if (any(is.na(subset))) { 
        stop("invalid subset indices specified")
    }
    return(unname(subset))
}

.get_feature_control_names <- function(object) {
    ## gets the names of the feature controls.
    object@featureControlInfo$name
}

.find_control_SF <- function(object) { 
    ## returns a list of indices and SFs for each control set.
    control_list <- list()
    for (fc in .get_feature_control_names(object)) {
        specific_sf <- suppressWarnings(sizeFactors(object, type=fc))
        if (!is.null(specific_sf)) {
            which.current <- fData(object)[[paste0("is_feature_control_", fc)]]
            control_list[[fc]] <- list(SF=specific_sf, ID=which.current)
        }
    }
    return(control_list)
}

.compute_exprs <- function(exprs_mat, size_factors, log = TRUE,
                           sum = FALSE, logExprsOffset = 1,
                           subset.row = NULL) {
    ## computes normalized expression values.

    size_factors <- size_factors / mean(size_factors)
    if (is.null(subset.row)) { 
        subset.row <- seq_len(nrow(exprs_mat))
    } else {
        subset.row <- .subset2index(subset.row, exprs_mat)
    }
    .checkedCall(cxx_calc_exprs, exprs_mat, as.double(size_factors), 
                 as.double(logExprsOffset), as.logical(log), 
                 as.logical(sum), subset.row - 1L)
}

