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
    for (fc in fcontrols) {
        specific_sf <- sizeFactors(object, type=fc)
        if (is.null(specific_sf)) {
            warning(sprintf("spike-in set '%s' should have its own size factors", fc))
        } else {
            counter <- counter+1L
            which.current <- isSpike(object, type=fc)
            to.use[which.current] <- counter # after increment, as 1 is the NULL sizeFactors.
            sf.list[[counter]] <- specific_sf
        }
    }

    # Returning the output.
    return(list(size.factors=sf.list[seq_len(counter)], index=to.use)) 
}

.compute_exprs <- function(exprs_mat, 
                           size_factor_val, size_factor_idx,
                           sum = FALSE, log = TRUE, logExprsOffset = 1,
                           subset_row = NULL) {

    ## Specify the rows to be subsetted.
    subset_row <- .subset2index(subset_row, exprs_mat, byrow=TRUE)
    
    ## computes normalized expression values.
    .Call(cxx_calc_exprs, exprs_mat, size_factor_val, size_factor_idx,
          as.double(logExprsOffset), as.logical(log),
          as.logical(sum), subset_row - 1L)
}

.qc_hunter <- function(object, qc_field, mode = "column", error = TRUE) 
# This function searches for QC fields in the various plotQC functions,
# accounting for potential compactness.
{
    if (mode=="column") {
        meta_data <- colData(object)
        protected <- c("feature_control", "endogenous") 
        setname <- "^feature_control"
    } else {
        meta_data <- rowData(object)
        protected <- c("cell_control", "non_control")
        setname <- "^cell_control"
    }

    # Simple is best.
    if (qc_field %in% colnames(meta_data)) { 
        return(qc_field)
    }

    # Looking inside.
    meta_data <- meta_data$scater_qc
    if (qc_field %in% colnames(meta_data)) {
        return(c("scater_qc", qc_field))
    }

    # Looking digging further for the fields.
    for (subfield in colnames(meta_data)) {
        sub_meta_data <- meta_data[[subfield]]
        if (!is(sub_meta_data, "DataFrame")) { 
            next
        }
        
        if (subfield == "all") {
            suffix <- ""
        } else if (subfield %in% protected) {
            suffix <- paste0("_", subfield)
        } else {
            suffix <- sub(setname, "", subfield)
        }
        
        renamed <- sprintf("%s%s", colnames(sub_meta_data), suffix)
        present <- qc_field==renamed
        if (any(present)) {
            return(c("scater_qc", subfield, colnames(sub_meta_data)[which(present)[1]]))
        }
    }

    if (error) {
        FUN <- stop
    } else {
        FUN <- warning
    }
    FUN(sprintf("failed to find '%s' in %s metadata", qc_field, mode))
    return(NULL)
}

########################################################
# matrixStats equivalents that are yet to have a home. #
########################################################

.realize_subsets <- function(x, rows, cols) {
    if (is.null(cols)) { 
        cols <- seq_len(ncol(x)) 
    } else {
        cols <- .subset2index(cols, x, byrow=FALSE)
    }
    if (is.null(rows)) { 
        rows <- seq_len(nrow(x))
    } else {
        rows <- .subset2index(rows, x, byrow=TRUE)
    }
    return(list(rows=rows, cols=cols))
}

.rowVars <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_vars, x, converted$rows - 1L, converted$cols - 1L)
}

.colVars <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_vars, x, converted$rows - 1L, converted$cols - 1L)
}

.rowSums <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_sums, x, converted$rows - 1L, converted$cols - 1L)
}

.colSums <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_sums, x, converted$rows - 1L, converted$cols - 1L)
}

.rowAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

.colAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

