# Other internal utilies that are placed here to make them easier to find.

.subset2index <- function(subset, target, byrow=TRUE) 
## Converts a subsetting vector into a integer equivalent.
## Requires some care to handle logical/character vectors.
{
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

#' @importFrom SummarizedExperiment colData rowData
#' @importFrom BiocGenerics colnames
.qc_hunter <- function(object, qc_field, mode = "column", error = TRUE) 
# This function searches for QC fields in the various plotQC functions,
# accounting for potential compactness and nesting of DataFrames.
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

