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

