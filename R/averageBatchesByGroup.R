#' Compute batch-aware averages from group-level statistics
#'
#' Compute an average statistic for each group in a batch-aware manner, by fitting a linear model and extracting the coefficients.
#' This handles statistics such as the average log-expression or the average number of cells with detected expression.
#'
#' @param x A numeric matrix containing statistics for each gene (row) and combination of group and block (column).
#' @param group A factor or vector specifying the group identity for each column of \code{x}, usually clusters or cell types.
#' @param block A factor or vector specifying the blocking level for each column of \code{x}, e.g., batch of origin.
#' @param transform String indicating how the differences between groups should be computed, for the batch adjustment.
#' @param offset Numeric scalar specifying the offset to use when \code{difference="log"} (default 1) or \code{difference="logit"} (default 0.01).
#' 
#' @return A numeric matrix with number of rows equal to \code{nrow(x)} and number of columns equal to the number of unique levels in \code{group}.
#' Each column corresponds to a group and contains the averaged statistic across batches.
#'
#' @details
#' The problem with directly averaging group-level statistics across batches is that some groups may not exist in particular batches,
#' e.g., due to the presence of unique cell types in different samples.
#' A direct average would be affected by batch effects, making it difficult to intepret the value.
#' 
#' To overcome this, we use any shared levels to correct for the batch effect.
#' For each gene, we fit a linear model to the (transformed) values containing both the group and block factors.
#' We then report the coefficient for each group as the batch-adjusted average, which is possible as the fitted model has no intercept.
#'
#' The default of \code{transform="raw"} will not transform the values, and is generally suitable for log-expression values.
#' Setting \code{transform="log"} will perform a log-transformation after adding \code{offset}, and is suitable for normalized counts.
#' Setting \code{transform="logit"} will perform a logit transformation after adding \code{offset} to the numerator and denominator (to shrink towards 0.5),
#' and is suitable for proportional data such as the proportion of detected cells.
#'
#' After the model is fitted to the transformed values, the reverse transformation is applied to the coefficients to obtain the batch-adjusted average.
#' For \code{transform="log"}, any negative values are coerced to zero,
#' while for \code{transform="logit"}, any values outside of [0, 1] are coerced to the closest boundary.
#'
#' @author Aaron Lun
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=1000)
#' group <- sample(10, ncol(y), replace=TRUE)
#' block <- sample(5, ncol(y), replace=TRUE)
#'
#' # Note to self: replace with scuttle::summarizeAssayByGroup
#' library(scater)
#' summed <- sumCountsAcrossCells(y, DataFrame(group=group, block=block), average=TRUE)
#'
#' # Computing batch-aware averages:
#' averaged <- averageBatchesByGroup(assay(summed), 
#'     group=summed$group, block=summed$block)
#' 
#' @export
#' @importFrom stats lm.fit model.matrix
#' @importFrom Matrix t
averageBatchesByGroup <- function(x, group, block, transform=c("raw", "log", "logit"), offset=NULL) { 
    transform <- match.arg(transform)
    if (transform=="log") {
        if (is.null(offset)) {
            offset <- 1
        }
        x <- log(x + offset)
    } else if (transform=="logit") {
        if (is.null(offset)) {
            offset <- 0.01
        }
        x <- (x + offset) / (1 + 2 * offset)
        x <- log(x/(1-x))
    }

    group <- factor(group)
    if (length(unique(block))==1L) {
        design <- model.matrix(~0 + group)
    } else {
        design <- model.matrix(~0 + group + factor(block))
    }

    # Replace with fit.
    fit <- lm.fit(y=t(x), x=design)
    averages <- t(fit$coefficients[seq_len(nlevels(group)),,drop=FALSE])
    colnames(averages) <- levels(group)
    rownames(averages) <- rownames(x)

    if (transform=="log") {
        averages <- exp(averages) - offset
        averages[averages < 0] <- 0
    } else if (transform=="logit") {
        averages <- exp(averages)
        averages <- averages/(averages + 1)
        averages <- averages * (1 + 2 * offset) - offset
        averages[averages < 0] <- 0
        averages[averages > 1] <- 1
    }

    averages
}
