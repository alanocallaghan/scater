#' Compute batch-corrected group-level averages
#'
#' Compute an average statistic for each group in a manner that corrects for batch effects, by fitting a linear model and extracting the coefficients.
#' This handles statistics such as the average log-expression or the proportion of cells with detected expression.
#'
#' @param x A numeric matrix containing statistics for each gene (row) and combination of group and block (column),
#' computed by functions such as \code{\link{summarizeAssayByGroup}} - see Examples.
#' @param group A factor or vector specifying the group identity for each column of \code{x}, usually clusters or cell types.
#' @param block A factor or vector specifying the blocking level for each column of \code{x}, e.g., batch of origin.
#' @param transform String indicating how the differences between groups should be computed, for the batch adjustment.
#' @param offset Numeric scalar specifying the offset to use when \code{difference="log"} (default 1) or \code{difference="logit"} (default 0.01).
#' 
#' @return A numeric matrix with number of rows equal to \code{nrow(x)} and number of columns equal to the number of unique levels in \code{group}.
#' Each column corresponds to a group and contains the averaged statistic across batches.
#'
#' @details
#' This function considers group-level statistics such as the average expression of all cells or the proportion with detectable expression.
#' These are helpful for any visualizations that operate on individual groups, e.g., \code{\link{plotGroupedHeatmap}}.
#' However, if groups are distributed across multiple batches, some manner of batch correction is required. 
#' The problem with directly averaging group-level statistics across batches is that some groups may not exist in particular batches,
#' e.g., due to the presence of unique cell types in different samples.
#' A direct average would be biased by variable contributions of the batch effect for each group.
#' 
#' To overcome this, we use groups that are present in multiple batches to correct for the batch effect.
#' (That is, any level of \code{groups} that occurs for multiple levels of \code{block}.)
#' For each gene, we fit a linear model to the (transformed) values containing both the group and block factors.
#' We then report the coefficient for each group as the batch-adjusted average for that group;
#' this is possible as the fitted model has no intercept.
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
#' @seealso
#' \code{\link{plotGroupedHeatmap}} and \code{\link{plotDots}}, where this function gets used.
#'
#' \code{regressBatches} from the \pkg{batchelor} package, to remove the batch effect from per-cell expression values.
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=1000)
#' group <- sample(10, ncol(y), replace=TRUE)
#' block <- sample(5, ncol(y), replace=TRUE)
#'
#' library(scuttle)
#' summaries <- summarizeAssayByGroup(y, DataFrame(group=group, block=block), 
#'     statistics=c("mean", "prop.detected"))
#'
#' # Computing batch-aware averages:
#' library(scater)
#' averaged <- batchCorrectedAverages(assay(summaries, "mean"), 
#'     group=summaries$group, block=summaries$block)
#' 
#' num <- batchCorrectedAverages(assay(summaries, "prop.detected"),
#'     group=summaries$group, block=summaries$block, transform="logit") 
#' 
#' @export
#' @importFrom stats lm.fit model.matrix
#' @importFrom Matrix t
batchCorrectedAverages <- function(x, group, block, transform=c("raw", "log", "logit"), offset=NULL) { 
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
