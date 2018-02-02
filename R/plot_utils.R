################################################
## Colour-related functions

.get_palette <- function(palette_name) 
# Function to define color palettes.
{
    switch(palette_name,
           tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                         "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                         "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                         "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
           tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                               "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                               "#CDCC5D", "#6DCCDA"),
           colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                            "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                            "#FFBC79", "#CFCFCF"),
           trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
                            "#FFC156", "#69B764", "#F26C64", "#FFDD71",
                            "#9FCD99"),
           purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
                            "#5F5A41", "#B4B19B", "#995688", "#D898BA",
                            "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
           bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
                         "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
                         "#BD0A36", "#F4737A"),
           greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                             "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                             "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
           cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
                      "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
                      "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
                      "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB")
    )
}

.resolve_plot_colours <- function(plot_out, colour_by, colour_by_name, fill = FALSE) 
# Get nice plotting colour schemes for very general colour variables
{
    if ( is.null(colour_by) ) {
        return(plot_out)
    }

    # Picking whether to fill or not.
    if ( fill ) {
        VIRIDFUN <- viridis::scale_fill_viridis
        SCALEFUN <- scale_fill_manual
    } else {
        VIRIDFUN <- viridis::scale_color_viridis
        SCALEFUN <- scale_color_manual
    }

    # Set a sensible colour scheme and return the plot_out object
    if ( is.numeric(colour_by) ) {
        plot_out <- plot_out + VIRIDFUN(name = colour_by_name)
    } else {
        nlevs_colour_by <- nlevels(as.factor(colour_by))
        if (nlevs_colour_by <= 10) {
            plot_out <- plot_out + SCALEFUN(
                values = .get_palette("tableau10medium"),
                name = colour_by_name)
        } else {
            if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                plot_out <- plot_out + SCALEFUN(
                    values = .get_palette("tableau20"),
                    name = colour_by_name)
            } else {
                plot_out <- plot_out + VIRIDFUN(
                    name = colour_by_name, discrete = TRUE)
            }
        }
    }
    plot_out
}

################################################
## Choosing visualization values.

.choose_vis_values <- function(x, by, mode=c("column", "row"), search=c("any", "metadata", "feature"),
                               exprs_values = "logcounts", coerce_factor = FALSE, level_limit = NA,
                               discard_solo = FALSE) 
# This function looks through the visualization data and returns the
# values to be visualized. Either 'by' itself, or a column of colData,
# or a column of rowData, or the expression values of a feature.
{
    vals <- NULL
    if (is.character(by)) {
        mode <- match.arg(mode)
        search <- match.arg(search)

        # Determining what to check, based on input 'by'.
        if (search=="any") { 
            check_metadata <- check_features <- TRUE

            if (length(by)==0) {
                check_metadata <- FALSE 
                check_features <- FALSE
            } else if (length(by)>1) {
                check_metadata <- TRUE
                check_features <- FALSE
            } else if (length(by)==1L) {
                cur_name <- names(by)
                if (!is.null(cur_name) && !is.na(cur_name)) { 
                    if (cur_name=="Metadata") {
                        check_features <- FALSE
                        check_metadata <- TRUE
                    } else if ((mode=="column" && cur_name=="Feature")
                               || (mode=="row" && cur_name=="Cell")) {
                        check_features <- TRUE
                        check_metadata <- FALSE
                    } 
                }
            }
        } else {
            check_metadata <- (search=="metadata")
            check_features <- !check_metadata            
        }
           
        # Checking the metadata; note the loop to account for nesting.
        if (check_metadata) { 
            if (mode=="column") {
                meta_data <- colData(x)
            } else {
                meta_data <- rowData(x)
            }
            
            for (field in by) {
                if (!field %in% colnames(meta_data)) {
                    break
                }
                vals <- meta_data[[field]]
                metadata <- vals
            }
            by <- paste(by, collapse=":") # collapsing to a single string for output.
        }

        # Metadata takes priority, so we don't bother searching if 'vals' is non-NULL.
        if (check_features) {
            if (is.null(vals)){
                exprs <- assay(x, i = exprs_values)
                if (mode=="column") {
                    vals <- exprs[by,] # coloring columns, so we take the row values.
                } else {
                    vals <- exprs[,by]
                }
            }
        }

        if (is.null(vals) && (check_metadata || check_features)) {
            stop("cannot find the supplied '*_by' in features or metadata")
        }
    } else if (is.data.frame(by)) {
        if (ncol(by) != 1L) {
            stop("'*_by' should be a data frame with one column")
        } else if (nrow(by) != ncol(x)) {
            stop("'nrow(*_by)' should be equal to number of columns in 'x'")
        }

        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)

    } else if (!is.null(by)) {
        stop("invalid value of 'by' supplied")
    }

    # Checking the level limit.
    if (coerce_factor && !is.null(vals)) {
        vals <- factor(vals)
        if (level_limit < nlevels(vals)) {
            stop(sprintf("number of unique levels exceeds %i", level_limit))
        }
    }
    
    # If only one level for the variable, set to NULL.
    if (length(unique(vals))<=1L && discard_solo) { 
        by <- NULL
        vals <- NULL
    }
    return(list(name = by, val = vals))
}

################################################
## Integrating shape/colour/size_by choices.

.get_point_args <- function(colour_by, shape_by, size_by, alpha=0.65, size=NULL) 
## Note the use of colour instead of fill when shape_by is set, as not all shapes have fill.
## (Fill is still the default as it looks nicer.)
{
    aes_args <- list()
    fill_colour <- TRUE
    if (!is.null(shape_by)) {
        aes_args$shape <- "shape_by"
        fill_colour <- FALSE
    }
    if (!is.null(colour_by)) {
        if (fill_colour) {
            aes_args$fill <- "colour_by"
        } else {
            aes_args$colour <- "colour_by"
        }
    }
    if (!is.null(size_by)) {
        aes_args$size <- "size_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    
    geom_args <- list(mapping=new_aes, alpha=alpha)
    if (is.null(colour_by) || fill_colour) {
        geom_args$colour <- "grey70"
    }
    if (is.null(colour_by) || !fill_colour) { # set fill when there is no fill colour, to distinguish between e.g., pch=16 and pch=21.
        geom_args$fill <- "grey20"
    }
    if (is.null(shape_by)) {
        geom_args$shape <- 21
    }
    if (is.null(size_by)) {
        geom_args$size <- size
    }
    return(list(args=geom_args, fill=fill_colour))
}

.add_extra_guide <- function(plot_out, shape_by, size_by) 
# Adding extra legend information on the shape and size.
{
    guide_args <- list()
    if (!is.null(shape_by)) {
        guide_args$shape <- guide_legend(title = shape_by)
    }
    if (!is.null(size_by)) { 
        guide_args$size <- guide_legend(title = size_by)
    }
    if (length(guide_args)) { 
        plot_out <- plot_out + do.call(guides, guide_args)
    }
    return(plot_out)
}

################################################
## Creating pair plots.

.makePairs <- function(data_matrix) 
# with thanks to Gaston Sanchez, who posted this code online
# https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
{
    if ( is.null(names(data_matrix)) )
        names(data_matrix) <- paste0("row", 1:nrow(data_matrix))
    exp_grid <- expand.grid(x = 1:ncol(data_matrix), y = 1:ncol(data_matrix))
    exp_grid <- exp_grid[exp_grid$x != exp_grid$y,]
    all_panels <- do.call("rbind", lapply(1:nrow(exp_grid), function(i) {
        xcol <- exp_grid[i, "x"]
        ycol <- exp_grid[i, "y"]
        data.frame(xvar = names(data_matrix)[ycol], yvar = names(data_matrix)[xcol],
                   x = data_matrix[, xcol], y = data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels = names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels = names(data_matrix))
    densities <- do.call("rbind", lapply(1:ncol(data_matrix), function(i) {
        data.frame(xvar = names(data_matrix)[i], yvar = names(data_matrix)[i],
                   x = data_matrix[, i])
    }))
    list(all = all_panels, densities = densities)
}
