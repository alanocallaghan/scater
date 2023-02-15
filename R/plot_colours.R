.get_palette <- function(palette_name) 
# Function to define colour palettes.
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
           colourblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
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

#' @importFrom ggplot2 scale_fill_manual scale_colour_manual
#' @importFrom viridis scale_fill_viridis scale_colour_viridis
.resolve_plot_colours <- function(plot_out, colour_by, colour_by_name, fill = FALSE, colour = FALSE) 
# Get nice plotting colour schemes for very general colour variables
{
    if (is.null(colour_by)) {
        return(plot_out)
    }

    # Picking whether to fill or not.
    aesthetics <- c("fill", "colour")[c(fill, colour)]

    if (fill) {
        VIRIDFUN <- scale_fill_viridis
        SCALEFUN <- scale_fill_manual
    } else if (colour) {
        VIRIDFUN <- scale_colour_viridis
        SCALEFUN <- scale_colour_manual
    }
    
    # Set a sensible colour scheme and return the plot_out object
    if (is.numeric(colour_by)) {
        plot_out <- plot_out + VIRIDFUN(
            name = colour_by_name, aesthetics = aesthetics
        )
    } else {
        nlevs_colour_by <- nlevels(as.factor(colour_by))
        if (nlevs_colour_by <= 10) {
            plot_out <- plot_out + SCALEFUN(
                values = .get_palette("tableau10medium"),
                name = colour_by_name,
                aesthetics = aesthetics)
        } else {
            if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                plot_out <- plot_out + SCALEFUN(
                    values = .get_palette("tableau20"),
                    name = colour_by_name,
                    aesthetics = aesthetics)
            } else {
                plot_out <- plot_out + VIRIDFUN(
                    name = colour_by_name, discrete = TRUE
                    # , aesthetics = aesthetics
                    )
            }
        }
    }
    plot_out
}
