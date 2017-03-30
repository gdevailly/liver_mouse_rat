library(cowplot)

plotMeanVariance <- function(tb, cutoff = 1, title = NULL) {
    tb <- mutate(tb, grouping = if_else(sd >= cutoff, TRUE, FALSE))
    mp <- ggplot(tb, aes(x = sd, y = mean)) +
        geom_point(mapping = aes(color = grouping), size = 0.1) +
        geom_density2d(color = "orange") +
        scale_color_manual(values = c("grey","black")) +
        geom_text(
            data = data_frame(
                cat = c("variable", "non variable"),
                eff = c(length(which(tb$grouping)), length(which(!tb$grouping))),
                xpo = c(4.5, 4.5),
                ypo = c(14.5, 13.5),
                gro = c(TRUE, FALSE)
            ),
            mapping = aes(x = xpo, y = ypo, label = eff, color = gro),
            hjust = "right"
        ) +
        theme(legend.position = "none") +
        labs(title = title, x = "Standard deviation", y = "Mean expression") +
        coord_cartesian(xlim = c(0, 4.5))

    return(mp)
}



