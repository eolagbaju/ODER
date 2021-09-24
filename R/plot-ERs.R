#' Plot Expressed regions
#'
#' Plots the median deltas and the number of ERs with a delta of 0 against the
#' MCCs on two separate graphs with a line for each of the various MRGs.
#'
#' @param ers_delta tibble/dataframe containing summarised delta values. One row
#'  per set of ERs.
#' @param opt_mcc_mrg vector containing the optimum mcc and mrg, in that order
#'
#' @return Plot of MCC against median delta and number of ERS with a delta of 0
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 element_text
#'
#' @examples
#' # gtex_SRP012682_SRX222703_lung_erdelta_1 is from the package data folder
#' eg_plots <- plot_ers(ers_delta = gtex_SRP012682_SRX222703_lung_erdelta_1, opt_mcc_mrg = c(
#'     "mcc_10",
#'     "mrg_20"
#' ))
#'
#' eg_plots
plot_ers <- function(ers_delta, opt_mcc_mrg) {
    if (missing(ers_delta)) {
        stop("No ER deltas were entered")
    } else if (missing(opt_mcc_mrg)) {
        stop("The optimum MCC and/or MRG were not entered")
    }
    opt_mcc <- as.double(stringr::str_remove(
        opt_mcc_mrg[[1]],
        stringr::fixed("mcc_")
    ))
    opt_mrg <- as.double(stringr::str_remove(
        opt_mcc_mrg[[2]],
        stringr::fixed("mrg_")
    ))
    opt_median <- ers_delta[["median"]][
        ers_delta[["mcc"]] == opt_mcc & ers_delta[["mrg"]] == opt_mrg
    ]
    opt_n_eq_0 <- ers_delta[["n_eq_0"]][
        ers_delta[["mcc"]] == opt_mcc & ers_delta[["mrg"]] == opt_mrg
    ]
    mrgs <- unique(ers_delta[["mrg"]])
    clr_num <- length(mrgs)
    maxgaps_colours <- data.frame(
        maxgap = mrgs,
        colours = ggpubr::get_palette("Blues", clr_num)
    ) %>%
        dplyr::mutate(colours = ifelse(maxgap == opt_mrg, "red", colours))

    exon_delta_median_plot <- plot_exon_delta_median(
        mrg, mcc, opt_mcc, ers_delta,
        median, opt_median,
        maxgaps_colours
    )
    num_exon_delta_eq_0_plot <- plot_num_exon_delta_eq_0(
        mrg, mcc, opt_mcc,
        ers_delta, n_eq_0,
        opt_n_eq_0,
        maxgaps_colours
    )
    optimised_exon_delta_plots <- ggpubr::ggarrange(exon_delta_median_plot,
        num_exon_delta_eq_0_plot,
        nrow = 2, ncol = 1, align = "v",
        common.legend = TRUE, legend = "right"
    )
    return(optimised_exon_delta_plots)
}

#' Plots the optimum exon deltas for the various MCCs and MRGs
#'
#' @param mrg Max Region Gap
#' @param mcc Mean Cutoff Coverages
#' @param opt_mcc the optimum mean cutoff coverage
#' @param ers_delta the er deltas
#' @param median the median exon deltas
#' @param opt_median the median exon delta produced by the optimum MCC and MRG
#' @param maxgaps_colours colours to use for the plot lines
#'
#' @return exon_delta_median_plot
#' @keywords internal
#' @noRd
plot_exon_delta_median <- function(mrg, mcc, opt_mcc, ers_delta,
    median, opt_median, maxgaps_colours) {
    exon_delta_median_plot <- ggplot(
        data = ers_delta,
        mapping = aes(x = mcc, y = median)
    ) +
        geom_line(aes(colour = as.factor(mrg))) +
        geom_vline(xintercept = opt_mcc, colour = "#177D87", linetype = 2) +
        geom_point(aes(x = opt_mcc, y = opt_median),
            colour = "black", shape = 4
        ) +
        ggrepel::geom_text_repel(
            data = dplyr::tibble(
                mcc = opt_mcc,
                median = opt_median
            ), aes(x = mcc, y = median), label = "optimum",
            min.segment.length = 0, force_pull = 0.5
        ) +
        scale_x_continuous(name = "MCC") +
        scale_y_continuous(name = expression("Median" ~ Delta)) +
        scale_colour_manual("MRG", values = maxgaps_colours$colours) +
        ggpubr::theme_pubr(legend = "right") +
        theme(
            legend.title = element_text(colour = "red"),
            axis.title.x = element_text(colour = "#177D87")
        )

    return(exon_delta_median_plot)
}

#' Plots the optimum exon deltas for the various MCCs and MRGs
#'
#' @param mrg Max Region Gap
#' @param mcc Mean Cutoff Coverages
#' @param opt_mcc the optimum mean cutoff coverage
#' @param ers_delta the er deltas
#' @param n_eq_0 the number of ers with an exon delta of 0
#' @param opt_n_eq_0 the number of ers with an exon delta of 0 for the optimum
#'  MCC and MRG
#' @param maxgaps_colours colours to use for the plot lines
#'
#' @return num_exon_delta_eq_0_plot
#' @keywords internal
#' @noRd
plot_num_exon_delta_eq_0 <- function(mrg, mcc, opt_mcc, ers_delta,
    n_eq_0, opt_n_eq_0, maxgaps_colours) {
    num_exon_delta_eq_0_plot <- ggplot(
        data = ers_delta,
        mapping = aes(x = mcc, y = n_eq_0)
    ) +
        geom_line(aes(colour = as.factor(mrg))) +
        geom_vline(xintercept = opt_mcc, colour = "#177D87", linetype = 2) +
        geom_point(aes(x = opt_mcc, y = opt_n_eq_0),
            colour = "black",
            shape = 4
        ) +
        scale_x_continuous(name = "MCC") +
        scale_y_continuous(
            name = expression(
                "Number of ERs with " ~ Delta ~ "= 0"
            )
        ) +
        scale_colour_manual("MRG", values = maxgaps_colours$colours) +
        ggpubr::theme_pubr(legend = "right") +
        theme(
            legend.title = element_text(colour = "red"),
            axis.title.x = element_text(colour = "#177D87")
        ) +
        ggrepel::geom_text_repel(
            data = dplyr::tibble(
                mcc = opt_mcc, n_eq_0 = opt_n_eq_0
            ), aes(x = mcc, y = n_eq_0),
            label = "optimum", min.segment.length = 0, force_pull = 0.5
        )

    return(num_exon_delta_eq_0_plot)
}
