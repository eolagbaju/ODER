test_plots <- plot_ers(
    ers_delta = gtex_lung_erdelta_1, # gtex_lung_erdelta_1 is from the data folder
    opt_mcc_mrg = c("mcc_10", "mrg_20")
)

test_that("plot_ers works", {
    expect_error(
        plot_ers(
            ers_delta = gtex_lung_erdelta_1
        ),
        "The optimum MCC and/or MRG were not entered"
    )
    expect_error(
        plot_ers(
            opt_mcc_mrg = c("mcc_10", "mrg_20")
        ),
        "No ER deltas were entered"
    )

    expect_true(methods::is(test_plots, "ggarrange"))
})
