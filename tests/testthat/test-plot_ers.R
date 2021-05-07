test_plots <- plot_ers(
    ers_delta = ers_delta_example, # ers_delta_example is from the data folder
    opt_mcc_mrg = c("mcc_10", "mrg_20")
)

test_that("plot_ers works", {
    expect_error(
        plot_ers(
            ers_delta = ers_delta_example
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
