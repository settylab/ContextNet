test_that("get_config returns expected structure", {
    config <- get_config()

    expect_type(config, "list")
    expect_true("cache_dir" %in% names(config))
    expect_true("parallel_backend" %in% names(config))
    expect_true("n_cores" %in% names(config))
    expect_true("version" %in% names(config))
})

test_that("reset_config resets to defaults", {
    # Set some custom values
    options(contextnet.log_level = "debug")

    reset_config()

    expect_equal(getOption("contextnet.log_level"), "info")
    expect_equal(getOption("contextnet.parallel_backend"), "serial")
})

test_that("setup_logging changes log level", {
    old_level <- setup_logging("debug")
    expect_equal(getOption("contextnet.log_level"), "debug")

    setup_logging("info")
    expect_equal(getOption("contextnet.log_level"), "info")
})

test_that("check_ecosystem returns status list", {
    status <- check_ecosystem(verbose = FALSE)

    expect_type(status, "list")
    expect_true("contextnet" %in% names(status))
    expect_true(status$contextnet$installed)
})
