test_that("validate_weights accepts valid weights", {
    weights <- c(grn = 0.5, nichenet = 0.5)
    result <- validate_weights(weights)
    expect_equal(sum(result), 1)
})

test_that("validate_weights rejects out-of-range weights", {
    expect_error(validate_weights(c(-0.1, 0.5)), "range")
    expect_error(validate_weights(c(0.5, 1.5)), "range")
})

test_that("validate_weights normalizes if needed", {
    weights <- c(grn = 0.3, nichenet = 0.3)
    expect_warning(result <- validate_weights(weights), "Normalizing")
    expect_equal(sum(result), 1)
})

test_that("validate_path checks file existence", {
    expect_error(
        ContextNet:::validate_path("/nonexistent/path", must_exist = TRUE),
        "does not exist"
    )
})
