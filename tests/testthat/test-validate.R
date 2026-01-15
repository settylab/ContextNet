test_that("validate_integration_mode accepts valid modes", {
    expect_invisible(validate_integration_mode("weighted"))
    expect_invisible(validate_integration_mode("replace"))
    expect_invisible(validate_integration_mode("augment"))
    expect_invisible(validate_integration_mode("hierarchical"))
})

test_that("validate_integration_mode rejects invalid modes", {
    expect_error(validate_integration_mode("invalid"), "Invalid integration mode")
})

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

test_that("validate_geneset_oi rejects missing geneset", {
    expect_error(validate_geneset_oi(NULL), "required")
    expect_error(validate_geneset_oi(character()), "empty")
})

test_that("validate_geneset_oi accepts valid geneset", {
    genes <- c("Gene1", "Gene2", "Gene3")
    expect_invisible(validate_geneset_oi(genes))
})

test_that("validate_path checks file existence", {
    expect_error(
        validate_path("/nonexistent/path", must_exist = TRUE),
        "does not exist"
    )
})

test_that("validate_positive_integer works correctly", {
    expect_invisible(validate_positive_integer(5, "test"))
    expect_error(validate_positive_integer(-1, "test"), "between")
    expect_error(validate_positive_integer(1.5, "test"), "integer")
    expect_error(validate_positive_integer("a", "test"), "number")
})

test_that("validate_grn_edges accepts valid edges", {
    edges <- data.frame(
        from = c("TF1", "TF2"),
        to = c("Gene1", "Gene2"),
        weight = c(0.5, 0.8)
    )
    expect_invisible(validate_grn_edges(edges))
})

test_that("validate_grn_edges rejects non-dataframe", {
    expect_error(validate_grn_edges("not a df"), "data.frame")
})

test_that("validate_grn_edges rejects empty dataframe", {
    expect_error(validate_grn_edges(data.frame()), "empty")
})

test_that("validate_context_grn accepts valid GRN", {
    edges <- create_test_edges(50)
    grn <- ContextGRN(edges)
    expect_invisible(validate_context_grn(grn))
})

test_that("validate_context_grn rejects non-GRN objects", {
    expect_error(validate_context_grn("not a grn"), "ContextGRN")
})
