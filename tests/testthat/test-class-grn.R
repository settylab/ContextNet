test_that("ContextGRN can be created from edges", {
    edges <- create_test_edges(100)
    grn <- ContextGRN(edges)

    expect_s4_class(grn, "ContextGRN")
    expect_equal(nrow(grn@edges), 100)
    expect_true(length(tfs(grn)) > 0)
    expect_true(length(targets(grn)) > 0)
})

test_that("ContextGRN validates edge weights", {
    edges <- data.frame(
        from = c("TF1", "TF2"),
        to = c("Gene1", "Gene2"),
        weight = c(0.5, 1.5)  # Invalid: weight > 1
    )

    expect_error(ContextGRN(edges), "range")
})

test_that("ContextGRN normalizes column names", {
    # Using alternative column names
    edges <- data.frame(
        tf = c("TF1", "TF2"),
        target = c("Gene1", "Gene2"),
        score = c(0.5, 0.8)
    )

    grn <- ContextGRN(edges)
    expect_true("from" %in% colnames(grn@edges))
    expect_true("to" %in% colnames(grn@edges))
    expect_true("weight" %in% colnames(grn@edges))
})

test_that("ContextGRN builds regulons correctly", {
    edges <- data.frame(
        from = c("TF1", "TF1", "TF2"),
        to = c("Gene1", "Gene2", "Gene1"),
        weight = c(0.5, 0.6, 0.7)
    )

    grn <- ContextGRN(edges)
    regs <- regulons(grn)

    expect_true("TF1" %in% names(regs))
    expect_true("TF2" %in% names(regs))
    expect_equal(length(regs$TF1), 2)
    expect_equal(length(regs$TF2), 1)
})

test_that("ContextGRN show method works", {
    edges <- create_test_edges(50)
    grn <- ContextGRN(edges)

    # Verify show() runs without error and returns invisible(NULL)
    # Note: cli functions write directly to console, not captured by capture.output()
    result <- expect_invisible(show(grn))
    expect_null(result)
})

test_that("ContextGRN can be subset by TF", {
    edges <- data.frame(
        from = c("TF1", "TF1", "TF2", "TF2"),
        to = c("Gene1", "Gene2", "Gene3", "Gene4"),
        weight = c(0.5, 0.6, 0.7, 0.8)
    )

    grn <- ContextGRN(edges)
    subset_grn <- grn["TF1"]

    expect_equal(nrow(subset_grn@edges), 2)
    expect_equal(unique(subset_grn@edges$from), "TF1")
})

test_that("ContextGRN length method works", {
    edges <- create_test_edges(75)
    grn <- ContextGRN(edges)

    expect_equal(length(grn), 75)
})

test_that("ContextGRN as.data.frame returns edges", {
    edges <- create_test_edges(50)
    grn <- ContextGRN(edges)
    df <- as.data.frame(grn)

    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 50)
})

test_that("ContextGRN provenance is set correctly", {
    edges <- create_test_edges(10)
    grn <- ContextGRN(edges)
    prov <- provenance(grn)

    expect_true("created" %in% names(prov))
    expect_true("package_version" %in% names(prov))
    expect_equal(prov$source, "manual")
})
