test_that("ContextNetwork can be created", {
    mat <- create_test_ligand_target_matrix(20, 50)
    network <- ContextNetwork(mat)

    expect_s4_class(network, "ContextNetwork")
    expect_equal(ncol(network@ligand_target_matrix), 20)
    expect_equal(nrow(network@ligand_target_matrix), 50)
})

test_that("ContextNetwork requires a matrix", {
    expect_error(
        ContextNetwork("not a matrix"),
        "must be a matrix"
    )
})

test_that("ContextNetwork ligands accessor works", {
    mat <- create_test_ligand_target_matrix(15, 30)
    network <- ContextNetwork(mat)

    ligs <- ligands(network)
    expect_equal(length(ligs), 15)
    expect_equal(ligs, colnames(mat))
})

test_that("ContextNetwork targets accessor works", {
    mat <- create_test_ligand_target_matrix(15, 30)
    network <- ContextNetwork(mat)

    tgts <- targets(network)
    expect_equal(length(tgts), 30)
    expect_equal(tgts, rownames(mat))
})

test_that("ContextNetwork integration_params accessor works", {
    mat <- create_test_ligand_target_matrix()
    params <- list(mode = "weighted", weights = c(grn = 0.6, nichenet = 0.4))
    network <- ContextNetwork(mat, integration_params = params)

    ip <- integration_params(network)
    expect_equal(ip$mode, "weighted")
    expect_equal(ip$weights, c(grn = 0.6, nichenet = 0.4))
})

test_that("ContextNetwork can store source GRN", {
    edges <- create_test_edges(50)
    grn <- ContextGRN(edges)

    mat <- create_test_ligand_target_matrix()
    network <- ContextNetwork(mat, source_grn = grn)

    expect_s4_class(network@source_grn, "ContextGRN")
})

test_that("ContextNetwork show method works", {
    mat <- create_test_ligand_target_matrix()
    network <- ContextNetwork(mat)

    # Verify show() runs without error and returns invisible(NULL)
    # Note: cli functions write directly to console, not captured by capture.output()
    result <- expect_invisible(show(network))
    expect_null(result)
})
