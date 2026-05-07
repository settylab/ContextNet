# Tests for network-integrate.R

test_that("integrate_networks requires ContextGRN object", {
    expect_error(
        integrate_networks("not a grn"),
        "ContextGRN"
    )
})

test_that("integrate_networks requires non-empty edges", {
    edges <- create_test_edges(10)
    grn <- ContextGRN(edges)
    grn@edges <- data.frame(from = character(), to = character(), weight = numeric())

    expect_error(
        integrate_networks(grn),
        "no edges"
    )
})

test_that("integrate_networks mode argument validation", {
    edges <- create_test_edges(10)
    grn <- ContextGRN(edges)

    expect_error(
        integrate_networks(grn, mode = "invalid_mode"),
        "arg"
    )
})

test_that(".validate_nichenet_networks validates lr_network", {
    # Not a data.frame
    expect_error(
        ContextNet:::.validate_nichenet_networks(list(), list(lr_sig = data.frame(), gr = data.frame())),
        "data.frame"
    )

    # Missing 'from' column
    expect_error(
        ContextNet:::.validate_nichenet_networks(
            data.frame(to = "A"),
            list(lr_sig = data.frame(), gr = data.frame())
        ),
        "from"
    )
})

test_that(".validate_nichenet_networks validates weighted_networks", {
    lr_net <- data.frame(from = "A", to = "B")

    # Not a list
    expect_error(
        ContextNet:::.validate_nichenet_networks(lr_net, "not a list"),
        "list"
    )

    # Missing components
    expect_error(
        ContextNet:::.validate_nichenet_networks(lr_net, list(lr_sig = data.frame())),
        "missing"
    )
})

test_that(".validate_nichenet_networks accepts valid inputs", {
    lr_net <- data.frame(from = "A", to = "B")
    wn <- list(
        lr_sig = data.frame(from = "A", to = "B", weight = 0.5),
        gr = data.frame(from = "TF1", to = "Gene1", weight = 0.8)
    )

    expect_invisible(ContextNet:::.validate_nichenet_networks(lr_net, wn))
})

test_that(".integrate_weighted_networks replace mode works", {
    weighted_networks <- list(
        lr_sig = data.frame(from = "L1", to = "R1", weight = 0.5),
        gr = data.frame(from = "TF_old", to = "Gene_old", weight = 0.3)
    )

    tf_target <- data.frame(
        from = c("TF1", "TF2"),
        to = c("Gene1", "Gene2"),
        weight = c(0.8, 0.6)
    )

    result <- ContextNet:::.integrate_weighted_networks(
        weighted_networks = weighted_networks,
        tf_target = tf_target,
        mode = "replace",
        weights = c(grn = 0.5, nichenet = 0.5),
        verbose = FALSE
    )

    expect_equal(result$gr, tf_target)
    expect_equal(result$lr_sig, weighted_networks$lr_sig)
})

test_that(".integrate_weighted_networks augment mode combines edges", {
    weighted_networks <- list(
        lr_sig = data.frame(from = "L1", to = "R1", weight = 0.5),
        gr = data.frame(
            from = c("TF_old", "TF_shared"),
            to = c("Gene_old", "Gene_shared"),
            weight = c(0.3, 0.4)
        )
    )

    tf_target <- data.frame(
        from = c("TF_new", "TF_shared"),
        to = c("Gene_new", "Gene_shared"),
        weight = c(0.8, 0.9)  # TF_shared->Gene_shared should override
    )

    result <- ContextNet:::.integrate_weighted_networks(
        weighted_networks = weighted_networks,
        tf_target = tf_target,
        mode = "augment",
        weights = c(grn = 0.5, nichenet = 0.5),
        verbose = FALSE
    )

    # Should have 3 unique edges (TF_new, TF_shared from custom, TF_old from original)
    expect_equal(nrow(result$gr), 3)

    # TF_shared->Gene_shared should have custom GRN weight (0.9), not original (0.4)
    shared_edge <- result$gr[result$gr$from == "TF_shared" & result$gr$to == "Gene_shared", ]
    expect_equal(shared_edge$weight, 0.9)
})

test_that(".integrate_weighted_networks weighted mode combines with weights", {
    weighted_networks <- list(
        lr_sig = data.frame(from = "L1", to = "R1", weight = 0.5),
        gr = data.frame(
            from = c("TF1", "TF2"),
            to = c("Gene1", "Gene2"),
            weight = c(0.4, 0.6)  # NicheNet weights
        )
    )

    tf_target <- data.frame(
        from = c("TF1", "TF3"),
        to = c("Gene1", "Gene3"),
        weight = c(0.8, 0.5)  # Custom GRN weights
    )

    weights <- c(grn = 0.6, nichenet = 0.4)

    result <- ContextNet:::.integrate_weighted_networks(
        weighted_networks = weighted_networks,
        tf_target = tf_target,
        mode = "weighted",
        weights = weights,
        verbose = FALSE
    )

    # TF1->Gene1 should be weighted combination: 0.6 * 0.8 + 0.4 * 0.4 = 0.64
    tf1_edge <- result$gr[result$gr$from == "TF1" & result$gr$to == "Gene1", ]
    expect_equal(tf1_edge$weight, 0.64, tolerance = 1e-6)

    # TF2->Gene2 only in NicheNet: 0.6 * 0 + 0.4 * 0.6 = 0.24
    tf2_edge <- result$gr[result$gr$from == "TF2" & result$gr$to == "Gene2", ]
    expect_equal(tf2_edge$weight, 0.24, tolerance = 1e-6)

    # TF3->Gene3 only in custom GRN: 0.6 * 0.5 + 0.4 * 0 = 0.3
    tf3_edge <- result$gr[result$gr$from == "TF3" & result$gr$to == "Gene3", ]
    expect_equal(tf3_edge$weight, 0.3, tolerance = 1e-6)
})

test_that("as_nichenet_input requires ContextNetwork object", {
    expect_error(
        as_nichenet_input("not a network"),
        "ContextNetwork"
    )
})

test_that("as_nichenet_input returns correct structure", {
    # Create a minimal ContextNetwork object for testing
    ltm <- create_test_ligand_target_matrix(10, 20)
    sig_net <- data.frame(from = "L1", to = "R1", weight = 0.5)
    reg_net <- data.frame(from = "TF1", to = "Gene1", weight = 0.8)

    network <- new("ContextNetwork",
        ligand_target_matrix = ltm,
        ligand_receptor_matrix = matrix(nrow = 0, ncol = 0),
        signaling_network = sig_net,
        regulatory_network = reg_net,
        weighted_networks = list(lr_sig = sig_net, gr = reg_net),
        integration_params = list(mode = "replace"),
        source_grn = new("ContextGRN"),
        metadata = list(),
        provenance = list(source = "test")
    )

    result <- as_nichenet_input(network)

    expect_type(result, "list")
    expect_true("ligand_target_matrix" %in% names(result))
    expect_true("lr_network" %in% names(result))
    expect_true("sig_network" %in% names(result))
    expect_true("weighted_networks" %in% names(result))

    expect_equal(result$ligand_target_matrix, ltm)
})

# Integration tests that require nichenetr
test_that("integrate_networks works with nichenetr", {
    skip_if_no_nichenetr()

    edges <- create_test_edges(100)
    grn <- ContextGRN(edges)

    # This test requires actual NicheNet data
    network <- integrate_networks(
        grn,
        nichenet_database = get_nichenet_database(),
        verbose = FALSE
    )

    expect_s4_class(network, "ContextNetwork")
    expect_true(nrow(network@ligand_target_matrix) > 0)
    expect_true(ncol(network@ligand_target_matrix) > 0)
})

test_that("integrate_networks respects integration mode", {
    skip_if_no_nichenetr()

    edges <- create_test_edges(100)
    grn <- ContextGRN(edges)

    network_replace <- integrate_networks(
        grn,
        mode = "replace",
        nichenet_database = get_nichenet_database(),
        verbose = FALSE
    )
    expect_equal(network_replace@integration_params$mode, "replace")

    network_augment <- integrate_networks(
        grn,
        mode = "augment",
        nichenet_database = get_nichenet_database(),
        verbose = FALSE
    )
    expect_equal(network_augment@integration_params$mode, "augment")
})

test_that("integrate_networks stores parameters correctly", {
    skip_if_no_nichenetr()

    edges <- create_test_edges(100)
    grn <- ContextGRN(edges)

    network <- integrate_networks(
        grn,
        mode = "weighted",
        weights = c(grn = 0.7, nichenet = 0.3),
        nichenet_database = get_nichenet_database(),
        verbose = FALSE
    )

    expect_equal(network@integration_params$mode, "weighted")
    expect_equal(network@integration_params$weights, c(grn = 0.7, nichenet = 0.3))
    expect_true(!is.null(network@integration_params$damping_factor))
    expect_true(!is.null(network@integration_params$ltf_cutoff))
    expect_equal(network@integration_params$algorithm, "PPR")
})
