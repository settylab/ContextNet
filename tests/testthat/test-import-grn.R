# Tests for import-grn.R

test_that("import_grn_pipeline reads CSV correctly", {
    # Create temporary test CSV
    tmp_file <- tempfile(fileext = ".csv")

    # Create a gene_x_tf matrix (genes as rows, TFs as columns)
    gene_x_tf <- data.frame(
        row.names = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
        TF1 = c(0.8, 0.5, 0.0, 0.3, 0.0),
        TF2 = c(0.0, 0.6, 0.9, 0.0, 0.4),
        TF3 = c(0.2, 0.0, 0.0, 0.7, 0.5)
    )
    write.csv(gene_x_tf, tmp_file)

    on.exit(unlink(tmp_file))

    grn <- import_grn_pipeline(tmp_file, verbose = FALSE)

    expect_s4_class(grn, "ContextGRN")
    expect_true(nrow(grn@edges) > 0)
    expect_true(all(c("from", "to", "weight") %in% colnames(grn@edges)))
})

test_that("import_grn_pipeline normalizes weights correctly", {
    tmp_file <- tempfile(fileext = ".csv")

    gene_x_tf <- data.frame(
        row.names = c("Gene1", "Gene2"),
        TF1 = c(10, 5),  # Max is 10, so should be normalized to 1, 0.5
        TF2 = c(2, 4)    # Max is 4, so should be normalized to 0.5, 1
    )
    write.csv(gene_x_tf, tmp_file)

    on.exit(unlink(tmp_file))

    grn <- import_grn_pipeline(tmp_file, normalize = TRUE, verbose = FALSE)

    # Check that max weight per TF is 1
    edges <- grn@edges
    tf1_edges <- edges[edges$from == "TF1", ]
    tf2_edges <- edges[edges$from == "TF2", ]

    expect_equal(max(tf1_edges$weight), 1, tolerance = 1e-6)
    expect_equal(max(tf2_edges$weight), 1, tolerance = 1e-6)
})

test_that("import_grn_pipeline removes zero-weight edges", {
    tmp_file <- tempfile(fileext = ".csv")

    gene_x_tf <- data.frame(
        row.names = c("Gene1", "Gene2", "Gene3"),
        TF1 = c(0.8, 0.0, 0.0),  # Two zeros
        TF2 = c(0.0, 0.6, 0.0)   # Two zeros
    )
    write.csv(gene_x_tf, tmp_file)

    on.exit(unlink(tmp_file))

    grn <- import_grn_pipeline(tmp_file, normalize = FALSE, verbose = FALSE)

    # Should only have non-zero edges
    expect_true(all(grn@edges$weight > 0))
    expect_equal(nrow(grn@edges), 2)  # Only Gene1-TF1 and Gene2-TF2
})

test_that("import_grn_pipeline applies min_weight filter", {
    tmp_file <- tempfile(fileext = ".csv")

    gene_x_tf <- data.frame(
        row.names = c("Gene1", "Gene2", "Gene3"),
        TF1 = c(0.8, 0.3, 0.1)
    )
    write.csv(gene_x_tf, tmp_file)

    on.exit(unlink(tmp_file))

    grn <- import_grn_pipeline(tmp_file, normalize = FALSE, min_weight = 0.25, verbose = FALSE)

    # Should only have edges >= 0.25
    expect_true(all(grn@edges$weight >= 0.25))
    expect_equal(nrow(grn@edges), 2)  # Gene1 (0.8) and Gene2 (0.3)
})

test_that("import_grn_pipeline errors on non-existent file", {
    expect_error(
        import_grn_pipeline("/nonexistent/file.csv", verbose = FALSE),
        "not found"
    )
})

test_that("import_grn_table accepts valid data.frame", {
    edges <- data.frame(
        from = c("TF1", "TF1", "TF2"),
        to = c("Gene1", "Gene2", "Gene1"),
        weight = c(0.8, 0.6, 0.7)
    )

    grn <- import_grn_table(edges, verbose = FALSE)

    expect_s4_class(grn, "ContextGRN")
    expect_equal(nrow(grn@edges), 3)
})

test_that("import_grn_table auto-detects column names", {
    edges <- data.frame(
        TF = c("TF1", "TF2"),
        Target = c("Gene1", "Gene2"),
        Score = c(0.8, 0.6)
    )

    grn <- import_grn_table(edges, verbose = FALSE)

    expect_s4_class(grn, "ContextGRN")
    expect_equal(nrow(grn@edges), 2)
})

test_that("import_grn_table accepts custom column names", {
    edges <- data.frame(
        regulator = c("TF1", "TF2"),
        downstream = c("Gene1", "Gene2"),
        importance = c(0.8, 0.6)
    )

    grn <- import_grn_table(
        edges,
        from_col = "regulator",
        to_col = "downstream",
        weight_col = "importance",
        verbose = FALSE
    )

    expect_s4_class(grn, "ContextGRN")
    expect_equal(nrow(grn@edges), 2)
})

test_that("import_grn_table rejects negative weights", {
    edges <- data.frame(
        from = c("TF1", "TF2"),
        to = c("Gene1", "Gene2"),
        weight = c(0.8, -0.5)  # Negative weight
    )

    expect_error(
        import_grn_table(edges, verbose = FALSE),
        "Negative"
    )
})

test_that("import_grn_table removes NA weights with warning", {
    edges <- data.frame(
        from = c("TF1", "TF2", "TF3"),
        to = c("Gene1", "Gene2", "Gene3"),
        weight = c(0.8, NA, 0.6)
    )

    expect_warning(
        grn <- import_grn_table(edges, verbose = FALSE),
        "NA"
    )

    expect_equal(nrow(grn@edges), 2)
})

test_that("import_grn_table normalizes weights when requested", {
    edges <- data.frame(
        from = c("TF1", "TF1", "TF2"),
        to = c("Gene1", "Gene2", "Gene1"),
        weight = c(10, 5, 8)  # Will be normalized per TF
    )

    grn <- import_grn_table(edges, normalize = TRUE, verbose = FALSE)

    # After normalization, max per TF should be 1
    tf1_edges <- grn@edges[grn@edges$from == "TF1", ]
    expect_equal(max(tf1_edges$weight), 1, tolerance = 1e-6)
})

test_that("import_grn_table scales weights > 1 to [0,1]", {
    edges <- data.frame(
        from = c("TF1", "TF2"),
        to = c("Gene1", "Gene2"),
        weight = c(100, 50)  # Max is 100, will be scaled
    )

    grn <- import_grn_table(edges, normalize = FALSE, verbose = FALSE)

    # Should be scaled to [0, 1]
    expect_true(all(grn@edges$weight <= 1))
    expect_true(all(grn@edges$weight >= 0))
})

test_that("import_grn_table rejects empty data.frame", {
    expect_error(
        import_grn_table(data.frame(), verbose = FALSE),
        "empty"
    )
})

test_that("import_grn auto-detects table format", {
    edges <- create_test_edges(50)

    grn <- import_grn(edges, verbose = FALSE)

    expect_s4_class(grn, "ContextGRN")
})

test_that("import_grn auto-detects pipeline format for CSV files", {
    tmp_file <- tempfile(fileext = ".csv")

    gene_x_tf <- data.frame(
        row.names = c("Gene1", "Gene2"),
        TF1 = c(0.8, 0.5)
    )
    write.csv(gene_x_tf, tmp_file)

    on.exit(unlink(tmp_file))

    grn <- import_grn(tmp_file, verbose = FALSE)

    expect_s4_class(grn, "ContextGRN")
})

test_that("import_grn_scenic is not yet implemented", {
    expect_error(
        import_grn_scenic("/path/to/scenic.h5ad"),
        "not yet implemented"
    )
})

test_that("import_grn_pando is not yet implemented", {
    expect_error(
        import_grn_pando("/path/to/pando"),
        "not yet implemented"
    )
})

test_that(".detect_edge_columns finds standard column names", {
    edges <- data.frame(
        from = c("TF1"),
        to = c("Gene1"),
        weight = c(0.8)
    )

    detected <- ContextNet:::.detect_edge_columns(edges)

    expect_equal(detected$from, "from")
    expect_equal(detected$to, "to")
    expect_equal(detected$weight, "weight")
})

test_that(".detect_edge_columns finds alternative column names", {
    edges <- data.frame(
        TF = c("TF1"),
        Target = c("Gene1"),
        Score = c(0.8)
    )

    detected <- ContextNet:::.detect_edge_columns(edges)

    expect_equal(detected$from, "TF")
    expect_equal(detected$to, "Target")
    expect_equal(detected$weight, "Score")
})

test_that(".detect_edge_columns errors on missing columns", {
    edges <- data.frame(
        unknown1 = c("TF1"),
        unknown2 = c("Gene1"),
        unknown3 = c(0.8)
    )

    expect_error(
        ContextNet:::.detect_edge_columns(edges),
        "Cannot detect"
    )
})
