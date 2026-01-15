#' Test helpers for skipping tests based on available dependencies
#'
#' @name test-helpers
#' @keywords internal
NULL

#' Skip if nichenetr is not available or data files are missing
skip_if_no_nichenetr <- function() {
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        testthat::skip("nichenetr not installed")
    }
    # Check if nichenetr data is available
    tryCatch({
        lr_net <- nichenetr::lr_network
        wn <- nichenetr::weighted_networks
    }, error = function(e) {
        testthat::skip("nichenetr data files not available")
    })
}

#' Skip if multinichenetr is not available
skip_if_no_multinichenetr <- function() {
    if (!requireNamespace("multinichenetr", quietly = TRUE)) {
        testthat::skip("multinichenetr not installed")
    }
}

#' Skip if Seurat is not available
skip_if_no_seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        testthat::skip("Seurat not installed")
    }
}

#' Skip if no AnnData backend is available
skip_if_no_anndata <- function() {
    backend <- ContextNet:::.detect_anndata_backend()
    if (backend == "none") {
        testthat::skip("No AnnData backend available")
    }
}

#' Skip if contextnet-pipeline is not available
skip_if_no_pipeline <- function() {
    if (!nzchar(Sys.which("snakemake"))) {
        testthat::skip("contextnet-pipeline (snakemake) not available")
    }
}

#' Skip if contextnet-plot is not available
skip_if_no_plotting <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        testthat::skip("reticulate not installed")
    }
    if (!reticulate::py_module_available("contextnet_plot")) {
        testthat::skip("contextnet-plot not available")
    }
}

#' Create example edge data.frame for testing
#' @param n_edges Number of edges to create
#' @return A data.frame with from, to, weight columns
create_test_edges <- function(n_edges = 100) {
    tfs <- paste0("TF", 1:10)
    targets <- paste0("Gene", 1:50)

    data.frame(
        from = sample(tfs, n_edges, replace = TRUE),
        to = sample(targets, n_edges, replace = TRUE),
        weight = runif(n_edges, 0.1, 1.0),
        stringsAsFactors = FALSE
    )
}

#' Create example ligand-target matrix for testing
#' @param n_ligands Number of ligands
#' @param n_targets Number of targets
#' @return A matrix
create_test_ligand_target_matrix <- function(n_ligands = 50, n_targets = 200) {
    mat <- matrix(
        runif(n_ligands * n_targets, 0, 1),
        nrow = n_targets,
        ncol = n_ligands
    )
    rownames(mat) <- paste0("Target", 1:n_targets)
    colnames(mat) <- paste0("Ligand", 1:n_ligands)
    mat
}
