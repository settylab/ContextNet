#' Test helpers for skipping tests based on available dependencies
#'
#' @name test-helpers
#' @keywords internal
NULL

#' Get the path to the local NicheNet data directory
#'
#' @return Path to data directory or NULL if not found
get_nichenet_database <- function() {
    # Check for data directory relative to package root
    pkg_root <- testthat::test_path("..", "..")
    nichenet_database <- file.path(pkg_root, "data")

    if (dir.exists(nichenet_database)) {
        # Check for required files
        required_files <- c(
            "lr_network_human.rds",
            "ligand_target_matrix_nsga2r_final.rds",
            "weighted_networks_nsga2r_final.rds"
        )
        if (all(file.exists(file.path(nichenet_database, required_files)))) {
            return(nichenet_database)
        }
    }
    return(NULL)
}

#' Load NicheNet data from local files
#'
#' @return List with lr_network, ligand_target_matrix, weighted_networks
load_local_nichenet_data <- function() {
    nichenet_database <- get_nichenet_database()
    if (is.null(nichenet_database)) {
        stop("Local NicheNet data files not found")
    }

    list(
        lr_network = readRDS(file.path(nichenet_database, "lr_network_human.rds")),
        ligand_target_matrix = readRDS(file.path(nichenet_database, "ligand_target_matrix_nsga2r_final.rds")),
        weighted_networks = readRDS(file.path(nichenet_database, "weighted_networks_nsga2r_final.rds"))
    )
}

#' Skip if nichenetr is not available or data files are missing
skip_if_no_nichenetr <- function() {
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        testthat::skip("nichenetr not installed")
    }

    # First check for local data files
    nichenet_database <- get_nichenet_database()
    if (!is.null(nichenet_database)) {
        return(invisible(TRUE))  # Local data available
    }

    # Fall back to checking nichenetr package data
    tryCatch({
        lr_net <- nichenetr::lr_network
        wn <- nichenetr::weighted_networks
    }, error = function(e) {
        testthat::skip("nichenetr data files not available")
    })
}

#' Skip if multinichenetr is not available or API is incompatible
skip_if_no_multinichenetr <- function() {
    if (!requireNamespace("multinichenetr", quietly = TRUE)) {
        testthat::skip("multinichenetr not installed")
    }
    # Check for required functions (API v2.1.0+)
    # multi_nichenet_analysis is the main entry point in v2.1.0
    if (!exists("multi_nichenet_analysis", envir = asNamespace("multinichenetr"))) {
        testthat::skip("multinichenetr API incompatible - multi_nichenet_analysis not found")
    }
}

#' Skip if Seurat is not available
skip_if_no_seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        testthat::skip("Seurat not installed")
    }
}

#' Skip full workflow tests if NicheNet data files are unavailable
#'
#' The simulated test data now uses REAL gene names from NicheNet networks.
#' This skip function checks if the NicheNet data files are available locally.
#' If available, the tests will run with real network subsets.
#'
#' @note As of v0.2.0, simulated data uses real gene names, so this skip
#'   is only triggered when NicheNet data files are missing, not by default.
skip_if_no_nichenet_data <- function() {
    nichenet_database <- get_nichenet_database()
    if (is.null(nichenet_database)) {
        testthat::skip("NicheNet data files not found locally")
    }
}

#' Legacy function - kept for backwards compatibility but no longer skips by default
#' @note Deprecated: simulated data now uses real gene names
skip_if_simulated_data <- function() {
    # Check if NicheNet data is available
    nichenet_database <- get_nichenet_database()
    if (is.null(nichenet_database)) {
        testthat::skip("NicheNet data files required for full workflow tests")
    }
    # Otherwise, proceed with the test using real gene names
    invisible(TRUE)
}

#' Skip if contextnet-pipeline is not available
skip_if_no_pipeline <- function() {
    if (!nzchar(Sys.which("snakemake"))) {
        testthat::skip("contextnet-pipeline (snakemake) not available")
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
