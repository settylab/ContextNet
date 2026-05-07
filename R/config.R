#' Configuration and Setup
#'
#' @description
#' Functions for configuring ContextNet package settings and checking
#' ecosystem installation status.
#'
#' @name config
NULL


#' Configure ContextNet package settings
#'
#' @description Setup cache directory and parallel processing.
#'
#' @param cache_dir Directory for caching NicheNet networks and other data.
#'   Defaults to the R user cache directory.
#' @param parallel_backend Parallelization backend: "auto", "serial",
#'   "multicore", or "snow"
#' @param n_cores Number of cores for parallel processing. NULL = auto-detect.
#'
#' @return A list of configuration settings (invisibly)
#'
#' @examples
#' \dontrun{
#' # Basic setup
#' setup_contextnet()
#'
#' # Custom cache directory
#' setup_contextnet(cache_dir = "~/.contextnet_cache")
#'
#' # Use 4 cores for parallel processing
#' setup_contextnet(n_cores = 4)
#' }
#'
#' @export
setup_contextnet <- function(cache_dir = NULL,
                              parallel_backend = c("auto", "serial", "multicore", "snow"),
                              n_cores = NULL) {
    parallel_backend <- match.arg(parallel_backend)

    # Cache directory
    if (is.null(cache_dir)) {
        cache_dir <- tools::R_user_dir("ContextNet", "cache")
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    options(contextnet.cache_dir = cache_dir)

    # Parallel backend
    if (parallel_backend == "auto") {
        parallel_backend <- if (.Platform$OS.type == "unix") "multicore" else "snow"
    }
    options(contextnet.parallel_backend = parallel_backend)

    # Number of cores
    if (is.null(n_cores)) {
        n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    options(contextnet.n_cores = as.integer(n_cores))

    cli::cli_alert_success("ContextNet configured")
    cli::cli_alert_info("Cache: {.path {cache_dir}}")
    cli::cli_alert_info("Parallel: {parallel_backend} ({n_cores} cores)")

    invisible(list(
        cache_dir = cache_dir,
        parallel_backend = parallel_backend,
        n_cores = n_cores
    ))
}


#' Get current ContextNet configuration
#'
#' @return A list of current configuration settings
#'
#' @examples
#' get_config()
#'
#' @export
get_config <- function() {
    list(
        cache_dir = getOption("contextnet.cache_dir",
                               tools::R_user_dir("ContextNet", "cache")),
        parallel_backend = getOption("contextnet.parallel_backend", "serial"),
        n_cores = getOption("contextnet.n_cores", 1L),
        version = as.character(utils::packageVersion("ContextNet"))
    )
}


#' Reset ContextNet to default configuration
#'
#' @return TRUE invisibly
#'
#' @export
reset_config <- function() {
    options(
        contextnet.cache_dir = NULL,
        contextnet.parallel_backend = "serial",
        contextnet.n_cores = 1L
    )
    cli::cli_alert_success("ContextNet configuration reset to defaults")
    invisible(TRUE)
}


#' Check ContextNet ecosystem installation
#'
#' @description Check which components of the ContextNet ecosystem are
#'   installed and available.
#'
#' @param verbose Logical, print status messages (default: TRUE)
#'
#' @return A list of installation status for each component (invisibly)
#'
#' @examples
#' check_ecosystem()
#'
#' @export
check_ecosystem <- function(verbose = TRUE) {
    status <- list()

    # ContextNet (always TRUE if this function is called)
    pkg_version <- as.character(utils::packageVersion("ContextNet"))
    status$contextnet <- list(installed = TRUE, version = pkg_version)
    if (verbose) {
        cli::cli_alert_success("ContextNet {pkg_version}")
    }

    # NicheNet (required for analysis)
    status$nichenetr <- list(
        installed = requireNamespace("nichenetr", quietly = TRUE)
    )
    if (verbose) {
        if (status$nichenetr$installed) {
            nn_version <- as.character(utils::packageVersion("nichenetr"))
            cli::cli_alert_success("nichenetr {nn_version}")
        } else {
            cli::cli_alert_warning("nichenetr not found (required for analysis)")
            cli::cli_alert_info("Install: {.code devtools::install_github('saeyslab/nichenetr')}")
        }
    }

    # MultiNicheNet (optional)
    status$multinichenetr <- list(
        installed = requireNamespace("multinichenetr", quietly = TRUE)
    )
    if (verbose) {
        if (status$multinichenetr$installed) {
            cli::cli_alert_success("multinichenetr installed")
        } else {
            cli::cli_alert_info("multinichenetr not found (optional, for multi-condition)")
        }
    }

    # Seurat (for expression data)
    status$seurat <- list(
        installed = requireNamespace("Seurat", quietly = TRUE)
    )
    if (verbose) {
        if (status$seurat$installed) {
            seurat_version <- as.character(utils::packageVersion("Seurat"))
            cli::cli_alert_success("Seurat {seurat_version}")
        } else {
            cli::cli_alert_warning("Seurat not found")
        }
    }

    # reticulate (optional, for Python integration)
    status$reticulate <- list(
        installed = requireNamespace("reticulate", quietly = TRUE)
    )
    if (verbose) {
        if (status$reticulate$installed) {
            cli::cli_alert_success("reticulate installed (Python integration available)")
        } else {
            cli::cli_alert_info("reticulate not found (optional, for Python plotting)")
        }
    }

    # Check for snakemake (pipeline)
    snakemake_path <- Sys.which("snakemake")
    status$pipeline <- list(installed = nzchar(snakemake_path))
    if (verbose) {
        if (status$pipeline$installed) {
            cli::cli_alert_success("snakemake found (GRN pipeline available)")
        } else {
            cli::cli_alert_info("contextnet-pipeline not found (optional, for GRN construction)")
        }
    }

    # ccflowR (optional, for sankey flow visualization)
    status$ccflowR <- list(
        installed = requireNamespace("ccflowR", quietly = TRUE)
    )
    if (verbose) {
        if (status$ccflowR$installed) {
            ccflowr_version <- as.character(utils::packageVersion("ccflowR"))
            cli::cli_alert_success("ccflowR {ccflowr_version} (sankey visualization available)")
        } else {
            cli::cli_alert_info("ccflowR not found (optional, for sankey flow visualization)")
            cli::cli_alert_info("Install: {.code devtools::install_github('username/ccflowR')}")
        }
    }

    invisible(status)
}


