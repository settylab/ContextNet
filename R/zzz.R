#' @importFrom rlang %||%
#' @importFrom methods setClass setMethod setGeneric setValidity new
#' @importFrom cli cli_abort cli_alert cli_alert_danger cli_alert_info
#'   cli_alert_success cli_alert_warning cli_warn cli_h1 cli_h2 cli_dl cli_text
#'   cli_progress_step
NULL


# Package-level environment for storing state
.contextnet_env <- new.env(parent = emptyenv())


#' @noRd
.onLoad <- function(libname, pkgname) {
    # Initialize default options
    op <- options()
    op_contextnet <- list(
        contextnet.cache_dir = NULL,
        contextnet.parallel_backend = "serial",
        contextnet.n_cores = 1L,
        contextnet.log_level = "info",
        contextnet.log_file = NULL
    )

    # Only set options that aren't already set
    toset <- !(names(op_contextnet) %in% names(op))
    if (any(toset)) {
        options(op_contextnet[toset])
    }

    invisible()
}


#' @noRd
.onAttach <- function(libname, pkgname) {
    version <- utils::packageVersion(pkgname)
    packageStartupMessage(sprintf("ContextNet version %s", version))

    # Check for crucial dependencies
    missing <- .check_crucial_deps()
    if (length(missing) > 0) {
        packageStartupMessage(
            sprintf(
                "Note: Missing recommended packages: %s",
                paste(missing, collapse = ", ")
            )
        )
        packageStartupMessage(
            "Install with: ContextNet::install_dependencies()"
        )
    } else {
        packageStartupMessage("For help, run: ?ContextNet or check_ecosystem()")
    }
}


#' Check for crucial dependencies
#' @noRd
.check_crucial_deps <- function() {
    crucial <- c("nichenetr", "multinichenetr")
    missing <- crucial[!vapply(crucial, requireNamespace, logical(1),
                               quietly = TRUE)]
    missing
}


#' Install ContextNet Dependencies
#'
#' Installs crucial dependencies for ContextNet that are not available on CRAN,
#' including nichenetr and multinichenetr from GitHub.
#'
#' @param dependencies Character vector of package names to install. Default
#'   installs all crucial dependencies: nichenetr and multinichenetr.
#' @param upgrade Whether to upgrade packages. Passed to
#'   \code{remotes::install_github}. Default is "never" to avoid breaking
#'   existing installations.
#' @param force Logical. If TRUE, force reinstallation even if packages are
#'   already installed. Default is FALSE.
#' @param ... Additional arguments passed to \code{remotes::install_github}.
#'
#' @return Invisibly returns a character vector of installed packages.
#'
#' @examples
#' \dontrun{
#' # Install all crucial dependencies
#' install_dependencies()
#'
#' # Install only nichenetr
#' install_dependencies("nichenetr")
#'
#' # Force reinstallation
#' install_dependencies(force = TRUE)
#' }
#'
#' @export
install_dependencies <- function(
    dependencies = c("nichenetr", "multinichenetr"),
    upgrade = "never",
    force = FALSE,
    ...
) {
    # Map package names to GitHub repos
    github_repos <- list(
        nichenetr = "saeyslab/nichenetr",
        multinichenetr = "saeyslab/multinichenetr"
    )

    # Validate dependencies
    unknown <- setdiff(dependencies, names(github_repos))
    if (length(unknown) > 0) {
        cli::cli_abort(
            "Unknown dependencies: {.val {unknown}}. ",
            "Available: {.val {names(github_repos)}}"
        )
    }

    # Check which need installation
    if (!force) {
        already_installed <- vapply(
            dependencies,
            requireNamespace,
            logical(1),
            quietly = TRUE
        )
        to_install <- dependencies[!already_installed]

        if (length(to_install) == 0) {
            cli::cli_alert_success("All dependencies are already installed.")
            return(invisible(character(0)))
        }
    } else {
        to_install <- dependencies
    }

    # Check for remotes package
    if (!requireNamespace("remotes", quietly = TRUE)) {
        cli::cli_alert_info("Installing {.pkg remotes} package...")
        utils::install.packages("remotes")
    }

    # Install each dependency
    installed <- character(0)
    for (pkg in to_install) {
        repo <- github_repos[[pkg]]
        cli::cli_alert_info("Installing {.pkg {pkg}} from GitHub ({repo})...")

        tryCatch(
            {
                remotes::install_github(repo, upgrade = upgrade, ...)
                cli::cli_alert_success("Successfully installed {.pkg {pkg}}")
                installed <- c(installed, pkg)
            },
            error = function(e) {
                cli::cli_alert_danger(
                    "Failed to install {.pkg {pkg}}: {e$message}"
                )
            }
        )
    }

    if (length(installed) > 0) {
        cli::cli_alert_success(
            "Installed {length(installed)} package{?s}: {.pkg {installed}}"
        )
    }

    invisible(installed)
}
