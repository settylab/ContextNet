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
    op <- options()
    op_contextnet <- list(
        contextnet.cache_dir = NULL,
        contextnet.parallel_backend = "serial",
        contextnet.n_cores = 1L
    )

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

    missing <- .check_crucial_deps()
    if (length(missing) > 0) {
        packageStartupMessage(
            sprintf(
                "Note: Missing recommended packages: %s",
                paste(missing, collapse = ", ")
            )
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
