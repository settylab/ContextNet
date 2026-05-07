#' Upstream Output Validation
#'
#' @description
#' Helper for validating outputs from upstream packages (nichenetr,
#' multinichenetr) so silent API changes are caught at the boundary.
#'
#' @name validate-data
#' @keywords internal
NULL


#' Assert upstream package output has expected structure
#'
#' @description Validates that output from nichenetr/multinichenetr functions
#'   has the expected columns and type. Catches silent API changes that would
#'   otherwise produce wrong results downstream.
#'
#' @param x The output object to validate
#' @param required_cols Character vector of required column names (for data.frames)
#' @param required_elements Character vector of required list element names (for lists)
#' @param caller Name of the upstream function that produced x (for error messages)
#' @param pkg Name of the upstream package (default: auto-detected from caller)
#'
#' @return x invisibly if valid, otherwise throws an informative error
#'
#' @keywords internal
.assert_upstream_output <- function(x, required_cols = NULL,
                                     required_elements = NULL,
                                     caller = "unknown",
                                     pkg = NULL) {
    .pkg_version_str <- function(pkg_name) {
        tryCatch(
            as.character(utils::packageVersion(pkg_name)),
            error = function(e) "unknown"
        )
    }

    version_info <- character()
    for (p in c("nichenetr", "multinichenetr")) {
        if (requireNamespace(p, quietly = TRUE)) {
            version_info <- c(version_info,
                              paste0(p, " ", .pkg_version_str(p)))
        }
    }
    version_hint <- if (length(version_info) > 0) {
        paste0("Installed: ", paste(version_info, collapse = ", "))
    } else {
        "No nichenetr/multinichenetr installed"
    }

    if (!is.null(required_cols)) {
        if (!is.data.frame(x)) {
            cli::cli_abort(c(
                "Unexpected output from {.fn {caller}}",
                "x" = "Expected data.frame, got {.cls {class(x)}}",
                "i" = "This may indicate a {.pkg nichenetr}/{.pkg multinichenetr} API change",
                "i" = version_hint
            ))
        }
        missing_cols <- setdiff(required_cols, colnames(x))
        if (length(missing_cols) > 0) {
            cli::cli_abort(c(
                "Unexpected columns from {.fn {caller}}",
                "x" = "Missing: {.val {missing_cols}}",
                "i" = "Got columns: {.val {colnames(x)}}",
                "i" = "This may indicate a {.pkg nichenetr}/{.pkg multinichenetr} API change",
                "i" = version_hint
            ))
        }
    }

    if (!is.null(required_elements)) {
        if (!is.list(x)) {
            cli::cli_abort(c(
                "Unexpected output from {.fn {caller}}",
                "x" = "Expected list, got {.cls {class(x)}}",
                "i" = "This may indicate a {.pkg nichenetr}/{.pkg multinichenetr} API change",
                "i" = version_hint
            ))
        }
        missing_els <- setdiff(required_elements, names(x))
        if (length(missing_els) > 0) {
            cli::cli_abort(c(
                "Unexpected structure from {.fn {caller}}",
                "x" = "Missing elements: {.val {missing_els}}",
                "i" = "Got elements: {.val {names(x)}}",
                "i" = "This may indicate a {.pkg nichenetr}/{.pkg multinichenetr} API change",
                "i" = version_hint
            ))
        }
    }

    invisible(x)
}
