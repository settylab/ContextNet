#' Parameter Validation Helpers
#'
#' @description
#' Validation helpers for user-facing function parameters.
#'
#' @name validate-params
#' @keywords internal
NULL


#' Validate weight parameters
#'
#' @param weights Numeric vector of weights
#' @param normalize Logical, whether to normalize weights to sum to 1
#' @param call The calling environment for error reporting
#'
#' @return The (possibly normalized) weights, or throws an error
#'
#' @export
validate_weights <- function(weights, normalize = TRUE, call = rlang::caller_env()) {
    if (!is.numeric(weights)) {
        cli::cli_abort(c(
            "Weights must be numeric",
            "x" = "Got: {.cls {class(weights)}}"
        ), call = call)
    }

    if (any(weights < 0) || any(weights > 1)) {
        cli::cli_abort(c(
            "Weights must be in range [0, 1]",
            "x" = "Got range: [{min(weights)}, {max(weights)}]"
        ), call = call)
    }

    if (normalize && !isTRUE(all.equal(sum(weights), 1))) {
        cli::cli_warn(c(
            "Weights do not sum to 1",
            "i" = "Normalizing automatically"
        ))
        weights <- weights / sum(weights)
    }

    invisible(weights)
}


#' Validate file path exists
#'
#' @param path File or directory path
#' @param must_exist Logical, whether path must exist
#' @param type Character: "file", "directory", or "any"
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_path <- function(path, must_exist = TRUE, type = "any",
                           call = rlang::caller_env()) {
    if (!is.character(path) || length(path) != 1) {
        cli::cli_abort(c(
            "Path must be a single character string",
            "x" = "Got: {.cls {class(path)}}"
        ), call = call)
    }

    if (must_exist && !file.exists(path)) {
        cli::cli_abort(c(
            "Path does not exist: {.path {path}}"
        ), call = call)
    }

    is_dir <- isTRUE(file.info(path)$isdir)

    if (must_exist && type == "file" && is_dir) {
        cli::cli_abort(c(
            "Expected a file but got a directory: {.path {path}}"
        ), call = call)
    }

    if (must_exist && type == "directory" && !is_dir) {
        cli::cli_abort(c(
            "Expected a directory but got a file: {.path {path}}"
        ), call = call)
    }

    invisible(TRUE)
}
