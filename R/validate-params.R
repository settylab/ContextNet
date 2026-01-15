#' Parameter Validation Functions
#'
#' @description
#' Structured validation functions for user-facing function parameters.
#' These functions provide informative error messages with suggestions.
#'
#' @name validate-params
#' @keywords internal
NULL


#' Validate cell type exists in object
#'
#' @param object A Seurat object or similar
#' @param cell_type Character string of cell type to validate
#' @param group_col Column name containing cell type annotations
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @export
validate_cell_type <- function(object, cell_type, group_col,
                                call = rlang::caller_env()) {
    # Check if group_col exists
    if (!group_col %in% colnames(object[[]])) {
        cli::cli_abort(c(
            "Column {.field {group_col}} not found in object metadata",
            "i" = "Available columns: {.val {colnames(object[[]])}}"
        ), call = call)
    }

    # Get available cell types
    available <- unique(object[[group_col, drop = TRUE]])

    # Check if cell type exists
    if (!cell_type %in% available) {
        cli::cli_abort(c(
            "Cell type {.val {cell_type}} not found in column {.field {group_col}}",
            "i" = "Available cell types: {.val {available}}"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate network integration mode
#'
#' @param mode Character string specifying integration mode
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @export
validate_integration_mode <- function(mode, call = rlang::caller_env()) {
    valid_modes <- c("weighted", "replace", "augment", "hierarchical")

    if (!mode %in% valid_modes) {
        cli::cli_abort(c(
            "Invalid integration mode: {.val {mode}}",
            "i" = "Valid modes: {.val {valid_modes}}"
        ), call = call)
    }

    invisible(TRUE)
}


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
    # Check type
    if (!is.numeric(weights)) {
        cli::cli_abort(c(
            "Weights must be numeric",
            "x" = "Got: {.cls {class(weights)}}"
        ), call = call)
    }

    # Check range
    if (any(weights < 0) || any(weights > 1)) {
        cli::cli_abort(c(
            "Weights must be in range [0, 1]",
            "x" = "Got range: [{min(weights)}, {max(weights)}]"
        ), call = call)
    }

    # Normalize if needed
    if (normalize && !isTRUE(all.equal(sum(weights), 1))) {
        cli::cli_warn(c(
            "Weights do not sum to 1",
            "i" = "Normalizing automatically"
        ))
        weights <- weights / sum(weights)
    }

    invisible(weights)
}


#' Validate condition parameters
#'
#' @param object A Seurat object or similar
#' @param condition_col Column name for condition grouping
#' @param condition_oi Condition of interest
#' @param condition_ref Reference condition
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_conditions <- function(object, condition_col, condition_oi, condition_ref,
                                 call = rlang::caller_env()) {
    # Check if condition_col exists
    if (!condition_col %in% colnames(object[[]])) {
        cli::cli_abort(c(
            "Column {.field {condition_col}} not found in object metadata",
            "i" = "Available columns: {.val {colnames(object[[]])}}"
        ), call = call)
    }

    # Get available conditions
    available <- unique(object[[condition_col, drop = TRUE]])

    # Check condition_oi
    if (!condition_oi %in% available) {
        cli::cli_abort(c(
            "Condition {.val {condition_oi}} not found in column {.field {condition_col}}",
            "i" = "Available conditions: {.val {available}}"
        ), call = call)
    }

    # Check condition_ref
    if (!condition_ref %in% available) {
        cli::cli_abort(c(
            "Reference condition {.val {condition_ref}} not found in column {.field {condition_col}}",
            "i" = "Available conditions: {.val {available}}"
        ), call = call)
    }

    # Check they're different
    if (condition_oi == condition_ref) {
        cli::cli_abort(c(
            "Condition of interest and reference condition must be different",
            "x" = "Both are: {.val {condition_oi}}"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate geneset_oi parameter
#'
#' @param geneset_oi Character vector of genes of interest
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_geneset_oi <- function(geneset_oi, call = rlang::caller_env()) {
    if (missing(geneset_oi) || is.null(geneset_oi)) {
        cli::cli_abort(c(
            "{.arg geneset_oi} is required but was not provided",
            "x" = "You must provide a character vector of genes of interest",
            "i" = "These are typically differentially expressed genes in the receiver",
            "i" = "Example: {.code geneset_oi = rownames(de_results)[de_results$padj < 0.05]}"
        ), call = call)
    }

    if (!is.character(geneset_oi)) {
        cli::cli_abort(c(
            "{.arg geneset_oi} must be a character vector",
            "x" = "Got: {.cls {class(geneset_oi)}}"
        ), call = call)
    }

    if (length(geneset_oi) == 0) {
        cli::cli_abort(c(
            "{.arg geneset_oi} is empty",
            "x" = "You must provide at least one gene of interest"
        ), call = call)
    }

    invisible(TRUE)
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

    if (must_exist && type == "file" && !file.info(path)$isdir == FALSE) {
        cli::cli_abort(c(
            "Expected a file but got a directory: {.path {path}}"
        ), call = call)
    }

    if (must_exist && type == "directory" && !file.info(path)$isdir == TRUE) {
        cli::cli_abort(c(
            "Expected a directory but got a file: {.path {path}}"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate positive integer
#'
#' @param x Value to validate
#' @param name Name of the parameter for error messages
#' @param min_value Minimum allowed value (default: 1)
#' @param max_value Maximum allowed value (default: Inf)
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_positive_integer <- function(x, name, min_value = 1, max_value = Inf,
                                       call = rlang::caller_env()) {
    if (!is.numeric(x) || length(x) != 1) {
        cli::cli_abort(c(
            "{.arg {name}} must be a single number",
            "x" = "Got: {.cls {class(x)}} of length {length(x)}"
        ), call = call)
    }

    if (x != as.integer(x)) {
        cli::cli_abort(c(
            "{.arg {name}} must be an integer",
            "x" = "Got: {x}"
        ), call = call)
    }

    if (x < min_value || x > max_value) {
        cli::cli_abort(c(
            "{.arg {name}} must be between {min_value} and {max_value}",
            "x" = "Got: {x}"
        ), call = call)
    }

    invisible(TRUE)
}
