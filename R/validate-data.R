#' Data Structure Validation Functions
#'
#' @description
#' Validation functions for data structures like GRN edges, AnnData keys,
#' and network matrices.
#'
#' @name validate-data
#' @keywords internal
NULL


#' Validate GRN edge data.frame structure
#'
#' @param edges A data.frame with edge information
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_grn_edges <- function(edges, call = rlang::caller_env()) {
    if (!is.data.frame(edges)) {
        cli::cli_abort(c(
            "Edges must be a data.frame",
            "x" = "Got: {.cls {class(edges)}}"
        ), call = call)
    }

    if (nrow(edges) == 0) {
        cli::cli_abort("Edges data.frame is empty", call = call)
    }

    required_cols <- c("from", "to", "weight")
    cols_lower <- tolower(colnames(edges))

    # Check for from column (or equivalent)
    from_patterns <- c("from", "tf", "source", "regulator")
    if (!any(cols_lower %in% from_patterns)) {
        cli::cli_abort(c(
            "Cannot find 'from' column in edges",
            "i" = "Expected one of: {.val {from_patterns}}",
            "i" = "Found columns: {.val {colnames(edges)}}"
        ), call = call)
    }

    # Check for to column (or equivalent)
    to_patterns <- c("to", "target", "gene", "dest")
    if (!any(cols_lower %in% to_patterns)) {
        cli::cli_abort(c(
            "Cannot find 'to' column in edges",
            "i" = "Expected one of: {.val {to_patterns}}",
            "i" = "Found columns: {.val {colnames(edges)}}"
        ), call = call)
    }

    # Check for weight column (or equivalent)
    weight_patterns <- c("weight", "score", "confidence", "correlation")
    if (!any(cols_lower %in% weight_patterns)) {
        cli::cli_abort(c(
            "Cannot find 'weight' column in edges",
            "i" = "Expected one of: {.val {weight_patterns}}",
            "i" = "Found columns: {.val {colnames(edges)}}"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate AnnData key exists
#'
#' @param adata_path Path to AnnData file
#' @param key Key to check for
#' @param slot Slot to check: "varm", "obsm", "uns"
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_anndata_key <- function(adata_path, key, slot = "varm",
                                  call = rlang::caller_env()) {
    # This function requires an h5ad reading backend
    available_keys <- .get_anndata_keys(adata_path, slot)

    if (!key %in% available_keys) {
        cli::cli_abort(c(
            "Key {.val {key}} not found in {.field {slot}}",
            "i" = "Available keys: {.val {available_keys}}"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate ligand-target matrix
#'
#' @param matrix A matrix of ligand-target regulatory potential
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_ligand_target_matrix <- function(matrix, call = rlang::caller_env()) {
    if (!is.matrix(matrix)) {
        cli::cli_abort(c(
            "Ligand-target matrix must be a matrix",
            "x" = "Got: {.cls {class(matrix)}}"
        ), call = call)
    }

    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
        cli::cli_abort("Ligand-target matrix is empty", call = call)
    }

    if (is.null(rownames(matrix))) {
        cli::cli_abort("Ligand-target matrix must have row names (target genes)", call = call)
    }

    if (is.null(colnames(matrix))) {
        cli::cli_abort("Ligand-target matrix must have column names (ligands)", call = call)
    }

    # Check for valid values
    if (any(is.na(matrix))) {
        n_na <- sum(is.na(matrix))
        cli::cli_warn(c(
            "Ligand-target matrix contains {n_na} NA values",
            "i" = "These will be treated as zero regulatory potential"
        ))
    }

    invisible(TRUE)
}


#' Validate Seurat object
#'
#' @param object An object to validate
#' @param require_assay Character vector of required assays
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_seurat <- function(object, require_assay = NULL,
                             call = rlang::caller_env()) {
    # Check if it's a Seurat object
    if (!inherits(object, "Seurat")) {
        cli::cli_abort(c(
            "Expected a Seurat object",
            "x" = "Got: {.cls {class(object)}}"
        ), call = call)
    }

    # Check required assays
    if (!is.null(require_assay)) {
        available_assays <- names(object@assays)
        missing <- setdiff(require_assay, available_assays)
        if (length(missing) > 0) {
            cli::cli_abort(c(
                "Required assay(s) not found in Seurat object",
                "x" = "Missing: {.val {missing}}",
                "i" = "Available: {.val {available_assays}}"
            ), call = call)
        }
    }

    invisible(TRUE)
}


#' Validate expression matrix
#'
#' @param expr Expression matrix (genes x cells or cells x genes)
#' @param genes_in_rows Logical, whether genes are in rows (default: TRUE)
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_expression_matrix <- function(expr, genes_in_rows = TRUE,
                                        call = rlang::caller_env()) {
    # Check type
    if (!is.matrix(expr) && !inherits(expr, "Matrix")) {
        cli::cli_abort(c(
            "Expression data must be a matrix",
            "x" = "Got: {.cls {class(expr)}}"
        ), call = call)
    }

    if (nrow(expr) == 0 || ncol(expr) == 0) {
        cli::cli_abort("Expression matrix is empty", call = call)
    }

    # Check for row/column names
    if (genes_in_rows) {
        if (is.null(rownames(expr))) {
            cli::cli_abort("Expression matrix must have row names (gene names)", call = call)
        }
    } else {
        if (is.null(colnames(expr))) {
            cli::cli_abort("Expression matrix must have column names (gene names)", call = call)
        }
    }

    invisible(TRUE)
}


#' Validate ContextGRN object
#'
#' @param grn A ContextGRN object
#' @param min_edges Minimum number of edges required
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_context_grn <- function(grn, min_edges = 1, call = rlang::caller_env()) {
    if (!inherits(grn, "ContextGRN")) {
        cli::cli_abort(c(
            "Expected a ContextGRN object",
            "x" = "Got: {.cls {class(grn)}}",
            "i" = "Create one with {.fn import_grn} or {.fn ContextGRN}"
        ), call = call)
    }

    n_edges <- nrow(grn@edges)
    if (n_edges < min_edges) {
        cli::cli_abort(c(
            "ContextGRN has too few edges",
            "x" = "Found {n_edges} edges, need at least {min_edges}",
            "i" = "Check that the GRN import completed successfully"
        ), call = call)
    }

    invisible(TRUE)
}


#' Validate ContextNetwork object
#'
#' @param network A ContextNetwork object
#' @param call The calling environment for error reporting
#'
#' @return TRUE invisibly if valid, otherwise throws an error
#'
#' @keywords internal
validate_context_network <- function(network, call = rlang::caller_env()) {
    if (!inherits(network, "ContextNetwork")) {
        cli::cli_abort(c(
            "Expected a ContextNetwork object",
            "x" = "Got: {.cls {class(network)}}",
            "i" = "Create one with {.fn integrate_networks}"
        ), call = call)
    }

    # Check ligand-target matrix is not empty
    if (ncol(network@ligand_target_matrix) == 0) {
        cli::cli_abort(c(
            "ContextNetwork has empty ligand-target matrix",
            "i" = "Check that network integration completed successfully"
        ), call = call)
    }

    invisible(TRUE)
}


#' Get keys from AnnData slot
#'
#' @param path Path to AnnData file
#' @param slot Slot to read: "varm", "obsm", "uns"
#'
#' @return Character vector of available keys
#'
#' @keywords internal
.get_anndata_keys <- function(path, slot = "varm") {
    # Try different backends
    if (requireNamespace("hdf5r", quietly = TRUE)) {
        h5 <- hdf5r::H5File$new(path, mode = "r")
        on.exit(h5$close_all())

        if (slot %in% names(h5)) {
            return(names(h5[[slot]]))
        }
    }

    # Fallback: return empty with warning
    cli::cli_warn(c(
        "Cannot read AnnData keys",
        "i" = "Install {.pkg hdf5r}: {.code install.packages('hdf5r')}"
    ))
    character()
}
