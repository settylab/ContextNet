#' Serialization Functions
#'
#' @description
#' Functions for saving and loading ContextNet objects with version tracking
#' and metadata preservation.
#'
#' @name serialize
NULL


#' Save a ContextNet object to file
#'
#' @description Save a ContextGRN, ContextNetwork, or ContextResult object
#'   to a file with version metadata for future compatibility checking.
#'
#' @param object A ContextBase-derived object (ContextGRN, ContextNetwork, or ContextResult)
#' @param file Path to output file (recommended extension: .cnet)
#' @param compress Compression type: TRUE (default, gzip), FALSE (none), or "xz", "bzip2"
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return The file path (invisibly)
#'
#' @details
#' The saved file includes:
#' - The complete object with all slots
#' - ContextNet package version at time of save
#' - Timestamp of save operation
#'
#' @examples
#' \dontrun{
#' grn <- import_grn("path/to/grn.csv")
#' save_contextnet(grn, "my_grn.cnet")
#'
#' network <- integrate_networks(grn)
#' save_contextnet(network, "my_network.cnet")
#'
#' result <- analyze_communication(seurat_obj, network, ...)
#' save_contextnet(result, "my_result.cnet")
#' }
#'
#' @seealso [load_contextnet()], [peek_contextnet()]
#'
#' @export
save_contextnet <- function(object,
                             file,
                             compress = TRUE,
                             verbose = TRUE) {
    # ── Validate inputs ──
    if (is.null(object)) {
        cli::cli_abort("{.arg object} is required and cannot be NULL")
    }

    if (!inherits(object, "ContextBase")) {
        cli::cli_abort(c(
            "{.arg object} must be a ContextNet object",
            "i" = "Expected: ContextGRN, ContextNetwork, or ContextResult",
            "x" = "Got: {.cls {class(object)}}"
        ))
    }

    if (is.null(file) || !nzchar(file)) {
        cli::cli_abort("{.arg file} is required and cannot be empty")
    }

    # ── Add version metadata ──
    attr(object, "contextnet_version") <- as.character(utils::packageVersion("ContextNet"))
    attr(object, "save_timestamp") <- Sys.time()
    attr(object, "r_version") <- as.character(getRversion())

    # ── Save object ──
    saveRDS(object, file, compress = compress)

    if (verbose) {
        obj_class <- class(object)[1]
        cli::cli_alert_success("Saved {.cls {obj_class}} to {.path {file}}")
    }

    invisible(file)
}


#' Load a ContextNet object from file
#'
#' @description Load a ContextGRN, ContextNetwork, or ContextResult object
#'   that was saved with save_contextnet().
#'
#' @param file Path to the saved file
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return The loaded ContextNet object
#'
#' @details
#' If the object was saved with a different version of ContextNet, a warning
#' is issued but the object is still loaded. The original version is preserved
#' in the object attributes.
#'
#' @examples
#' \dontrun{
#' # Load a saved GRN
#' grn <- load_contextnet("my_grn.cnet")
#'
#' # Load a saved network
#' network <- load_contextnet("my_network.cnet")
#' }
#'
#' @seealso [save_contextnet()], [peek_contextnet()]
#'
#' @export
load_contextnet <- function(file, verbose = TRUE) {
    # ── Validate inputs ──
    if (is.null(file) || !nzchar(file)) {
        cli::cli_abort("{.arg file} is required and cannot be empty")
    }

    if (!file.exists(file)) {
        cli::cli_abort(c(
            "File not found",
            "x" = "{.path {file}} does not exist"
        ))
    }

    # ── Load object ──
    object <- readRDS(file)

    # ── Check version compatibility ──
    saved_version <- attr(object, "contextnet_version")
    current_version <- as.character(utils::packageVersion("ContextNet"))

    if (!is.null(saved_version) && saved_version != current_version) {
        if (verbose) {
            cli::cli_warn(c(
                "Version mismatch detected",
                "i" = "Object was saved with ContextNet {saved_version}",
                "i" = "Current version is {current_version}",
                "!" = "Object may behave differently"
            ))
        }
    }

    # ── Validate object class ──
    if (!inherits(object, "ContextBase")) {
        cli::cli_abort(c(
            "Invalid ContextNet file",
            "x" = "Object is not a ContextNet class",
            "i" = "Got: {.cls {class(object)}}"
        ))
    }

    if (verbose) {
        obj_class <- class(object)[1]
        cli::cli_alert_success("Loaded {.cls {obj_class}} from {.path {file}}")
    }

    object
}


#' Peek at ContextNet file metadata
#'
#' @description Inspect metadata of a saved ContextNet file without loading
#'   the full object. Useful for checking file contents before loading large objects.
#'
#' @param file Path to the saved file
#'
#' @return A list containing metadata about the saved object:
#'   - class: The object class (ContextGRN, ContextNetwork, ContextResult)
#'   - version: ContextNet version used to save
#'   - save_timestamp: When the object was saved
#'   - Additional class-specific information
#'
#' @examples
#' \dontrun{
#' # Check what's in a file before loading
#' info <- peek_contextnet("my_result.cnet")
#' print(info$class)
#' print(info$version)
#' }
#'
#' @seealso [save_contextnet()], [load_contextnet()]
#'
#' @export
peek_contextnet <- function(file) {
    # ── Validate inputs ──
    if (is.null(file) || !nzchar(file)) {
        cli::cli_abort("{.arg file} is required and cannot be empty")
    }

    if (!file.exists(file)) {
        cli::cli_abort(c(
            "File not found",
            "x" = "{.path {file}} does not exist"
        ))
    }

    # ── Load object ──
    object <- readRDS(file)

    # ── Extract metadata ──
    info <- list(
        class = class(object)[1],
        version = attr(object, "contextnet_version") %||% "unknown",
        save_timestamp = attr(object, "save_timestamp"),
        r_version = attr(object, "r_version"),
        file_size = file.size(file)
    )

    # ── Add class-specific info ──
    if (inherits(object, "ContextGRN")) {
        info$n_edges <- nrow(object@edges)
        info$n_tfs <- length(unique(object@edges$from))
        info$n_targets <- length(unique(object@edges$to))
        info$has_regulons <- length(object@regulons) > 0
    } else if (inherits(object, "ContextNetwork")) {
        info$ltm_dims <- dim(object@ligand_target_matrix)
        info$n_ligands <- ncol(object@ligand_target_matrix)
        info$n_targets <- nrow(object@ligand_target_matrix)
        info$integration_mode <- object@integration_params$mode
        info$n_signaling_edges <- nrow(object@signaling_network)
        info$n_regulatory_edges <- nrow(object@regulatory_network)
    } else if (inherits(object, "ContextResult")) {
        info$analysis_type <- object@analysis_type
        info$n_ligand_activities <- nrow(object@ligand_activities)
        info$n_ligand_target_links <- nrow(object@ligand_target_links)
        info$top_ligands <- object@top_ligands
        info$n_top_ligands <- length(object@top_ligands)
        info$n_top_targets <- length(object@top_targets)
        info$sender <- object@parameters$sender
        info$receiver <- object@parameters$receiver
    }

    info
}
