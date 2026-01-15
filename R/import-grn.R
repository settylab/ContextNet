#' GRN Import Functions
#'
#' @description
#' Functions for importing gene regulatory networks from various sources
#' including the contextnet-pipeline, SCENIC+, Pando, and table formats.
#'
#' @importFrom utils read.csv
#' @importFrom rlang .data
#'
#' @name import-grn
NULL


#' Import GRN from various sources
#'
#' @description Auto-detect format and import GRN from file or object.
#'
#' @param x Path to file/directory, or a data.frame/matrix
#' @param format Format of input: "auto", "pipeline", "table", "scenic", "pando"
#' @param ... Additional arguments passed to format-specific import functions
#'
#' @return A ContextGRN object
#'
#' @examples
#' \dontrun{
#' # Import from pipeline output
#' grn <- import_grn("path/to/pipeline/output/")
#'
#' # Import from table
#' grn <- import_grn(my_edges_df, format = "table")
#' }
#'
#' @export
import_grn <- function(x, format = "auto", ...) {
    # Auto-detect format
    if (format == "auto") {
        format <- .detect_grn_format(x)
        .log_debug("Auto-detected format: {format}")
    }

    # Dispatch to appropriate importer
    result <- switch(format,
        pipeline = import_grn_pipeline(x, ...),
        table = import_grn_table(x, ...),
        scenic = import_grn_scenic(x, ...),
        pando = import_grn_pando(x, ...),
        cli::cli_abort(c(
            "Unknown GRN format: {.val {format}}",
            "i" = "Valid formats: pipeline, table, scenic, pando"
        ))
    )

    result
}


#' Import GRN from gene x TF matrix (CSV format)
#'
#' @description Import a GRN from a gene x TF matrix stored as CSV.
#'   This is the format produced by the contextnet-pipeline.
#'   The matrix has genes as rows and TFs as columns, with regulatory
#'   weights as values.
#'
#' @param path Path to the gene_x_tf.csv file
#' @param normalize Logical, whether to normalize weights per TF (default: TRUE)
#' @param min_weight Minimum weight threshold to include edge (default: 0)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A ContextGRN object
#'
#' @details
#' The normalization step scales weights so that the maximum value for each

#' TF (column) is 1. This ensures comparable regulatory strengths across TFs.
#'
#' @examples
#' \dontrun{
#' grn <- import_grn_pipeline("path/to/gene_x_tf.csv")
#' }
#'
#' @export
import_grn_pipeline <- function(path,
                                 normalize = TRUE,
                                 min_weight = 0,
                                 verbose = TRUE) {
    # ── Validate path ──
    if (!file.exists(path)) {
        cli::cli_abort(c(
            "File not found: {.path {path}}",
            "i" = "Expected a gene_x_tf.csv file from contextnet-pipeline"
        ))
    }

    if (verbose) cli::cli_alert("Loading GRN from: {.path {path}}")

    # ── Load gene_x_target matrix ──
    gene_x_tf <- tryCatch({
        read.csv(path, row.names = 1, check.names = FALSE)
    }, error = function(e) {
        cli::cli_abort(c(
            "Failed to read CSV file: {.path {path}}",
            "x" = conditionMessage(e)
        ))
    })

    # ── Validate matrix structure ──
    if (nrow(gene_x_tf) == 0 || ncol(gene_x_tf) == 0) {
        cli::cli_abort("gene_x_tf matrix is empty")
    }

    if (verbose) {
        cli::cli_alert_info("Matrix dimensions: {nrow(gene_x_tf)} genes x {ncol(gene_x_tf)} TFs")
    }

    # ── Normalize weights per TF (column) ──
    # Logic preserved from notebook: scale so max value in each column is 1
    if (normalize) {
        if (verbose) cli::cli_alert("Normalizing weights per TF...")

        gene_x_tf_norm <- gene_x_tf
        for (colname in colnames(gene_x_tf_norm)) {
            max_val <- max(gene_x_tf_norm[[colname]], na.rm = TRUE)
            if (!is.na(max_val) && max_val != 0) {
                gene_x_tf_norm[[colname]] <- gene_x_tf_norm[[colname]] / max_val
            }
        }
        gene_x_tf <- gene_x_tf_norm
    }

    # ── Convert matrix to edge list ──
    # Logic preserved from notebook: melt to from/to/weight format
    if (verbose) cli::cli_alert("Converting to edge list...")

    gene_x_tf$to <- rownames(gene_x_tf)
    tf_target <- reshape2::melt(
        gene_x_tf,
        id.vars = "to",
        variable.name = "from",
        value.name = "weight"
    )

    # Reorder columns and ensure from is character
    tf_target <- tf_target[, c("from", "to", "weight")]
    tf_target$from <- as.character(tf_target$from)

    # ── Apply weight threshold ──
    if (min_weight > 0) {
        n_before <- nrow(tf_target)
        tf_target <- tf_target[tf_target$weight >= min_weight, ]
        n_after <- nrow(tf_target)
        if (verbose) {
            cli::cli_alert_info("Filtered {n_before - n_after} edges below threshold {min_weight}")
        }
    }

    # ── Remove zero-weight edges ──
    tf_target <- tf_target[tf_target$weight > 0, ]

    # ── Report statistics ──
    n_tfs <- length(unique(tf_target$from))
    n_targets <- length(unique(tf_target$to))
    n_edges <- nrow(tf_target)

    if (verbose) {
        cli::cli_alert_success("Imported GRN:")
        cli::cli_dl(c(
            "TFs" = "{.val {n_tfs}}",
            "Targets" = "{.val {n_targets}}",
            "Edges" = "{.val {format(n_edges, big.mark = ',')}}"
        ))
    }

    # ── Create ContextGRN object ──
    grn <- new("ContextGRN",
        edges = tf_target,
        nodes = data.frame(
            gene = unique(c(tf_target$from, tf_target$to)),
            stringsAsFactors = FALSE
        ),
        regulons = .build_regulons(tf_target),
        peak_gene_links = data.frame(),
        tf_binding = data.frame(),
        metadata = list(
            source_file = path,
            normalized = normalize,
            min_weight = min_weight
        ),
        provenance = .init_provenance(
            source = "contextnet-pipeline",
            source_version = NA_character_
        )
    )

    grn
}


#' Import GRN from edge table
#'
#' @description Import a GRN from a data.frame or matrix containing edges.
#'
#' @param x A data.frame with columns for TF, target, and weight
#' @param from_col Name of column containing TF/regulator names
#' @param to_col Name of column containing target gene names
#' @param weight_col Name of column containing edge weights
#' @param normalize Logical, whether to normalize weights (default: FALSE)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A ContextGRN object
#'
#' @examples
#' \dontrun{
#' edges <- data.frame(
#'     TF = c("TF1", "TF1", "TF2"),
#'     target = c("Gene1", "Gene2", "Gene1"),
#'     score = c(0.8, 0.6, 0.7)
#' )
#' grn <- import_grn_table(edges, from_col = "TF", to_col = "target", weight_col = "score")
#' }
#'
#' @export
import_grn_table <- function(x,
                              from_col = NULL,
                              to_col = NULL,
                              weight_col = NULL,
                              normalize = FALSE,
                              verbose = TRUE) {
    # ── Validate input ──
    if (!is.data.frame(x)) {
        cli::cli_abort("{.arg x} must be a data.frame")
    }

    if (nrow(x) == 0) {
        cli::cli_abort("Input data.frame is empty")
    }

    # ── Auto-detect column names if not provided ──
    if (is.null(from_col) || is.null(to_col) || is.null(weight_col)) {
        detected <- .detect_edge_columns(x)
        from_col <- from_col %||% detected$from
        to_col <- to_col %||% detected$to
        weight_col <- weight_col %||% detected$weight

        if (verbose) {
            cli::cli_alert_info("Using columns: from={.val {from_col}}, to={.val {to_col}}, weight={.val {weight_col}}")
        }
    }

    # ── Validate columns exist ──
    missing_cols <- setdiff(c(from_col, to_col, weight_col), colnames(x))
    if (length(missing_cols) > 0) {
        cli::cli_abort(c(
            "Column(s) not found: {.val {missing_cols}}",
            "i" = "Available columns: {.val {colnames(x)}}"
        ))
    }

    # ── Extract and rename columns ──
    edges <- data.frame(
        from = as.character(x[[from_col]]),
        to = as.character(x[[to_col]]),
        weight = as.numeric(x[[weight_col]]),
        stringsAsFactors = FALSE
    )

    # ── Validate weights ──
    if (any(is.na(edges$weight))) {
        n_na <- sum(is.na(edges$weight))
        cli::cli_warn("Removing {n_na} edges with NA weights")
        edges <- edges[!is.na(edges$weight), ]
    }

    # ── Check weight range ──
    weight_range <- range(edges$weight)
    if (weight_range[1] < 0) {
        cli::cli_abort(c(
            "Negative weights found",
            "x" = "Weight range: [{weight_range[1]}, {weight_range[2]}]",
            "i" = "Weights must be non-negative"
        ))
    }

    # ── Normalize if requested ──
    if (normalize) {
        if (verbose) cli::cli_alert("Normalizing weights...")

        # Normalize per TF
        edges <- edges %>%
            dplyr::group_by(.data$from) %>%
            dplyr::mutate(weight = .data$weight / max(.data$weight, na.rm = TRUE)) %>%
            dplyr::ungroup() %>%
            as.data.frame()
    } else if (weight_range[2] > 1) {
        # Scale to [0, 1] if max > 1
        if (verbose) cli::cli_alert("Scaling weights to [0, 1]...")
        edges$weight <- edges$weight / weight_range[2]
    }

    # ── Create ContextGRN ──
    grn <- ContextGRN(
        edges = edges,
        metadata = list(
            from_col = from_col,
            to_col = to_col,
            weight_col = weight_col,
            normalized = normalize
        )
    )

    # Update provenance
    grn@provenance$source <- "table"

    if (verbose) {
        cli::cli_alert_success("Imported {nrow(edges)} edges from table")
    }

    grn
}


#' Import GRN from SCENIC+ output
#'
#' @description Import a GRN from SCENIC+ AnnData output.
#'
#' @param path Path to SCENIC+ output (h5ad file or directory)
#' @param ... Additional arguments (currently unused)
#'
#' @return A ContextGRN object
#'
#' @export
import_grn_scenic <- function(path, ...) {
    cli::cli_abort(c(
        "SCENIC+ import not yet implemented",
        "i" = "Use {.fn import_grn_table} with exported edge list"
    ))
}


#' Import GRN from Pando output
#'
#' @description Import a GRN from Pando output.
#'
#' @param path Path to Pando output
#' @param ... Additional arguments (currently unused)
#'
#' @return A ContextGRN object
#'
#' @export
import_grn_pando <- function(path, ...) {
    cli::cli_abort(c(
        "Pando import not yet implemented",
        "i" = "Use {.fn import_grn_table} with exported edge list"
    ))
}


#' Detect GRN format from input
#'
#' @param x Input path or object
#' @return Character string indicating format
#'
#' @keywords internal
.detect_grn_format <- function(x) {
    # Data.frame -> table format
    if (is.data.frame(x)) {
        return("table")
    }

    # Path-based detection
    if (is.character(x) && length(x) == 1) {
        if (dir.exists(x)) {
            # Directory - check for pipeline output
            if (file.exists(file.path(x, "gene_x_tf.csv"))) {
                return("pipeline")
            }
            if (file.exists(file.path(x, "rna_sc_ad.h5ad"))) {
                return("pipeline")
            }
        } else if (file.exists(x)) {
            # File - check extension
            ext <- tools::file_ext(x)
            if (ext == "csv") {
                return("pipeline")
            }
            if (ext == "h5ad") {
                # Could be SCENIC+ or pipeline
                return("scenic")
            }
            if (ext %in% c("rds", "RDS")) {
                return("table")
            }
        }
    }

    cli::cli_abort(c(
        "Cannot auto-detect GRN format",
        "i" = "Specify format explicitly: {.code import_grn(x, format = 'table')}"
    ))
}


#' Detect edge column names in a data.frame
#'
#' @param x A data.frame
#' @return A list with from, to, weight column names
#'
#' @keywords internal
.detect_edge_columns <- function(x) {
    cols <- colnames(x)
    cols_lower <- tolower(cols)

    # Patterns for each column type
    from_patterns <- c("from", "tf", "source", "regulator", "transcription_factor")
    to_patterns <- c("to", "target", "gene", "dest", "destination")
    weight_patterns <- c("weight", "score", "confidence", "correlation", "cor", "value")

    # Find matches
    from_match <- cols[cols_lower %in% from_patterns]
    to_match <- cols[cols_lower %in% to_patterns]
    weight_match <- cols[cols_lower %in% weight_patterns]

    if (length(from_match) == 0) {
        cli::cli_abort(c(
            "Cannot detect 'from' column",
            "i" = "Expected one of: {.val {from_patterns}}",
            "i" = "Found columns: {.val {cols}}"
        ))
    }

    if (length(to_match) == 0) {
        cli::cli_abort(c(
            "Cannot detect 'to' column",
            "i" = "Expected one of: {.val {to_patterns}}",
            "i" = "Found columns: {.val {cols}}"
        ))
    }

    if (length(weight_match) == 0) {
        cli::cli_abort(c(
            "Cannot detect 'weight' column",
            "i" = "Expected one of: {.val {weight_patterns}}",
            "i" = "Found columns: {.val {cols}}"
        ))
    }

    list(
        from = from_match[1],
        to = to_match[1],
        weight = weight_match[1]
    )
}
